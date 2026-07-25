/**
# Chemistry solver
This header contains the implementation of the chemistry solver.
The general form of the equation is we solve an equation of the form:
$$\frac{dMs_{i}}{dt} = \sum^{NR}_j R_{ji} \nu_{ji}(1-\epsilon)$$
where y is the mass (or mass fractions) of species i, R is the reaction rate and nu is the 
stoichiometric coefficient of species i in reaction j.

We use the OpenSMOKE++ library to solve the ODE system as this offers
a wide variety of stiff solvers optimized for chemical kinetics.

We also have the option to use an explicit solver for non-stiff problems.
*/

#include "reactors.h"

extern scalar zeta;
extern scalar T;
extern scalar porosity;

#ifdef STORE_SOURCES
/**
Optional diagnostic: when `STORE_SOURCES` is defined, the full per-cell
`data.sources` vector (length NEQ) computed by the reactor is copied into the
user-provided field list `sourcesList`. */
extern scalar * sourcesList;
#endif

#ifndef TURN_OFF_REACTIONS

/**
## Explicit ODE solvers
These are simple implementations of explicit ODE solvers: Euler and
Runge-Kutta 4(5). They can be used for non-stiff problems.
*/
void ODESolverEXP (odefunction ode, unsigned int neq, double dt, double* y, void* args) {

  double dy[neq];
  ode(y, dt, dy, args);

  for (int jj=0; jj<neq; jj++)
    y[jj] += dt*dy[jj];
}

void RungeKutta45EXP(odefunction ode, unsigned int neq, double dt, double *y, void *args) {

  // Allocate arrays for the k values and temporary y values
  double k1[neq], k2[neq], k3[neq], k4[neq], k5[neq], k6[neq];
  double ytmp[neq];

  // Coefficients for the RK45 method
  const double a2 = 1.0 / 4.0;
  const double a3 = 3.0 / 8.0;
  const double a4 = 12.0 / 13.0;
  const double a6 = 1.0 / 2.0;

  const double b31 = 3.0 / 32.0;
  const double b32 = 9.0 / 32.0;

  const double b41 = 1932.0 / 2197.0;
  const double b42 = -7200.0 / 2197.0;
  const double b43 = 7296.0 / 2197.0;

  const double b51 = 439.0 / 216.0;
  const double b52 = -8.0;
  const double b53 = 3680.0 / 513.0;
  const double b54 = -845.0 / 4104.0;

  const double b61 = -8.0 / 27.0;
  const double b62 = 2.0;
  const double b63 = -3544.0 / 2565.0;
  const double b64 = 1859.0 / 4104.0;
  const double b65 = -11.0 / 40.0;

  // Coefficients for the 5th order solution
  const double c1 = 16.0 / 135.0;
  const double c3 = 6656.0 / 12825.0;
  const double c4 = 28561.0 / 56430.0;
  const double c5 = -9.0 / 50.0;
  const double c6 = 2.0 / 55.0;

  // Step 1: Calculate k1 = f(t, y)
  ode(y, 0, k1, args);

  // Step 2: Calculate k2 = f(t + a2*dt, y + a2*k1*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * a2 * k1[j];
  ode(ytmp, a2 * dt, k2, args);

  // Step 3: Calculate k3 = f(t + a3*dt, y + b31*k1*dt + b32*k2*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b31 * k1[j] + b32 * k2[j]);
  ode(ytmp, a3 * dt, k3, args);

  // Step 4: Calculate k4 = f(t + a4*dt, y + b41*k1*dt + b42*k2*dt + b43*k3*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b41 * k1[j] + b42 * k2[j] + b43 * k3[j]);
  ode(ytmp, a4 * dt, k4, args);

  // Step 5: Calculate k5 = f(t + a5*dt, y + b51*k1*dt + b52*k2*dt + b53*k3*dt + b54*k4*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b51 * k1[j] + b52 * k2[j] + b53 * k3[j] + b54 * k4[j]);
  ode(ytmp, dt, k5, args);

  // Step 6: Calculate k6 = f(t + a6*dt, y + b61*k1*dt + b62*k2*dt + b63*k3*dt + b64*k4*dt + b65*k5*dt)
  for (int j = 0; j < neq; j++)
    ytmp[j] = y[j] + dt * (b61 * k1[j] + b62 * k2[j] + b63 * k3[j] + b64 * k4[j] + b65 * k5[j]);
  ode(ytmp, a6 * dt, k6, args);

  // Update y using the 5th order solution
  for (int j = 0; j < neq; j++) {
    y[j] += dt * (c1 * k1[j] + c3 * k3[j] + c4 * k4[j] + c5 * k5[j] + c6 * k6[j]);
    y[j] = y[j] < 0 ? 0 : y[j]; // Ensure non-negativity
  }

  y[NGS+NSS] = clamp (y[NGS+NSS], 0., 1.); // Ensure boundness for porosity
}

event init (i = 0) {
  OpenSMOKE_InitODESolver ();
}

event cleanup (t = end) {
  OpenSMOKE_CleanODESolver ();
}

event reset_sources (i++) {
  foreach()
    omega[] = 0.;
}

#ifdef CHEMISTRY_LOG
scalar t_solid[], t_gas[];
#endif

#ifdef BINNING
/**
Scale the gas tracers (species + temperature) of the current cell by a common
factor. `1/(1-f)` maps the VOF-tracer form (`Y*(1-f)`) to the actual mass
fractions the reactor expects; `(1-f)` is the inverse. */

static void scale_gas_tracers (Point point, double factor) {
  for (int jj = 0; jj < NGS; jj++) {
    scalar YG = YGList_G[jj];
    YG[] *= factor;
  }
  TG[] *= factor;
}
#endif

#ifdef VARPROP
/**
Predict the gas-phase reaction source at the (post-integration) state `y0ode`
(species mass fractions in `[0..NGS-1]`, temperature in `[NGS]`) and accumulate
it into the divergence/energy source fields of the current cell. Shared by the
binning and per-cell paths. */

static void accumulate_gas_sources (Point point, double * y0ode) {
  UserDataODE data;
  data.P = Pref + p[];
  data.T = y0ode[NGS];
  double sources[NGS + 1];
  data.sources = sources;

  double dy_tmp[NGS + 1];
  gas_batch_nonisothermal_constantpressure (y0ode, dt, dy_tmp, &data);

  for (int jj = 0; jj < NGS; jj++) {
    scalar DYDtGjj = DYDtG_G[jj];
    DYDtGjj[] += sources[jj]*cm[];
  }
  DTDtG[] += sources[NGS]*cm[];
}
#endif

event chemistry (i++) {

#ifdef CHEMISTRY_LOG
  reset ({t_solid, t_gas}, 0.);
  double time_mpi[npe()];

  for (int pe = 0; pe < npe(); pe++)
    time_mpi[pe] = 0.;

  struct timespec start, end;
  clock_gettime (CLOCK_MONOTONIC, &start);
#endif

#ifdef SOLVE_TEMPERATURE
  odefunction batch = &solid_batch_nonisothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1 + 1; //NGS + NSS + porosity + T
#else
  odefunction batch = &solid_batch_isothermal_constantpressure;
  unsigned int NEQ = NGS + NSS + 1;
#endif
  /**
  ## Solid-gas reactions
  We solve the solid-gas reaction system in each cell where there is
  solid present (i.e. f > F_ERR). The system is solved in terms of mass
  because the volume of the solid phase is variable due to porosity changes.
  */
  foreach ()
    if (f[] > F_ERR) {
      double temperature = TS[]/f[];
      // Reject two FPE triggers before mutating state, both of which make the
      // gas-species mole-fraction conversion in the RHS divide by sum(y/MW)==0:
      //  - sliver-garbage temperature (TS/f outside a physical window);
      //  - an empty reactor seed (gasmass = YG/f * rhoGvh * porosity, built below).
      double rhoGvh_seed;
      #ifdef VARPROP
      rhoGvh_seed = rhoGv_S[];
      #else
      rhoGvh_seed = rhoG;
      #endif
      double ygsum_seed = 0.;
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        ygsum_seed += YG[];
      }
      if (!(temperature > 273.) || !(temperature < 3500.) ||
          !(ygsum_seed > 0.) || !(porosity[] > 0.) || !(rhoGvh_seed > 0.))
        continue;

      porosity[] /= f[];

      double y0ode[NEQ];
      UserDataODE data;
      data.P = Pref + p[];
#ifdef VARPROP
      data.rhos = rhoSv[];
      data.rhog = rhoGv_S[];
#else
      data.rhos = rhoS;
      data.rhog = rhoG;
#endif
      data.zeta = zeta[];
#ifdef SOLVE_TEMPERATURE
# ifdef VARPROP
      data.cps = cpSv[];
      data.cpg = cpGv_S[];
# else
      data.cps = cpS;
      data.cpg = cpG;
# endif
#endif
      double sources[NEQ];
#ifdef STORE_SOURCES
      for (int jj = 0; jj < NEQ; jj++) 
        sources[jj] = 0.; // solid-species slots are never written
#endif
      data.sources = NULL; // do not fill sources during integration; predict after the solve

      double gasmass[NGS];
      double rhoGvh;
      #ifdef VARPROP
      rhoGvh = rhoGv_S[];
      #else
      rhoGvh = rhoG;
      #endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        gasmass[jj] = YG[]/f[]*rhoGvh*porosity[];
        y0ode[jj] = gasmass[jj];
      }

      double solidmass[NSS];
      for (int jj = 0; jj < NSS; jj++) {
        scalar YS = YSList[jj];
        solidmass[jj] = YS[]/f[]*rhoSv[]*(1. - porosity[]);
        y0ode[jj+NGS] = solidmass[jj];
      }

      y0ode[NGS+NSS] = porosity[];

#ifdef SOLVE_TEMPERATURE
      y0ode[NGS+NSS+1] = TS[]/f[];
#endif

#ifdef EXPLICIT_REACTIONS
    // ODESolverEXP (batch, NEQ, dt, y0ode, &data);
      RungeKutta45EXP (batch, NEQ, dt, y0ode, &data);
#else //default
      OpenSMOKE_ODESolver (batch, NEQ, dt, y0ode, &data);
#endif

      /**
      The source term is predicted once, at the converged end-of-step state
      (exact as dt -> 0), rather than being accumulated from the solver's
      internal RHS evaluations. */

      data.sources = sources;
      double dy_tmp[NEQ];
      batch (y0ode, dt, dy_tmp, &data);

      double totgasmass = 0;
      for (int jj = 0; jj < NGS; jj++)
        totgasmass += fmax (0., y0ode[jj]);

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        YG[] = (totgasmass < 1e-8) ? 0. : fmax (0., y0ode[jj])/totgasmass*f[];
      }

      double totsolidmass = 0;
      for (int jj = 0; jj < NSS; jj++)
        totsolidmass += fmax (0., y0ode[jj+NGS]);

      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        YS[] = (totsolidmass < 1e-8) ? 0. : fmax (0., y0ode[jj+NGS])/totsolidmass*f[];
      }

      porosity[] = y0ode[NGS+NSS]*f[];

#ifdef VARPROP
      for (int jj = 0; jj < NGS; jj++) {
        scalar DYDtGjj = DYDtG_S[jj];
        DYDtGjj[] += sources[jj]*cm[];
      }
#endif

#ifdef SOLVE_TEMPERATURE
      TS[] = y0ode[NGS+NSS+1]*f[];
# ifdef VARPROP
      DTDtS[] += sources[NGS+NSS+1]*cm[];
# endif
#endif
      omega[] = sources[NGS+NSS];
#ifdef STORE_SOURCES
      for (int jj = 0; jj < NEQ; jj++) {
        scalar src = sourcesList[jj];
        src[] = sources[jj];
      }
#endif
    }

  /**
  ## Gas-phase reactions
  We solve the gas-phase reaction system in every cell that contains gas
  (i.e. f < 1 - F_ERR). The system is solved in terms of mass fraction.
  */

#ifdef BINNING
# ifndef VARPROP
#   error "BINNING requires VARPROP (it uses rhoGv_G/cpGv_G and the DYDtG_G/DTDtG sources)"
# endif

  /**
  The bin partitioning is driven by a case-provided list of thermochemical
  `targets` and a per-target tolerance `eps` (mixed-radix bin id, see
  binning.h). */

  extern scalar * targets;
  extern double * eps;

  /**
  Flag the pure-gas cells (the same set integrated by the non-binning branch)
  and convert their gas fields from VOF-tracer form (`Y*(1-f)`) to the actual
  mass fractions the reactor expects. */

  scalar gasmask[];
  foreach() {
    gasmask[] = (f[] < 1. - F_ERR && TG[] > 0.) ? 1. : 0.;
    if (gasmask[])
      scale_gas_tracers (point, 1./(1. - f[]));
  }

  /**
  Agglomerate the flagged cells into bins of similar thermochemical state and
  integrate the stiff chemistry ODE once per bin. `bin->phi[j]` holds the
  mass-averaged value of `fields[j]`: entries `[0..NGS-1]` are the gas species
  and entry `[NGS]` is the temperature. */

  scalar * fields = list_concat (YGList_G, {TG});

  BinTable * table = binning (fields, targets, eps, rhoGv_G, cpGv_G, gasmask);

#ifdef CHEMISTRY_LOG
  static FILE * fp = NULL;
  if (!fp) {
    char name[20];
    sprintf (name, "bin-%d", pid());
    fp = fopen (name, "w");
  }
  fprintf (fp, "%g %ld %ld\n", t, grid->n, binning_stats(table).nactive);
  fflush (fp);
#endif

  foreach_bin (table) {
    double y0ode[NGS + 1];
    for (size_t j = 0; j < bin->nfields; j++)
      y0ode[j] = bin->phi[j];

    UserDataODE data;
    data.P = Pref + bin_average (bin, p);
    data.sources = NULL;
    data.rhog = bin->rho;
    data.cpg = bin->cp;

    OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure,
        NGS + 1, dt, y0ode, &data);

    for (size_t j = 0; j < bin->nfields; j++)
      bin->phi[j] = (j < (size_t)NGS) ? fmax (0., y0ode[j]) : y0ode[j];

    bin->rho = data.rhog;
    bin->cp = data.cpg;
  }

  binning_remap (table, fields, rhoGv_G, cpGv_G);
  binning_cleanup (table);
  free (fields), fields = NULL;

  /**
  Predict the gas-phase source terms at the post-reaction state (one local
  evaluation per cell), then restore the VOF-tracer form of the gas fields. */

  foreach() {
    if (gasmask[]) {
      double y0ode[NGS + 1]; // NGS + T
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        y0ode[jj] = YG[];
      }
      y0ode[NGS] = TG[];

      accumulate_gas_sources (point, y0ode);
      scale_gas_tracers (point, 1. - f[]); // restore VOF-tracer form
    }
  }
#else // !BINNING

  foreach() {
    if (f[] < 1. - F_ERR) {
      double temperature = TG[]/(1. - f[]);
      if (!(temperature > 273.) || !(temperature < 3500.))
        continue;

      double y0ode[NGS + 1]; // NGS + T
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        y0ode[jj] = YG[]/(1. - f[]);
      }
      y0ode[NGS] = temperature;

      UserDataODE data;
      data.P = Pref + p[];
      data.T = y0ode[NGS];
      data.sources = NULL; // do not fill sources during integration; predict after the solve
# ifdef VARPROP
      data.rhog = rhoGv_G[];
      data.cpg = cpGv_G[];
# else
      data.rhog = rhoG;
      data.cpg = cpG;
# endif
      /**
        Using an explicit solver for gas-phase reactions is not
        recommended as they are usually stiff.
        */
      OpenSMOKE_ODESolver (&gas_batch_nonisothermal_constantpressure, NGS + 1, dt, y0ode, &data);

      /**
        The source term is predicted once, at the converged end-of-step state
        (exact as dt -> 0), rather than accumulated during integration. */

# ifdef VARPROP
      accumulate_gas_sources (point, y0ode);
# endif

      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        YG[] = (y0ode[jj] > 0.) ? y0ode[jj]*(1. - f[]) : 0.;
      }
      TG[] = y0ode[NGS]*(1. - f[]);
    }
  }
#endif // BINNING

#ifdef CHEMISTRY_LOG
  clock_gettime (CLOCK_MONOTONIC, &end);
  time_mpi[pid()] = (end.tv_sec - start.tv_sec) +
                    (end.tv_nsec - start.tv_nsec)*1e-9;
@if _MPI
  if (pid() == 0) {
    MPI_Reduce(MPI_IN_PLACE, time_mpi, npe(), MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(time_mpi, NULL, npe(), MPI_DOUBLE,
        MPI_SUM, 0, MPI_COMM_WORLD);
  }
@endif

  fprintf (stderr, "%g ", t);

  for (int pe = 0; pe < npe(); pe++)
    fprintf (stderr, "%g ", time_mpi[pe]);

  fprintf (stderr, "\n");
#endif
}
#endif // TURN_OFF_REACTIONS
