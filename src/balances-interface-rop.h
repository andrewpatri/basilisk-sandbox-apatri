#ifndef BALANCES
# define BALANCES
#endif

#include "fracface.h"
#include "intgrad.h"

extern int maxlevel;
extern scalar f, porosity;

double Ptest;
vector u_i[];

// ---- NOTES ----
// 1) use padding=1 in adapt_wavelet_leave_interface to avoid problems with the interface
// 2) adapt seems to cause some problems with the diffusive fluxes at the first iterations
// 3) total mass conservation is better than it seems, sum over the non inert species to 
//    reconstruct the total gas mass produced

struct MassBalances 
{
  //total mass variables
  double tot_sol_mass_start;  // Total initial solid mass
  double tot_gas_mass_start;  // Total initial gas mass
  double tot_sol_mass;        // Total current solid mass
  double tot_gas_mass;        // Total current gas mass
  double tot_mass_start;      // Total initial mass 
  double tot_mass;            // Total current mass
  double tot_gas_mass_intnow; // Total gas mass that has crossed the interface in the last iteration
  double tot_gas_mass_int;    // Total gas mass that has crossed the interface

  //gas mass variables
  double* gas_mass_intnow;
  double* gas_mass;
  double* gas_mass_int;
  double* gas_mass_start;

  //solid mass variables
  double* sol_mass;
  double* sol_mass_start;

  //file variables
  FILE* fb;
  char name[80];
  int print_iter;
  int inert_index; // index of the inert species in the gas phase
};

struct MassBalances mb = {0};
scalar rhot[];

static void interface_fluxes (Point point) {

  // Calculate interfacial area and normal vector
  coord n = facet_normal (point, fS, fsS), p;
  double alpha = plane_alpha (fS[], n);
  double area = plane_area_center (n, alpha, &p);
  normalize (&n);

  // Calculate interfacial density
  double rhoGInt;
  #ifdef VARPROP
  rhoGInt = interpolate(rhot, x+p.x*Delta, y+p.y*Delta);
  #else
  rhoGInt = rhoG;
  #endif

  // CONVECTIVE FLUXES
  // Calculate velocity at the interface
  coord ui = {interpolate (uf.x, x+p.x*Delta, y+p.y*Delta),
              interpolate (uf.y, x+p.y*Delta, y+p.y*Delta)};
  coord us = {interpolate (ubf.x, x+p.x*Delta, y+p.y*Delta),
              interpolate (ubf.y, x+p.x*Delta, y+p.y*Delta)};

  foreach_dimension() {
    mb.tot_gas_mass_intnow += rhoGInt*ui.x*n.x*area*Delta;
    mb.tot_gas_mass_intnow += -rhoGInt*us.x*n.x*area*Delta*cm[]; // gas phase left out due to the interface movement
    // mb.tot_gas_mass_intnow += rhoGInt*u.x[]*n.x*area*Delta*cm[]*dt;
  }

#ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++) {
    scalar YGInt = YGList_Int[jj];
    foreach_dimension() {
      mb.gas_mass_intnow[jj] += rhoGInt*ui.x*n.x*area*Delta*YGInt[];
      mb.gas_mass_intnow[jj] += -rhoGInt*us.x*n.x*area*Delta*cm[]*YGInt[]; // gas phase left out due to the interface movement
      // mb.gas_mass_intnow[jj] += rhoGInt*u.x[]*n.x*area*Delta*cm[]*YG[]/f[]*dt;
    }
  }
  
  // DIFFUSIVE FLUXES
  // We need to recalculate the fluxes since sGexp gets modified
  // by the diffusion solver
  double jG[NGS];
  for (int jj = 0; jj < NGS; jj++) {
    scalar DmixG = DmixGList_G[jj];
# ifdef MOLAR_DIFFUSION
    scalar XG = XGList_G[jj];
    scalar XGInt = XGList_Int[jj];
    double Gtrgrad = ebmgrad(point, XG, fS, fG, fsS, fsG, true, XGInt[], &success);
    jG[jj] = rhoGInt * DmixG[] * Gtrgrad * gas_MWs[jj] / MWmixG_G[];
# else
    scalar YG = YGList_G[jj];
    scalar YGInt = YGList_Int[jj];
    double Gtrgrad = ebmgrad(point, YG, fS, fG, fsS, fsG, true, YGInt[], &success);
    jG[jj] = rhoGInt * DmixG[] * Gtrgrad;
# endif
  }

  double jGtot = 0.;
# ifdef FICK_CORRECTED
  for (int jj = 0; jj < NGS; jj++)
    jGtot += jG[jj];
# endif

  for (int jj = 0; jj < NGS; jj++) {
    scalar YGInt = YGList_Int[jj];
    mb.gas_mass_intnow[jj] += (jG[jj] - jGtot*YGInt[])*area*Delta*cm[]; 
  }
#endif
}

static void compute_initial_state (void) {
  mb.tot_sol_mass_start = 0.;
  mb.tot_gas_mass_start = 0.;
  
  //compute initial solid mass
  foreach() {
    if (f[] > F_ERR) {
      mb.tot_sol_mass_start += (f[] -porosity[])*rhoS*dv();

      double rhoGh;
      #ifdef VARPROP
      rhoGh = rhoGv_S[];
      #else
      rhoGh = rhoG;
      #endif
      mb.tot_gas_mass_start += porosity[]*rhoGh*dv(); //ef
    }
  }

  #ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++)
    mb.gas_mass_start[jj] = mb.tot_gas_mass_start*gas_start[jj];
  for (int jj=0; jj<NSS; jj++)
    mb.sol_mass_start[jj] = mb.tot_sol_mass_start*sol_start[jj]; 
  #endif
}

static void compute_balances(void) {
  mb.tot_sol_mass = 0.;
  mb.tot_gas_mass = 0.;
  mb.tot_mass = 0.;
  mb.tot_gas_mass_intnow = 0.;

#ifdef MULTICOMPONENT
  for (int jj=0; jj<NGS; jj++) {
    mb.gas_mass_intnow[jj] = 0.;
    mb.gas_mass[jj] = 0.;
  }

  for (int jj=0; jj<NSS; jj++)
    mb.sol_mass[jj] = 0.;
#endif

  foreach()
    fS[] = f[];

  face_fraction (fS, fsS);

  #ifdef MULTICOMPONENT
  // We need to lose tracer form for YGList_G as it is used for the diffusive fluxes
  foreach() {
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG_G = YGList_G[jj];
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[]/(1. - f[]) : 0.; // gas phase
    }
  }
  #endif

  #ifdef VARPROP
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  #endif

    foreach()
    if (f[] > F_ERR && f[] < 1.-F_ERR)
        interface_fluxes (point);

    foreach() {
      if (f[] > F_ERR) {
        mb.tot_sol_mass += (f[] - porosity[]) * rhoS * dv(); //(f-ef)

#ifdef MULTICOMPONENT
        for (int jj = 0; jj < NSS; jj++) {
          scalar YS = YSList[jj];
          mb.sol_mass[jj] += (f[] - porosity[]) * rhoS * (YS[] / f[]) * dv(); // ef/f*(YG*f) = ef*YG
        }
#endif
      }
    }

#ifdef MULTICOMPONENT
  // Recover tracer form for YGList_G
  foreach() {
    for (int jj = 0; jj < NGS; jj++) {
      scalar YG_G = YGList_G[jj];
      YG_G[] = YG_G[]*(1. - f[]);
    }
  }
#endif
}

static void write_balances(void) {
  //print total mass
  fprintf(mb.fb, "%-18.12f %-18.12f %-18.12f %-18.12f ", t,
                                        mb.tot_sol_mass/mb.tot_sol_mass_start,
                                        mb.tot_gas_mass_intnow/mb.tot_sol_mass_start,
                                        mb.tot_mass);

  #ifdef MULTICOMPONENT
  //print individual gas species mass
  for (int jj=0; jj<NGS; jj++)
    fprintf(mb.fb, "%-18.12f", mb.gas_mass_intnow[jj]/mb.tot_sol_mass_start);
  //print individual solid species mass
  for (int jj=0; jj<NSS; jj++)
    fprintf(mb.fb, "%-18.12f", mb.sol_mass[jj]/mb.tot_sol_mass_start);
  #endif
  fprintf(mb.fb, "\n");
  fflush(mb.fb);
}

event defaults (i = 0) {
  mb.tot_gas_mass_intnow = 0.;
  mb.tot_gas_mass = 0.;
  mb.tot_gas_mass_int = 0.;
  mb.tot_sol_mass = 0.;
  mb.tot_sol_mass_start = 0.;
  mb.tot_gas_mass_start = 0.;
  mb.tot_mass_start = 0.;
  mb.tot_mass = 0.;

  mb.print_iter = 10;

  mb.gas_mass_start = NULL;
  mb.gas_mass_intnow = NULL;
  mb.gas_mass_int = NULL;
  mb.gas_mass = NULL;

  mb.sol_mass_start = NULL;
  mb.sol_mass = NULL;

  sprintf(mb.name, "balances-%d-rop", maxlevel);
}

event init (i = 0) { 
  #ifdef MULTICOMPONENT
  // mb.inert_index = OpenSMOKE_InertSpeciesIndex();
  mb.gas_mass_start =   (double*) malloc(NGS*sizeof(double));
  mb.gas_mass_int =     (double*) malloc(NGS*sizeof(double));
  mb.gas_mass_intnow =  (double*) malloc(NGS*sizeof(double));
  mb.gas_mass =         (double*) malloc(NGS*sizeof(double));

  for (int jj=0; jj<NGS; jj++) {
    mb.gas_mass_start[jj] = 0.;
    mb.gas_mass_int[jj] = 0.;
    mb.gas_mass_intnow[jj] = 0.;
    mb.gas_mass[jj] = 0.;
  }

  mb.sol_mass_start = (double*) malloc(NSS*sizeof(double));
  mb.sol_mass =       (double*) malloc(NSS*sizeof(double));

  for (int jj=0; jj<NSS; jj++) {
    mb.sol_mass_start[jj] = 0.;
    mb.sol_mass[jj] = 0.;
  }
  #endif

  mb.fb = fopen(mb.name, "w");
  fprintf(mb.fb, "%-18s %-18s %-18s %-18s ", "t(1)", "tot_sol_mass(2)", 
                                        "tot_gas_mass(3)", "tot_mass(4)");
#ifdef MULTICOMPONENT
  int counter = 5;
  //individual gas species mass header
  for (int jj=0; jj<NGS; jj++) {
    char buffer[80];
    sprintf(buffer, "mG_%s(%d)", OpenSMOKE_NamesOfSpecies(jj), counter++);
    fprintf(mb.fb, "%-18s", buffer);
  }

  //solid species mass header
  for (int jj=0; jj<NSS; jj++) {
    char buffer[80];
    sprintf(buffer, "mS_%s(%d)",  OpenSMOKE_NamesOfSolidSpecies(jj), counter++);
    fprintf(mb.fb, "%-18s", buffer);
  }
#endif
  
  fprintf(mb.fb, "\n");
  fflush(mb.fb);
}

event cleanup (t = end) {
  write_balances();
  fclose(mb.fb);
  free(mb.gas_mass_start);
  free(mb.gas_mass_int);
  free(mb.gas_mass_intnow);
  free(mb.gas_mass);
  free(mb.sol_mass_start);
  free(mb.sol_mass);
}

event reset_sources (i++) {
  // Compute initial state if not already done
  if (mb.tot_gas_mass_start == 0. && mb.tot_sol_mass_start == 0.)
    compute_initial_state();

  // Compute balances
  compute_balances();

  // Write balances to file
  if (i % mb.print_iter == 0 && pid() == 0)
    write_balances();
}