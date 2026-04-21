/**
## Update Properties

Update the thermodynamic properties for the multicomponent phase
change model, and compute the lagrangian derivative of the density,
which is used as a sorce term for the velocity divergence, to
describe low Mach compressibility effects. */

#ifdef VARPROP
#include "solid-thermal-conductivity.h"

scalar rhoGv_G0[], rhoGv_S0[];
extern scalar porosity;
scalar DTDtS[], DTDtG[];
scalar * DYDtG_G = NULL;
scalar * DYDtG_S = NULL;

trace
void update_properties (void) {

  foreach()
    rhoGv_S0[] = rhoGv_S[]*f[] + (1.-f[])*rhoGv_G[]; //field looks nicer done in one field

  // Reset all the properties fields
  reset ({rhoGv_S, rhoGv_G, rhoSv,
          muGv_S, muGv_G,
          lambdaGv_S, lambdaGv_G, lambdaSv,
          cpGv_S, cpGv_G, cpSv}, 0.);
  reset (DmixGList_S, 0.);
  reset (DmixGList_G, 0.);
  reset ({MWmixG_S, MWmixG_G}, 0.);

  foreach() {
    ThermoState tsGh, tsSh;
    double Diff_coeff[NGS];
    if (f[] > F_ERR && TS[] > 0.) {
      double xG[NGS], yG[NGS];
      double MWmixG;
      // Update internal gas properties
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[]/f[];
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, yG);
      MWmixG_S[] = MWmixG;

      tsGh.T = TS[]/f[];
      tsGh.P = Pref+p[];
      tsGh.x = xG;

      rhoGv_S[] = tpG.rhov (&tsGh);
      cpGv_S[] = tpG.cpv (&tsGh);
      lambdaGv_S[] = tpG.lambdav (&tsGh);
      muGv_S[] = tpG.muv (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);
#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_S[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for(int jj=0; jj<NGS; jj++) {
        scalar DmixGv = DmixGList_S[jj];
        #ifdef CONST_DIFF
        DmixGv[] = CONST_DIFF*pow(porosity[]/f[], 4./3.);
        #else
        DmixGv[] = Diff_coeff[jj]*pow(porosity[]/f[], 4./3.);
        #endif
      }
      
      // Update internal solid properties
      double xS[NSS], yS[NSS];
      double MWmixS;
      for (int jj=0; jj<NSS; jj++) {
        scalar YS = YSList[jj];
        yS[jj] = YS[]/f[];
      }
      OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (xS, &MWmixS, yS);

      tsSh.T = TS[]/f[];
      tsSh.P = Pref+p[];
      tsSh.x = xS;

      rhoSv[] = tpS.rhov (&tsSh);
      cpSv[] = tpS.cpv (&tsSh);
    
      coord lambda_pf = pseudo_phase_thermal_conductivity(point, lambdaGv_S[], 
                                                          porosity[]/f[], 
                                                          tsSh.T, 
                                                          f);
      foreach_dimension()
        lambda1v.x[] = lambda_pf.x;
    }

    if (f[] < 1. - F_ERR && TG[] > 0.) {
      // Update external gas properties
      double xG[NGS], yG[NGS];
      double MWmixG;
      for (int jj=0; jj<NGS; jj++) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[]/(1.-f[]);
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixG, yG);
      MWmixG_G[] = MWmixG;

      tsGh.T = TG[]/(1.-f[]);
      tsGh.P = Pref+p[];
      tsGh.x = xG;

      rhoGv_G[] = tpG.rhov (&tsGh);
      muGv_G[] = tpG.muv (&tsGh);
      cpGv_G[] = tpG.cpv (&tsGh);
      lambdaGv_G[] = tpG.lambdav (&tsGh);
      tpG.diff (&tsGh, Diff_coeff);

#ifdef MASS_DIFFUSION_ENTHALPY
      double cpG[NGS];
      tpG.cpvs (&tsGh, cpG);
      for(int jj=0; jj<NGS; jj++) {
        scalar cpGv = cpGList_G[jj];
        cpGv[] = cpG[jj];
      }
#endif // MASS_DIFFUSION_ENTHALPY

      for (int jj=0; jj<NGS; jj++) {
        scalar Dmix2v = DmixGList_G[jj];
# ifdef CONST_DIFF
        Dmix2v[] = CONST_DIFF;
# else
        Dmix2v[] = Diff_coeff[jj];
# endif
      }

      foreach_dimension()
        lambda2v.x[] = lambdaGv_G[];
    }
  }
}

//Should  done in the default event but is executed before OS++ initialization otherwise
event init (i = 0) {
  DYDtG_G = NULL;
  DYDtG_S = NULL;

  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "DYDtG_%s_G", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_G = list_append (DYDtG_G, a);
  }
  reset (DYDtG_G, 0.);
  
  for (int jj=0; jj<NGS; jj++) {
    scalar a = new scalar;
    free (a.name);
    char name[64];
    snprintf (name, sizeof(name), "DYDtG_%s_S", OpenSMOKE_NamesOfSpecies(jj));
    a.name = strdup (name);
    a.nodump = true;
    DYDtG_S = list_append (DYDtG_S, a);
  }
  reset (DYDtG_S, 0.);

  MWmixG_G.dirty = true;
  MWmixG_S.dirty = true;

#if TREE
  for (scalar s in {drhodt}) {
#if EMBED
    s.refine = s.prolongation = refine_embed_linear;
#else
    s.refine  = refine_linear;
#endif
    s.restriction = restriction_volume_average;
    s.dirty = true; // boundary conditions need to be updated
  }
#endif
}

event properties (i = 0) {
  update_properties();
}

event cleanup (t = end)
{
  delete (DYDtG_G), free (DYDtG_G), DYDtG_G = NULL;
  delete (DYDtG_S), free (DYDtG_S), DYDtG_S = NULL;
}

event reset_sources (i++) {
  foreach() {
    DTDtG[] = 0.;
    DTDtS[] = 0.;
  }

  reset (DYDtG_G, 0.);
  reset (DYDtG_S, 0.);
}

void update_divergence (void) {

//   // ENSURE THAT THE TRACER FORM IS LOST
//   /**
//   We define the variables used to compute the lagrangian derivative
//   on each level. */

  restriction ({T,TS,TG});
  restriction (YSList);
  restriction (YGList_G);
  restriction (YGList_S);
#ifdef MOLAR_DIFFUSION
  restriction (XGList_G);
  restriction (XGList_S);
#endif

//   /**
//   We calculate the Lagrangian derivative of the temperature fields. */

  face vector lambdagradTS[], lambdagradTG[];
  foreach_face() {
    lambdagradTS.x[] = face_value(lambda1v.x, 0)*face_gradient_x (TS, 0)*fm.x[]*fsS.x[];
    lambdagradTG.x[] = face_value(lambda2v.x, 0)*face_gradient_x (TG, 0)*fm.x[]*fsG.x[];
  }

  foreach() {
    foreach_dimension()
      DTDtS[] += (lambdagradTS.x[1] - lambdagradTS.x[])/Delta;
    DTDtS[] += sST[];

    foreach_dimension()
      DTDtG[] += (lambdagradTG.x[1] - lambdagradTG.x[])/Delta;
    DTDtG[] += sGT[];
  }

  // EXTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_G[jj];
    scalar DmixGv = DmixGList_G[jj];
    scalar DYDtGjj = DYDtG_G[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
    }

    scalar sgexp = sGexpList[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDtGjj[] += sgexp[];
    }
  }

  /**
  We add diffusion correction contributions to the chemical species
  mass fraction derivatives. */

  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList_G[jj];

      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
# ifdef MOLAR_DIFFUSION
      double MWmixGf = face_value (MWmixG_G, 0);

      scalar XG = XGList_G[jj];
      phicGtot.x[] += (MWmixGf > 0.) ?
        rhoGf*DmixGf*gas_MWs[jj]/MWmixGf*face_gradient_x (XG, 0)*fm.x[]*fsG.x[] : 0.;
# else
      scalar YG = YGList_G[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsG.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicGjj[];
    foreach_face() {
      phicGjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixGv = DmixGList_G[jj];

      double rhoGf = face_value (rhoGv_G, 0);
      double DmixGf = face_value (DmixGv, 0);
      double MWmixGf = face_value (MWmixG_G, 0);

      phicGjj.x[] -= (MWmixGf > 0.) ?
        rhoGf*DmixGf/MWmixGf*face_gradient_x (MWmixG_G, 0)*fm.x[]*fsG.x[] : 0.;
#endif

      scalar YG = YGList_G[jj];
      phicGjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_G[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicGjj.x[1] - phicGjj.x[])/Delta;
  }

  // INTERNAL GAS PHASE
  /**
  We calculate the Lagrangian derivative for the chemical species mass
  fractions. */ 

  for (int jj=0; jj<NGS; jj++) {
    scalar YG = YGList_S[jj];
    scalar DmixGv = DmixGList_S[jj];
    scalar DYDtGjj = DYDtG_S[jj];

    face vector rhoDmixYGjj[];
    foreach_face() {
      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
      rhoDmixYGjj.x[] = rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
    }

    scalar ssexp = sSexpList[jj];

    foreach() {
      foreach_dimension()
        DYDtGjj[] += (rhoDmixYGjj.x[1] - rhoDmixYGjj.x[])/Delta;
      DYDtGjj[] += ssexp[];
    }
  }

  face vector phicStot[];
  foreach_face() {
    phicStot.x[] = 0.;
#ifdef FICK_CORRECTED
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixGv = DmixGList_S[jj];

      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
# ifdef MOLAR_DIFFUSION
      double MWmixGf = face_value (MWmixG_S, 0);

      scalar XG = XGList_S[jj];
      phicStot.x[] += (MWmixGf > 0.) ?
        rhoGf*DmixGf*gas_MWs[jj]/MWmixGf*face_gradient_x (XG, 0)*fm.x[]*fsS.x[] : 0.;
# else
      scalar YG = YGList_S[jj];
      phicStot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0)*fm.x[]*fsS.x[];
# endif // MOLAR_DIFFUSION
    }
#endif  // FICK_CORRECTED
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicSjj[];
    foreach_face() {
      phicSjj.x[] = phicStot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixGv = DmixGList_S[jj];

      double rhoGf = face_value (rhoGv_S, 0);
      double DmixGf = face_value (DmixGv, 0);
      double MWmixGf = face_value (MWmixG_S, 0);

      phicSjj.x[] -= (MWmixGf > 0.) ?
        rhoGf*DmixGf/MWmixGf*face_gradient_x (MWmixG_S, 0)*fm.x[]*fsS.x[] : 0.;
#endif

      scalar YG = YGList_S[jj];
      phicSjj.x[] *= face_value (YG, 0);
    }

    scalar DYDtGjj = DYDtG_S[jj];

    foreach()
      foreach_dimension()
        DYDtGjj[] -= (phicSjj.x[1] - phicSjj.x[])/Delta;
  }

  // We calculate the one-field divergence by volume-averaging the liquid and the
  // gas-phase contributions.

  foreach() {
    double divu1 = 0., divu2 = 0.;

    // Add internal gas temperature contribution
    divu1 += (TS[]*rhoGv_S[]*cpGv_S[] > 0.) ?
      1./(TS[]*(rhoGv_S[]*cpGv_S[]*porosity[]/f[] + rhoSv[]*cpSv[]*(1-porosity[]/f[])))*DTDtS[] : 0.;

    // Add external gas temperature contribution
    divu2 += (TG[]*rhoGv_G[]*cpGv_G[] > 0.) ?
      1./(TG[]*rhoGv_G[]*cpGv_G[])*DTDtG[] : 0.;

    // Add internal gas chemical species contribution
    double divu1species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_S[jj];
      divu1species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu1 += (rhoGv_S[] > 0.) ? MWmixG_S[]/rhoGv_S[]*divu1species : 0.;

    // Add external gas chemical species contribution
    double divu2species = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DYDtGjj = DYDtG_G[jj];
      divu2species += 1./gas_MWs[jj]*DYDtGjj[];
    }
    divu2 += (rhoGv_G[] > 0.) ? MWmixG_G[]/rhoGv_G[]*divu2species : 0.;
    
    // Volume averaged contributions
    drhodt[] = divu1*f[] + divu2*(1. - f[]);

    // Adjust sign for internal convention
    drhodt[] *= -1.;
  }
}

void update_divergence_density_u (void) {

  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  
  scalar DrhoDt[];
  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  vector grho[];
  gradients ({rhot}, {grho});

  foreach()
    foreach_dimension()
      DrhoDt[] += u.x[]*grho.x[];


  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density_uf (void) {
  
  scalar rhot[];
  foreach()
    rhot[] = rhoGv_S[]*f[] + rhoGv_G[]*(1. - f[]);
  
  scalar eps[];
  foreach()
    eps[] = porosity[] + (1. - f[]);
  scalar DrhoDt[];

  foreach()
    DrhoDt[] = (rhot[] - rhoGv_S0[])*eps[]/dt;
  
  face vector rhoGflux[];
  tracer_fluxes (rhot, uf, rhoGflux, dt, zeroc);

  foreach() 
    foreach_dimension()
      DrhoDt[] += (rhoGflux.x[1] - rhoGflux.x[] - rhot[]*(uf.x[1] - uf.x[]))/Delta;

  foreach(){
    DrhoDt[] = DrhoDt[]*cm[];

    double one_over_rho = (rhot[] > 0.) ? 1./rhot[] : 0.;

    if (iter > 1) {
      drhodt[] = DrhoDt[]*one_over_rho;
    }
  }
}

void update_divergence_density (void) {
  // Choose one of the two methods to compute drhodt
  update_divergence_density_u();
  // update_divergence_density_uf();
}
#endif
