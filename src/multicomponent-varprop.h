#ifndef MULTICOMPONENT
  #define MULTICOMPONENT 1
#endif

#include "intgrad.h"

#ifdef EXPLICIT_DIFFUSION
  #include "diffusion-explicit.h"
#else
  #include "diffusion.h"
#endif

#include "common-phasechange.h"
#include "memoryallocation-varprop.h"
#include "int-temperature-mod.h"
#include "int-concentration.h"
#include "multicomponent-properties.h"
#include "chemistry.h"

event reset_sources (i++) {
#ifdef SOLVE_TEMPERATURE
  foreach() {
    sST[] = 0.;
    sGT[] = 0.;
  }
#endif

  reset (sGexpList, 0.);
  reset (sSexpList, 0.);
}

extern face vector ufsave;
face vector u_prime[];
#ifndef STOP_TRACER_ADVECTION
event tracer_advection (i++) {

foreach() {
    f[] = clamp (f[], 0., 1.);
    f[] = (f[] > F_ERR) ? f[] : 0.;
    f[] = (f[] < 1.-F_ERR) ? f[] : 1.;
    fS[] = f[]; fG[] = 1. - f[];
  }

  //Compute face gradients
  face_fraction (fS, fsS);
  face_fraction (fG, fsG);

  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
  check_and_correct_fractions (YSList,   NSS, false);

#ifdef VARPROP
  update_properties();
#else
  update_properties_constant();
#endif

  // lose tracer form and extrapolate fields
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1.-f[]) > F_ERR) ? TG[]/(1.-f[]) : 0.;

    TS[] = (f[] > F_ERR) ? TS[] : TG[];
    TG[] = (f[] < 1.-F_ERR) ? TG[] : TS[];
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]/f[] : 0.;
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[]/(1.-f[]) : 0.;

      YG_S[] = (f[] > F_ERR) ? YG_S[] : YG_G[];
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[] : YG_S[];
    }
  }

  advection_div(YGList_S, ufsave, dt);
  advection_div(YGList_G, ufsave, dt);

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    double ef = clamp(face_value(porosity, 0), 0., 1.);

    double rhoGvh_S, rhoSvh;
    double cpGvh_S, cpSvh;

    #ifdef VARPROP
    rhoGvh_S = face_value(rhoGv_S, 0); rhoSvh = face_value(rhoSv, 0);
    cpGvh_S = face_value(cpGv_S, 0); cpSvh = face_value(cpSv, 0);
    #else
    rhoGvh_S = rhoG; rhoSvh = rhoS;
    cpGvh_S = cpG; cpSvh = cpS;
    #endif

    u_prime.x[] = (fsS.x[] > F_ERR) ? 
                  fsS.x[]*ufsave.x[]*(rhoGvh_S*cpGvh_S)/
                  (rhoGvh_S*cpGvh_S*ef + rhoSvh*cpSvh*(1. - ef))
                  : 0.;
  }

  advection_div({TS}, u_prime, dt);
# ifndef TEMPERATURE_PROFILE
  advection_div({TG}, ufsave, dt);
# endif
#endif

  // recover tracer form
  foreach() {
    porosity[] = (f[] > F_ERR) ? porosity[]*f[] : 0.;
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1.-f[]) > F_ERR) ? TG[]*(1.-f[]) : 0.;
#endif

    for (int jj=0; jj<NGS; jj++) { 
      scalar YG_S = YGList_S[jj];
      scalar YG_G = YGList_G[jj];

      YG_S[] = (f[] > F_ERR) ? YG_S[]*f[] : 0.;
      YG_G[] = (f[] < 1.-F_ERR) ? YG_G[]*(1. - f[]) : 0.;
    }
  }
}
#endif

void update_mole_fields() {
  #ifdef MOLAR_DIFFUSION
  foreach() {
    double xG[NGS], yG[NGS];
    double MWmix;
    if (f[] > F_ERR) { // Internal gas phase
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_S[jj];
        yG[jj] = YG[];
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmix, yG);
      // MWmixG_S[] = MWmix; // Already done in update_properties()
      for (int jj = 0; jj < NGS; jj++) {
        scalar XG = XGList_S[jj];
        XG[] = xG[jj];
      }
    }

    if (f[] < 1. - F_ERR) { // External gas phase
      for (int jj = 0; jj < NGS; jj++) {
        scalar YG = YGList_G[jj];
        yG[jj] = YG[];
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmix, yG);
      // MWmixG_G[] = MWmix; // Already done in update_properties()
      for (int jj = 0; jj < NGS; jj++) {
        scalar XG = XGList_G[jj];
        XG[] = xG[jj];
      }
    }
  }
  boundary (XGList_S);
  boundary (XGList_G); // Ensure boundary conditions are applied
  #endif
}

event tracer_diffusion (i++) {

  //Check the mass fractions Can be removed for performance
  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
  check_and_correct_fractions (YSList,   NSS, false);

  foreach() {
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]/f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]/(1. - f[]) : 0.;
#endif

    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_S[jj];
      YG[] = (f[] > F_ERR) ? YG[]/f[] : 0.;
    }
    
    for (int jj=0; jj<NGS; jj++) {
      scalar YG = YGList_G[jj];
      YG[] = ((1. - f[]) > F_ERR) ? YG[]/(1. - f[]) : 0.;
    }
  }

#ifdef MOLAR_DIFFUSION
  update_mole_fields();
#endif


#ifdef SOLVE_TEMPERATURE
  //interface temperature first guess
  foreach() {
    TInt[] = 0.;
    if (f[] > F_ERR && f[] < 1.-F_ERR)
      TInt[] = (TS[] + TG[])/2;
  }

  #ifdef FIXED_INT_TEMP //Force interface temperature = TG0
  foreach()
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      TInt[] = TG0;

  #elif TEMPERATURE_PROFILE
  double tv = TemperatureProfile_GetT(t);
  foreach() {
    if (f[] > F_ERR && f[] < 1. - F_ERR)
      TInt[] = tv;
    
    if (f[] < 1. - F_ERR)
      TG[] = tv;
  }

  #else //default: solve for interface temperature
  ijc_CoupledTemperature();
  #endif
#endif

  // first guess for species interface concentration
  foreach() {
    for (int jj=0; jj<NGS; jj++) {
      scalar YGInt = YGList_Int[jj];
      YGInt[] = 0.;
      if (f[] > F_ERR && f[] < 1. - F_ERR) {
        scalar YG_S = YGList_S[jj];
        scalar YG_G = YGList_G[jj];
        YGInt[] = (YG_G[] + YG_S[])/2;
        YGInt[] = clamp (YGInt[], 0., 1.);
      }
    }
  }

  //find the interface concentration for each species
  intConcentration();

#ifdef MOLAR_DIFFUSION // calculate the mole fractions at the interface
  foreach()
    if (f[] > F_ERR && f[] < 1. - F_ERR) {
      double xG[NGS], yG[NGS], MWmixInt;
      for (int jj=0; jj<NGS; jj++) {
        scalar YGInt = YGList_Int[jj];
        yG[jj] = YGInt[];
      }
      OpenSMOKE_MoleFractions_From_MassFractions (xG, &MWmixInt, yG);
      for (int jj=0; jj<NGS; jj++) {
        scalar XGInt = XGList_Int[jj];
        XGInt[] = xG[jj];
      }
    }
#endif

  //Calculate the source therm
  foreach() {
    if (f[] > F_ERR && f[] < 1.-F_ERR) {
      coord n = facet_normal (point, fS, fsS), p;
      double alpha = plane_alpha (fS[], n);
      double area = plane_area_center (n, alpha, &p);
      normalize (&n);

      //Solid side
      double jS[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar DmixG = DmixGList_S[jj];

        double rhoGvh_S;
        #ifdef VARPROP
        rhoGvh_S = rhoGv_S[];
        #else
        rhoGvh_S = rhoG;
        #endif
        
        #ifdef MOLAR_DIFFUSION
        scalar XG = XGList_S[jj];
        scalar XGInt = XGList_Int[jj];
        double Strgrad = ebmgrad (point, XG, fS, fG, fsS, fsG, false, XGInt[], &success);
        jS[jj] = rhoGvh_S*DmixG[]*Strgrad*gas_MWs[jj]/MWmixG_S[];
        #else
        scalar YG    = YGList_S[jj];
        scalar YGInt = YGList_Int[jj];
        double Strgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, false, YGInt[], &success);
        jS[jj] = rhoGvh_S*DmixG[]*Strgrad; 
        #endif
      }

      double jStot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++)
        jStot += jS[jj];
#endif

      for (int jj=0; jj<NGS; jj++) {
        scalar sSexp = sSexpList[jj];
        scalar YGInt = YGList_Int[jj];
        jS[jj] -= jStot*YGInt[];
#ifdef AXI
        sSexp[] += jS[jj]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sSexp[] += jS[jj]*area/Delta*cm[];
#endif
      }

      //Gas side
      double jG[NGS];
      for (int jj=0; jj<NGS; jj++) {
        scalar DmixG = DmixGList_G[jj];

        double rhoGvh_G;
#ifdef VARPROP
        rhoGvh_G = rhoGv_G[];
#else
        rhoGvh_G = rhoG;
#endif

#ifdef MOLAR_DIFFUSION
        scalar XG = XGList_G[jj];
        scalar XGInt = XGList_Int[jj];
        double Gtrgrad = ebmgrad (point, XG, fS, fG, fsS, fsG, true, XGInt[], &success);
        jG[jj] = rhoGvh_G*DmixG[]*Gtrgrad*gas_MWs[jj]/MWmixG_G[];
#else
        scalar YG    = YGList_G[jj];
        scalar YGInt = YGList_Int[jj];
        double Gtrgrad = ebmgrad (point, YG, fS, fG, fsS, fsG, true, YGInt[], &success);
        jG[jj] = rhoGvh_G*DmixG[]*Gtrgrad;
#endif
      }

      double jGtot = 0.;
#ifdef FICK_CORRECTED
      for (int jj=0; jj<NGS; jj++)
        jGtot += jG[jj];
#endif

      for (int jj=0; jj<NGS; jj++) {
        scalar sGexp = sGexpList[jj];
        scalar YGInt = YGList_Int[jj];
        jG[jj] -= jGtot*YGInt[];
#ifdef AXI
        sGexp[] += jG[jj]*area*(y + p.y*Delta)/(Delta*y)*cm[];
#else
        sGexp[] += jG[jj]*area/Delta*cm[];
#endif
      }

#ifdef SOLVE_TEMPERATURE
      double bc = TInt[];
      double Strgrad = ebmgrad (point, TS, fS, fG, fsS, fsG, false, bc, &success);
      double Gtrgrad = ebmgrad (point, TG, fS, fG, fsS, fsG, true , bc, &success);

      n.x = fabs(n.x); n.y = fabs(n.y);

      double lambda1vh = n.x/(n.x+n.y)*lambda1v.x[] + n.y/(n.x+n.y)*lambda1v.y[];
      double lambda2vh = n.x/(n.x+n.y)*lambda2v.x[] + n.y/(n.x+n.y)*lambda2v.y[];

      double Sheatflux = lambda1vh*Strgrad;
      double Gheatflux = lambda2vh*Gtrgrad;

# ifdef AXI
      sST[] += Sheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
      sGT[] += Gheatflux*area*(y + p.y*Delta)/(Delta*y)*cm[];
# else
      sST[] += Sheatflux*area/Delta*cm[];
      sGT[] += Gheatflux*area/Delta*cm[];
# endif
#endif
    }
  }

#ifdef MASS_DIFFUSION_ENTHALPY
  foreach() {
    if (f[] > 1. - F_ERR) { //Internal gas phase
      double mdeGS = 0.;
      coord gTS = {0., 0., 0.};
      coord gYGj_S = {0., 0., 0.};
      coord gYGsum_S = {0., 0., 0.};

      foreach_dimension()
        gTS.x = (TS[1] - TS[-1])/(2.*Delta);
      
      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGsum_S.x -= (MWmixG_S[] > 0.) ?
            rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList_S[jj];
          gYGsum_S.x -= rhoGv_S[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_S[jj];
          scalar cpGv = cpGList_S[jj];
          scalar Dmixv = DmixGList_S[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_S[jj];
          gYGj_S.x = (MWmixG_S[] > 0.) ?
            -rhoGv_S[]*Dmixv[]*gas_MWs[jj]/MWmixG_S[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_S.x = -rhoGv_S[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGS += cpGv[]*(gYGj_S.x - YG[]*gYGsum_S.x)*gTS.x;
        }
      }
    sST[] -= mdeGS*cm[];
    }
   
    if (f[] < F_ERR) { //Internal gas phase
      double mdeGG = 0.;
      coord gTG = {0., 0., 0.};
      coord gYGj_G = {0., 0., 0.};
      coord gYGsum_G = {0., 0., 0.};

      foreach_dimension()
        gTG.x = (TG[1] - TG[-1])/(2.*Delta);
      
      foreach_dimension() {
        for (int jj=0; jj<NGS; jj++) {
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGsum_G.x -= (MWmixG_G[] > 0.) ?
            rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          scalar YG = YGList_G[jj];
          gYGsum_G.x -= rhoGv_G[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
        }

        for (int jj=0; jj<NGS; jj++) {
          scalar YG = YGList_G[jj];
          scalar cpGv = cpGList_G[jj];
          scalar Dmixv = DmixGList_G[jj];
  # ifdef MOLAR_DIFFUSION
          scalar XG = XGList_G[jj];
          gYGj_G.x = (MWmixG_G[] > 0.) ?
            -rhoGv_G[]*Dmixv[]*gas_MWs[jj]/MWmixG_G[]*(XG[1] - XG[-1])/(2.*Delta) : 0.;
  # else
          gYGj_G.x = -rhoGv_G[]*Dmixv[]*(YG[1] - YG[-1])/(2.*Delta);
  # endif
          mdeGG += cpGv[]*(gYGj_G.x - YG[]*gYGsum_G.x)*gTG.x;
        }
      }
    sGT[] -= mdeGG*cm[];
    }
  }
#endif //MASS_DIFFUSION_ENTHALPY

#if defined VARPROP && !defined NO_EXPANSION
  update_divergence();
  // update_divergence_density();
#endif

#ifdef FICK_CORRECTED
  face vector phicGtot[];
  foreach_face() {
    phicGtot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixG = DmixGList_G[jj];
      double DmixGf = face_value(DmixG, 0);
      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_G, 0);
# else
      rhoGf = rhoG;
# endif

# ifdef MOLAR_DIFFUSION
      scalar XG = XGList_G[jj];
      double MWmixf = face_value(MWmixG_G, 0);
      phicGtot.x[] += (MWmixf > 0.) ? rhoGf*DmixGf*face_gradient_x (XG, 0)*gas_MWs[jj]/MWmixf : 0.;
# else
      scalar YG = YGList_G[jj];
      phicGtot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0);
# endif
    }
    phicGtot.x[] *= fsG.x[]*fm.x[];
  }

  face vector phicStot[];
  foreach_face() {
    phicStot.x[] = 0.;
    for (int jj=0; jj<NGS; jj++) {
      scalar DmixG = DmixGList_S[jj];
      double DmixGf = face_value(DmixG, 0);
      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_S, 0);
# else
      rhoGf = rhoG;
# endif

# ifdef MOLAR_DIFFUSION
      scalar XG = XGList_S[jj];
      double MWmixf = face_value(MWmixG_S, 0);
      phicStot.x[] += (MWmixf > 0.) ? rhoGf*DmixGf*face_gradient_x (XG, 0)*gas_MWs[jj]/MWmixf : 0.;
# else
      scalar YG = YGList_S[jj];
      phicStot.x[] += rhoGf*DmixGf*face_gradient_x (YG, 0);
# endif
    }
    phicStot.x[] *= fsS.x[]*fm.x[];
  }

  //Apply the Fick's law correction
  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicGtot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixG = DmixGList_G[jj];
      double DmixGf = face_value(DmixG, 0);
      double MWmixf = face_value(MWmixG_G, 0);

      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_G, 0);
# else
      rhoGf = rhoG;
# endif
      phicjj.x[] -= (MWmixf > 0.) ? rhoGf*DmixGf/MWmixf*face_gradient_x (MWmixG_G, 0)*fsG.x[]*fm.x[] : 0.;
#endif
    }

    scalar YG = YGList_G[jj];
    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = NULL; //reset the gradient
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function
    
    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_G[] > 0.) ? dt/(rhoGv_G[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.; 
  }

  for (int jj=0; jj<NGS; jj++) {
    face vector phicjj[];
    foreach_face() {
      phicjj.x[] = phicStot.x[];
#ifdef MOLAR_DIFFUSION
      scalar DmixG = DmixGList_S[jj];
      double DmixGf = face_value(DmixG, 0);
      double MWmixf = face_value(MWmixG_S, 0);

      double rhoGf;
# ifdef VARPROP
      rhoGf = face_value(rhoGv_S, 0);
# else
      rhoGf = rhoG;
# endif
      phicjj.x[] -= (MWmixf > 0.) ? rhoGf*DmixGf/MWmixf*face_gradient_x (MWmixG_S, 0)*fsS.x[]*fm.x[] : 0.;
#endif
  }

    scalar YG = YGList_S[jj];
    double (* gradient_backup)(double, double, double) = YG.gradient; // we need to backup the gradient function
    YG.gradient = NULL; //reset the gradient
    face vector flux[];
    tracer_fluxes (YG, phicjj, flux, dt, zeroc); //calculate the fluxes using the corrective velocity
    YG.gradient = gradient_backup; // restore the gradient function
    
    // apply the corrective fluxes
    foreach()
      foreach_dimension()
        YG[] += (rhoGv_S[] > 0.) ? dt/(rhoGv_S[])*(flux.x[] - flux.x[1])/(Delta*cm[]) : 0.; 
  }
  #endif //FICK_CORRECTED

  scalar theta1[], theta2[];

#if TREE
  theta1.refine = fraction_refine;
  set_prolongation (theta1, fraction_refine);
  theta2.refine = fraction_refine;
  set_prolongation (theta2, fraction_refine);
#endif

  // Internal gas diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_S[jj];
    foreach_face() {
      double rhoGvh_S;
#ifdef VARPROP
      rhoGvh_S = face_value(rhoGv_S, 0);
#else
      rhoGvh_S = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_S*fsS.x[]*fm.x[];
    }

    foreach() {
#ifdef VARPROP
      theta1[] = cm[]*max(rhoGv_S[]*porosity[], F_ERR); // porosity is already multiplied by fS
#else
      theta1[] = cm[]*max(rhoG*porosity[], F_ERR); // porosity is already multiplied by fS
#endif
    }

    scalar YG = YGList_S[jj];
    scalar sSexp = sSexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta1);
#else
    diffusion (YG, dt, D=DmixGf, r=sSexp, theta=theta1);
#endif
  }

  //external diffusion
  for (int jj=0; jj<NGS; jj++) {
    face vector DmixGf[];
    scalar DmixG = DmixGList_G[jj];
    foreach_face() {
      double rhoGvh_G;
#ifdef VARPROP
      rhoGvh_G = face_value(rhoGv_G, 0);
#else
      rhoGvh_G = rhoG;
#endif
      DmixGf.x[] = face_value(DmixG, 0)*rhoGvh_G*fsG.x[]*fm.x[];
    }
    foreach() {
#ifdef VARPROP
      theta2[] = cm[]*max(fG[]*rhoGv_G[], F_ERR);
#else
      theta2[] = cm[]*max(fG[]*rhoG, F_ERR);
#endif
    }

    scalar YG = YGList_G[jj];
    scalar sGexp = sGexpList[jj];

#ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (YG, dt, D=DmixGf, theta=theta2);
#else
    diffusion (YG, dt, D=DmixGf, r=sGexp, theta=theta2);
#endif
  }

#ifdef SOLVE_TEMPERATURE
  foreach_face() {
    lambda1f.x[] = face_value(lambda1v.x, 0)*fsS.x[]*fm.x[];
    lambda2f.x[] = face_value(lambda2v.x, 0)*fsG.x[]*fm.x[];
  }

  foreach() {
    double theta1vh, theta2vh;
# ifdef VARPROP
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoGv_S[]*cpGv_S[] + (1. - porosity[]/fS[])*rhoSv[]*cpSv[] : 0.;
    theta2vh = rhoGv_G[]*cpGv_G[];
# else
    theta1vh = fS[] > F_ERR ? porosity[]/fS[]*rhoG*cpG + (1. - porosity[]/fS[])*rhoS*cpS : 0.;
    theta2vh = rhoG*cpG;
# endif

    theta1[] = cm[]*max(fS[]*theta1vh, F_ERR);
    theta2[] = cm[]*max(fG[]*theta2vh, F_ERR);
  }

#ifdef VARCOEFF
  foreach()
    porosity[] = (f[] > F_ERR) ? porosity[]/f[] : 0;

  foreach_face() {
    double ef = face_value(porosity, 0);
    lambda1f.x[] = (ef > F_ERR) ? lambda1f.x[] / (rhoG*cpG*ef + rhoS*cpS*(1. - ef)) : 0.;
    lambda2f.x[] = lambda2f.x[] / (rhoG*cpG);
  }

  foreach() {
    theta1[] = cm[] * max(fS[], F_ERR);
    theta2[] = cm[] * max(fG[], F_ERR);
  }

  foreach() {
    sST[] = (f[] > F_ERR) ? sST[] / (rhoG*cpG*porosity[] + rhoS*cpS*(1. - porosity[])) : 0.;
    sGT[] = sGT[] / (rhoG*cpG);
  }

  foreach()
    porosity[] *= f[];
#endif

# ifdef EXPLICIT_DIFFUSION
    diffusion_explicit (TS, dt, D=lambda1f, r=sST, theta=theta1);
    diffusion_explicit (TG, dt, D=lambda2f, r=sGT, theta=theta2);
# else
    diffusion (TS, dt, D=lambda1f, r=sST, theta=theta1);
#  ifndef TEMPERATURE_PROFILE
    diffusion (TG, dt, D=lambda2f, r=sGT, theta=theta2);
#  endif
# endif
#endif

  //recover tracer form
  foreach() {
    for (scalar YG in YGList_S)
      YG[] = (f[] > F_ERR) ? YG[]*f[] : 0.;

    for (scalar YG in YGList_G)
      YG[] = ((1. - f[]) > F_ERR) ? YG[]*(1. - f[]) : 0.;
    
#ifdef SOLVE_TEMPERATURE
    TS[] = (f[] > F_ERR) ? TS[]*f[] : 0.;
    TG[] = ((1. - f[]) > F_ERR) ? TG[]*(1. - f[]) : 0.;
    T[] = TS[] + TG[];
#endif
  }

  check_and_correct_fractions (YGList_S, NGS, false);
  check_and_correct_fractions (YGList_G, NGS, true);
}
