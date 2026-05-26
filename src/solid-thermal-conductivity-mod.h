/**
# Solid thermal conductivity
This file defines the pseudophase thermal conductivity field and its update function. 
It defines various models for the solid thermal conductivity, which can be selected
through the 'lambdaSmodel' variable. The solid thermal conductivity is updated based
on the local temperature and composition of the solid phase, as well as on the
porosity and volume fraction fields.
 */

extern scalar* YSList; // List of solid species mass fraction fields

enum solid_thermal_conductivity_model {
  L_CONST,
  L_CORBETTA,
  L_HUANG,
  L_HUANG2,
  L_HUANG05,
  L_ANCACOUCE,
  L_KK,
  L_LU,
  L_GALGANO,
  L_LU_CONV,
};

enum solid_thermal_conductivity_model lambdaSmodel = L_CONST;

// We return a coord because some models can be anisotropic
// Constant solid thermal conductivity model
coord lambda_const (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
  double lambda = (1. - porosity) * lambdaS + porosity * lambdaG;
  return (coord){lambda, lambda};
}

// Corbetta model for solid thermal conductivity
coord lambda_corbetta (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    NOT_UNUSED(Temperature);
    double char_cond = 0.1405;
    double bio_cond = 0.1937;
    double char_fraction = calculate_char_fraction(point, YSList, f);
    double lambda = char_cond * char_fraction + bio_cond * (1. - char_fraction) + porosity * lambdaG;
    return (coord){lambda, lambda};
}

// Huang model for solid thermal conductivity
coord lambda_huang (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    double char_cond = 0.071;
    double bio_cond = 0.21;
    double char_fraction = calculate_char_fraction(point, YSList, f);

    int ash_index = OpenSMOKE_IndexOfSolidSpeciesWithoutError("ASH");
    scalar YASH;
    if (ash_index >= 0)
        YASH = YSList[ash_index];

    double ash_fraction = (ash_index >= 0) ? YASH[]/f[] : 0.;
    double lambda = (char_cond*char_fraction + 
                     bio_cond*(1. - char_fraction))*(1. - porosity) 
                    + 13.5*5.67e-8*pow(Temperature, 3)*80e-06/emissivity (char_fraction, ash_fraction)
                    + porosity*lambdaG;
    return (coord){lambda, lambda};
}

// Ancacouce model for solid thermal conductivity
coord lambda_ancacouce (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    NOT_UNUSED(porosity);
    double char_cond = 0.125;
    double bio_cond = 0.056 + 2.6e-4*Temperature;
    double char_fraction = calculate_char_fraction(point, YSList, f);
    double lambda = char_cond*char_fraction + bio_cond*(1. - char_fraction);
    return (coord){lambda, lambda};
}

// KK model for solid thermal conductivity
coord lambda_kk (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    NOT_UNUSED(Temperature);
    double lS_per = 0.430;
    double lS_par = 0.766;
    double leff_per = 1./((1. - porosity)/lS_per + porosity/lambdaG);
    double leff_par = (1. - porosity)*lS_par + porosity*lambdaG;
    // longitudinal direction theta = 1.0
    double lambda_par = leff_par;
    // transversal direction theta = 0.58
    double lambda_per = 0.58*leff_par + (1. - 0.58)*leff_per;
    // This model is anisotropic, we return the longitudinal conductivity as the effective one for simplicity
    return (coord){lambda_par, lambda_per};
}

// Lu model for solid thermal conductivity
coord lambda_lu (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    double Cw = calculate_moisture_fraction(point, YSList, f);
    double char_fraction = calculate_char_fraction(point, YSList, f);
    scalar YASH = YSList[OpenSMOKE_IndexOfSolidSpecies("ASH")];
    double ash_fraction = YASH[]/f[];

    double char_cond = 0.071;
    double ash_cond = 1.2;
    double wood_cond = (0.129 - 4.9e-2*Cw)*(1. + (2.05 + 4.*Cw)*1e-3*(Temperature - 273.15))*(0.986 + 2.695*Cw);
    double rad_cont = 5.67e-8*pow(Temperature, 3)*3.2e-6/emissivity(char_fraction, ash_fraction)*porosity;
    double lambda = (char_cond*char_fraction + wood_cond*(1. - char_fraction - ash_fraction) + ash_cond*ash_fraction)*\
                    (1. - porosity) + rad_cont + porosity*lambdaG;
    return (coord){lambda, lambda};
}

coord lambda_galgano (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
  double char_fraction = calculate_char_fraction(point, YSList, f);
  double moisture_fraction = calculate_moisture_fraction(point, YSList, f);
  double wood_fraction = 1. - char_fraction - moisture_fraction;
  scalar YASH = YSList[OpenSMOKE_IndexOfSolidSpecies("ASH")];
  double ash_fraction = YASH[]/f[];

  // Pore diameter
  double dp_wood = 4.0e-5;
  double dp_char = 2.0e-4;
  double dp = wood_fraction*dp_wood + char_fraction*dp_char;

  // Char conductivity with correction
  double char_cond = 0.071;
  double wood_cond = 0.15;
  double moist_cond = 0.23;

  // Effective thermal conductivity
  double lambda = (wood_fraction*wood_cond
                + char_fraction*char_cond
                + moisture_fraction*moist_cond)*(1. - porosity)
                + porosity*lambdaG
                + 4.*dp*porosity/(1. - porosity)
                  *pow(Temperature, 3)*5.67e-8
                  *emissivity(char_fraction, ash_fraction);

  return (coord){lambda, lambda};
}

double wood_fraction_0;
coord lambda_lu_conv (Point point, double lambdaG, double porosity, double Temperature, scalar f) {
    double Cw = calculate_moisture_fraction(point, YSList, f);
    double char_fraction = calculate_char_fraction(point, YSList, f);
    scalar YASH = YSList[OpenSMOKE_IndexOfSolidSpecies("ASH")];
    double ash_fraction = YASH[]/f[];
    double wood_fraction = 1. - char_fraction - Cw;

    double conversion_wood = 1. - (1. - porosity)/(1. - eps0)*(1. - wood_fraction)/(1. - wood_fraction_0);
    conversion_wood = clamp(conversion_wood, 0., 1.);

    double char_cond = 0.071;
    double ash_cond = 1.2;
    double wood_cond = (0.129 - 4.9e-2*Cw)*(1. + (2.05 + 4.*Cw)*1e-3*(Temperature - 273.15))*(0.986 + 2.695*Cw);
    double rad_cont = 5.67e-8*pow(Temperature, 3)*3.2e-6/emissivity(char_fraction, ash_fraction)*porosity;

    double lambda = (char_cond*conversion_wood 
                     + wood_cond*(1. - conversion_wood)
                     + ash_cond*ash_fraction)*(1. - porosity) 
                     + rad_cont 
                     + porosity*lambdaG;
    return (coord){lambda, lambda};
}

// Function pointer for the pseudophase thermal conductivity model
coord (*pseudo_phase_thermal_conductivity) (Point point, double lambdaG, double porosity, double Temperature, scalar f) = lambda_const;

event defaults (i = 0) {
  switch (lambdaSmodel) {
    case L_CONST:
      pseudo_phase_thermal_conductivity = lambda_const;
      break;
    case L_CORBETTA:
      pseudo_phase_thermal_conductivity = lambda_corbetta;
      break;
    case L_HUANG:
      pseudo_phase_thermal_conductivity = lambda_huang;
      break;
    case L_HUANG2:
      pseudo_phase_thermal_conductivity = 2.*lambda_huang;
      break;
    case L_HUANG05:
      pseudo_phase_thermal_conductivity = 0.5*lambda_huang;
      break;
    case L_ANCACOUCE:
      pseudo_phase_thermal_conductivity = lambda_ancacouce;
      break;
    case L_KK:
      pseudo_phase_thermal_conductivity = lambda_kk;
      break;
    case L_LU:
      pseudo_phase_thermal_conductivity = lambda_lu;
      break;
    case L_GALGANO:
      pseudo_phase_thermal_conductivity = lambda_galgano;
      break;
    case L_LU_CONV:
      pseudo_phase_thermal_conductivity = lambda_lu_conv;
      break;
    default:
      fprintf(stderr, "ERROR: Unknown solid thermal conductivity model, unkown lambdaSmodel = %d \n", lambdaSmodel);
      abort();
      break;
  }
}

// event properties (i = 0) {
//   if (lambdaSmodel == L_LU_CONV)
//     wood_fraction_0 = 1. - calculate_char_fraction(point, YSList, f) - calculate_moisture_fraction(point, YSList, f);
// }