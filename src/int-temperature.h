/**
# Interface temperature solver
This header contains the implementation of the solver for the interface temperature.
We exploit the GSL library to solve the non-linear equation for the interface temperature.

This a modifid version of Riccardo Caraccio's one in order to accomodate a q_sorg [in W/m2] at the interface of solid anche gas here it is constant at any point so it is a double not a scalar, if not uncomment the comment close to it
https://github.com/Riccaraccio/basilisk-sandbox-rcaraccio/blob/master/src/int-temperature.h
*/

#include "intgrad.h"
#ifndef USE_GSL
# define USE_GSL 1
#endif
#include "fsolve-gsl.h"

#ifndef RADIATION_TEMP
#define RADIATION_TEMP TG0
#endif

//Extern variables
extern scalar fS, fG;
extern face vector fsS, fsG;
extern vector lambda1v, lambda2v;
extern double TG0;
extern scalar TInt, TS, TG;
extern double q_sorg; // extern scalar q_sorg

typedef struct {
  coord c;
  double emissivity;
} UserDataNls;

#ifndef RADIATION_INTERFACE
# define RADIATION_INTERFACE 0.85 // Gronli thesis 1996
#endif

// Radiative heat flux at the interface
double divq_rad_int (double TInti, double Tbulk = 300., double alphacorr = 1.) {
  return 5.670373e-8*alphacorr*(pow(Tbulk, 4.) - pow(TInti, 4.));
}

// From Galgano & Diblasi 2003
double emissivity_diblasi (const double char_fraction, const double ash_fraction) {
  NOT_UNUSED(ash_fraction);
  const double char_emissivity = 1;
  const double wood_emissivity = 0.6;
  return char_fraction * char_emissivity + (1. - char_fraction) * wood_emissivity;
}

// From Lu et al 2008
double emissivity_lu (const double char_fraction, const double ash_fraction){
  const double char_emissivity = 0.95;
  const double wood_emissivity = 0.85;
  const double ash_emissivity = 0.7;
  return char_fraction * char_emissivity + (1. - char_fraction - ash_fraction) * wood_emissivity + ash_fraction * ash_emissivity;
}

double emissivity_constant (const double char_fraction, const double ash_fraction) {
  NOT_UNUSED(char_fraction);
  NOT_UNUSED(ash_fraction);
  return RADIATION_INTERFACE;
}

double (*emissivity) (const double char_fraction, const double ash_fraction) = NULL;

event defaults (i = 0) {
  if (emissivity == NULL)
    emissivity = emissivity_constant;
}

// Non-linear equation for the interface temperature to be zeroed
int EqTemperature (const gsl_vector * xdata, void * params, gsl_vector * fdata) {
  UserDataNls * data = (UserDataNls *)params;
  double TInti = gsl_vector_get(xdata, 0);
  
  Point point = locate(data->c.x, data->c.y, data->c.z);
  // foreach_point(data->c.x, data->c.y, data->c.z, serial) {

  bool success = false;

  double gradTGn = ebmgrad(point, TG, fS, fG, fsS, fsG, true, TInti, &success);
  double gradTSn = ebmgrad(point, TS, fS, fG, fsS, fsG, false, TInti, &success);

  coord n = facet_normal(point, fS, fsS);
  normalize(&n);
  n.x = fabs(n.x); n.y = fabs(n.y);

  double lambda1vh = n.x / (n.x + n.y) * lambda1v.x[] + n.y / (n.x + n.y) * lambda1v.y[];
  double lambda2vh = n.x / (n.x + n.y) * lambda2v.x[] + n.y / (n.x + n.y) * lambda2v.y[];
 
  double q_sorgi = q_sorg; // double q_sorgi = q_sorg[]; // here i stands for the value in the the point
  gsl_vector_set(fdata, 0,
                 -divq_rad_int(TInti, RADIATION_TEMP, data->emissivity)
                 + lambda1vh * gradTSn 
                 + lambda2vh * gradTGn
                 + q_sorgi * emissivity);
  // }
  return GSL_SUCCESS;
}

void ijc_CoupledTemperature() {

  scalar YASH;
  int ash_index = OpenSMOKE_IndexOfSolidSpeciesWithoutError("ASH");
  if (ash_index >= 0) {
    YASH = YSList[ash_index];
  }
  foreach() {
    if (f[]>F_ERR && f[] < 1.-F_ERR) {
      gsl_vector *unk = gsl_vector_alloc(1);
      gsl_vector_set(unk, 0, TInt[]);

      UserDataNls data;
      coord o = {x,y,z};
      foreach_dimension()
        data.c.x = o.x;

      double char_fraction = calculate_char_fraction(point, YSList, f);
      if (ash_index >= 0) {
        data.emissivity = emissivity(char_fraction, YASH[]/f[]);
      } else {
        data.emissivity = emissivity(char_fraction, 0.);
      }

      fsolve_gsl (EqTemperature, unk, &data, "EqTemperature");

      TInt[] = gsl_vector_get(unk, 0);
      gsl_vector_free(unk);
    }
  }
}