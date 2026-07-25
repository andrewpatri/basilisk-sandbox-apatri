/**
## Default OpenSMOKE++ header file
Requires the OpenSMOKE++ interface to link the library functions. See 
[OpenSMOKE++ Interface](https://github.com/edocipriano/OpenSMOKEppInterface) 
and [OpenSMOKE++](https://www.opensmokepp.polimi.it/).

This library provides functions to read thermodynamic, transport properties and reaction rates, 
all key necessity for computing reacting flows.
*/

#include "OpenSMOKE_Interface.h"

#define OPENSMOKE 1

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

char* kinfolder;

/**
## Safe composition conversion

`OpenSMOKE_MoleFractions_From_MassFractions` computes `MW = 1/sum(y_i/MW_i)`, which is a
`1/0` divide-by-zero (SIGFPE where FPE trapping is on) whenever the mass-fraction vector is
all zeros — a fresh / just-uncovered / fully-consumed cell whose composition is unfilled.
These wrappers own that empty-vector policy in one place: they return `false` for an empty
input (mole fractions zeroed, `MW` set to a nonzero fallback so callers that divide by `MW`
stay finite) and `true` otherwise, so a caller can also choose to skip a dependent fill. */

static inline bool mole_from_mass (double * x, double * MW, double * y, int n) {
  double ytot = 0.;
  for (int jj = 0; jj < n; jj++) 
    ytot += y[jj];

  if (ytot > 0.) {
    OpenSMOKE_MoleFractions_From_MassFractions (x, MW, y);
    return true;
  }

  for (int jj = 0; jj < n; jj++) 
    x[jj] = 0.;
  *MW = 1.;
  return false;
}

static inline bool solid_mole_from_mass (double * x, double * MW, double * y, int n) {
  double ytot = 0.;
  for (int jj = 0; jj < n; jj++) 
    ytot += y[jj];

  if (ytot > 0.) {
    OpenSMOKE_SolidMoleFractions_From_SolidMassFractions (x, MW, y);
    return true;
  }

  for (int jj = 0; jj < n; jj++) 
    x[jj] = 0.;
  *MW = 1.;
  return false;
}

event cleanup (t = end)
{
  OpenSMOKE_Clean ();
}
