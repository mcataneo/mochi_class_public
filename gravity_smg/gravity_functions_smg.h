#ifndef __GRAVITY_FUNCTIONS_SMG__
#define __GRAVITY_FUNCTIONS_SMG__

#include "common.h"
#include "background.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int gravity_functions_As_from_alphas_smg(struct background *pba,
                                         double * pvecback,
                                         double * pvecback_derivs);

int gravity_functions_Cs_from_Bs_smg(struct background *pba,
                                     double * pvecback,
                                     double * pvecback_derivs);


#ifdef __cplusplus
}
#endif

#endif
