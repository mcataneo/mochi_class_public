#ifndef __NONLINEAR_SMG__
#define __NONLINEAR_SMG__

#include "common.h"
#include "nonlinear.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int nonlinear_hmcode_correct_Delta_v_0_smg(
  struct background *pba,
  double z_at_tau,
  double * Delta_v_0
);

#ifdef __cplusplus
}
#endif

#endif
