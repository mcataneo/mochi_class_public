#ifndef __BACKGROUND_SMG__
#define __BACKGROUND_SMG__

#include "background.h"
#include "common.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int background_gravity_functions_smg(
        struct background *pba,
        double * pvecback_B,
        short return_format,
        double * pvecback
        );

#ifdef __cplusplus
}
#endif

#endif
