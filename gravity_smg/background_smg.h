#ifndef __BACKGROUND_SMG__
#define __BACKGROUND_SMG__

#include "common.h"
#include "background.h"

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
