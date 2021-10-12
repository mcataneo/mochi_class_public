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

int hi_class_define_indices_bg(
				struct background *pba,
				int * index_bg
			  );

int hi_class_define_indices_bi(
				struct background *pba,
				int * index_bg
			  );

int background_derivs_smg(
				struct background *pba,
        double * pvecback,
        double * pvecback_derivs,
        int i
			  );

int background_gravity_functions_A_C_smg(
        struct background *pba,
        double * pvecback,
        double * pvecback_derivs,
        int i
        );


#ifdef __cplusplus
}
#endif

#endif
