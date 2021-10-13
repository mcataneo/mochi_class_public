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
				int * index_bi
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

int background_stability_tests_smg(
        struct background *pba,
        double * pvecback,
        double * pvecback_integration
        );

int background_hi_class_second_loop(
        struct background *pba,
        double * pvecback
        );

int background_hi_class_third_loop(
        struct background *pba,
        double * pvecback,
        double * pvecback_integration
        );

int background_initial_conditions_smg(
        struct background *pba,
        double * pvecback,
        double * pvecback_integration,
        double rho_rad
        );

int hi_class_store_columntitles(
				struct background *pba,
				char titles[_MAXTITLESTRINGLENGTH_]
			  );

int hi_class_store_doubles(
				struct background *pba,
        double * pvecback,
        double * dataptr,
        int * storeidx
			  );


#ifdef __cplusplus
}
#endif

#endif
