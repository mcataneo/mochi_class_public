#ifndef __BACKGROUND_SMG__
#define __BACKGROUND_SMG__

#include "common.h"
#include "background.h"
#include "rootfinder.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int background_gravity_functions_smg(
        struct background *pba,
        double a,
        double * pvecback_B,
        short return_format,
        double * pvecback,
        double * ptr_rho_tot,
        double * ptr_p_tot,
        double * ptr_rho_de
        );

int hi_class_define_indices_bg(
				struct background *pba,
				int * index_bg
			  );

int hi_class_define_indices_bi(
				struct background *pba,
				int * index_bi
			  );

int background_derivs_alphas_smg(
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
        double a,
        double * pvecback,
        double * pvecback_integration,
        double * ptr_rho_rad
        );

int background_store_columntitles_smg(
				struct background *pba,
				char titles[_MAXTITLESTRINGLENGTH_]
			  );

int background_store_doubles_smg(
				struct background *pba,
        double * pvecback,
        double * dataptr,
        int * ptr_storeidx
			  );

int background_gravity_parameters(
			  struct background *pba
			  );

int background_free_smg(
			  struct background *pba
			  );

int background_print_smg(
			  struct background *pba,
				double * pvecback,
        double * pvecback_integration
			);

int background_derivs_smg(
			  struct background *pba,
        double a,
				double * pvecback,
				double * y,
				double * dy
			);

#ifdef __cplusplus
}
#endif

#endif
