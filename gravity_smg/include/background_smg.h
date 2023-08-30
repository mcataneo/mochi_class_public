#ifndef __BACKGROUND_SMG__
#define __BACKGROUND_SMG__

#include "common.h"
#include "background.h"
#include "rootfinder.h"
#include "gravity_functions_smg.h"
#include "gravity_models_smg.h"


/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif


int background_gravity_functions_smg(struct background *pba,
                                     double a,
                                     double * pvecback_B,
                                     double * pvecback,
                                     double * ptr_rho_tot,
                                     double * ptr_p_tot,
                                     double * ptr_rho_de);

int background_free_smg(struct background *pba);

int background_define_indices_bg_smg(struct background *pba,
				                             int * index_bg);

int background_define_indices_bi_smg(struct background *pba,
				                             int * index_bi);

int background_define_indices_bibw_smg(struct background *pba,
				                             int * index_bibw);

int background_solve_smg(struct precision *ppr,
                         struct background *pba,
                         double * pvecback,
                         double * pvecback_integration,
                         double * pvecback_bw_integration);

int background_solve_rho_smg(struct precision *ppr,
                             struct background *pba,
                             double * pback_rho_smg_bw_integration);

int interpolate_rho_smg_p_smg(struct background *pba,
                        double loga,
                        double loga_transition,
                        double * pvecback
                        );

int background_print_stdout_smg(struct background *pba,
				                        double * pvecback,
                                double * pvecback_integration);

int background_initial_conditions_smg(struct background *pba,
                                      double a,
                                      double * pvecback,
                                      double * pvecback_integration,
                                      double * pvecback_bw_integration,
                                      double * ptr_rho_rad);

int background_ic_rho_smg(struct background *pba,
                          double * pback_rho_smg_bw_integration);

int background_derivs_bw_smg(
                        double loga,
                        double * y,
                        double * dy,
                        void * parameters_and_workspace,
                        ErrorMsg error_message
                        );

int background_sources_bw_smg(
                       double loga,
                       double * y,
                       double * dy,
                       int index_loga,
                       void * parameters_and_workspace,
                       ErrorMsg error_message
                       );

int background_derivs_bw_rho_smg(
                        double loga,
                        double * y,
                        double * dy,
                        void * parameters_and_workspace,
                        ErrorMsg error_message
                        );

int background_sources_bw_rho_smg(
                       double loga,
                       double * y,
                       double * dy,
                       int index_loga,
                       void * parameters_and_workspace,
                       ErrorMsg error_message
                       );

int background_store_columntitles_smg(struct background *pba,
				                              char titles[_MAXTITLESTRINGLENGTH_]);

int background_output_data_smg(struct background *pba,
                               double * pvecback,
                               double * dataptr,
                               int * ptr_storeidx);

int background_derivs_smg(struct background *pba,
                          double a,
				                  double * pvecback,
				                  double * y,
				                  double * dy);

int stability_tests_smg(struct background *pba,
                        double * pvecback,
                        double * pvecback_integration);

int derivatives_alphas_smg(struct background *pba,
                           double * pvecback,
                           double * pvecback_derivs,
                           int i);

int copy_to_background_table_smg(struct background *pba,
                                 int row,
                                 int column,
                                 double value);


#ifdef __cplusplus
}
#endif

#endif
