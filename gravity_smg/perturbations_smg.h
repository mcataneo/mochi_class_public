#ifndef __PERTURBATIONS_SMG__
#define __PERTURBATIONS_SMG__

#include "common.h"
#include "perturbations.h"

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int perturb_tests_smg(struct precision * ppr,
                     struct background * pba,
                     struct perturbs * ppt);

int perturb_qs_functions_at_tau_and_k_qs_smg(struct background * pba,
                                            struct perturbs * ppt,
                                            double k,
                                            double tau,
                                            double *mass2,
                                            double *mass2_p,
                                            double *rad2,
                                            double *friction,
                                            double *slope);

int perturb_test_at_k_qs_smg(struct precision * ppr,
                            struct background * pba,
                            struct perturbs * ppt,
                            double k,
                            double tau,
                            int *approx);

int perturb_test_ini_qs_smg(struct precision * ppr,
                           struct background * pba,
                           struct perturbs * ppt,
                           double k_min,
                           double k_max,
                           double a_ini);

int perturb_find_scheme_qs_smg(struct precision * ppr,
                              struct background * pba,
                              struct perturbs * ppt,
                              struct perturb_workspace * ppw,
                              double k,
                              double tau_ini,
                              double tau_end);

int sample_functions_qs_smg(struct precision * ppr,
                           struct background * pba,
                           struct perturbs * ppt,
                           double k,
                           double tau_ini,
                           double tau_end,
                           double * tau_sample,
                           double * mass2_sample,
                           double * rad2_sample,
                           double * slope_sample,
                           int *size_sample);

int functions_to_approx_qs_smg(struct precision * ppr,
                              struct background * pba,
                              struct perturbs * ppt,
                              double tau_ini,
                              double tau_end,
                              double * tau_sample,
                              double * mass_sample,
                              double * rad_sample,
                              int * approx_sample,
                              int size_sample);

int shorten_first_qs_smg(double * tau_sample,
                        double * slope_sample,
                        int * approx_sample,
                        int size_sample,
                        double * tau_array,
                        double * slope_array,
                        int * approx_array,
                        int *size_array,
                        double tau_end);

int correct_with_slope_qs_smg(struct precision * ppr,
                             struct background * pba,
                             struct perturbs * ppt,
                             double tau_ini,
                             double tau_end,
                             double * tau_array,
                             double * slope_array,
                             int * approx_array,
                             int size_array);

int shorten_second_qs_smg(double * tau_array,
                         int * approx_array,
                         int size_array,
                         double * tau_scheme,
                         int * approx_scheme,
                         int *size_scheme);

int fit_real_scheme_qs_smg(double tau_end,
                          int * approx_scheme,
                          double * tau_scheme,
                          int size_scheme,
                          double * tau_export);

int perturb_test_ini_grav_ic_smg(struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt);

int perturb_test_ini_extfld_ic_smg(struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt);

int calc_extfld_ampl(int n,  double kin, double bra, double dbra, double run, double ten, double DelM2,
                      double Omx, double wx, double l1, double l2, double l3, double l4,
                      double l5, double l6,double l7,double l8, double cs2num, double Dd, double ic_regulator_smg,
                      double * amplitude);

int get_gravity_coefficients_smg(
                               struct perturbs * ppt,
                               struct background * pba,
                               double * pvecback,
                               double * delM2, double * M2, double * kin, double * bra,
                               double * ten, double * run, double * beh, double * res,
                               double * cD, double * cK, double * cB, double * cH, double * c0,
                               double * c1, double * c2, double * c3, double * c4,
                               double * c5, double * c6, double * c7, double * c8,
                               double * c9, double * c10, double * c11, double * c12,
                               double * c13, double * c14, double * c15, double * c16,
                               double * res_p, double *  cD_p, double *  cB_p, double *  cH_p,
                               double * c9_p, double * c10_p, double * c12_p, double * c13_p
                             );

int get_x_x_prime_qs_smg(
                        struct precision * ppr,
                        struct background * pba,
                        struct perturbs * ppt,
                        struct perturb_workspace * ppw,
                        double k, double * x_qs_smg, double * x_prime_qs_smg
                        );

int hi_class_define_indices_tp(
         struct perturbs * ppt,
				 int * index_type
			 );

int hi_class_define_indices_mt(
        struct perturb_workspace * ppw,
				int * index_mt
			 );

int perturb_hi_class_qs(struct precision * ppr,
                       struct background * pba,
                       struct perturbs * ppt,
                       struct perturb_workspace * ppw,
                       double k,
                       double * tau_ini,
                       double tau_end);

int perturb_store_columntitles_smg(
				struct perturbs * ppt
			  );

// int perturb_store_doubles_smg(
// 				struct background *pba,
//        double * pvecback,
//        double * dataptr,
//        int * ptr_storeidx
// 			  );

int perturb_verbose_qs_smg(
				struct perturbs * ppt,
        struct perturb_workspace * ppw,
        double k,
        double tau_switch,
        int * ap_ini,
        int * ap_end
			  );

int hi_class_define_indices_pt(
      struct perturb_workspace * ppw,
      struct perturb_vector * ppv,
			int * index_pt
			);

int perturb_vector_init_smg(
      struct perturb_workspace * ppw,
      struct perturb_vector * ppv,
      int * pa_old
			);

int perturb_adiabatic_ic_smg(
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct perturb_workspace * ppw,
      double * ptr_eta,
      double * ptr_delta_ur,
      double * ptr_theta_ur,
      double * ptr_shear_ur,
      double * ptr_l3_ur,
      double * ptr_delta_dr,
      double tau,
      double k,
      double fracnu,
      double om,
      double rho_r
			 );

int perturb_isocurvature_cdm_ic_smg(
     struct precision * ppr,
     struct background * pba,
     struct perturbs * ppt,
     struct perturb_workspace * ppw,
     double tau,
     double k,
     double fraccdm,
     double om
			 );

int perturb_isocurvature_b_ic_smg(
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct perturb_workspace * ppw,
    double tau,
    double k,
    double fracb,
    double om
			 );

int perturb_isocurvature_urd_ic_smg(
   struct precision * ppr,
   struct background * pba,
   struct perturbs * ppt,
   struct perturb_workspace * ppw,
   double tau,
   double k,
   double fracnu,
   double fracg,
   double fracb,
   double om
			 );

int perturb_isocurvature_urv_ic_smg(
  struct precision * ppr,
  struct background * pba,
  struct perturbs * ppt,
  struct perturb_workspace * ppw,
  double tau,
  double k,
  double fracnu,
  double fracg,
  double fracb,
  double om
			 );

int perturb_get_x_x_prime_newtonian(
 struct perturb_workspace * ppw
			 );

int perturb_get_h_prime_ic_from_00(
  struct background * pba,
  struct perturb_workspace * ppw,
  double k,
  double eta,
  double delta_rho_tot
			 );

int perturb_approximations_smg(
  struct perturbs * ppt,
  struct perturb_workspace * ppw,
  double tau
);

int perturb_einstein_smg(
  struct precision * ppr,
  struct background * pba,
  struct thermo * pth,
  struct perturbs * ppt,
  struct perturb_workspace * ppw,
  double k,
  double tau,
  double * y
);

int perturb_einstein_tensor_smg(
  struct background * pba,
  struct perturb_workspace * ppw,
  double k,
  double tau,
  double * y
);

int perturb_print_variables_smg(
  struct background * pba,
  struct perturbs * ppt,
  struct perturb_workspace * ppw,
  double k,
  double tau,
  double * dataptr,
  int * ptr_storeidx
);

int perturb_derivs_smg(
  struct perturbs * ppt,
  struct perturb_workspace * ppw,
  struct perturb_vector * pv,
  double * dy,
  double * pvecmetric
);

#ifdef __cplusplus
}
#endif

#endif
