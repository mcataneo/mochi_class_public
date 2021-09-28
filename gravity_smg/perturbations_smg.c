#include "perturbations_smg.h"


// //Here we calculate adiabatic initial conditions.
// int perturb_tests_smg(struct precision * ppr,
//                       struct background * pba,
//                       struct perturbs * ppt) {
//
//         return _SUCCESS_;
// }

//Here we do the smg tests.
int perturb_tests_smg(struct precision * ppr,
                      struct background * pba,
                      struct perturbs * ppt) {

        class_test(ppt->gauge == newtonian,
                   ppt->error_message,
                   "Asked for scalar modified gravity AND Newtonian gauge. Not yet implemented");

        class_test(ppr->a_ini_test_qs_smg < ppr->a_ini_over_a_today_default,
                   ppt->error_message,
                   "The initial time for testing the QS approximation (qs_smg) must be bigger than the background initial time (a_ini_test_qs_smg>=a_ini_over_a_today_default).");

        if ( ppt->pert_initial_conditions_smg == gravitating_attr ) {
                class_test_except((ppt->has_cdi == _TRUE_) || (ppt->has_bi == _TRUE_) || (ppt->has_nid == _TRUE_) || (ppt->has_niv == _TRUE_),
                                  ppt->error_message,
                                  perturb_free_nosource(ppt),
                                  "Isocurvature initial conditions for early modified gravity (Gravitating Attractor) not implemented.");
        }


        // TODO: think of some suitable tests for the scalar field

        if (ppt->method_qs_smg == automatic) {
                //Check if at the initial time all the k modes start with the same kind of qs_smg approximation
                class_call_except(perturb_test_ini_qs_smg(ppr,
                                                          pba,
                                                          ppt,
                                                          ppt->k[ppt->index_md_scalars][0],
                                                          ppt->k[ppt->index_md_scalars][ppt->k_size[ppt->index_md_scalars]-1],
                                                          ppr->a_ini_test_qs_smg),
                                  ppt->error_message,
                                  ppt->error_message,
                                  perturb_free_nosource(ppt));
        }

        if (!((ppt->method_qs_smg == automatic) && (ppt->initial_approx_qs_smg==_TRUE_))) {

                // if scalar is dynamical or always quasi-static, test for stability at the initial time.
                // Only in the case it is QS because of a trigger test (through "automatic" method_qs),
                // we already know mass is positive and therefore can assume it is stable, so skip this.

                if( ppt->pert_initial_conditions_smg == gravitating_attr) {
                        // If we are in gravitating_attr ICs, make sure the standard solution is dominant at some early redshift.
                        // If it is not, curvature is not conserved and we have lost the connection between the amplitude from inflation and
                        // the initial amplitude supplied to hi_class.
                        class_call_except(perturb_test_ini_grav_ic_smg(ppr,
                                                                       pba,
                                                                       ppt),
                                          ppt->error_message,
                                          ppt->error_message,
                                          perturb_free_nosource(ppt));
                }
                else if( ppt->pert_initial_conditions_smg == ext_field_attr) {
                        //If we have the ext_field_attr, test for tachyon instability in RD before pert initialisation
                        // If have it, fail, because we can't set the ICs properly

                        class_call_except(perturb_test_ini_extfld_ic_smg(ppr,
                                                                         pba,
                                                                         ppt),
                                          ppt->error_message,
                                          ppt->error_message,
                                          perturb_free_nosource(ppt));
                }
        }

        return _SUCCESS_;

}

int perturb_qs_functions_at_tau_and_k_qs_smg(struct background * pba,
                                             struct perturbs * ppt,
                                             double k,
                                             double tau,
                                             double *mass2,
                                             double *mass2_p,
                                             double *rad2,
                                             double *friction,
                                             double *slope) {

        /* Definition of local variables */
        double mass2_qs, mass2_qs_p, rad2_qs, friction_qs, slope_qs;
        double * pvecback;
        int first_index_back;

        class_alloc(pvecback,pba->bg_size*sizeof(double),ppt->error_message);
        class_call(background_at_tau(pba,
                                     tau,
                                     pba->normal_info,
                                     pba->inter_normal,
                                     &first_index_back,
                                     pvecback),
                   pba->error_message,
                   ppt->error_message);

        double delM2, M2, kin, bra, ten, run, beh;
        double res, cD, cK, cB, cH;
        double c0, c1, c2, c3, c4, c5, c6, c7, c8;
        double c9, c10, c11, c12, c13, c14, c15, c16;
        double c9_p, c10_p, c12_p, c13_p;
        double res_p, cD_p, cB_p, cH_p;
        double x_prime_qs_smg_num, x_prime_qs_smg_den;
        double a, H, rho_tot, p_tot, rho_smg, p_smg, rho_r;
        double k2 = k*k;

        a = pvecback[pba->index_bg_a];
        H = pvecback[pba->index_bg_H];
        rho_r = pvecback[pba->index_bg_rho_g] + pvecback[pba->index_bg_rho_ur];
        rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
        p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
        rho_smg = pvecback[pba->index_bg_rho_smg];
        p_smg = pvecback[pba->index_bg_p_smg];

        class_call(
                get_gravity_coefficients_smg(
                        ppt, pba, pvecback,
                        &delM2, &M2, &kin, &bra, &ten, &run, &beh, &res,
                        &cD, &cK, &cB, &cH, &c0, &c1, &c2, &c3,
                        &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11,
                        &c12, &c13, &c14, &c15, &c16, &res_p, &cD_p, &cB_p,
                        &cH_p, &c9_p, &c10_p, &c12_p, &c13_p
                        ),
                ppt->error_message,
                ppt->error_message);


        mass2_qs = -(c12 + c13*k2*pow(a*H,-2))/cD;

        mass2_qs_p =
                -(
                        +c12_p - c12*cD_p/cD
                        + (c13_p - c13*cD_p/cD + (rho_tot + rho_smg + 3.*p_tot + 3.*p_smg)*c13*a/H)*pow(a*H,-2)*k2
                        )/cD;

        rad2_qs = 3.*mass2_qs*pow(H,4)*pow(rho_r,-2)*pow(a*H,2)/k2;

        friction_qs = -(c11 - c3*k2*pow(a*H,-2))/cD;

        slope_qs = -1./4.*(1. - 2.*friction_qs + 3.*(p_tot + p_smg)/(rho_tot + rho_smg) - mass2_qs_p/mass2_qs/a/H);

        *mass2 = mass2_qs;
        *mass2_p = mass2_qs_p;
        *rad2 = rad2_qs;
        *friction = friction_qs;
        *slope = slope_qs;

        free(pvecback);

        return _SUCCESS_;

}

int perturb_test_at_k_qs_smg(struct precision * ppr,
                             struct background * pba,
                             struct perturbs * ppt,
                             double k,
                             double tau,
                             int *approx) {

        //Define local variables
        double mass2_qs, mass2_qs_p, rad2_qs, friction_qs, slope_qs;

        perturb_qs_functions_at_tau_and_k_qs_smg(
                pba,
                ppt,
                k,
                tau,
                &mass2_qs,
                &mass2_qs_p,
                &rad2_qs,
                &friction_qs,
                &slope_qs);

        double tau_fd;
        short proposal;

        class_call(background_tau_of_z(pba,
                                       ppr->z_fd_qs_smg,
                                       &tau_fd),
                   pba->error_message,
                   ppt->error_message);
        //Approximation
        if ((mass2_qs > pow(ppr->trigger_mass_qs_smg,2)) && (rad2_qs > pow(ppr->trigger_rad_qs_smg,2))) {
                proposal = _TRUE_;
        }
        else {
                proposal = _FALSE_;
        }
        if (tau <= tau_fd) {
                *approx = proposal;
        }
        else {
                *approx = _FALSE_;
        }

        return _SUCCESS_;

}

int perturb_test_ini_qs_smg(struct precision * ppr,
                            struct background * pba,
                            struct perturbs * ppt,
                            double k_min,
                            double k_max,
                            double a_ini) {
        //Define local variables
        double * pvecback;
        int first_index_back;
        double tau;
        int approx_k_min, approx_k_max;

        //Get background quantities at a_ini
        class_call(background_tau_of_z(pba,
                                       1./a_ini-1.,
                                       &tau),
                   pba->error_message,
                   ppt->error_message);

        //Approximation for k_min
        perturb_test_at_k_qs_smg(
                ppr,
                pba,
                ppt,
                k_min,
                tau,
                &approx_k_min
                );

        //Approximation for k_max
        perturb_test_at_k_qs_smg(
                ppr,
                pba,
                ppt,
                k_max,
                tau,
                &approx_k_max
                );

        class_test_except(approx_k_min != approx_k_max,
                          ppt->error_message,
                          free(pvecback),
                          "\n All the k modes should start evolving with the same type of initial conditions (either fully_dynamic or quasi_static).\n This is not the case at a = %e. Try to decrease a_ini_over_a_today_default.\n", ppr->a_ini_over_a_today_default);

        ppt->initial_approx_qs_smg = approx_k_min;

        free(pvecback);

        return _SUCCESS_;

}

int perturb_find_scheme_qs_smg(struct precision * ppr,
                               struct background * pba,
                               struct perturbs * ppt,
                               double k,
                               double tau_ini,
                               double tau_end,
                               double * tau_export) {

        int size_sample = ppr->n_max_qs_smg;

        double * tau_sample;
        double * mass2_sample;
        double * rad2_sample;
        double * slope_sample;

        class_alloc(tau_sample,size_sample*sizeof(double),ppt->error_message);
        class_alloc(mass2_sample,size_sample*sizeof(double),ppt->error_message);
        class_alloc(rad2_sample,size_sample*sizeof(double),ppt->error_message);
        class_alloc(slope_sample,size_sample*sizeof(double),ppt->error_message);

        /**
         * Input: background table
         * Output: sample of the time, mass, decaying rate of the oscillations (slope)
         *   and radiation density.
         **/
        sample_functions_qs_smg(
                ppr,
                pba,
                ppt,
                k,
                tau_ini,
                tau_end,
                tau_sample,
                mass2_sample,
                rad2_sample,
                slope_sample,
                &size_sample
                );


        int * approx_sample;

        class_alloc(approx_sample,size_sample*sizeof(int),ppt->error_message);

        /**
         * Input: sample of the time, mass and radiation density
         * Output: sample of the approx scheme
         **/
        functions_to_approx_qs_smg(
                ppr,
                pba,
                ppt,
                tau_ini,
                tau_end,
                tau_sample,
                mass2_sample,
                rad2_sample,
                approx_sample,
                size_sample
                );

        free(mass2_sample);
        free(rad2_sample);


        double * tau_array;
        double * slope_array;
        int * approx_array;
        int size_array = size_sample;

        class_alloc(tau_array,size_array*sizeof(double),ppt->error_message);
        class_alloc(slope_array,size_array*sizeof(double),ppt->error_message);
        class_alloc(approx_array,size_array*sizeof(int),ppt->error_message);

        /**
         * Input: sample of the time, slope and approximation scheme
         *   at small time interval
         * Output: arrays containing the time, the slope and the approximation
         * scheme only when it changes
         **/
        shorten_first_qs_smg(
                tau_sample,
                slope_sample,
                approx_sample,
                size_sample,
                tau_array,
                slope_array,
                approx_array,
                &size_array,
                tau_end
                );

        free(tau_sample);
        free(slope_sample);
        free(approx_sample);

        /**
         * Input: arrays with time, slope and approximation schemes
         * Output: arrays with time and approximation scheme corrected with the slope
         **/
        correct_with_slope_qs_smg(
                ppr,
                pba,
                ppt,
                tau_ini,
                tau_end,
                tau_array,
                slope_array,
                approx_array,
                size_array
                );

        free(slope_array);

        double * tau_scheme;
        int * approx_scheme;
        int size_scheme = size_array;

        class_alloc(tau_scheme,size_scheme*sizeof(double),ppt->error_message);
        class_alloc(approx_scheme,size_scheme*sizeof(int),ppt->error_message);

        /**
         * Input: arrays of time and approximation after correcting with the slope
         *   (there is the possibility that the same number in approx_array is repeated)
         * Output: shortened arrays of time and approximation
         **/
        shorten_second_qs_smg(
                tau_array,
                approx_array,
                size_array,
                tau_scheme,
                approx_scheme,
                &size_scheme
                );

        free(tau_array);
        free(approx_array);

        /**
         * Input: real approx_scheme and tau_scheme
         * Output: approx scheme (tau_export) adjusted to fit the implemented one
         **/
        fit_real_scheme_qs_smg(
                tau_end,
                approx_scheme,
                tau_scheme,
                size_scheme,
                tau_export
                );

        free(tau_scheme);
        free(approx_scheme);

//   // DEBUG: Initial and final times
//   printf("6 - Interval tau       = {%.1e, %.1e}\n", tau_ini, tau_end);
//   printf("7 - k mode             = {%.1e}\n", k);


        return _SUCCESS_;

}


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
                            int *size_sample) {

        /* Definition of local variables */
        double mass2_qs, mass2_qs_p, rad2_qs, friction_qs, slope_qs;
        double tau = tau_ini;
        double delta_tau = (tau_end - tau_ini)/ppr->n_max_qs_smg;
        int count = 0;


        /* Scan the time evolution and build several arrays containing
        * interesting quantities for the quasi-static approximation */
        while (tau < tau_end) {

                perturb_qs_functions_at_tau_and_k_qs_smg(
                        pba,
                        ppt,
                        k,
                        tau,
                        &mass2_qs,
                        &mass2_qs_p,
                        &rad2_qs,
                        &friction_qs,
                        &slope_qs);

//     DEBUG: To debug uncomment this and define a convenient function of time for each of these quantities
//     double x = (tau - tau_ini)/(tau_end - tau_ini);
//     mass2_qs = 1.5 + cos(10*_PI_*x);
//     rad2_qs = 1.;
//     slope_qs = 1.;

                tau_sample[count] = tau;
                mass2_sample[count] = mass2_qs;
                rad2_sample[count] = rad2_qs;
                slope_sample[count] = slope_qs;

                delta_tau = fabs(2.*mass2_qs/mass2_qs_p)/sqrt(ppr->n_min_qs_smg*ppr->n_max_qs_smg);
                delta_tau = MIN(delta_tau, (tau_end - tau_ini)/ppr->n_min_qs_smg);
                delta_tau = MAX(delta_tau, (tau_end - tau_ini)/ppr->n_max_qs_smg);

                tau += delta_tau;
                count += 1;

        }

        *size_sample = count;

        return _SUCCESS_;

}


int functions_to_approx_qs_smg(struct precision * ppr,
                               struct background * pba,
                               struct perturbs * ppt,
                               double tau_ini,
                               double tau_end,
                               double * tau_sample,
                               double * mass2_sample,
                               double * rad2_sample,
                               int * approx_sample,
                               int size_sample
                               ) {


        // Convert the input parameter z_fd into the corresponding conformal time
        double tau_fd;
        short proposal;

        class_call(background_tau_of_z(pba,
                                       ppr->z_fd_qs_smg,
                                       &tau_fd),
                   pba->error_message,
                   ppt->error_message);

        int i;
        for (i = 0; i < size_sample; i++) {

                if ((mass2_sample[i] > pow(ppr->trigger_mass_qs_smg,2)) && (rad2_sample[i] > pow(ppr->trigger_rad_qs_smg,2))) {
                        proposal = 1;
                }
                else {
                        proposal = 0;
                }

                if (tau_sample[i] <= tau_fd) {
                        approx_sample[i] = proposal;
                }
                else {
                        approx_sample[i] = 0;
                }

        }

        return _SUCCESS_;

}


int shorten_first_qs_smg(double * tau_sample,
                         double * slope_sample,
                         int * approx_sample,
                         int size_sample,
                         double * tau_array,
                         double * slope_array,
                         int * approx_array,
                         int *size_array,
                         double tau_end) {

        int i, j, count = 0;
        int last_switch = 0;
        int this_switch = 0;
        double slope_weighted;

        tau_array[0] = tau_sample[0];
        approx_array[0] = approx_sample[0];

        for (i = 1; i < size_sample; i++) {
                if (approx_sample[i] != approx_array[count]) {

                        count += 1;
                        // Shorten approximation scheme
                        approx_array[count] = approx_sample[i];

                        // Shorten time
                        if (approx_array[count-1] < approx_array[count]) {
                                this_switch = i;
                        }
                        else {
                                this_switch = i-1;
                        }
                        tau_array[count] = tau_sample[this_switch];

                        // Shorten slope
                        slope_weighted = 0.;
                        for (j = last_switch; j < this_switch; j++) {
                                slope_weighted += slope_sample[j]*(tau_sample[j+1] - tau_sample[j]);
                        }
                        slope_array[count -1] = slope_weighted/(tau_sample[this_switch] - tau_sample[last_switch]);
                        last_switch = this_switch;
                }
        }

        // Shorten slope last element
        slope_weighted = 0.;
        for (i = last_switch; i < size_sample-1; i++) {
                slope_weighted += slope_sample[i]*(tau_sample[i+1] - tau_sample[i]);
        }
        slope_array[count] = (slope_weighted + slope_sample[size_sample-1]*(tau_end - tau_sample[size_sample-1]))/(tau_end - tau_sample[last_switch]);

        *size_array = count + 1;

        return _SUCCESS_;

}


int correct_with_slope_qs_smg(struct precision * ppr,
                              struct background * pba,
                              struct perturbs * ppt,
                              double tau_ini,
                              double tau_end,
                              double * tau_array,
                              double * slope_array,
                              int * approx_array,
                              int size_array) {


        double * pvecback;
        int first_index_back;
        int i, j, count;
        for (i = 1; i < size_array; i++) {
                if ((approx_array[i-1] == 0) && (approx_array[i] == 1)) {

                        // Routine to calculate the time interval necessary to relax the oscillations
                        class_alloc(pvecback,pba->bg_size*sizeof(double),ppt->error_message);
                        class_call(background_at_tau(pba,
                                                     tau_array[i],
                                                     pba->short_info,
                                                     pba->inter_normal,
                                                     &first_index_back,
                                                     pvecback),
                                   pba->error_message,
                                   ppt->error_message);

                        double a_final = pvecback[pba->index_bg_a] * pow(ppr->eps_s_qs_smg, -1./slope_array[i]);
                        double tau_final;

                        class_call(background_tau_of_z(pba,
                                                       1./a_final-1.,
                                                       &tau_final),
                                   pba->error_message,
                                   ppt->error_message);

                        double delta_tau = tau_final - tau_array[i];

                        // Adjust time and approx to take into account the oscillations
                        double next_tau;
                        if (i+1<size_array) {
                                next_tau = tau_array[i+1];
                        }
                        else {
                                next_tau = tau_end;
                        }

                        if (tau_array[i] + delta_tau < next_tau) {
                                tau_array[i] += delta_tau;
                        }
                        else {
                                approx_array[i] = 0;
                        }

                        free(pvecback);
                }
        }

        return _SUCCESS_;

}


int shorten_second_qs_smg(double * tau_array,
                          int * approx_array,
                          int size_array,
                          double * tau_scheme,
                          int * approx_scheme,
                          int *size_scheme) {

        tau_scheme[0] = tau_array[0];
        approx_scheme[0] = approx_array[0];
        int i, j = 0;

        for (i = 0; i < size_array; i++) {
                if (approx_array[i] != approx_scheme[j]) {
                        j += 1;
                        approx_scheme[j] = approx_array[i];
                        tau_scheme[j] = tau_array[i];
                }
        }

        *size_scheme = j + 1;

        return _SUCCESS_;

}


int fit_real_scheme_qs_smg(double tau_end,
                           int * approx_scheme,
                           double * tau_scheme,
                           int size_scheme,
                           double * tau_export
                           ) {

        /* Definition of local variables */
        int implemented_scheme[] = _VALUES_QS_SMG_FLAGS_;
        int size_implemented_scheme = sizeof(implemented_scheme)/sizeof(int);

        int i, j;
        int start_position = 0;
        short scheme_fits = _FALSE_;

//   // DEBUG: print the implemented scheme
//   int count;
//   printf("1 - Implemented scheme = {");
//   for (count = 0; count < size_implemented_scheme; count++) {
//     printf("%d", implemented_scheme[count]);
//     if (count < size_implemented_scheme - 1) {
//       printf(", ");
//     }
//   }
//   printf("}\n");
//   // DEBUG: print the real scheme
//   printf("2 - Real scheme        = {");
//   for (count = 0; count < size_scheme; count++) {
//     printf("%d", approx_scheme[count]);
//     if (count < size_scheme - 1) {
//       printf(", ");
//     }
//   }
//   printf("}\n");


        while (scheme_fits == _FALSE_) {

                /* Check if the real approximation scheme fits the implemented one */
                for (i = 0; i < size_implemented_scheme - size_scheme + 1; i++) {
                        j = 0;
                        while (j < size_scheme - 1) {
                                if (approx_scheme[j] == implemented_scheme[i + j]) {
                                        j += 1;
                                }
                                else {
                                        break;
                                }
                        }
                        if ((j == size_scheme - 1) && (approx_scheme[j] == implemented_scheme[i + j])) {
                                start_position = i;
                                scheme_fits = _TRUE_;
                                break;
                        }
                }

                /* Shorten the real approximation scheme */
                if (scheme_fits == _FALSE_) {
                        if ((approx_scheme[size_scheme - 2]==0) && (approx_scheme[size_scheme - 1]==1)) {
                                size_scheme += -1;
                        }
                        else if ((approx_scheme[size_scheme - 3]==0) && (approx_scheme[size_scheme - 2]==1) && (approx_scheme[size_scheme - 1]==0)) {
                                size_scheme += -2;
                        }
                }
        }

        /* Generate the vector of times at which the approximation switches */
        for (i = 0; i < size_implemented_scheme; i++) {
                tau_export[i] = -1.;
        }

        for (i = 0; i < size_scheme; i++) {
                tau_export[start_position + i] = tau_scheme[i];
        }

        for (i = start_position + size_scheme; i < size_implemented_scheme; i++) {
                tau_export[i] = tau_end + 1.; // The +1 is here to make the final elements larger than everything else
        }

//   // DEBUG: print the fitted scheme
//   printf("3 - Fitted scheme      = {");
//   for (count = 0; count < size_scheme; count++) {
//     printf("%d", approx_scheme[count]);
//     if (count < size_scheme - 1) {
//       printf(", ");
//     }
//   }
//   printf("}\n");
//   // DEBUG: print the real tau switches
//   printf("4 - Real tau           = {");
//   for (count = 0; count < size_scheme; count++) {
//     printf("%.1e", tau_scheme[count]);
//     if (count < size_scheme - 1) {
//       printf(", ");
//     }
//   }
//   printf("}\n");
//   // DEBUG: print the tau switches after the fitting
//   printf("5 - Fitted tau         = {");
//   for (count = 0; count < size_implemented_scheme; count++) {
//     printf("%.1e", tau_export[count]);
//     if (count < size_implemented_scheme - 1) {
//       printf(", ");
//     }
//   }
//   printf("}\n");

        return _SUCCESS_;
}

/*
 * Test for stability of solutions in RD before initialisation of
 * perturbations: if standard solution not stable, cannot set ICs properly.
 */
int perturb_test_ini_grav_ic_smg(struct precision * ppr,
                                 struct background * pba,
                                 struct perturbs * ppt){
// test stability of gravitating_attr ICs

        double kin, bra, run, ten, DelM2, Omx, wx;
        double c3, c2, c1,c0, den1, den2, ic_regulator_smg;
        double tau_ini, z_ref;
        int i;
        double fastest_growth, wouldbe_adiab;
        double * pvecback;
        int first_index_back;
        double sols[3];
        int complex;

        class_alloc(pvecback,pba->bg_size*sizeof(double),ppt->error_message);

        z_ref = ppr->pert_ic_ini_z_ref_smg;

        class_call(background_tau_of_z(pba, z_ref,&tau_ini),
                   pba->error_message,
                   ppt->error_message);

        class_call(background_at_tau(pba,
                                     tau_ini,
                                     pba->long_info,
                                     pba->inter_normal,
                                     &first_index_back,
                                     pvecback),
                   pba->error_message,
                   ppt->error_message);

        // define alphas
        wx = pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg];
        Omx = pvecback[pba->index_bg_rho_smg]/pow(pvecback[pba->index_bg_H],2);
        kin = pvecback[pba->index_bg_kineticity_smg];
        bra = pvecback[pba->index_bg_braiding_smg];
        run = pvecback[pba->index_bg_mpl_running_smg];
        ten = pvecback[pba->index_bg_tensor_excess_smg];
        DelM2 = pvecback[pba->index_bg_delta_M2_smg];//M2-1

        /* Determine the solutions
         *
         *   h = C * tau^2+x
         *
         * where n is the solution of
         *
         *   c3 x^3 + c2 x^2 + c1 x + c0 = 0
         *
         * Note: if complex solutions then take the real part for the test
         * These solutions are exact when run=0. If not, then they are approximate.
         * These coefficients were obtain by solving the radiation+smg system
         * to obtain a 5th-order ODE for h. Removing two gauge modes, leaves
         * a cubic with three solutions relevant in the k->0 limit.
         * Note that we approximate the ci to O(run).
         */

        // Note: The denominators in the expressions below can be zero. We try to trap this and regulate.
        // We assume that M*^2>0 and D>0 which are tested for in the background routine.
        // Doing this gives wrong ICs, but it's better than segmentation faults.

        ic_regulator_smg =  ppr->pert_ic_regulator_smg;//  read in the minimum size that will get regulated
        ic_regulator_smg *= fabs(kin)+fabs(bra)+fabs(ten); //  scale it relative to the alphas

        c3  =   1.;

        c2  =   5. + 2.*run;

        den1 = (3.*bra*ten + kin*(2. + ten));

        if(ic_regulator_smg>0 &&(fabs(den1)<ic_regulator_smg)) {
                den1 = copysign(ic_regulator_smg,den1);
        }

        den2 =  4.*(9.*bra*(1. + DelM2) + (1. + DelM2)*kin - 12.*(DelM2 + Omx))*(3.*pow(bra,2.)*
                                                                                 (1. + DelM2) + 2.*kin*(DelM2 + Omx))*(-6.*(DelM2 + Omx)*(-2. + ten) + 9.*bra*(1. + DelM2)*(-1. + ten) + 2.*(1. + DelM2)*
                                                                                                                       kin*(1. + ten));

        if(ic_regulator_smg>0 &&(fabs(den2)<ic_regulator_smg)) {
                den2 = copysign(ic_regulator_smg,den2);
        }

        c2  +=  ((-1. + Omx)*run*(27.*pow(bra,4)*pow(1. + DelM2,2.)*ten*(432. - 373.*ten + 6.*pow(ten,2.)) + 9.*pow(bra,3.)*(1. + DelM2)*
                                  (864.*(DelM2 + Omx)*(-2. + ten)*ten + (1. + DelM2)*kin*(864. - 698.*ten - 329.*pow(ten,2.) + 6.*pow(ten,3.))) -
                                  3.*pow(bra,2.)*(1. + DelM2)*kin*(3456.*(DelM2 + Omx) - 4320.*(DelM2 + Omx)*ten + 6.*(1. + 446.*DelM2 + 445.*Omx)*
                                                                   pow(ten,2.) - 36.*(DelM2 + Omx)*pow(ten,3.) + (1. + DelM2)*kin*(768. + 227.*ten - 259.*pow(ten,2.) + 12.*pow(ten,3.))) -
                                  2.*pow(kin,2.)*(-6.*(DelM2 + Omx)*(-768.*(DelM2 + Omx) + (-1. + 191.*DelM2 + 192.*Omx)*pow(ten,2.)) + pow(1. + DelM2,2.)*
                                                  pow(kin,2.)*(-14. - 19.*ten - 4.*pow(ten,2.) + pow(ten,3.)) - (1. + DelM2)*kin*(-384.*(DelM2 + Omx) +
                                                                                                                                  (1. - 851.*DelM2 - 852.*Omx)*ten + (1. - 317.*DelM2 - 318.*Omx)*pow(ten,2.) + 6.*(DelM2 + Omx)*pow(ten,3.))) -
                                  6.*bra*kin*(-1152.*pow(DelM2 + Omx,2.)*(-2. + ten)*ten + pow(1. + DelM2,2.)*pow(kin,2.)*(-32. - 99.*ten - 40.*pow(ten,2.) +
                                                                                                                           3.*pow(ten,3.)) - (1. + DelM2)*kin*(1440.*(DelM2 + Omx) - 2.*(1. + 325.*DelM2 + 324.*Omx)*ten + (1. - 905.*DelM2 - 906.*Omx)*
                                                                                                                                                               pow(ten,2.) + 12.*(DelM2 + Omx)*pow(ten,3.)))))/(den2*den1);

        c1  =   (9*pow(bra,3)*pow(1 + DelM2,2)*(6 + 5*run)*ten + 3*pow(bra,2)*(1 + DelM2)*(-12*(-1 + Omx)*(-3 + run)*ten +
                                                                                           (1 + DelM2)*kin*(6 + 5*run)*(2 + ten)) + 6*bra*(-24*(-1 + Omx)*(DelM2 + Omx)*ten + (1 + DelM2)*kin*
                                                                                                                                           (12*(-1 + Omx) + (-2 + 6*DelM2 + 8*Omx + 5*(1 + DelM2)*run)*ten)) + 2*kin*((1 + DelM2)*kin*
                                                                                                                                                                                                                      (2*(2 + 3*DelM2 + Omx) + 5*(1 + DelM2)*run)*(2 + ten) - 12*(-1 + Omx)*((1 + DelM2)*run*ten + 2*(DelM2 + Omx)
                                                                                                                                                                                                                                                                                             *(2 + ten))))/(pow(1 + DelM2,2)*(3*pow(bra,2) + 2*kin)*den1);

        den2 = 4.*(1 + DelM2)*(9*bra*(1 + DelM2) + (1 + DelM2)*kin - 12*(DelM2 + Omx))*
               (3*pow(bra,2)*(1 + DelM2) + 2*kin*(DelM2 + Omx))*(-6*(DelM2 + Omx)*(-2 + ten) + 9*bra*(1 + DelM2)*(-1 + ten) +
                                                                 2*(1 + DelM2)*kin*(1 + ten));

        if(ic_regulator_smg>0 && (fabs(den2)<ic_regulator_smg)) {
                den2 = copysign(ic_regulator_smg,den2);
        }


        c1  +=  ((-1 + Omx)*run*(135*pow(bra,4)*pow(1 + DelM2,3)*ten*(288 - 229*ten + 6*pow(ten,2)) + 9*pow(bra,3)*
                                 pow(1 + DelM2,2)*(2880*(DelM2 + Omx)*(-2 + ten)*ten + (1 + DelM2)*kin*(3744 - 1780*ten - 1855*pow(ten,2) +
                                                                                                        66*pow(ten,3))) + 2*kin*(3456*pow(DelM2 + Omx,3)*(-2 + ten)*ten + 6*(1 + DelM2)*kin*(DelM2 + Omx)*(-2112*(DelM2 + Omx) -
                                                                                                                                                                                                           4*(1 + 25*DelM2 + 24*Omx)*ten + 3*(-1 + 95*DelM2 + 96*Omx)*pow(ten,2)) - pow(1 + DelM2,3)*pow(kin,3)*
                                                                                                                                 (-14 - 19*ten - 4*pow(ten,2) + pow(ten,3)) + pow(1 + DelM2,2)*pow(kin,2)*(-528*(DelM2 + Omx) +
                                                                                                                                                                                                           (1 - 1523*DelM2 - 1524*Omx)*ten + (1 - 545*DelM2 - 546*Omx)*pow(ten,2) + 18*(DelM2 + Omx)*pow(ten,3))) +
                                 3*pow(bra,2)*pow(1 + DelM2,2)*kin*((1 + DelM2)*kin*(-1296 - 2087*ten - 449*pow(ten,2) + 36*pow(ten,3)) +
                                                                    6*(-3072*(DelM2 + Omx) + 1532*(DelM2 + Omx)*ten + (-5 + 28*DelM2 + 33*Omx)*pow(ten,2) + 18*(DelM2 + Omx)*pow(ten,3))) -
                                 6*bra*(1 + DelM2)*kin*(576*pow(DelM2 + Omx,2)*(-4 + 5*ten) + pow(1 + DelM2,2)*pow(kin,2)*(-4 - 61*ten - 32*pow(ten,2) +
                                                                                                                           pow(ten,3)) - (1 + DelM2)*kin*(3552*(DelM2 + Omx) - 4*(1 + 121*DelM2 + 120*Omx)*ten - (1 + 1279*DelM2 + 1278*Omx)*
                                                                                                                                                          pow(ten,2) + 36*(DelM2 + Omx)*pow(ten,3)))))/(den2*den1);


        c0  =   (24*(-1 + Omx)*run*(4*kin*Omx - 3*pow(bra,2)*(-2 + ten) - DelM2*(3*pow(bra,2) + 2*kin)*(-2 + ten) +
                                    2*kin*(-2 + Omx)*ten + 6*bra*(-1 + Omx)*ten))/(pow(1 + DelM2,2)*(3*pow(bra,2) + 2*kin)*den1);

        den2 = (9*bra*(-1 + ten) + 2*(kin + 6*Omx + kin*ten - 3*Omx*ten) + DelM2*(9*bra*(-1 + ten) +
                                                                                  2*(6 + kin - 3*ten + kin*ten)));

        if(ic_regulator_smg>0 && (fabs(den2)<ic_regulator_smg)) {
                den2 = copysign(ic_regulator_smg,den2);
        }

        c0  +=  -((-1 + Omx)*run*(9*pow(bra,3)*(-288 + 98*ten + 119*pow(ten,2) + 6*pow(ten,3)) + 6*bra*(288*Omx*(-2 + ten)*
                                                                                                        ten + pow(kin,2)*(16 + 85*ten + 48*pow(ten,2) + 3*pow(ten,3)) + kin*(288*Omx + (314 - 216*Omx)*ten -
                                                                                                                                                                             (163 + 246*Omx)*pow(ten,2) - 12*(-1 + Omx)*pow(ten,3))) + 3*pow(bra,2)*(kin*(-480 + 335*ten + 383*pow(ten,2) +
                                                                                                                                                                                                                                                          18*pow(ten,3)) - 6*(-192*Omx + 48*(-3 + Omx)*ten + (85 + 83*Omx)*pow(ten,2) + 6*(-1 + Omx)*pow(ten,3))) +
                                  2*kin*(6*ten*(-192*Omx + (-1 + 97*Omx)*ten) + pow(kin,2)*(18 + 29*ten + 12*pow(ten,2) + pow(ten,3)) -
                                         kin*(192*Omx + (-107 + 300*Omx)*ten + (31 + 114*Omx)*pow(ten,2) + 6*(-1 + Omx)*pow(ten,3))) +
                                  DelM2*(9*pow(bra,3)*(-288 + 98*ten + 119*pow(ten,2) + 6*pow(ten,3)) + 2*kin*(576*(-2 + ten)*ten -
                                                                                                               kin*(192 + 193*ten + 145*pow(ten,2)) + pow(kin,2)*(18 + 29*ten + 12*pow(ten,2) + pow(ten,3))) +
                                         6*bra*(288*(-2 + ten)*ten + kin*(288 + 98*ten - 409*pow(ten,2)) + pow(kin,2)*(16 + 85*ten +
                                                                                                                       48*pow(ten,2) + 3*pow(ten,3))) + 3*pow(bra,2)*(-144*(-8 - 4*ten + 7*pow(ten,2)) + kin*(-480 + 335*ten +
                                                                                                                                                                                                              383*pow(ten,2) + 18*pow(ten,3))))))/(2.*(1 + DelM2)*(3*pow(bra,2) + 2*kin)*den1*den2);


        // Solve cubic to find the three solutions
        rf_solve_poly_3(c3,c2,c1,c0,sols,&complex);

        if (ppt->perturbations_verbose > 1) {
                printf("\nGravitating attractor ICs give growing modes at z=%e: \n (Approximate) polynomial",z_ref);
                printf(" solutions h ~ (k_tau)^n (complex = %i) with exponents: \n",complex);
        }

        fastest_growth = sols[0]; //want fastest
        wouldbe_adiab = sols[0]; //want closest to zero
        for (i=0; i<3; i+=1) {
                if (sols[i]  > fastest_growth) {
                        fastest_growth = sols[i];
                }
                if (fabs(sols[i]) < fabs(wouldbe_adiab)) {
                        wouldbe_adiab = sols[i];
                }
                if (ppt->perturbations_verbose > 1) {
                        printf("   n_%i = %f\n",i, 2+sols[i]);
                }
        }
        if (ppt->perturbations_verbose > 1) {
                printf("  fastest growing mode n = %f\n",2+fastest_growth);
                printf("  mode closest to standard adiabatic solution, n=%f\n",2+wouldbe_adiab);
                printf("  omx = %e, dM* = %e\n",Omx,DelM2);
        }

        // Check that would-be adiabatic mode is actually the fastest mode, otherwise
        // the would-be adiabatic attractor destabilises to the fastest mode, i.e. we cannot assume that the curvature was
        // conserved between inflation and the beginning of hi_class and therefore there is no
        // relation between the inflational amplitude A_S and the parameter we use for normalisation of curvature.

        /* We don't need this: te closest to zero mode actually conserves eta/zeta in any case
           class_test_except(ppr->pert_ic_tolerance_smg>0 && (fabs(wouldbe_adiab) > ppr->pert_ic_tolerance_smg),
                ppt->error_message,
                free(pvecback),
                "\n   Cannot set initial conditions for early_smg: adiabatic mode h ~ tau^2 lost, h ~ tau^n with n = %f",2+wouldbe_adiab);
         */

        if (fabs(fastest_growth)>fabs(wouldbe_adiab)) {
                class_test_except(ppr->pert_ic_tolerance_smg>0 && (fabs(fastest_growth) > ppr->pert_ic_tolerance_smg),
                                  ppt->error_message,
                                  free(pvecback),
                                  "\n   Cannot set initial conditions for early_smg:\n    There does exist a mode where curvature is conserved n=%f, but solution destabilises to a faster-growing non-conserving mode with n=%f.",2+wouldbe_adiab,2+fastest_growth);
        }

        free(pvecback);

        // If we get here, then initialise modes and evolve them!

        return _SUCCESS_;

}

/*
 * Test for tachyonic instability of x_smg in RD before initialisation of
 * perturbations: if not stable, cannot set ICs properly.
 */
int perturb_test_ini_extfld_ic_smg(struct precision * ppr,
                                   struct background * pba,
                                   struct perturbs * ppt){


        double kin, bra, run, ten, DelM2, Omx, wx;
        double l1,l2, l3, l4,l5,l6,l7,l8, cs2num, Dd;
        double B1_smg, B2_smg;
        double tau_ini, z_ref;
        double x_growth_smg;
        double * pvecback;
        int first_index_back;

        class_alloc(pvecback,pba->bg_size*sizeof(double),ppt->error_message);

        z_ref = ppr->pert_ic_ini_z_ref_smg;

        class_call(background_tau_of_z(pba, z_ref,&tau_ini),
                   pba->error_message,
                   ppt->error_message);

        class_call(background_at_tau(pba,
                                     tau_ini,
                                     pba->long_info,
                                     pba->inter_normal,
                                     &first_index_back,
                                     pvecback),
                   pba->error_message,
                   ppt->error_message);

        // look up alphas etc. at z_ref
        wx = pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg];
        Omx = pvecback[pba->index_bg_rho_smg]/pow(pvecback[pba->index_bg_H],2);
        kin = pvecback[pba->index_bg_kineticity_smg];
        bra = pvecback[pba->index_bg_braiding_smg];
        run = pvecback[pba->index_bg_mpl_running_smg];
        ten = pvecback[pba->index_bg_tensor_excess_smg];
        DelM2 = pvecback[pba->index_bg_delta_M2_smg];//M2-1
        l1 = pvecback[pba->index_bg_lambda_1_smg];
        l2 = pvecback[pba->index_bg_lambda_2_smg];
        l3 = pvecback[pba->index_bg_lambda_3_smg];
        l4 = pvecback[pba->index_bg_lambda_4_smg];
        l5 = pvecback[pba->index_bg_lambda_5_smg];
        l6 = pvecback[pba->index_bg_lambda_6_smg];
        l7 = pvecback[pba->index_bg_lambda_7_smg];
        l8 = pvecback[pba->index_bg_lambda_8_smg];
        cs2num = pvecback[pba->index_bg_cs2num_smg];
        Dd = pvecback[pba->index_bg_kinetic_D_smg];

        B1_smg = (bra/Dd)*(bra/(2.*(-2 + bra)*(kin + l1)))*((-6 + kin)*l1 + 3*l4);
        B1_smg +=  (3*pow(bra,3))*(l1/Dd)/(2.*(-2 + bra)*(kin + l1));
        B1_smg += 2*(cs2num/Dd)*(3*bra*kin + pow(kin,2) - 3*l4)/(2.*(-2. + bra)*(kin + l1));
        B1_smg += 2*(3*l2*l4/Dd + (kin/Dd)*(l1*l2 - 8*l7) - 8*l1/Dd*l7)/(2.*(-2 + bra)*(kin + l1));
        B1_smg -= 2*(bra/Dd)*((kin*l1/(kin + l1) - 3*l1*l2/(kin + l1) + 3*l4/(kin + l1))/(2.*(-2 + bra)));

        B2_smg =  8*(1 + DelM2)*(3*l2*l6/Dd + 4*kin*l8/Dd);
        B2_smg += 4*(l1/Dd)*(8*(1 + DelM2)*l8 + l2*(12 - 12*Omx + (1 + DelM2)*(-12 + kin + Omx*(3 - 9*wx))));
        B2_smg += 2*(bra/Dd)*bra*(6*(1 + DelM2)*l6 + l1*(12 - 12*Omx + (1 + DelM2)*(-30 + kin + 6*Omx*(1 - 3*wx))));
        B2_smg += 3*pow(bra,3)*(1 + DelM2)*(l1/Dd)*(6 + Omx*(-1 + 3*wx));
        B2_smg += 2*(cs2num/Dd)*(2*(1 + DelM2)*pow(kin,2) - 12*(1 + DelM2)*l6 + 3*kin*(8 - 8*Omx + (1 + DelM2)*(-8 + Omx*(2 - 6*wx) + bra*(6 + Omx*(-1 + 3*wx)))));
        B2_smg -= 2*(bra/Dd)*(12*(1 + DelM2)*l6 + l1*(24 - 24*Omx + (1 + DelM2)*(2*kin - 3*(8 + 2*Omx*(-1 + 3*wx) + l2*(6 + Omx*(-1 + 3*wx))))));
        B2_smg /= (4.*(-2 + bra)*(1 + DelM2)*(kin + l1));

        x_growth_smg = 0.5*(1.-B1_smg);

        if (1.-2.*B1_smg + B1_smg*B1_smg -4.*B2_smg >=0) {
                x_growth_smg += 0.5*sqrt(1. -2.*B1_smg + B1_smg*B1_smg -4.*B2_smg);
        }

        if (ppt->perturbations_verbose > 1) {
                printf("\nExternal field attractor ICs at z=%e. Standard solution for grav. field, h = (k tau)^2.\n",z_ref);
                if(x_growth_smg<=3) {
                        printf("  smg evolves on standard attractor in external field with x_smg = k^2 tau^3;\n\n");
                }
                else{
                        printf("  tachyonic instability in smg dominates, x_smg = k^2 tau^n with n=%f.\n",x_growth_smg);
                        printf("  smg is sensitive to its initial conditions at end of inflation.\n");
                }
        }

        class_test_except(ppr->pert_ic_tolerance_smg>0 && (x_growth_smg > 3.+ppr->pert_ic_tolerance_smg),
                          ppt->error_message,
                          free(pvecback),
                          "\n   Cannot set initial conditions for smg: tachyonic instability dominates superhorizon attractor.\n");

        free(pvecback);

        // If we get here, then initialise modes and evolve them!

        return _SUCCESS_;
}

int calc_extfld_ampl(int nexpo,  double kin, double bra, double dbra, double run, double ten, double DelM2,
                     double Omx, double wx, double l1, double l2, double l3, double l4,
                     double l5, double l6,double l7,double l8, double cs2num, double Dd,
                     double ic_regulator_smg, double * amplitude){

        /* Solutions assuming the alphas are small, i.e. x_smg does not gravitate but moves
         * on an attractor provided by collapsing radiation. (w!=1/3 terms included properly here!)
           // We have already tested for an RD tachyon at z=pert_ic_ini_z_ref_smg and it wasn't there.
         * We can thus assume that h has the standard solution (tau^2 for adiabatic)
         * and solve the x_smg e.o.m. assuming C1=C2=0.
         *
         *   x_smg = C1 tau^n1 + C2 tau^n2 + A k^2 tau^n
         *
         * This requires that if the tachyon has appeared at some later time, the system will be moving into it slowly.
         *
         * We do not correct any other fields, since it would be inconsistent to include them
         * here, but not in the calculation of the exponent. If this is importnant, use gravitating_attr ICs.
         *
         *
         * The on-attractor solution for the scalar velocity x_smg is x_smg = amplitude * k^2 tau^n * ppr->curvature_ini
         * with amplitude = -B3/(6 + 3*B1 + B2).
         */





        // Calculate the amplitude of x_smg in ext_field_attr ICs, both for adiabatic and isocurvature
        // Since the scalar does not backreact, the different Ad and ISOcurv solutions differ
        // only by the exponent in h, h = C*tau^n. The only n-dependent terms are in B3 and amplitude


        double B1_smg, B2_smg, B3_smg, B3num_smg, B3denom_smg;
        double den1, den2, den3, den4, reg_rescaled;

        reg_rescaled = ic_regulator_smg*(fabs(bra)+fabs(kin)+fabs(l1)); //rescale the regulator to be proportional to the alphas


        den1 = (2.*(-2 + bra)*(kin + l1));
        if(reg_rescaled>0 && (fabs(den1)<reg_rescaled)) {
                den1 = copysign(reg_rescaled,den1);
        }

        B1_smg = (bra/Dd)*(bra/den1)*((-6 + kin)*l1 + 3*l4);
        B1_smg +=  (3*pow(bra,3))*(l1/Dd)/den1;
        B1_smg += 2*(cs2num/Dd)*(3*bra*kin + pow(kin,2) - 3*l4)/(2.*(-2. + bra)*(kin + l1));
        B1_smg += 2*(3*l2*l4/Dd + (kin/Dd)*(l1*l2 - 8*l7) - 8*l1/Dd*l7)/den1;
        B1_smg -= 2*(bra/Dd)*((kin*l1/(kin + l1) - 3*l1*l2/(kin + l1) + 3*l4/(kin + l1))/(2.*(-2 + bra)));

        den2 = (4.*(-2 + bra)*(1 + DelM2)*(kin + l1));
        if(reg_rescaled>0 && (fabs(den2)<reg_rescaled)) {
                den2 = copysign(reg_rescaled,den2);
        }

        B2_smg =  8*(1 + DelM2)*(3*l2*l6/Dd + 4*kin*l8/Dd);
        B2_smg += 4*(l1/Dd)*(8*(1 + DelM2)*l8 + l2*(12 - 12*Omx + (1 + DelM2)*(-12 + kin + Omx*(3 - 9*wx))));
        B2_smg += 2*(bra/Dd)*bra*(6*(1 + DelM2)*l6 + l1*(12 - 12*Omx + (1 + DelM2)*(-30 + kin + 6*Omx*(1 - 3*wx))));
        B2_smg += 3*pow(bra,3)*(1 + DelM2)*(l1/Dd)*(6 + Omx*(-1 + 3*wx));
        B2_smg += 2*(cs2num/Dd)*(2*(1 + DelM2)*pow(kin,2) - 12*(1 + DelM2)*l6 + 3*kin*(8 - 8*Omx + (1 + DelM2)*(-8 + Omx*(2 - 6*wx) + bra*(6 + Omx*(-1 + 3*wx)))));
        B2_smg -= 2*(bra/Dd)*(12*(1 + DelM2)*l6 + l1*(24 - 24*Omx + (1 + DelM2)*(2*kin - 3*(8 + 2*Omx*(-1 + 3*wx) + l2*(6 + Omx*(-1 + 3*wx))))));
        B2_smg /= den2;


        den3 = ((2. * Omx)*(kin + l1));
        reg_rescaled *=Omx;
        if(reg_rescaled>0 && (fabs(den3)<reg_rescaled)) {
                den3 = copysign(reg_rescaled,den3);
        }
        B3num_smg = ((-(((-2. + bra) * bra + 2 * l2) *
                        ((-2. + bra) * l1 - 4 * l3 + 2 * Dd * (-1. + nexpo))) +
                      cs2num * (-2 * (-2. + bra) * kin - 8 * l3 + 4 * Dd * (-1. + nexpo))) *
                     nexpo) / den3;

        B3denom_smg = 4*(Dd/Omx)*(-2 + bra);

        B3_smg = B3num_smg/B3denom_smg;

        reg_rescaled = ic_regulator_smg*(fabs(B1_smg)+fabs(B2_smg));
        den4 = B1_smg + B2_smg + nexpo + B1_smg*nexpo + pow(nexpo,2);

        if(reg_rescaled>0 && (fabs(den4)<reg_rescaled)) {
                den4 = copysign(reg_rescaled,den4);
        }


        *amplitude = -B3_smg/den4;

        return _SUCCESS_;
}

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
        double * res_p, double * cD_p, double * cB_p, double * cH_p,
        double * c9_p, double * c10_p, double * c12_p, double * c13_p
        ){
        /* It returns the alphas and the coefficients of the Einstein equations
           that will be used to evaluate the perturbations and their initial
           conditions. This function uses use_pert_var_deltaphi_smg to decide which
           coefficients to output.
         */

        double a = pvecback[pba->index_bg_a];
        double H = pvecback[pba->index_bg_H];
        double Hp = pvecback[pba->index_bg_H_prime];

        *delM2 = pvecback[pba->index_bg_delta_M2_smg];
        *M2 = pvecback[pba->index_bg_M2_smg];
        *kin = pvecback[pba->index_bg_kineticity_smg];
        *bra = pvecback[pba->index_bg_braiding_smg];
        *ten = pvecback[pba->index_bg_tensor_excess_smg];
        *run = pvecback[pba->index_bg_mpl_running_smg];
        *beh = pvecback[pba->index_bg_beyond_horndeski_smg];

        if (ppt->use_pert_var_deltaphi_smg == _TRUE_) {
                *res = 1.;
                *cD  = pvecback[pba->index_bg_kinetic_D_over_phiphi_smg];
                *cK  = pvecback[pba->index_bg_kineticity_over_phiphi_smg];
                *cB  = pvecback[pba->index_bg_braiding_over_phi_smg];
                *cH  = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg];
                *c0  = pvecback[pba->index_bg_C0_smg];
                *c1  = pvecback[pba->index_bg_C1_smg];
                *c2  = pvecback[pba->index_bg_C2_smg];
                *c3  = pvecback[pba->index_bg_C3_smg];
                *c4  = pvecback[pba->index_bg_C4_smg];
                *c5  = pvecback[pba->index_bg_C5_smg];
                *c6  = pvecback[pba->index_bg_C6_smg];
                *c7  = pvecback[pba->index_bg_C7_smg];
                *c8  = pvecback[pba->index_bg_C8_smg];
                *c9  = pvecback[pba->index_bg_C9_smg];
                *c10 = pvecback[pba->index_bg_C10_smg];
                *c11 = pvecback[pba->index_bg_C11_smg];
                *c12 = pvecback[pba->index_bg_C12_smg];
                *c13 = pvecback[pba->index_bg_C13_smg];
                *c14 = pvecback[pba->index_bg_C14_smg];
                *c15 = pvecback[pba->index_bg_C15_smg];
                *c16 = pvecback[pba->index_bg_C16_smg];
                *res_p = 0.;
                *cD_p  = pvecback[pba->index_bg_kinetic_D_over_phiphi_prime_smg];
                *cB_p  = pvecback[pba->index_bg_braiding_over_phi_prime_smg];
                *cH_p  = pvecback[pba->index_bg_beyond_horndeski_over_phi_prime_smg];
                *c9_p  = pvecback[pba->index_bg_C9_prime_smg];
                *c10_p = pvecback[pba->index_bg_C10_prime_smg];
                *c12_p = pvecback[pba->index_bg_C12_prime_smg];
                *c13_p = pvecback[pba->index_bg_C13_prime_smg];
        }
        else if (ppt->use_pert_var_deltaphi_smg == _FALSE_) {
                *res = -a*H;
                *cD  = pvecback[pba->index_bg_kinetic_D_smg];
                *cK  = pvecback[pba->index_bg_kineticity_smg];
                *cB  = pvecback[pba->index_bg_braiding_smg];
                *cH  = pvecback[pba->index_bg_beyond_horndeski_smg];
                *c0  = pvecback[pba->index_bg_A0_smg];
                *c1  = pvecback[pba->index_bg_A1_smg];
                *c2  = pvecback[pba->index_bg_A2_smg];
                *c3  = pvecback[pba->index_bg_A3_smg];
                *c4  = pvecback[pba->index_bg_A4_smg];
                *c5  = pvecback[pba->index_bg_A5_smg];
                *c6  = pvecback[pba->index_bg_A6_smg];
                *c7  = pvecback[pba->index_bg_A7_smg];
                *c8  = pvecback[pba->index_bg_A8_smg];
                *c9  = pvecback[pba->index_bg_A9_smg];
                *c10 = pvecback[pba->index_bg_A10_smg];
                *c11 = pvecback[pba->index_bg_A11_smg];
                *c12 = pvecback[pba->index_bg_A12_smg];
                *c13 = pvecback[pba->index_bg_A13_smg];
                *c14 = pvecback[pba->index_bg_A14_smg];
                *c15 = pvecback[pba->index_bg_A15_smg];
                *c16 = pvecback[pba->index_bg_A16_smg];
                *res_p = -a*(Hp + a*H);
                *cD_p  = pvecback[pba->index_bg_kinetic_D_prime_smg];
                *cB_p  = pvecback[pba->index_bg_braiding_prime_smg];
                *cH_p  = pvecback[pba->index_bg_beyond_horndeski_prime_smg];
                *c9_p  = pvecback[pba->index_bg_A9_prime_smg];
                *c10_p = pvecback[pba->index_bg_A10_prime_smg];
                *c12_p = pvecback[pba->index_bg_A12_prime_smg];
                *c13_p = pvecback[pba->index_bg_A13_prime_smg];
        }
        else {
                printf("It was not possible to determine if oscillations of the background scalar field should be allowed or not.\n");
                return _FAILURE_;
        }

        return _SUCCESS_;
}

int get_x_x_prime_qs_smg(
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct perturb_workspace * ppw,
        double k, double * x_qs_smg, double * x_prime_qs_smg
        ){

        double k2 = k*k;
        double rho_r, p_tot, p_smg;
        double a, H, delM2, M2, kin, bra, ten, run, beh;
        double res, cD, cK, cB, cH;
        double c0, c1, c2, c3, c4, c5, c6, c7, c8;
        double c9, c10, c11, c12, c13, c14, c15, c16;
        double c9_p, c10_p, c12_p, c13_p;
        double res_p, cD_p, cB_p, cH_p;
        double x_prime_qs_smg_num, x_prime_qs_smg_den;

        a = ppw->pvecback[pba->index_bg_a];
        H = ppw->pvecback[pba->index_bg_H];
        rho_r = ppw->pvecback[pba->index_bg_rho_g] + ppw->pvecback[pba->index_bg_rho_ur];
        p_tot = ppw->pvecback[pba->index_bg_p_tot_wo_smg];
        p_smg = ppw->pvecback[pba->index_bg_p_smg];

        class_call(
                get_gravity_coefficients_smg(
                        ppt, pba, ppw->pvecback,
                        &delM2, &M2, &kin, &bra, &ten, &run, &beh, &res,
                        &cD, &cK, &cB, &cH, &c0, &c1, &c2, &c3,
                        &c4, &c5, &c6, &c7, &c8, &c9, &c10, &c11,
                        &c12, &c13, &c14, &c15, &c16, &res_p, &cD_p, &cB_p,
                        &cH_p, &c9_p, &c10_p, &c12_p, &c13_p
                        ),
                ppt->error_message,
                ppt->error_message);

        /* This is the expression for the scalar field in the quasi static approximation */
        if (ppt->get_h_from_trace == _TRUE_) {
                /* Scalar field in QS with h' */
                *x_qs_smg =
                        +1./res*(
                                -9./2.*cB*pow(a,2)*ppw->delta_p/M2
                                + c10*k2*ppw->pvecmetric[ppw->index_mt_eta]
                                + (c9*pow(a*H,2) - 1./3.*cH*k2)*ppw->pvecmetric[ppw->index_mt_h_prime]/a/H
                                + 2./3.*cH*pow(k2,2)*ppw->pvecmetric[ppw->index_mt_alpha]/a/H
                                )/(c13*k2 + c12*pow(a*H,2));
        }
        else {
                /* Scalar field in QS without h' */
                *x_qs_smg =
                        1./res*(
                                +k2*(
                                        +4.*(1. + beh)*cH*k2
                                        - 3.*pow(a*H,2)*((2. - bra)*c10 + 4.*(1. + beh)*c9)
                                        )*ppw->pvecmetric[ppw->index_mt_eta]
                                - 2.*(2. - bra)*cH*H*pow(k2,2)*a*ppw->pvecmetric[ppw->index_mt_alpha]
                                + 6.*(
                                        +(cH*k2 - 3.*c9*pow(a*H,2))*ppw->delta_rho
                                        + 9./4.*cB*(2. - bra)*pow(a*H,2)*ppw->delta_p
                                        )*pow(a,2)/M2
                                )/(
                                +4.*c15*cH*pow(k2,2)
                                - k2*pow(a*H,2)*(3.*c13*(2. - bra) + 12.*c15*c9 - 4.*c16*cH)
                                - 3.*pow(a*H,4)*(c12*(2. - bra) + 4.*c16*c9)
                                );
        }


        /* scalar field derivative equation
         * In order to estimate it we followed this procedure:
         * - we calculated analytically the time derivative of the x_smg equation
         * - we used delta_p' = delta_rho_r'/3 (radiation is the only component that contributes to delta_p')
         * - we used the conservation equation for radiation to get rid of delta_rho_r'
         * - we used the Einstein equations to get rid of eta', h'', alpha'
         * The result is approximated when rsa is on since the velocity of radiation gets updated only after
         * this call in perturb_einstein */

        if (ppt->get_h_from_trace == _TRUE_) {
                /* Numerator of the scalar field derivative in QS with h' */
                x_prime_qs_smg_num =
                        +3.*(
                                +3.*(2.*c9*cK - cB*cD*(2. + run) - cD*(cB*res_p/res - cB_p)/a/H)
                                - 2.*cH*cK*pow(a*H,-2)*k2
                                )*ppw->delta_rho_r*a/H/M2
                        + 9.*(
                                +2.*cD*(cH*res_p/res - cH_p)/a/H + 6.*c3*c9
                                - cD*(cB + c10) + 3.*cD*cH*(1. + 2./3.*run - (p_tot + p_smg)*pow(H,-2))
                                - 2.*c3*cH*k2*pow(a*H,-2)
                                )*pow(H,-2)*ppw->rho_plus_p_theta_r/M2
                        + 18.*cD*cH*pow(a*H,-2)*k2*ppw->rho_plus_p_shear*a/H/M2
                        + 4.*k2*pow(a*H,-2)*(
                                +cH*(c1 - cD - ten*cD)*k2*pow(a*H,-2)
                                - 3./2.*(2.*c1*c9 - (c10*cD*res_p/res - cD*c10_p)/a/H)
                                )*a*H*ppw->pvecmetric[ppw->index_mt_eta]
                        + 3.*(
                                +2.*cD*(c9*res_p/res - c9_p)/a/H - 2.*c2*c9 + c9*cD
                                - cD*(2.*cB*rho_r/M2 - 3.*c9*(p_tot + p_smg))*pow(H,-2)
                                + 2./3.*cH*(c2 + 2.*cD + run*cD)*k2*pow(a*H,-2)
                                )*ppw->pvecmetric[ppw->index_mt_h_prime]
                        + 6.*a*H*res*(
                                +c6*c9 + cD*(c12_p/a/H - c12 - 3.*c12*(p_tot + p_smg)*pow(H,-2))
                                - (
                                        +cD*(2.*c0*cH*res_p/res - c13_p - 2.*c0*cH_p)/a/H
                                        + c9*(6.*c0*c3 - c7) - c0*c10*cD + c6*cH/3.
                                        + 3.*c0*cD*cH*(1. + 2./3.*run - (p_tot + p_smg)*pow(H,-2))
                                        )*k2*pow(a*H,-2)
                                + 1./3.*cH*(6.*c0*c3 - c7 + 2.*c8*cD)*pow(k2,2)*pow(a*H,-4)
                                )*ppw->pvecmetric[ppw->index_mt_x_smg];

                /* Denominator of the scalar field derivative in QS with h' */
                x_prime_qs_smg_den =
                        -6.*res*(
                                +c4*c9 + c12*cD
                                - k2*(
                                        +6.*cB*cD*(cH*res_p/res - cH_p)/a/H
                                        - 12.*c9*(c5 - 3./2.*c3*cB)
                                        - 3.*cD*(c10*cB + 2.*c13)
                                        + 2.*c4*cH
                                        + 3.*cB*cD*cH*(3. + 2.*run)
                                        - 9.*cB*cD*cH*(p_tot + p_smg)*pow(H,-2)
                                        )/6.*pow(a*H,-2)
                                - cH*pow(k2,2)*(2.*c5 - 3.*c3*cB + 2.*cD*cH)/3.*pow(a*H,-4)
                                );
        }
        else {
                /* Numerator of the scalar field derivative in QS without h' */
                x_prime_qs_smg_num =
                        -18.*(2. - bra)*cD*cH*pow(H,-3)*k2*ppw->rho_plus_p_shear/a/M2
                        + (2. - bra)*(
                                +6.*cH*cK*pow(H,-3)*k2/a
                                + 9.*(
                                        +cD*(cB*res_p - cB_p*res)
                                        + ((2. + run)*cB*cD - 2.*c9*cK)*H*res*a
                                        )*pow(H,-2)/res
                                )*ppw->delta_rho_r/M2
                        + 9.*(2. - bra)*(
                                +2.*c3*cH*pow(H,-4)*pow(a,-2)*k2
                                - (
                                        +2.*cD*H*(cH*res_p - cH_p*res)/a/res
                                        + (6.*c3*c9 - c10*cD - cB*cD + 3.*cD*cH + 2.*run*cD*cH)*pow(H,2)
                                        - 3.*cD*cH*(p_tot + p_smg)
                                        )*pow(H,-4)
                                )/M2*ppw->rho_plus_p_theta_r
                        - (
                                +12.*(c2 + 2.*cD + run*cD)*cH*pow(H,-3)*k2/a
                                + 18.*(
                                        +2.*cD*H*(c9*res_p/res - c9_p)
                                        - 2.*cB*cD*rho_r*a/M2
                                        + c9*(cD - 2.*c2)*pow(H,2)*a
                                        + 3.*c9*cD*(p_tot + p_smg)*a
                                        )*pow(H,-3)
                                )/M2*ppw->delta_rho
                        + (
                                +4.*cH*(
                                        +(2. - bra)*(cD + ten*cD - c1)
                                        - 2.*(c2 + 2.*cD + run*cD)*(1. + beh)
                                        )*pow(k2,2)*pow(a*H,-3)
                                - 6.*(
                                        +(2. - bra)*(cD*(c10*res_p/res - c10_p)/a/H - 2.*c1*c9)
                                        + 2.*(1. + beh)*(
                                                +2.*cD*(c9*res_p/res - c9_p)/a/H
                                                + c9*(cD - 2.*c2)
                                                - 2.*cB*cD*rho_r*pow(H,-2)/M2
                                                + 3.*c9*cD*(p_tot + p_smg)*pow(H,-2)
                                                )
                                        )*k2/a/H
                                )*ppw->pvecmetric[ppw->index_mt_eta]
                        + (
                                +(
                                        -(2. - bra)*(6.*c0*c3 - c7 + 2.*c8*cD)
                                        + 4.*c15*(c2 + 2.*cD + run*cD)
                                        )*2.*cH*pow(k2,2)*res*pow(a*H,-3)
                                + 6.*(
                                        +cD*H*(
                                                +4.*c16*c9*res_p/res
                                                - c12_p*(2. - bra)
                                                - 4.*c16*c9_p
                                                )/a
                                        - 4.*c16*cB*cD*rho_r/M2
                                        - 4.*c16*c2*c9*pow(H,2)
                                        - c6*c9*(2. - bra)*pow(H,2)
                                        + cD*((2. - bra)*c12 + 2.*c16*c9)*(pow(H,2) + 3.*(p_tot + p_smg))
                                        )*res*a/H
                                + 2.*(
                                        +12.*cD*c15*(c9*res_p/res - c9_p)/H/a
                                        - 12.*c15*cB*cD*rho_r/M2*pow(H,-2)
                                        + 2.*(
                                                +3.*c15*c9*(cD - 2.*c2)
                                                + 2.*cH*c16*(c2 + 2.*cD + run*cD)
                                                )
                                        + (2. - bra)*(
                                                -3.*cD*(2.*c0*cH_p - 2.*c0*cH*res_p/res + c13_p)/a/H
                                                + 18.*c0*c3*c9 - 3.*c7*c9
                                                - 3.*c0*c10*cD + c6*cH
                                                + 6.*(3./2. + run)*c0*cD*cH
                                                )
                                        - 9.*cD*((2. - bra)*c0*cH - 2.*c15*c9)*(p_tot + p_smg)*pow(H,-2)
                                        )*res*k2/a/H
                                )*ppw->pvecmetric[ppw->index_mt_x_smg];


                /* Denominator of the scalar field derivative in QS without h' */
                x_prime_qs_smg_den =
                        -2.*cH*res*(2. - bra)*(2.*c5 - 3.*c3*cB + 2.*cD*cH)*pow(a*H,-4)*pow(k2,2)
                        + a*(
                                -6.*cB*cD*(2. - bra)*(cH*res_p - cH_p*res)*H
                                + (2. - bra)*(
                                        +12.*c5*c9 - 18.*c3*c9*cB
                                        + 6.*c13*cD + 3.*c10*cB*cD - 2.*c4*cH
                                        - 3.*(3. + 2.*run)*cB*cD*cH
                                        + 9.*cB*cD*cH*(p_tot + p_smg)*pow(H,-2)
                                        )*res*a*pow(H,2)
                                - 8.*(c2 + 2.*cD + run*cD)*c14*cH*pow(H,2)*res*a
                                )*k2*pow(a*H,-4)
                        + 6.*(
                                -4.*c14*cD*H*(c9*res_p - c9_p*res)
                                + 4.*c14*(cB*cD*rho_r/M2 + c2*c9*pow(H,2))*res*a
                                + (2. - bra)*(c4*c9 + c12*cD)*pow(H,2)*res*a
                                - 2.*c14*c9*cD*(pow(H,2) + 3.*(p_tot + p_smg))*res*a
                                )*pow(H,-2)/a;
        }

        *x_prime_qs_smg = x_prime_qs_smg_num/x_prime_qs_smg_den;

        return _SUCCESS_;
}
