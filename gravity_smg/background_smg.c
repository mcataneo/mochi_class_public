/** @file background_smg.c Documented background_smg module
 *
 * Emilio Bellini, Ignacy Sawicki, Miguel Zumalacarregui, TODO_EB: date here xx.xx.xxxx
 *
 * Additional functions for the background module.
 * It contains all the hi_class related functions (_smg)
 * that are used by background.c. In this way the main hi_class
 * modifications are stored here and the standard Class modules
 * remain cleaner.
 *
 * The following nomenclature has been adopted:
 *
 * -# all the functions end with "_smg" to make them easily
 *    recognizable
 * -# all the functions starting with "background" are
 *    directly called by background.c or the classy wrapper
 * -# all the functions that do not start with "background"
 *    are only used internally in background_smg.c
 */

#include "background_smg.h"


/**
* Fill the modified gravity part of background_functions.
* First all the Horndeski functions G_i(X,phi) are computed.
* A loop is used to allow different implementations without erasing previous ones
* Note that in CLASS units a canonical field has G2 = X ~ [Mpc]^-2
* This implies that phi_prime ~ [Mpc]^-1
* and phi is dimensionless (we can think of phi as given in units of the Planck mass
* - A default module to numerically compute the derivatives when no analytic functions are given should be added.
* Numerical derivatives may further serve as a consistency check.
* TODO: Add a background_write_alpha_primes
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pvecback      Output: vector of background quantities (assumed to be already allocated)
* @param ptr_rho_tot   Output: pointer to the total (without smg) energy density
* @param ptr_p_tot     Output: pointer to the total (without smg) pressure
* @param ptr_rho_de    Output: pointer to the dark energy density (smg)
* @return the error status
*/
int background_gravity_functions_smg(
				 														 struct background *pba,
				 													 	 double a,
				 													   double * pvecback_B,
				 													   double * pvecback,
         									 					 double * ptr_rho_tot,
				 													   double * ptr_p_tot,
				 													   double * ptr_rho_de
				 												 		 ){


	/** - save rho_tot and p_tot without smg (it is important that smg is the last thing to be computed) */
	pvecback[pba->index_bg_rho_tot_wo_smg] = *ptr_rho_tot;
	pvecback[pba->index_bg_p_tot_wo_smg] = *ptr_p_tot;

	// scalar field + curvature not yet implemented
	// TODO_EB: look for a better place to put this test
	class_test(pba->K !=0 ,
	     pba->error_message,
	     "has_smg with curvature K = %e not yet implemented",pba->K);

  if (pba->field_evolution_smg == _TRUE_) {

		/* Define local variables */
		double x,f,df;
    int n, n_max=100;
		double H;

		/* Get phi and phi_prime from the integrator */
		pvecback[pba->index_bg_phi_smg] = pvecback_B[pba->index_bi_phi_smg];
    pvecback[pba->index_bg_phi_prime_smg] = pvecback_B[pba->index_bi_phi_prime_smg];

		/* declare G functions and set defaults */
		struct G_functions_and_derivs gf = DEFAULT_G_FUNCTIONS_AND_DERIVS;

		/* update G functions and derivatives */
 	  class_call(gravity_models_get_Gs_smg(pba, a, pvecback_B, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* get Es functions. Notation:
		* E0 + E1 H + E3 H^3 = E2 H^2
		*/
		class_call(gravity_functions_Es_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
			pba->error_message,
			pba->error_message
		);

		double E0 = pvecback[pba->index_bg_E0_smg];
		double E1 = pvecback[pba->index_bg_E1_smg];
		double E2 = pvecback[pba->index_bg_E2_smg];
		double E3 = pvecback[pba->index_bg_E3_smg];

		class_test(E3*pow(E0,1./2.) > 1e-10 && pba->initial_conditions_set_smg == _FALSE_,
							 pba->error_message,
							 "E3=%e is large in Friedmann constraint when setting ICs ",  E3);

		/* get Hubble, either solving the cubic Friedmann
		* equation or from the integrator.
		*/
		if ((pba->hubble_evolution == _FALSE_) || (pba->initial_conditions_set_smg == _FALSE_)) {
      /* Use Newton's method to solve the cubic Friedmann equation */
      x = sqrt(E0);
      f = E3*x*x*x -E2*x*x + E1*x + E0;
      n = 0;
      while (fabs(f/E0)> 1e-8 && n < 100){
				f = E3*x*x*x - E2*x*x + E1*x + E0;
				df = 3.*E3*x*x - 2.*E2*x + E1;
				x -= f/df;
				n++;
      }
      H=x;
      if ((pba->background_verbose > 5) && (pba->initial_conditions_set_smg == _FALSE_ ))
				printf(" Initial H = %e, sqrt(rho) = %e, ratio = %e, n=%i \n", H, sqrt(pvecback[pba->index_bg_rho_tot_wo_smg]),sqrt(pvecback[pba->index_bg_rho_tot_wo_smg])/H,n);

 	    class_test(isnan(H),
 		       			 pba->error_message,
 	               "H=%e is not a number at a = %e. phi = %e, phi_prime = %e, E0=%e, E1=%e, E2=%e, E3=%e ",
 		       			 H,a,pvecback[pba->index_bg_phi_smg],pvecback[pba->index_bg_phi_prime_smg],E0,E1,E2,E3);

 		}
		else {
			H = exp(pvecback_B[pba->index_bi_logH]);
		}
		pvecback[pba->index_bg_H] = H;

		/* get Ps and Rs functions. Notation:
		* P0 phi'' + P1 H' + P2 = 0
		* R0 phi'' + R1 H' + R2 = 0
		*/
 	  class_call(gravity_functions_Ps_and_Rs_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		double P0 = pvecback[pba->index_bg_P0_smg];
    double P1 = pvecback[pba->index_bg_P1_smg];
    double P2 = pvecback[pba->index_bg_P2_smg];
    double R0 = pvecback[pba->index_bg_R0_smg];
    double R1 = pvecback[pba->index_bg_R1_smg];
    double R2 = pvecback[pba->index_bg_R2_smg];

    class_test_except((P1*R0 - P0*R1) == 0 ,
	       pba->error_message,
         free(pvecback);free(pvecback_B),
               "scalar field mixing with metric has degenerate denominator at a = %e, phi = %e, phi_prime = %e \n with P1 = %e, R0 =%e, P0=%e, R1=%e, \n H=%e, E0=%e, E1=%e, E2=%e, E3=%e",
	       a,pvecback[pba->index_bg_phi_smg],pvecback[pba->index_bg_phi_prime_smg], P1, R0, P0, R1,
	       pvecback[pba->index_bg_H],E0,E1,E2,E3);

    /* Friedmann space-space equation with friction added */
    pvecback[pba->index_bg_H_prime] = (R0*P2 - P0*R2)/(P0*R1 - P1*R0);
    /* choose sign for friction depending on the derivative */
    if ((2.*E2*H - E1 - 3*E3*H*H)>=0)
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(E2*H*H - E0 - E1*H - E3*H*H*H);
    else{
      pvecback[pba->index_bg_H_prime] += a*pba->hubble_friction*(E2*H*H - E0 - E1*H - E3*H*H*H);
    }

    /* Field equation */
    pvecback[pba->index_bg_phi_prime_prime_smg] = (P1*R2 - R1*P2)/(P0*R1 - P1*R0);

		/* get alphas, background density and pressure, shift and current. */
		class_call(gravity_functions_building_blocks_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* get Bs, intermediate step for perturbations. */
		class_call(gravity_functions_Bs_from_Gs_smg(pba, a, pvecback_B, pvecback, &gf),
 	    pba->error_message,
 	    pba->error_message
 	  );

	}
	// end of if pba->field_evolution_smg
  else{

		double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
	  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

		/* get background parametrizations. */
		class_call(gravity_models_get_back_par_smg(pba, a, pvecback, pvecback_B),
 	    pba->error_message,
 	    pba->error_message
 	  );

		/* add _smg to rho_tot */
		rho_tot += pvecback[pba->index_bg_rho_smg];
    p_tot += pvecback[pba->index_bg_p_smg];

    pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

    //add friction term
    if (pba->hubble_evolution == _TRUE_ && pba->initial_conditions_set_smg == _TRUE_){
      pvecback[pba->index_bg_H] = exp(pvecback_B[pba->index_bi_logH]);
      /** - compute derivative of H with respect to conformal time */
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H] - rho_tot - pba->K/a/a);
    }

    // Compute time derivative of rho_smg
    if (pba->rho_evolution_smg == _TRUE_){
        pvecback[pba->index_bg_rho_prime_smg] = -3.*a*pvecback[pba->index_bg_H]*(1.+pvecback[pba->index_bg_w_smg])*pvecback[pba->index_bg_rho_smg];
    }


		/* initialize the values to the defaults */
    pvecback[pba->index_bg_kineticity_smg] = 0;
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
    pvecback[pba->index_bg_delta_M2_smg] = 0.;
    pvecback[pba->index_bg_mpl_running_smg] = 0.;

	if(pba->gravity_model_smg != stable_params){
		/* get background parametrizations. */
		class_call(gravity_models_get_alphas_par_smg(pba, a, pvecback, pvecback_B),
 	    pba->error_message,
 	    pba->error_message
 	  );
  	}

	}
	//end of parameterized mode

  // add a value to the kineticity to avoid problems with perturbations in certain models.
  // NOTE: this needs to be done here to avoid interfering with the equations
  pvecback[pba->index_bg_kineticity_smg] += pba->kineticity_safe_smg;

  //Derivatives of the BS functions and others. Set to zero here and computed numerically once the background is integrated (needed so that debuggers don't complain).

  pvecback[pba->index_bg_kineticity_prime_smg] = 0.;
  pvecback[pba->index_bg_braiding_prime_smg] = 0.;
  pvecback[pba->index_bg_mpl_running_prime_smg] = 0.;
  pvecback[pba->index_bg_tensor_excess_prime_smg] = 0.;
  pvecback[pba->index_bg_beyond_horndeski_prime_smg] = 0.;
  pvecback[pba->index_bg_H_prime_prime] = 0.;
  pvecback[pba->index_bg_p_tot_wo_prime_smg] = 0.;
  pvecback[pba->index_bg_p_tot_wo_prime_prime_smg] = 0.;
  pvecback[pba->index_bg_p_prime_smg] = 0.;
  pvecback[pba->index_bg_p_prime_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_prime_smg] = 0.;
	if (pba->field_evolution_smg == _TRUE_) {
	  pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] = 0.;
	  pvecback[pba->index_bg_kinetic_D_over_phiphi_prime_smg] = 0.;
	}
  pvecback[pba->index_bg_cs2num_smg] = 0.;
  pvecback[pba->index_bg_cs2num_prime_smg] = 0.;
  pvecback[pba->index_bg_mu_p_smg] = 0.;
  pvecback[pba->index_bg_mu_inf_smg] = 0.;
  pvecback[pba->index_bg_muZ_inf_smg] = 0.;
  pvecback[pba->index_bg_mu_p_prime_smg] = 0.;
  pvecback[pba->index_bg_mu_inf_prime_smg] = 0.;
  pvecback[pba->index_bg_muZ_inf_prime_smg] = 0.;
  pvecback[pba->index_bg_A0_smg] = 0.;
  pvecback[pba->index_bg_A1_smg] = 0.;
  pvecback[pba->index_bg_A2_smg] = 0.;
  pvecback[pba->index_bg_A3_smg] = 0.;
  pvecback[pba->index_bg_A4_smg] = 0.;
  pvecback[pba->index_bg_A5_smg] = 0.;
  pvecback[pba->index_bg_A6_smg] = 0.;
  pvecback[pba->index_bg_A7_smg] = 0.;
  pvecback[pba->index_bg_A8_smg] = 0.;
  pvecback[pba->index_bg_A9_smg] = 0.;
  pvecback[pba->index_bg_A10_smg] = 0.;
  pvecback[pba->index_bg_A11_smg] = 0.;
  pvecback[pba->index_bg_A12_smg] = 0.;
  pvecback[pba->index_bg_A13_smg] = 0.;
  pvecback[pba->index_bg_A14_smg] = 0.;
  pvecback[pba->index_bg_A15_smg] = 0.;
  pvecback[pba->index_bg_A16_smg] = 0.;
  pvecback[pba->index_bg_A9_prime_smg] = 0.;
  pvecback[pba->index_bg_A10_prime_smg] = 0.;
  pvecback[pba->index_bg_A12_prime_smg] = 0.;
  pvecback[pba->index_bg_A13_prime_smg] = 0.;
	if (pba->field_evolution_smg == _TRUE_) {
		pvecback[pba->index_bg_C0_smg] = 0.;
	  pvecback[pba->index_bg_C1_smg] = 0.;
	  pvecback[pba->index_bg_C2_smg] = 0.;
	  pvecback[pba->index_bg_C3_smg] = 0.;
	  pvecback[pba->index_bg_C4_smg] = 0.;
	  pvecback[pba->index_bg_C5_smg] = 0.;
	  pvecback[pba->index_bg_C6_smg] = 0.;
	  pvecback[pba->index_bg_C7_smg] = 0.;
	  pvecback[pba->index_bg_C8_smg] = 0.;
	  pvecback[pba->index_bg_C9_smg] = 0.;
	  pvecback[pba->index_bg_C10_smg] = 0.;
	  pvecback[pba->index_bg_C11_smg] = 0.;
	  pvecback[pba->index_bg_C12_smg] = 0.;
	  pvecback[pba->index_bg_C13_smg] = 0.;
	  pvecback[pba->index_bg_C14_smg] = 0.;
	  pvecback[pba->index_bg_C15_smg] = 0.;
	  pvecback[pba->index_bg_C16_smg] = 0.;
	  pvecback[pba->index_bg_C9_prime_smg] = 0.;
	  pvecback[pba->index_bg_C10_prime_smg] = 0.;
	  pvecback[pba->index_bg_C12_prime_smg] = 0.;
	  pvecback[pba->index_bg_C13_prime_smg] = 0.;
	}
  pvecback[pba->index_bg_lambda_1_smg] = 0.;
  pvecback[pba->index_bg_lambda_2_smg] = 0.;
  pvecback[pba->index_bg_lambda_3_smg] = 0.;
  pvecback[pba->index_bg_lambda_4_smg] = 0.;
  pvecback[pba->index_bg_lambda_5_smg] = 0.;
  pvecback[pba->index_bg_lambda_6_smg] = 0.;
  pvecback[pba->index_bg_lambda_7_smg] = 0.;
  pvecback[pba->index_bg_lambda_8_smg] = 0.;
  pvecback[pba->index_bg_lambda_9_smg] = 0.;
  pvecback[pba->index_bg_lambda_10_smg] = 0.;
  pvecback[pba->index_bg_lambda_11_smg] = 0.;
  pvecback[pba->index_bg_lambda_2_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_8_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_9_prime_smg] = 0.;
  pvecback[pba->index_bg_lambda_11_prime_smg] = 0.;
  pvecback[pba->index_bg_G_eff_smg] = 0.;
  pvecback[pba->index_bg_slip_eff_smg] = 0.;
  if (pba->field_evolution_smg == _FALSE_ && pba->M_pl_evolution_smg == _FALSE_){
    pvecback[pba->index_bg_mpl_running_smg] = 0.;
  }

    /* Check required conditions for the gravity_models. */
  if ( (pba->skip_stability_tests_smg == _FALSE_) && (pba->parameters_tuned_smg == _TRUE_) && (pba->Omega_smg_debug == 0) ) {
		if (pba->is_quintessence_smg == _TRUE_){
		    /*  Check that w is not lower than w < -1 for quintessence */
		   class_test( (pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg] < -(1 + pba->quintessence_w_safe_smg)),
		       pba->error_message,
		       "Dark Energy equation of state at a = %e is w = %e, lower than w = -(1 + th), with threshold, th = %e. Quintessence models can only have  w > -1. Aborting.\n", a, pvecback[pba->index_bg_p_smg]/pvecback[pba->index_bg_rho_smg], pba->quintessence_w_safe_smg);

		   class_test(( (pvecback[pba->index_bg_rho_smg] / (pba->Omega0_smg*pow(pba->H0,2))  + 1e-3)  < 1),
		           pba->error_message,
		           "Dark Energy density, rho_smg = %e, at a = %e is lower than its given present value, Omega0_smg*H0^2= %e, since quintessence models have  w > -1, it will never be reached. Aborting.\n", pvecback[pba->index_bg_rho_smg], a, pba->Omega0_smg*pow(pba->H0,2) );
		}
  }

	*ptr_rho_tot += pvecback[pba->index_bg_rho_smg];
	*ptr_p_tot += pvecback[pba->index_bg_p_smg];
	//divide relativistic & nonrelativistic (not very meaningful for oscillatory models)

	//TODO: need to define menaingfully -> separate early universe (IC, BBN...) from late (Halofit...)
	//BUG: causes problem with halofit!, if not, causes bug with Brans-Dicke
	*ptr_rho_de += pvecback[pba->index_bg_rho_smg];

	/** - compute w_smg */
	if (pba->rho_evolution_smg == _FALSE_ && pba->expansion_model_smg != wext && pba->expansion_model_smg != rho_de) {
		pvecback[pba->index_bg_w_smg] = pvecback[pba->index_bg_p_smg] / pvecback[pba->index_bg_rho_smg];
	}

  return _SUCCESS_;

}


/**
 * Free all memory space allocated by hi_class in background.
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */
int background_free_smg(
			  								struct background *pba
											  ) {

	if (pba->parameters_smg != NULL)
	  free(pba->parameters_smg);
	//dealocate parameters_2_smg only for parameterizations
	if (pba->field_evolution_smg == _FALSE_ && pba->parameters_2_smg != NULL)
	  free(pba->parameters_2_smg);
	if(pba->gravity_model_smg == stable_params && pba->stable_params_lna_smg != NULL)
		free(pba->stable_params_lna_smg);
	if(pba->gravity_model_smg == stable_params && pba->stable_params_smg != NULL)
		free(pba->stable_params_smg);
	if(pba->gravity_model_smg == stable_params && pba->ddstable_params_smg != NULL)
		free(pba->ddstable_params_smg);
	if(pba->gravity_model_smg == stable_params && pba->stable_params_derived_smg != NULL)
		free(pba->stable_params_derived_smg);
	if(pba->gravity_model_smg == stable_params && pba->ddstable_params_derived_smg != NULL)
		free(pba->ddstable_params_derived_smg);
	if(pba->gravity_model_smg == stable_params && pba->stable_params_aux_smg != NULL)
		free(pba->stable_params_aux_smg);
	if(pba->gravity_model_smg == stable_params && (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) && pba->stable_wext_lna_smg != NULL)
		free(pba->stable_wext_lna_smg);
	if(pba->gravity_model_smg == stable_params && (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) && pba->stable_rho_smg != NULL)
		free(pba->stable_rho_smg);
	if(pba->gravity_model_smg == stable_params && (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) && pba->ddstable_rho_smg != NULL)
		free(pba->ddstable_rho_smg);
	if(pba->gravity_model_smg == stable_params && pba->expansion_model_smg == wext && pba->stable_wext_smg != NULL)
		free(pba->stable_wext_smg);
	if(pba->gravity_model_smg == stable_params && pba->expansion_model_smg == wext && pba->ddstable_wext_smg != NULL)
		free(pba->ddstable_wext_smg);

	return _SUCCESS_;
}


/**
 * Define hi_class bg indices.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */
int background_define_indices_bg_smg(
				 														 struct background *pba,
				 													   int * index_bg
			 													 		 ) {

	class_define_index(pba->index_bg_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_rho_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_rho_prime_smg,pba->rho_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_w_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_current_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_shift_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_delta_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kineticity_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_tensor_excess_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mpl_running_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_kineticity_over_phiphi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_prime_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_kineticity_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mpl_running_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_tensor_excess_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_cs2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_cs2num_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_cs2num_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kinetic_D_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kinetic_D_prime_smg,_TRUE_,*index_bg,1);
	if (pba->field_evolution_smg == _TRUE_) {
		class_define_index(pba->index_bg_kinetic_D_over_phiphi_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_kinetic_D_over_phiphi_prime_smg,_TRUE_,*index_bg,1);
	}
	class_define_index(pba->index_bg_A0_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A1_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A3_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A4_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A5_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A6_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A7_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A8_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A9_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A10_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A11_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A12_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A13_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A14_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A15_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A16_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A9_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A10_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A12_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_A13_prime_smg,_TRUE_,*index_bg,1);

	if (pba->field_evolution_smg == _TRUE_) {
		class_define_index(pba->index_bg_B0_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B1_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B2_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B3_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B4_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B5_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B6_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B7_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B8_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B9_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B10_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B11_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_B12_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C0_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C1_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C2_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C3_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C4_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C5_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C6_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C7_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C8_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C9_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C10_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C11_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C12_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C13_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C14_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C15_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C16_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C9_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C10_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C12_prime_smg,_TRUE_,*index_bg,1);
		class_define_index(pba->index_bg_C13_prime_smg,_TRUE_,*index_bg,1);
	}

	class_define_index(pba->index_bg_lambda_1_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_3_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_4_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_5_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_6_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_7_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_8_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_9_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_10_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_11_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_2_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_8_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_9_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_lambda_11_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_E0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_E3_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_P0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_P1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_P2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R0_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R1_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_R2_smg,pba->field_evolution_smg == _TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_rho_tot_wo_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_H_prime_prime,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_prime_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_prime_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mu_p_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mu_inf_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_muZ_inf_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mu_p_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mu_inf_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_muZ_inf_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_G_eff_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_slip_eff_smg,_TRUE_,*index_bg,1);

  return _SUCCESS_;
}


/**
 * Define hi_class bi indices.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */
int background_define_indices_bi_smg(
				 														 struct background *pba,
				 													   int * index_bi
			 													     ) {

	/* -> scalar field and its derivative wrt conformal time (only if needs to evolve the field)
	* plus other parameters that might be integrated in certain parameterizations
	*/
	class_define_index(pba->index_bi_phi_smg, pba->field_evolution_smg,*index_bi,1);
	class_define_index(pba->index_bi_phi_prime_smg, pba->field_evolution_smg,*index_bi,1);

	//if model needs to integrate M_pl from alpha_M, declare an index
	class_define_index(pba->index_bi_delta_M_pl_smg, pba->M_pl_evolution_smg,*index_bi,1);

	/* index for the smg energy density */
	class_define_index(pba->index_bi_rho_smg, pba->rho_evolution_smg,*index_bi,1);

  return _SUCCESS_;
}

/**
 * Define hi_class bi indices for backward integration.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */
int background_define_indices_bibw_smg(
				 														 struct background *pba,
				 													   int * index_bibw
			 													     ) {


	//index for \tilde{B} and d\tilde{B}/dlna used to compute alpha_B
	// class_define_index(pba->index_bibw_B_tilde_smg,_TRUE_,*index_bibw,1);
	// class_define_index(pba->index_bibw_dB_tilde_smg,_TRUE_,*index_bibw,1);

	//index for \tilde{B}=alpha_B -- 1st order ODE
	class_define_index(pba->index_bibw_B_tilde_smg,_TRUE_,*index_bibw,1);

  return _SUCCESS_;
}

/**
* In background_solve, after quantities {B} and {C} have been
* inegrated, this routine computes derived quantities, such as
* numerical derivatives of the alphas and functions that depend
* on combinations of alphas and derivatives. This is also the
* place where stability tests are done (in stability_tests_smg).
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @param pvecback_bw_integration Output: vector of background quantities to be backward integrated, returned with proper initial values
* @return the error status
*/
int background_solve_smg(
						 struct precision *ppr,
						 struct background *pba,
                         double * pvecback,
                         double * pvecback_integration,// in original hi_class only passed to be freed in case stability tests fail
						 double * pvecback_bw_integration
											   ) {

	/* parameters and workspace for the background_derivs function */
  	struct background_parameters_and_workspace bpaw;
	double * pvec_stable_params_smg;
	/* evolvers */
	extern int evolver_rk(EVOLVER_PROTOTYPE);
	extern int evolver_ndf15(EVOLVER_PROTOTYPE);
	int (*generic_evolver)(EVOLVER_PROTOTYPE) = evolver_ndf15;
	/* initial and final time for backward integration in stable_params */
	double loga_final = log((1. - ppr->eps_bw_integration_braid)/(1+pba->z_gr_smg)); // deep in the matter-dominated era, all MG models considered are effectively GR there
	double loga_ini = 0.; // time used for initial conditions in backward integration
	/* indices for the different arrays */
	int index_loga, index_out;
	/* necessary for calling array_interpolate(), but never used */
	int last_index;
	/* what times are used in the output? */
	int * used_in_output;
	/* re-ordered time array for alpha_B and alpha_K interpolations */
	double * loga_fw_table;
	double loga;
	/* workspace for interpolation of derived parameters */
	double * pvec_stable_params_derived_smg;
	/* flag used to find first loga >= loga_final in loga_table */
	short find_min = _TRUE_;
	/* value of braiding at min_loga */
	double braid_min = 0.;
	/* index of braiding element in background_table with value closest to braid_activation_threshold */
	int braid_target_idx = 0;

    /* When gravity model is "stable parametrization" perform backward integration and overwrite alpha's as well as M_pl */ 
    if(pba->gravity_model_smg == stable_params){

		/** - setup background workspace */
		bpaw.pba = pba;
		bpaw.pvecback = pvecback;
		bpaw.pvecback_B = pvecback_integration;
		class_alloc(pvec_stable_params_smg,pba->num_stable_params*sizeof(double),pba->error_message);
		bpaw.pvec_stable_params_smg = pvec_stable_params_smg;
		class_alloc(pvec_stable_params_derived_smg,pba->num_stable_params_derived*sizeof(double),pba->error_message);
		/* allocate time and output flag arrays */
		pba->bt_bw_size = pba->stable_params_size_smg;
		class_alloc(pba->loga_bw_table,pba->bt_bw_size*sizeof(double),pba->error_message);
		class_alloc(used_in_output, pba->bt_bw_size*sizeof(int), pba->error_message);
		class_alloc(loga_fw_table,pba->bt_bw_size*sizeof(double),pba->error_message);

		/** - define values of loga at which results will be stored from loga_ini to loga_final*/
		for (index_loga=0; index_loga<pba->bt_bw_size; index_loga++) {
			pba->loga_bw_table[index_loga] = loga_ini + index_loga*(loga_final-loga_ini)/(pba->bt_bw_size-1);
			loga_fw_table[index_loga] = loga_final + index_loga*(loga_ini-loga_final)/(pba->bt_bw_size-1);
		}
		/** - define which parameters will be used in the output */
		for (index_out=0; index_out<pba->bi_bw_B_size; index_out++) {
			used_in_output[index_out] = 1;
		}
		
		class_call_except(generic_evolver(background_derivs_bw_smg,
									loga_ini,
									loga_final,
									pvecback_bw_integration,
									used_in_output,
									pba->bi_bw_B_size,
									&bpaw,
									ppr->tol_background_bw_integration,
									ppr->smallest_allowed_variation,
									background_timescale, //'evaluate_timescale', required by evolver_rk but not by ndf15
									ppr->background_integration_stepsize,
									pba->loga_bw_table,
									pba->bt_bw_size,
									background_sources_bw_smg,
									NULL, //'print_variables' in evolver_rk could be set, but, not required
									pba->error_message),
					pba->error_message,
					pba->error_message,
					free(pba->loga_bw_table);
					free(pvecback);
					free(pvecback_integration);
					free(pvecback_bw_integration);
					free(used_in_output);
					free(loga_fw_table);
					free(pvec_stable_params_smg);
					free(pvec_stable_params_derived_smg);
					background_free(pba);
					);

		// Spline new alpha_B and alpha_K
		class_call(array_spline_table_lines(loga_fw_table,
                                        	pba->bt_bw_size,
											pba->stable_params_derived_smg,
											pba->num_stable_params_derived,
											pba->ddstable_params_derived_smg,
											_SPLINE_EST_DERIV_,
											pba->error_message),
               pba->error_message,
               pba->error_message);

		// Update background table
		for (index_loga=0; index_loga<pba->bt_size; index_loga++) {

			loga = pba->loga_table[index_loga];

			if(loga >= loga_final){
				// interpolate Delta_M_pl^2, D_kin, cs2 and alpha_M
				class_call(array_interpolate_spline(
										pba->stable_params_lna_smg,
										pba->stable_params_size_smg,
										pba->stable_params_smg,
										pba->ddstable_params_smg,
										pba->num_stable_params,
										loga,
										&last_index, // not used
										pvec_stable_params_smg,
										pba->num_stable_params,
										pba->error_message),
				pba->error_message,
				pba->error_message);
				// interpolate alpha_B, alpha_B' and alpha_K
					class_call(array_interpolate_spline(
										loga_fw_table,
										pba->bt_bw_size,
										pba->stable_params_derived_smg,
										pba->ddstable_params_derived_smg,
										pba->num_stable_params_derived,
										loga,
										&last_index, // not used
										pvec_stable_params_derived_smg,
										pba->num_stable_params_derived,
										pba->error_message),
					pba->error_message,
					pba->error_message);

				/** store value of braiding at loga_final */
				if(find_min == _TRUE_){
					braid_min = pvec_stable_params_derived_smg[pba->index_derived_braiding_smg];
					find_min = _FALSE_;
				}

				copy_to_background_table_smg(pba, index_loga, pba->index_bg_delta_M2_smg, pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg]);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_M2_smg, 1. + pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg]);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_mpl_running_smg, pvec_stable_params_smg[pba->index_stable_Mpl_running_smg]);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_kinetic_D_smg, pvec_stable_params_smg[pba->index_stable_Dkin_smg] + pba->kineticity_safe_smg);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_cs2_smg, pvec_stable_params_smg[pba->index_stable_cs2_smg]);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_cs2num_smg, pvec_stable_params_smg[pba->index_stable_cs2_smg]*(pvec_stable_params_smg[pba->index_stable_Dkin_smg] + pba->kineticity_safe_smg));
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_kineticity_smg, pvec_stable_params_derived_smg[pba->index_derived_kineticity_smg] + pba->kineticity_safe_smg);				
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_braiding_smg, pvec_stable_params_derived_smg[pba->index_derived_braiding_smg]);
				/* stable_params only works for scalar Horndeski */
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_tensor_excess_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_beyond_horndeski_smg, 0.);

				if(pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) {
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_delta_M2_smg] = pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg];
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_M2_smg] = 1. + pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg];
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_mpl_running_smg] = pvec_stable_params_smg[pba->index_stable_Mpl_running_smg];
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_kinetic_D_smg] = pvec_stable_params_smg[pba->index_stable_Dkin_smg] + pba->kineticity_safe_smg;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_cs2_smg] = pvec_stable_params_smg[pba->index_stable_cs2_smg];
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_cs2num_smg] = pvec_stable_params_smg[pba->index_stable_cs2_smg]*(pvec_stable_params_smg[pba->index_stable_Dkin_smg] + pba->kineticity_safe_smg);
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_kineticity_smg] = pvec_stable_params_derived_smg[pba->index_derived_kineticity_smg] + pba->kineticity_safe_smg;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_braiding_smg] = pvec_stable_params_derived_smg[pba->index_derived_braiding_smg];
					/* stable_params only works for scalar Horndeski */
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_tensor_excess_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_beyond_horndeski_smg] = 0.;
				}
			}
			else {
				//Set Horndeski parameters to their GR-LCDM limit. Not used anyway				
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_delta_M2_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_M2_smg, 1.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_mpl_running_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_kinetic_D_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_cs2_smg, 1.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_braiding_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_kineticity_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_tensor_excess_smg, 0.);
				copy_to_background_table_smg(pba, index_loga, pba->index_bg_beyond_horndeski_smg, 0.);

				if(pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) {
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_delta_M2_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_M2_smg] = 1.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_mpl_running_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_kinetic_D_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_cs2_smg] = 1.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_kineticity_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_braiding_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_tensor_excess_smg] = 0.;
					pba->background_table_late[index_loga*pba->bg_size + pba->index_bg_beyond_horndeski_smg] = 0.;
				}
			}

		}

		// Re-spline background table
		class_call(array_spline_table_lines(
									pba->loga_table,
									pba->bt_size,
									pba->background_table,
									pba->bg_size,
									pba->d2background_dloga2_table,
									_SPLINE_EST_DERIV_,
									pba->error_message),
			pba->error_message,
			pba->error_message);

		if(pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de){
			class_call(array_spline_table_lines(pba->loga_table,
											pba->bt_size,
											pba->background_table_late,
											pba->bg_size,
											pba->d2background_dloga2_table_late,
											_SPLINE_EST_DERIV_,
											pba->error_message),
					pba->error_message,
					pba->error_message);
		}

		/** if braiding smaller than braid_activation_threshold and braid(0) larger that this threshold z_gr_smg such that abs(bra(z_gr_smg))~1e-12 */
		if(fabs(braid_min) < ppr->braid_activation_threshold && fabs(pba->parameters_2_smg[0]) > ppr->braid_activation_threshold){
			class_call(array_find_closest_background_table(pba->background_table,
														pba->bg_size,
														pba->index_bg_braiding_smg, 
														pba->bt_size, 
														ppr->braid_activation_threshold,
														&braid_target_idx,
														pba->error_message),
					pba->error_message,
					pba->error_message);
			/* update GR -> SMG transition redshift used by the perturbations */
			pba->z_gr_smg = 1/exp(pba->loga_table[braid_target_idx]) - 1.;
			
			if (pba->background_verbose > 0) {
				printf(" -> transition to smg delayed to z_gr_smg = %e\n",pba->z_gr_smg);
			}
		}

    } // End stable parametrization



	/* define local variables */
	double a;
	int i;

	/** - second loop over lines, overwrite derivatives that can't be analytically computed from background_functions
	 * Fill the derivatives of the Bellini-Sawicki functions in pvecback
	 * This is done just by overwriting the pvecback entries corresponding to the relevant indice
	 */
	double * pvecback_derivs;
	class_alloc(pvecback_derivs,pba->bg_size*sizeof(double),pba->error_message);
	
 	for (i=0; i < pba->bt_size; i++) {

		if (pba->gravity_model_smg == stable_params && (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de)) {
			// write the derivatives in the structure
			class_call(array_derivate_spline(pba->loga_table, // x_array
						pba->bt_size, // int n_lines
						pba->background_table_late, // array
						pba->d2background_dloga2_table_late, // double * array_splined
						pba->bg_size, // n_columns
						pba->loga_table[i], // double x -> tau
						&last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
						pvecback_derivs, // double * result
						pba->bg_size, //result_size, from 1 to n_columns
						pba->error_message),
			pba->error_message,
			pba->error_message);	

			/* - indices for scalar field (modified gravity) */
			a = pba->background_table_late[i*pba->bg_size + pba->index_bg_a];

			class_call(derivatives_alphas_smg(pba, pba->background_table_late + i*pba->bg_size, pvecback_derivs, i),
				pba->error_message,
				pba->error_message
			);

			class_call(gravity_functions_As_from_alphas_smg(pba, pba->background_table_late + i*pba->bg_size, pvecback_derivs),
								pba->error_message,
								pba->error_message);

			copy_to_background_table_smg(pba, i, pba->index_bg_A0_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A0_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A1_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A1_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A2_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A2_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A3_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A3_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A4_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A4_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A5_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A5_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A6_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A6_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A7_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A7_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A8_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A8_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A9_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A9_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A10_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A10_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A11_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A11_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A12_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A12_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A13_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A13_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A14_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A14_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A15_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A15_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A16_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_A16_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_1_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_1_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_2_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_2_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_3_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_3_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_4_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_4_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_5_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_5_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_6_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_6_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_7_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_7_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_8_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_8_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_9_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_9_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_10_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_10_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_11_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_lambda_11_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_G_eff_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_G_eff_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_slip_eff_smg, pba->background_table_late[i*pba->bg_size + pba->index_bg_slip_eff_smg]);

			/* Here we update the minimum values of the stability quantities
			* even though not required for stable parametrisation 
			*/
			if (a > pba->a_min_stability_test_smg && pba->parameters_tuned_smg == _TRUE_){
				if (pba->background_table[pba->index_bg_kinetic_D_smg + i*pba->bg_size] < pba->min_D_smg){
					pba->min_D_smg = pba->background_table[pba->index_bg_kinetic_D_smg + i*pba->bg_size];
					pba->a_min_D_smg = a;
				}
				if (pba->background_table[pba->index_bg_cs2_smg + i*pba->bg_size] <= pba->min_cs2_smg){
					pba->min_cs2_smg = pba->background_table[pba->index_bg_cs2_smg + i*pba->bg_size];
					pba->a_min_cs2_smg = a;
				}
				if (pba->background_table[pba->index_bg_M2_smg + i*pba->bg_size] < pba->min_M2_smg){
					pba->min_M2_smg = pba->background_table[pba->index_bg_M2_smg + i*pba->bg_size];
					pba->a_min_M2_smg = a;
				}
				if (pba->background_table[pba->index_bg_tensor_excess_smg + i*pba->bg_size] + 1. < pba->min_ct2_smg){
					pba->min_ct2_smg = 1. + pba->background_table[pba->index_bg_tensor_excess_smg + i*pba->bg_size];
					pba->a_min_ct2_smg = a;
				}
				if (pba->background_table[pba->index_bg_braiding_smg + i*pba->bg_size] < pba->min_bra_smg){
					pba->min_bra_smg = pba->background_table[pba->index_bg_braiding_smg + i*pba->bg_size];
				}
				if (pba->background_table[pba->index_bg_braiding_smg + i*pba->bg_size] > pba->max_bra_smg){
					pba->max_bra_smg = pba->background_table[pba->index_bg_braiding_smg + i*pba->bg_size];
				}
			}

		}
		else {
			// write the derivatives in the structure
			class_call(array_derivate_spline(pba->loga_table, // x_array
						pba->bt_size, // int n_lines
						pba->background_table, // array
						pba->d2background_dloga2_table, // double * array_splined
						pba->bg_size, // n_columns
						pba->loga_table[i], // double x -> tau
						&last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
						pvecback_derivs, // double * result
						pba->bg_size, //result_size, from 1 to n_columns
						pba->error_message),
			pba->error_message,
			pba->error_message);

			//Need to update pvecback
			class_call(background_at_tau(pba,
					pba->tau_table[i],
					long_info,
					inter_normal,
					&last_index, //should be no problem to use the same one as for the derivatives
					pvecback),
			pba->error_message,
			pba->error_message);

			a = pvecback[pba->index_bg_a];

			/* - indices for scalar field (modified gravity) */
			class_call(derivatives_alphas_smg(pba, pvecback, pvecback_derivs, i),
				pba->error_message,
				pba->error_message
			);

			class_call(gravity_functions_As_from_alphas_smg(pba, pvecback, pvecback_derivs),
								pba->error_message,
								pba->error_message);

			if(pba->gravity_model_smg != stable_params){
				copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_smg, pvecback[pba->index_bg_kinetic_D_smg]);
			}
			copy_to_background_table_smg(pba, i, pba->index_bg_A0_smg, pvecback[pba->index_bg_A0_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A1_smg, pvecback[pba->index_bg_A1_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A2_smg, pvecback[pba->index_bg_A2_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A3_smg, pvecback[pba->index_bg_A3_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A4_smg, pvecback[pba->index_bg_A4_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A5_smg, pvecback[pba->index_bg_A5_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A6_smg, pvecback[pba->index_bg_A6_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A7_smg, pvecback[pba->index_bg_A7_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A8_smg, pvecback[pba->index_bg_A8_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A9_smg, pvecback[pba->index_bg_A9_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A10_smg, pvecback[pba->index_bg_A10_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A11_smg, pvecback[pba->index_bg_A11_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A12_smg, pvecback[pba->index_bg_A12_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A13_smg, pvecback[pba->index_bg_A13_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A14_smg, pvecback[pba->index_bg_A14_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A15_smg, pvecback[pba->index_bg_A15_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_A16_smg, pvecback[pba->index_bg_A16_smg]);

			if (pba->field_evolution_smg == _TRUE_) {

				class_call(gravity_functions_Cs_from_Bs_smg(pba, pvecback, pvecback_derivs),
									pba->error_message,
									pba->error_message);

				copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_over_phiphi_smg, pvecback[pba->index_bg_kinetic_D_over_phiphi_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C0_smg, pvecback[pba->index_bg_C0_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C1_smg, pvecback[pba->index_bg_C1_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C2_smg, pvecback[pba->index_bg_C2_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C3_smg, pvecback[pba->index_bg_C3_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C4_smg, pvecback[pba->index_bg_C4_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C5_smg, pvecback[pba->index_bg_C5_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C6_smg, pvecback[pba->index_bg_C6_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C7_smg, pvecback[pba->index_bg_C7_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C8_smg, pvecback[pba->index_bg_C8_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C9_smg, pvecback[pba->index_bg_C9_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C10_smg, pvecback[pba->index_bg_C10_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C11_smg, pvecback[pba->index_bg_C11_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C12_smg, pvecback[pba->index_bg_C12_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C13_smg, pvecback[pba->index_bg_C13_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C14_smg, pvecback[pba->index_bg_C14_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C15_smg, pvecback[pba->index_bg_C15_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_C16_smg, pvecback[pba->index_bg_C16_smg]);
			}

			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_1_smg, pvecback[pba->index_bg_lambda_1_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_2_smg, pvecback[pba->index_bg_lambda_2_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_3_smg, pvecback[pba->index_bg_lambda_3_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_4_smg, pvecback[pba->index_bg_lambda_4_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_5_smg, pvecback[pba->index_bg_lambda_5_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_6_smg, pvecback[pba->index_bg_lambda_6_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_7_smg, pvecback[pba->index_bg_lambda_7_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_8_smg, pvecback[pba->index_bg_lambda_8_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_9_smg, pvecback[pba->index_bg_lambda_9_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_10_smg, pvecback[pba->index_bg_lambda_10_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_lambda_11_smg, pvecback[pba->index_bg_lambda_11_smg]);
			if(pba->gravity_model_smg != stable_params){
				copy_to_background_table_smg(pba, i, pba->index_bg_cs2num_smg, pvecback[pba->index_bg_cs2num_smg]);
				copy_to_background_table_smg(pba, i, pba->index_bg_cs2_smg, pvecback[pba->index_bg_cs2_smg]);
			}
			copy_to_background_table_smg(pba, i, pba->index_bg_G_eff_smg, pvecback[pba->index_bg_G_eff_smg]);
			copy_to_background_table_smg(pba, i, pba->index_bg_slip_eff_smg, pvecback[pba->index_bg_slip_eff_smg]);

			/* Here we update the minimum values of the stability quantities
			* test will be performed based on the lowest values
			*/
			if (a > pba->a_min_stability_test_smg && pba->parameters_tuned_smg == _TRUE_){
				if (pvecback[pba->index_bg_kinetic_D_smg] < pba->min_D_smg){
					pba->min_D_smg = pvecback[pba->index_bg_kinetic_D_smg];
					pba->a_min_D_smg = a;
				}
				if (pvecback[pba->index_bg_cs2_smg] <= pba->min_cs2_smg){
					pba->min_cs2_smg = pvecback[pba->index_bg_cs2_smg];
					pba->a_min_cs2_smg = a;
				}
				if (pvecback[pba->index_bg_M2_smg] < pba->min_M2_smg){
					pba->min_M2_smg = pvecback[pba->index_bg_M2_smg];
					pba->a_min_M2_smg = a;
				}
				if (pvecback[pba->index_bg_tensor_excess_smg] + 1. < pba->min_ct2_smg){
					pba->min_ct2_smg = 1. + pvecback[pba->index_bg_tensor_excess_smg];
					pba->a_min_ct2_smg = a;
				}
				if (pvecback[pba->index_bg_braiding_smg] < pba->min_bra_smg){
					pba->min_bra_smg = pvecback[pba->index_bg_braiding_smg];
				}
				if (pvecback[pba->index_bg_braiding_smg] > pba->max_bra_smg){
					pba->max_bra_smg = pvecback[pba->index_bg_braiding_smg];
				}
			}

		}
	}
	
	class_call(
		stability_tests_smg(pba, pvecback, pvecback_integration),
		pba->error_message,
		pba->error_message
	);

	/* Check that braiding doesn't cross 2., else we have instability at some point during evolution */
	class_test_except(pba->min_bra_smg < 2. && pba->max_bra_smg > 2.,
	      pba->error_message,
	      free(pvecback);free(pvecback_integration);background_free(pba),
	      "Instability for scalar field perturbations with alpha_B crossing 2 during evolution: alpha_B_min = %e and alpha_B_max = %e.\n", pba->min_bra_smg, pba->max_bra_smg);

	 /* Yet another (third!) loop to make sure the background table makes sense
	 */
	double factor;
	/* These below are used to compute QSA quantities */
	double cB, cM, M2, cs2num, H, rho_tot, p_tot, rho_tot_wo_smg, p_tot_wo_smg, rho_smg, p_smg;
	double p_tot_wo_smg_prime, p_smg_prime, p_tot_wo_smg_prime_prime, p_smg_prime_prime, cs2num_p, cM_p;
	for (i=0; i < pba->bt_size; i++) {
		//Need to update pvecback
		class_call(background_at_tau(pba,
					pba->tau_table[i],
					long_info,
					inter_normal,
					&last_index, //should be no problem to use the same one as for the derivatives
					pvecback),
			pba->error_message,
			pba->error_message);

		if (pba->gravity_model_smg == stable_params && (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de)) {
			// For wext or rho_de expansion we want to make sure we calculate the appropriate scale factor a*H 
			// when extending smg to slightly earlier times than z_gr_smg.
			factor = pba->background_table_late[i*pba->bg_size + pba->index_bg_a]*pba->background_table_late[i*pba->bg_size + pba->index_bg_H];
		}
		else {
			// TODO_EB: note that the derivative is now calculated w.r.t. loga, while our _prime are w.r.t. tau
			factor = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H];
		}

		double d_over_dtau;

		//write the derivatives in the structure
		class_call(array_derivate_spline(pba->loga_table, // x_array
					pba->bt_size, // int n_lines
					pba->background_table, // array
					pba->d2background_dloga2_table, // double * array_splined
					pba->bg_size, // n_columns
					pba->loga_table[i], // double x -> tau
					&last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
					pvecback_derivs, // double * result
					pba->bg_size, //result_size, from 1 to n_columns
					pba->error_message),
		pba->error_message,
		pba->error_message);

		//cs2num'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_cs2num_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_cs2num_prime_smg, d_over_dtau);

		//D'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_kinetic_D_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_prime_smg, d_over_dtau);

		if (pba->field_evolution_smg == _TRUE_) {
			//D_over_phiphi'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_kinetic_D_over_phiphi_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_kinetic_D_over_phiphi_prime_smg, d_over_dtau);
		}

		//A9'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A9_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A9_prime_smg, d_over_dtau);

		//A10'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A10_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A10_prime_smg, d_over_dtau);

		//A12'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A12_smg];

		//A13'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_A13_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_A13_prime_smg, d_over_dtau);

		if (pba->field_evolution_smg == _TRUE_) {
		//C9'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C9_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C9_prime_smg, d_over_dtau);

		//C10'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C10_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C10_prime_smg, d_over_dtau);

		//C12'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C12_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C12_prime_smg, d_over_dtau);

		//C13'
			d_over_dtau = factor*pvecback_derivs[pba->index_bg_C13_smg];
			copy_to_background_table_smg(pba, i, pba->index_bg_C13_prime_smg, d_over_dtau);
		}

		//lambda_2'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_2_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_2_prime_smg, d_over_dtau);

		//lambda_8'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_8_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_8_prime_smg, d_over_dtau);

		//lambda_9'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_9_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_9_prime_smg, d_over_dtau);

		//lambda_11'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_lambda_11_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_lambda_11_prime_smg, d_over_dtau);

		//p_tot_wo_smg''
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_tot_wo_prime_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_p_tot_wo_prime_prime_smg, d_over_dtau);

		//p_smg''
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_prime_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_p_prime_prime_smg, d_over_dtau);

		//calculations for mu_p, mu_inf, muZ_inf and derivatives when using stable parametrisation
		if(pba->gravity_model_smg == stable_params){
			/** QSA expressions */
			a = pba->background_table[i*pba->bg_size + pba->index_bg_a];
			H = pba->background_table[i*pba->bg_size + pba->index_bg_H];
			cB = pba->background_table[i*pba->bg_size + pba->index_bg_braiding_smg];
			cM = pba->background_table[i*pba->bg_size + pba->index_bg_mpl_running_smg];
			M2 = pba->background_table[i*pba->bg_size + pba->index_bg_M2_smg];
			cs2num = pba->background_table[i*pba->bg_size + pba->index_bg_cs2num_smg];
			cs2num_p = pba->background_table[i*pba->bg_size + pba->index_bg_cs2num_prime_smg];
			cM_p = pba->background_table[i*pba->bg_size + pba->index_bg_mpl_running_prime_smg];
			rho_tot = pba->background_table[i*pba->bg_size + pba->index_bg_rho_tot];
			rho_tot_wo_smg = pba->background_table[i*pba->bg_size + pba->index_bg_rho_tot_wo_smg];
			p_tot = pba->background_table[i*pba->bg_size + pba->index_bg_p_tot];
			p_tot_wo_smg = pba->background_table[i*pba->bg_size + pba->index_bg_p_tot_wo_smg];
			p_tot_wo_smg_prime = pba->background_table[i*pba->bg_size + pba->index_bg_p_tot_wo_prime_smg];
			p_smg_prime = pba->background_table[i*pba->bg_size + pba->index_bg_p_prime_smg];
			p_tot_wo_smg_prime_prime = pba->background_table[i*pba->bg_size + pba->index_bg_p_tot_wo_prime_prime_smg];
			p_smg_prime_prime = pba->background_table[i*pba->bg_size + pba->index_bg_p_prime_prime_smg];

			/** mu_p */
			pba->background_table[i*pba->bg_size + pba->index_bg_mu_p_smg] = 9./(4*a*pow(H,3.)) * (a*H*(2.*cs2num + (cB - 2.)*cB + 4.*(cB - 1.)*cM)*(rho_tot + p_tot) + 2.*cB*(p_tot_wo_smg_prime + p_smg_prime));
			/** mu_inf */
			pba->background_table[i*pba->bg_size + pba->index_bg_mu_inf_smg] = (2.*cs2num + pow(cB + 2.*cM, 2.))/(2.*cs2num*M2);
			/** muZ_inf */
			pba->background_table[i*pba->bg_size + pba->index_bg_muZ_inf_smg] = (2.*cs2num + cB*(cB + 2.*cM))/(2.*cs2num*M2);
			/** mu_p' */
			pba->background_table[i*pba->bg_size + pba->index_bg_mu_p_prime_smg] = (9*(2*a*pow(H,3)*(2*cs2num + (-2 + cB)*cB + 4*(-1 + cB)*cM)*(p_tot + rho_tot) - 3*a*H*(2*cs2num + (-2 + cB)*cB 
																				+ 4*(-1 + cB)*cM)*pow(p_tot + rho_tot,2) + 4*(cs2num*pow(H,2) + (3*(p_tot_wo_smg + rho_tot_wo_smg))/M2 + ((-2 + cB)*(3*p_tot 
																				+ pow(H,2)*(cB + 2*cM) + 3*rho_tot))/2.)*(p_tot_wo_smg_prime + p_smg_prime) + 2*pow(H,2)*(2*cs2num + (-2 + cB)*cB 
																				+ 4*(-1 + cB)*cM)*(-3*a*H*(p_tot + rho_tot) + p_tot_wo_smg_prime + p_smg_prime) - 2*pow(H,2)*(a*H*(2*cs2num 
																				+ (-2 + cB)*cB + 4*(-1 + cB)*cM)*(p_tot + rho_tot) + 2*cB*(p_tot_wo_smg_prime + p_smg_prime)) + 9*(p_tot 
																				+ rho_tot)*(a*H*(2*cs2num + (-2 + cB)*cB + 4*(-1 + cB)*cM)*(p_tot + rho_tot) + 2*cB*(p_tot_wo_smg_prime + p_smg_prime)) 
																				+ (2*H*(p_tot + rho_tot)*(a*(-1 + cB + 2*cM)*(6*(p_tot_wo_smg + rho_tot_wo_smg) + M2*(2*cs2num*pow(H,2) + (-2 + cB)*(3*p_tot
																				+ pow(H,2)*(cB + 2*cM) + 3*rho_tot))) + 2*H*M2*(cs2num_p + 2*(-1 + cB)*cM_p)))/M2 + (4*H*cB*(p_tot_wo_smg_prime_prime + p_smg_prime_prime))/a))/(8.*pow(H,4));
			/** mu_inf' */
			pba->background_table[i*pba->bg_size + pba->index_bg_mu_inf_prime_smg] = -0.5*(a*cs2num*H*M2*cM*(2*cs2num + pow(cB + 2*cM,2)) - (a*cs2num*(cB + 2*cM)*(6*(p_tot_wo_smg + rho_tot_wo_smg) 
                																	+ M2*(2*cs2num*pow(H,2) + (-2 + cB)*(3*p_tot + pow(H,2)*(cB + 2*cM) + 3*rho_tot))))/H 
                																	+ M2*(cB + 2*cM)*((cB + 2*cM)*cs2num_p - 4*cs2num*cM_p))/(pow(cs2num,2)*pow(M2,2));
			/** muZ_inf' */
			pba->background_table[i*pba->bg_size + pba->index_bg_muZ_inf_prime_smg] = (-2*a*pow(cs2num,2)*H*cM - a*cs2num*H*cB*cM*(cB + 2*cM) + (a*cs2num*(cB + 2*cM)*(cs2num*pow(H,2) + (3*(p_tot_wo_smg + rho_tot_wo_smg))/M2 + 
                  																	((-2 + cB)*(3*p_tot + pow(H,2)*(cB + 2*cM) + 3*rho_tot))/2.))/H - cB*(cB + 2*cM)*cs2num_p 
                  																	+ cs2num*cB*((a*(cs2num*pow(H,2) + (3*(p_tot_wo_smg + rho_tot_wo_smg))/M2 + ((-2 + cB)*(3*p_tot + pow(H,2)*(cB + 2*cM) 
                  																	+ 3*rho_tot))/2.))/H + 2*cM_p))/(2.*pow(cs2num,2)*M2);

		}


		int j = 0;
		while (j < pba->bg_size){
			class_test_except(isnan(pvecback[j]) && (pba->parameters_tuned_smg == _TRUE_),
					pba->error_message,
					free(pvecback_derivs);free(pvecback);free(pvecback_integration);background_free(pba),
					"pvecback[%i] = %e at a = %e in background!",j,pvecback[j],pvecback[pba->index_bg_a]);
				j++;
		}

	}
	
	free(pvecback_derivs);  //free the structure
	if(pba->gravity_model_smg == stable_params){
		free(pba->loga_bw_table);
		free(loga_fw_table);
		free(pvec_stable_params_smg);
		free(pvec_stable_params_derived_smg);
	}

	return _SUCCESS_;

}

/**
* In background_solve_rho_smg we integrate the scalar field background
* energy-density when expansion_model = wext
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @param pvecback_bw_integration Output: vector of background quantities to be backward integrated, returned with proper initial values
* @return the error status
*/
int background_solve_rho_smg(
						 struct precision *ppr,
						 struct background *pba,
						 double * pback_rho_smg_bw_integration
											   ) {

	/* parameters and workspace for the background_derivs function */
  	struct background_parameters_and_workspace bpaw;
	/* evolvers */
	extern int evolver_rk(EVOLVER_PROTOTYPE);
	extern int evolver_ndf15(EVOLVER_PROTOTYPE);
	int (*generic_evolver)(EVOLVER_PROTOTYPE) = evolver_ndf15;
	/* initial and final time for backward integration in stable_params */
	double loga_final = log((1. - ppr->eps_bw_integration_rho_smg)/(1+pba->z_gr_smg)); // deep in the matter-dominated era, all MG models considered are effectively GR there
	double loga_ini = 0.; // time used for initial conditions in backward integration
	/* indices for the different arrays */
	int index_loga; // index_out;
	/* necessary for calling array_interpolate(), but never used */
	int last_index;
	/* what times are used in the output? */
	int * used_in_output;
	
	pba->loga_final_bw_integration = loga_final;
	class_test(pba->a_file_lcdm_smg > exp(loga_final),
				pba->error_message,
				"minimum scale factor in input file for w is %f. For wext parametrization it must be < %f for accurate numerical derivatives of smg parameters at late times", pba->a_file_lcdm_smg, exp(loga_final));

	/** - setup background workspace */
	bpaw.pba = pba;
	/* allocate time and output flag arrays */
	pba->bt_bw_rho_smg_size = pba->stable_wext_size_smg;
	class_alloc(pba->loga_bw_table_rho_smg,pba->bt_bw_rho_smg_size*sizeof(double),pba->error_message);
	// class_alloc(pba->loga_bw_table,pba->bt_bw_rho_smg_size*sizeof(double),pba->error_message);
	class_alloc(pba->loga_fw_table_rho_smg,pba->bt_bw_rho_smg_size*sizeof(double),pba->error_message);
	// class_alloc(used_in_output, pba->bt_bw_rho_smg_size*sizeof(int), pba->error_message);
	class_alloc(used_in_output, 1*sizeof(int), pba->error_message);

	/** - define values of loga at which results will be stored from loga_ini to loga_final*/
	for (index_loga=0; index_loga<pba->bt_bw_rho_smg_size; index_loga++) {
		pba->loga_bw_table_rho_smg[index_loga] = loga_ini + index_loga*(loga_final-loga_ini)/(pba->bt_bw_rho_smg_size-1);
		// pba->loga_bw_table[index_loga] = loga_ini + index_loga*(loga_final-loga_ini)/(pba->bt_bw_rho_smg_size-1);
		pba->loga_fw_table_rho_smg[index_loga] = loga_final + index_loga*(loga_ini-loga_final)/(pba->bt_bw_rho_smg_size-1);
		// used_in_output[index_loga] = 1;
	}
	used_in_output[0] = 1;

	/* - store final integration time */
	pba->loga_final_rho_smg = loga_final;
	/** - perform the backward integration for rho_smg*/
	class_call(generic_evolver(background_derivs_bw_rho_smg,
								loga_ini,
								loga_final,
								pback_rho_smg_bw_integration,
								used_in_output,
								1,
								&bpaw,
								ppr->tol_background_bw_integration,
								ppr->smallest_allowed_variation,
								background_timescale, //'evaluate_timescale', required by evolver_rk but not by ndf15
								ppr->background_integration_stepsize,
								pba->loga_bw_table_rho_smg, // pba->loga_bw_table,
								pba->bt_bw_rho_smg_size,
								background_sources_bw_rho_smg,
								NULL, //'print_variables' in evolver_rk could be set, but, not required
								pba->error_message),
				pba->error_message,
				pba->error_message);
	// Spline integrated rho_smg
	class_call(array_spline_table_lines(pba->loga_fw_table_rho_smg,
										pba->bt_bw_rho_smg_size,
										pba->stable_rho_smg,
										1,
										pba->ddstable_rho_smg,
										_SPLINE_EST_DERIV_,
										pba->error_message),
			pba->error_message,
			pba->error_message);

	free(pba->loga_bw_table_rho_smg);
	free(used_in_output);

	return _SUCCESS_;

}

/**
* Interpolation function for the integrated rho_smg and derived p_smg.
*
* @param pba                  Input: pointer to background structure
* @param loga				  Input: scale factor ln(a)
* @param loga_transition	  Input: for loga < loga_transition rho_smg behaves like a cosmological constant
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @return the error status
*/
int interpolate_rho_smg_p_smg(struct background *pba,
						double loga,
                        double loga_transition,
                        double * pvecback
                        ){

double w = 0.;
double drho_de = 0.;
int last_index; // necessary for calling array_interpolate(), but never used 

// TODO_MC: test if loga_transition is larger than loga_final_rho_smg
if(pba->expansion_model_smg == wext) {
	if (loga >= loga_transition) {
	// interpolate integrated rho_smg
	class_call(array_interpolate_spline(
										pba->loga_fw_table_rho_smg,
										pba->bt_bw_rho_smg_size,
										pba->stable_rho_smg,
										pba->ddstable_rho_smg,
										1,
										loga,
										&last_index, // not used
										&pvecback[pba->index_bg_rho_smg],
										1,
										pba->error_message),
				pba->error_message,
				pba->error_message);    
	// interpolate w
	class_call(array_interpolate_spline(
										pba->stable_wext_lna_smg,
										pba->stable_wext_size_smg,
										pba->stable_wext_smg,
										pba->ddstable_wext_smg,
										1,
										loga,
										&last_index,
										&w,
										1,
										pba->error_message),
			pba->error_message,
			pba->error_message);

	pvecback[pba->index_bg_p_smg] = w*pvecback[pba->index_bg_rho_smg];

	} 
	else if (loga < loga_transition) {

		// interpolate integrated rho_smg at transition time
		class_call(array_interpolate_spline(
											pba->loga_fw_table_rho_smg,
											pba->bt_bw_rho_smg_size,
											pba->stable_rho_smg,
											pba->ddstable_rho_smg,
											1,
											loga_transition,
											&last_index, // not used
											&pvecback[pba->index_bg_rho_smg],
											1,
											pba->error_message),
					pba->error_message,
					pba->error_message); 

		pvecback[pba->index_bg_p_smg] = -pvecback[pba->index_bg_rho_smg]; // cosmological constant
	}
}
else if (pba->expansion_model_smg == rho_de) {
	if (loga >= loga_transition) {
		// interpolate rho_smg
		class_call(array_interpolate_spline(
											pba->stable_wext_lna_smg,
											pba->stable_wext_size_smg,
											pba->stable_rho_smg,
											pba->ddstable_rho_smg,
											1,
											loga,
											&last_index, // not used
											&pvecback[pba->index_bg_rho_smg],
											1,
											pba->error_message),
					pba->error_message,
					pba->error_message);    
		// differentiate rho_de
		class_call(array_derivate_spline(
						pba->stable_wext_lna_smg, // x_array
						pba->stable_wext_size_smg, // int n_lines
						pba->stable_rho_smg, // array
						pba->ddstable_rho_smg, // double * array_splined
						1, // n_columns
						loga, // double x -> loga
						&last_index, // int* last_index // this is something for the interpolation to talk to each other when using a loop
						&drho_de, // double * result
						1, //result_size, from 1 to n_columns
						pba->error_message),
			pba->error_message,
			pba->error_message);
		// derive p_smg from continuity equation
		pvecback[pba->index_bg_p_smg] = -(pvecback[pba->index_bg_rho_smg] + drho_de/3.);

	} 
	else if (loga < loga_transition) {

		// interpolate rho_smg at transition time
		class_call(array_interpolate_spline(
											pba->stable_wext_lna_smg,
											pba->stable_wext_size_smg,
											pba->stable_rho_smg,
											pba->ddstable_rho_smg,
											1,
											loga_transition,
											&last_index, // not used
											&pvecback[pba->index_bg_rho_smg],
											1,
											pba->error_message),
					pba->error_message,
					pba->error_message); 

		pvecback[pba->index_bg_p_smg] = -pvecback[pba->index_bg_rho_smg]; // cosmological constant
	}
}

return _SUCCESS_;

}

/**
* Send _smg information to the standard output.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @return the error status
*/
int background_print_stdout_smg(
			  											  struct background *pba,
																double * pvecback,
																double * pvecback_integration
															  ) {

	printf(" -> Omega_smg = %f, wanted %f ",pvecback[pba->index_bg_rho_smg]/pvecback[pba->index_bg_rho_crit], pba->Omega0_smg);
	if(pba->has_lambda == _TRUE_)
		printf(", Omega_Lambda = %f", pba->Omega0_lambda);
	printf("\n");
	if (pba->background_verbose > 3) {
		if (pba->gravity_model_smg != stable_params)
			printf("Minimal stability values: cs2 = %g, ct2 = %g, D = %g, M2 = %g \n",pba->min_cs2_smg,pba->min_ct2_smg,pba->min_D_smg,pba->min_M2_smg);
	}

	if (pba->field_evolution_smg == _TRUE_){
	    pba->xi_0_smg = pvecback_integration[pba->index_bi_phi_prime_smg]*pvecback[pba->index_bg_H]/pow(pba->H0,2);
	    pba->phi_0_smg = pvecback_integration[pba->index_bi_phi_smg];
	}

	gravity_models_print_stdout_smg(pba);

	return _SUCCESS_;
}


/**
* For each integrated variable fix the initial condition.
*
* @param pba                  Input: pointer to background structure
* @param a                    Input: scale factor
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @param pvecback_bw_integration Output: vector of background quantities to be backward integrated, returned with proper initial values
* @param pback_rho_smg_bw_integration Output: scalar quantity rho_smg to be backward integrated, returned with proper initial values
* @param ptr_rho_rad          Input: pointer to the density of radiation
* @return the error status
*/
int background_initial_conditions_smg(
        															struct background *pba,
																			double a,
																			double * pvecback,
        															double * pvecback_integration,
																	double * pvecback_bw_integration,
																			double * ptr_rho_rad
														  			  ) {

	double rho_rad = *ptr_rho_rad;
	double phi_scale, V_scale,p1,p2,p3; //smg related variables
	int i = 0;

	/** - fix initial value of modified gravity
	* run over all possible model cases
	*/

	pba->initial_conditions_set_smg = _FALSE_;

	//default value, can override later
	if (pba->M_pl_evolution_smg ==_TRUE_){
	  pvecback_integration[pba->index_bi_delta_M_pl_smg] = 0.;
	}

	class_call(gravity_models_initial_conditions_smg(pba, a, pvecback, pvecback_integration, pvecback_bw_integration, &rho_rad),
						 pba->error_message,
						 pba->error_message);

	if (pba->field_evolution_smg == _TRUE_){
	  if (pba->background_verbose>3)
	    printf(" -> Initial conditions: phi = %e, phi' = %e \n",pvecback_integration[pba->index_bi_phi_smg],pvecback_integration[pba->index_bi_phi_prime_smg]);

	  class_test_except(!isfinite(pvecback_integration[pba->index_bi_phi_smg]) || !isfinite(pvecback_integration[pba->index_bi_phi_smg]),
		 								 pba->error_message,
	   						 	   free(pvecback);free(pvecback_integration);background_free(pba),
		 							   "initial phi = %e phi_prime = %e -> check initial conditions",
		 							 	 pvecback_integration[pba->index_bi_phi_smg],pvecback_integration[pba->index_bi_phi_smg]);
	 }
	if (pba->M_pl_evolution_smg == _TRUE_)
	  if (pba->background_verbose>3)
	    printf(" -> Initial conditions: delta_M_pl = %e \n",pvecback_integration[pba->index_bi_delta_M_pl_smg]);

	if(pba->gravity_model_smg == stable_params) 
		if(pba->background_verbose>3)
			printf(" -> Initial (a=1) B_tilde = %e dB_tilde = %.10e \n",
				pvecback_bw_integration[pba->index_bibw_B_tilde_smg],pvecback_bw_integration[pba->index_bibw_dB_tilde_smg]);

	return _SUCCESS_;
}

/**
* For each integrated variable fix the initial condition.
*
* @param pba                  Input: pointer to background structure
* @param pback_rho_smg_bw_integration Output: scalar quantity rho_smg to be backward integrated, returned with proper initial values
* @return the error status
*/
int background_ic_rho_smg(
        								struct background *pba,
										double * pback_rho_smg_bw_integration
									) {

	double Omega_const_smg = pba->parameters_smg[0];
    pback_rho_smg_bw_integration[0] = Omega_const_smg * pow(pba->H0,2);	

	if(pba->background_verbose>3)
		printf(" -> Initial (a=1) rho_smg = %e \n",pback_rho_smg_bw_integration[0]);


	return _SUCCESS_;
}

/**
 * Subroutine for formatting background _smg output
 *
 * @param pba                  Input: pointer to background structure
 * @param titles               Ouput: name of columns when printing the background table
 * @return the error status
 */
int background_store_columntitles_smg(
																		  struct background *pba,
																			char titles[_MAXTITLESTRINGLENGTH_]
																		  ) {

	class_store_columntitle(titles,"(.)rho_smg",_TRUE_);
  	class_store_columntitle(titles,"(.)p_smg",_TRUE_);
	class_store_columntitle(titles,"(.)p_smg_prime",_TRUE_);
	class_store_columntitle(titles,"(.)p_smg_prime_prime",_TRUE_);

	if (pba->output_background_smg >= 1){
    class_store_columntitle(titles,"M*^2_smg",_TRUE_);
    class_store_columntitle(titles,"D_M*^2_smg",_TRUE_);
    class_store_columntitle(titles,"kineticity_smg",_TRUE_);
    class_store_columntitle(titles,"braiding_smg",_TRUE_);
    class_store_columntitle(titles,"tensor_excess_smg",_TRUE_);
    class_store_columntitle(titles,"Mpl_running_smg",_TRUE_);
    class_store_columntitle(titles,"beyond_horndeski_smg",_TRUE_);
    class_store_columntitle(titles,"c_s^2",_TRUE_);
    class_store_columntitle(titles,"kin (D)",_TRUE_);
    class_store_columntitle(titles,"Current",pba->field_evolution_smg);
    class_store_columntitle(titles,"Shift",pba->field_evolution_smg);
  }

  if (pba->output_background_smg >= 2){
    class_store_columntitle(titles,"phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"phi'",pba->field_evolution_smg);
    class_store_columntitle(titles,"phi''",pba->field_evolution_smg);
    class_store_columntitle(titles,"E0",pba->field_evolution_smg);
    class_store_columntitle(titles,"E1",pba->field_evolution_smg);
    class_store_columntitle(titles,"E2",pba->field_evolution_smg);
    class_store_columntitle(titles,"E3",pba->field_evolution_smg);
    class_store_columntitle(titles,"P0",pba->field_evolution_smg);
    class_store_columntitle(titles,"P1",pba->field_evolution_smg);
    class_store_columntitle(titles,"P2",pba->field_evolution_smg);
    class_store_columntitle(titles,"R0",pba->field_evolution_smg);
    class_store_columntitle(titles,"R1",pba->field_evolution_smg);
    class_store_columntitle(titles,"R2",pba->field_evolution_smg);
    class_store_columntitle(titles,"G_eff_smg",_TRUE_);
    class_store_columntitle(titles,"slip_eff_smg",_TRUE_);
  }

  //TODO: add in output background trigger
  if (pba->output_background_smg >= 3){
    class_store_columntitle(titles,"kineticity_prime_smg",_TRUE_);
    class_store_columntitle(titles,"braiding_prime_smg",_TRUE_);
	class_store_columntitle(titles,"Mpl_running_prime_smg",_TRUE_);
    class_store_columntitle(titles,"kineticity_over_phiphi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"braiding_over_phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"beyond_horndeski_over_phi_smg",pba->field_evolution_smg);
    class_store_columntitle(titles,"kin_over_phiphi (D)",pba->field_evolution_smg);
    class_store_columntitle(titles,"A0",_TRUE_);
    class_store_columntitle(titles,"A1",_TRUE_);
    class_store_columntitle(titles,"A2",_TRUE_);
    class_store_columntitle(titles,"A3",_TRUE_);
    class_store_columntitle(titles,"A4",_TRUE_);
    class_store_columntitle(titles,"A5",_TRUE_);
    class_store_columntitle(titles,"A6",_TRUE_);
    class_store_columntitle(titles,"A7",_TRUE_);
    class_store_columntitle(titles,"A8",_TRUE_);
    class_store_columntitle(titles,"A9",_TRUE_);
    class_store_columntitle(titles,"A10",_TRUE_);
    class_store_columntitle(titles,"A11",_TRUE_);
    class_store_columntitle(titles,"A12",_TRUE_);
    class_store_columntitle(titles,"A13",_TRUE_);
    class_store_columntitle(titles,"A14",_TRUE_);
    class_store_columntitle(titles,"A15",_TRUE_);
    class_store_columntitle(titles,"A16",_TRUE_);
    class_store_columntitle(titles,"B0",pba->field_evolution_smg);
    class_store_columntitle(titles,"B1",pba->field_evolution_smg);
    class_store_columntitle(titles,"B2",pba->field_evolution_smg);
    class_store_columntitle(titles,"B3",pba->field_evolution_smg);
    class_store_columntitle(titles,"B4",pba->field_evolution_smg);
    class_store_columntitle(titles,"B5",pba->field_evolution_smg);
    class_store_columntitle(titles,"B6",pba->field_evolution_smg);
    class_store_columntitle(titles,"B7",pba->field_evolution_smg);
    class_store_columntitle(titles,"B8",pba->field_evolution_smg);
    class_store_columntitle(titles,"B9",pba->field_evolution_smg);
    class_store_columntitle(titles,"B10",pba->field_evolution_smg);
    class_store_columntitle(titles,"B11",pba->field_evolution_smg);
    class_store_columntitle(titles,"B12",pba->field_evolution_smg);
    class_store_columntitle(titles,"C0",pba->field_evolution_smg);
    class_store_columntitle(titles,"C1",pba->field_evolution_smg);
    class_store_columntitle(titles,"C2",pba->field_evolution_smg);
    class_store_columntitle(titles,"C3",pba->field_evolution_smg);
    class_store_columntitle(titles,"C4",pba->field_evolution_smg);
    class_store_columntitle(titles,"C5",pba->field_evolution_smg);
    class_store_columntitle(titles,"C6",pba->field_evolution_smg);
    class_store_columntitle(titles,"C7",pba->field_evolution_smg);
    class_store_columntitle(titles,"C8",pba->field_evolution_smg);
    class_store_columntitle(titles,"C9",pba->field_evolution_smg);
    class_store_columntitle(titles,"C10",pba->field_evolution_smg);
    class_store_columntitle(titles,"C11",pba->field_evolution_smg);
    class_store_columntitle(titles,"C12",pba->field_evolution_smg);
    class_store_columntitle(titles,"C13",pba->field_evolution_smg);
    class_store_columntitle(titles,"C14",pba->field_evolution_smg);
    class_store_columntitle(titles,"C15",pba->field_evolution_smg);
    class_store_columntitle(titles,"C16",pba->field_evolution_smg);
    class_store_columntitle(titles,"lambda_1",_TRUE_);
    class_store_columntitle(titles,"lambda_2",_TRUE_);
    class_store_columntitle(titles,"lambda_3",_TRUE_);
    class_store_columntitle(titles,"lambda_4",_TRUE_);
    class_store_columntitle(titles,"lambda_5",_TRUE_);
    class_store_columntitle(titles,"lambda_6",_TRUE_);
    class_store_columntitle(titles,"lambda_7",_TRUE_);
    class_store_columntitle(titles,"lambda_8",_TRUE_);
    class_store_columntitle(titles,"lambda_9",_TRUE_);
    class_store_columntitle(titles,"lambda_10",_TRUE_);
    class_store_columntitle(titles,"lambda_11",_TRUE_);
    class_store_columntitle(titles,"lambda_2_p",_TRUE_);
    class_store_columntitle(titles,"lambda_8_p",_TRUE_);
    class_store_columntitle(titles,"lambda_9_p",_TRUE_);
    class_store_columntitle(titles,"lambda_11_p",_TRUE_);
    class_store_columntitle(titles,"cs2num",_TRUE_);
    class_store_columntitle(titles,"cs2num_p",_TRUE_);
	class_store_columntitle(titles,"mu_p",_TRUE_);
	class_store_columntitle(titles,"mu_inf",_TRUE_);
	class_store_columntitle(titles,"muZ_inf",_TRUE_);
	class_store_columntitle(titles,"mu_p_prime",_TRUE_);
	class_store_columntitle(titles,"mu_inf_prime",_TRUE_);
	class_store_columntitle(titles,"muZ_inf_prime",_TRUE_);
  }

	return _SUCCESS_;
}


/**
 * Subroutine for writing the background _smg output
 *
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities
 * @param dataptr              Ouput: pointer to 1d array storing all the background table
 * @param ptr_storeidx         Ouput: pointer to index with number of columns
 * @return the error status
 */
int background_output_data_smg(
															 struct background *pba,
        							 				 double * pvecback,
															 double * dataptr,
															 int * ptr_storeidx
														   ) {

	int storeidx = *ptr_storeidx;

	class_store_double(dataptr,pvecback[pba->index_bg_rho_smg],_TRUE_,storeidx);
  	class_store_double(dataptr,pvecback[pba->index_bg_p_smg],_TRUE_,storeidx);
	class_store_double(dataptr,pvecback[pba->index_bg_p_prime_smg],_TRUE_,storeidx);
	class_store_double(dataptr,pvecback[pba->index_bg_p_prime_prime_smg],_TRUE_,storeidx);

	if (pba->output_background_smg >= 1){
		class_store_double(dataptr,pvecback[pba->index_bg_M2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_delta_M2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_tensor_excess_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mpl_running_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_beyond_horndeski_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kinetic_D_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_current_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_shift_smg],pba->field_evolution_smg,storeidx);
	}

	if (pba->output_background_smg >= 2){
		class_store_double(dataptr,pvecback[pba->index_bg_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_prime_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_E3_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_P2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_R2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_G_eff_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_slip_eff_smg],_TRUE_,storeidx);
	}

	if (pba->output_background_smg >= 3){
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mpl_running_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kineticity_over_phiphi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_braiding_over_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_beyond_horndeski_over_phi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_kinetic_D_over_phiphi_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A0_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A1_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A3_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A4_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A5_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A6_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A7_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A8_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A9_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A10_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A11_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A12_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A13_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A14_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A15_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_A16_smg],_TRUE_,storeidx);

		if (pba->field_evolution_smg == _TRUE_) {
			class_store_double(dataptr,pvecback[pba->index_bg_B0_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B1_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B3_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B4_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B5_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B6_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B7_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B8_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B9_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B10_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B11_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_B12_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C0_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C1_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C3_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C4_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C5_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C6_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C7_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C8_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C9_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C10_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C11_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C12_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C13_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C15_smg],_TRUE_,storeidx);
			class_store_double(dataptr,pvecback[pba->index_bg_C16_smg],_TRUE_,storeidx);
		}

		class_store_double(dataptr,pvecback[pba->index_bg_lambda_1_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_2_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_3_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_4_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_5_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_6_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_7_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_8_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_9_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_10_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_11_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_2_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_8_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_9_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_lambda_11_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2num_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_cs2num_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mu_p_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mu_inf_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_muZ_inf_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mu_p_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_mu_inf_prime_smg],_TRUE_,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_muZ_inf_prime_smg],_TRUE_,storeidx);
	}

	*ptr_storeidx = storeidx;

	return _SUCCESS_;
}


/**
 * Derivatives of {B} _smg quantities
 *
 * @param pba                  Input: pointer to background structure
 * @param a                    Input: scale factor
 * @param pvecback             Input: vector of background quantities
 * @param y                    Input: current vector of integrated quantities (with index_bi)
 * @param dy                   Output: current derivative of y w.r.t log(a)
 * @return the error status
 */
int background_derivs_smg(
			  									struct background *pba,
													double a,
													double * pvecback,
													double * y,
													double * dy
												  ) {

	double H = pvecback[pba->index_bg_H];

	/** - calculate /f$ \rho'(\tau)= -3aH (1+w) \rho /f$ written as \f$ d\rho/dloga = \rho' / (aH) \f$ */
	if (pba->rho_evolution_smg == _TRUE_){
	  dy[pba->index_bi_rho_smg] = pvecback[pba->index_bg_rho_prime_smg]/a/H;
	}

	/** - Scalar field equation: \f$ \phi''(t) + 2 a H \phi'(t) + a^2 dV = 0 \f$
				(note H is wrt cosmic time) **/
	if(pba->field_evolution_smg){
	  dy[pba->index_bi_phi_smg] = y[pba->index_bi_phi_prime_smg]/a/H;
	  dy[pba->index_bi_phi_prime_smg] = pvecback[pba->index_bg_phi_prime_prime_smg]/a/H;
	}
	/** - Planck mass equation (if parameterization in terms of alpha_m **/
	if (pba->M_pl_evolution_smg == _TRUE_)
	  dy[pba->index_bi_delta_M_pl_smg] = pvecback[pba->index_bg_mpl_running_smg]*(y[pba->index_bi_delta_M_pl_smg] + 1);   //in this case the running has to be integrated (eq 3.3 of 1404.3713 yields M2' = aH\alpha_M)

	return _SUCCESS_;
}

/**
 * Derivatives of {bw_B} _smg quantities
 *
 * @param pba                  Input: pointer to background structure
 * @param loga                    Input: scale factor
 * @param y                    Input: current vector of integrated quantities (with index_bi)
 * @param dy                   Output: current derivative of y w.r.t log(a)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_derivs_bw_smg(
                        double loga,
                    	double * y, /* vector with argument y[index_bi_bw] (must be already allocated with size pba->bi_bw_B_size) */
                        double * dy, /* vector with argument dy[index_bi_bw]
                                     (must be already allocated with
                                     size pba->bi_bw_B_size) */
                        void * parameters_and_workspace,
                        ErrorMsg error_message
                        ){

	// fill in workspace and interpolate all relavant quantities. Also find rho_tot and P_tot including all spieces except dark energy component
	/** - define local variables */
	struct background_parameters_and_workspace * pbpaw;
	struct background * pba;
	double * pvecback;
	double * pvecback_B; // used to infer non-MG integrated quantities from interpolation and fed to background_functions
	double * pvec_stable_params_smg;
	double a, z, H, dH, rho_tot, P_tot, rho_smg, p_smg;
	double Delta_Mpl, alpha_M, D_kin, cs2;
	int last_index; // necessary for calling array_interpolate(), but never used
	int pvecback_size;
	int pvec_stable_params_size; // size of output vector, controlled by input parameter return_format

	pbpaw = parameters_and_workspace;
  	pba =  pbpaw->pba;
	pvecback = pbpaw->pvecback;
	pvecback_B = pbpaw->pvecback_B;
	pvecback_size = pba->bg_size;
	pvec_stable_params_smg = pbpaw->pvec_stable_params_smg;
	pvec_stable_params_size = pba->num_stable_params;

	a = exp(loga);
	// Interpolate dH/dlna from background_table_late (GR->MG transition happening at loga_final_bw_integration) if expansion parametrized by rho_de
	// TODO_MC: do we still need to evaluate dH here? We do that below anyway (!!!)
	// if (pba->expansion_model_smg == rho_de) {
	// 	class_call(array_interpolate_spline(
    //                                     pba->loga_table,
    //                                     pba->bt_size,
    //                                     pba->background_table_late,
    //                                     pba->d2background_dloga2_table_late,
    //                                     pba->bg_size,
    //                                     loga,
    //                                     &last_index,
    //                                     pvecback,
    //                                     pvecback_size,
    //                                     pba->error_message),
    //            pba->error_message,
    //            pba->error_message);

	// 	dH = pvecback[pba->index_bg_H_prime]; // this is dH/dtau
	// }

	// Assign pvecback_B elements based on interpolation background table above
	if (pba->has_dcdm == _TRUE_ || pba->has_dr == _TRUE_) {
		class_call(array_interpolate_spline(
                                        pba->loga_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dloga2_table,
                                        pba->bg_size,
                                        loga,
                                        &last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
		/* dcdm */
		if (pba->has_dcdm == _TRUE_) {
			pvecback_B[pba->index_bi_rho_dcdm] = pvecback[pba->index_bg_rho_dcdm];
		}
		/* dr */
		if (pba->has_dr == _TRUE_) {
			pvecback_B[pba->index_bi_rho_dr] = pvecback[pba->index_bg_rho_dr];
		}
	}

	// call directly background_functions to speed up entire code and increase accuracy
	 class_call(background_functions(pba, a, pvecback_B, normal_info, pvecback),
             pba->error_message,
             error_message);

	if (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) {
		class_call(interpolate_rho_smg_p_smg(pba, log(a), pba->loga_final_bw_integration, pvecback),
                 pba->error_message,
                 pba->error_message
      	);
		// Update quantities depending on rho_smg and p_smg
		H = sqrt(pvecback[pba->index_bg_rho_tot_wo_smg] + pvecback[pba->index_bg_rho_smg] - pba->K/a/a);
		dH = (-1.5*a*(pvecback[pba->index_bg_rho_tot_wo_smg] + pvecback[pba->index_bg_rho_smg]
													+ pvecback[pba->index_bg_p_tot_wo_smg] + pvecback[pba->index_bg_p_smg]) + pba->K/a)/a/H; // dH/dloga = 1/aH * dH/dtau
	}
	else {
		H = pvecback[pba->index_bg_H];
		dH = pvecback[pba->index_bg_H_prime]/a/H; // dH/dloga = 1/aH * dH/dtau
	}
	
	rho_smg = pvecback[pba->index_bg_rho_smg];
	p_smg = pvecback[pba->index_bg_p_smg];
	rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg]; // all matter exluding scalar field
	P_tot = pvecback[pba->index_bg_p_tot_wo_smg]; // all matter excluding scalar field

	// interpolate Delta_M_pl^2, D_kin, cs2 and alpha_M
	class_call(array_interpolate_spline(
                                        pba->stable_params_lna_smg,
                                        pba->stable_params_size_smg,
                                        pba->stable_params_smg,
                                        pba->ddstable_params_smg,
                                        pba->num_stable_params,
                                        loga,
                                        &last_index,
                                        pvec_stable_params_smg,
                                        pvec_stable_params_size,
                                        error_message),
               error_message,
               error_message);

	Delta_Mpl = pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg];
	D_kin = pvec_stable_params_smg[pba->index_stable_Dkin_smg];
	cs2 = pvec_stable_params_smg[pba->index_stable_cs2_smg];
	alpha_M = pvec_stable_params_smg[pba->index_stable_Mpl_running_smg];

	/* Derivatives here are w.r.t. loga */

	/* system of two first-order linear ODE */
	// dy[pba->index_bibw_B_tilde_smg] = y[pba->index_bibw_dB_tilde_smg];
	// dy[pba->index_bibw_dB_tilde_smg] = (1. + alpha_M - dH/H)*y[pba->index_bibw_dB_tilde_smg] - (1.5 * (rho_tot + P_tot)/(H * H * (1.+Delta_Mpl)) + 0.5 * D_kin * cs2)*y[pba->index_bibw_B_tilde_smg];

	/* single first-order non-linear ODE */
	// dy[pba->index_bibw_B_tilde_smg] = (y[pba->index_bibw_B_tilde_smg] - 2.) * (0.5 * y[pba->index_bibw_B_tilde_smg] + alpha_M - dH/H) + D_kin * cs2 + 3. * (rho_tot + P_tot)/(H * H * (1.+Delta_Mpl));
	/* re-written for increased numerical accuracy */
	dy[pba->index_bibw_B_tilde_smg] = 0.5*y[pba->index_bibw_B_tilde_smg]*y[pba->index_bibw_B_tilde_smg] -(1. - alpha_M + dH/H)*y[pba->index_bibw_B_tilde_smg] - (2.*alpha_M - D_kin*cs2 - (1. - 1./(1.+Delta_Mpl))*2.*dH/H + 3*(rho_smg + p_smg)/(1.+ Delta_Mpl)/H/H);

	return _SUCCESS_;

}

/**
 * Derivative of rho_smg
 *
 * @param pba                  Input: pointer to background structure
 * @param loga                    Input: scale factor
 * @param y                    Input: current vector of integrated quantities (with index_bi)
 * @param dy                   Output: current derivative of y w.r.t log(a)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_derivs_bw_rho_smg(
                        double loga,
                        double * y, /* vector with argument y[index_bi_bw] (must be already allocated with size pba->bi_bw_B_size) */
                        double * dy, /* vector with argument dy[index_bi_bw]
                                     (must be already allocated with
                                     size pba->bi_bw_B_size) */
                        void * parameters_and_workspace,
                        ErrorMsg error_message
                        ){

	// fill in workspace and interpolate all relavant quantities. Also find rho_tot and P_tot including all spieces except dark energy component
	/** - define local variables */
	struct background_parameters_and_workspace * pbpaw;
	struct background * pba;
	double w;
	int last_index; // necessary for calling array_interpolate(), but never used
	
	pbpaw = parameters_and_workspace;
  	pba =  pbpaw->pba;
	// interpolate w
	class_call(array_interpolate_spline(
                                        pba->stable_wext_lna_smg,
                                        pba->stable_wext_size_smg,
                                        pba->stable_wext_smg,
                                        pba->ddstable_wext_smg,
                                        1,
                                        loga,
                                        &last_index,
                                        &w,
                                        1,
                                        error_message),
               error_message,
               error_message);
	/* Derivative here are w.r.t. loga */
	dy[0] = -3. * y[0] * (1. + w);

	return _SUCCESS_;

}

/**
 * At some step during the backward integraton of the background equations,
 * this function extracts the quantities that we want to keep memory
 * of, and stores them in a row of the stable_params_smg table.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pba->error_message, but to a generic
 *   error_message passed in the list of arguments.
 *
 * @param loga                     Input: current value of log(a)
 * @param y                        Input: current vector of backward integrated quantities (with index_bibw)
 * @param dy                       Input: current derivative of y w.r.t log(a)
 * @param index_loga               Input: index of the log(a) value within loga_bw_table
 * @param parameters_and_workspace Input/output: fixed parameters (e.g. indices), workspace, background structure used to derive relevant quantities
 * @param error_message            Output: error message
 */

int background_sources_bw_smg(
                       double loga,
                       double * y,
                       double * dy,
                       int index_loga,
                       void * parameters_and_workspace,
                       ErrorMsg error_message
                       ) {

	struct background_parameters_and_workspace * pbpaw;
	struct background * pba;
	double * pvec_stable_params_smg;
	/** experimental bit*/
	double * pvecback;
	double * pvecback_B; // used to infer non-MG integrated quantities from interpolation and fed to background_functions
	double a, H, dH, rho_tot, P_tot;
	double Delta_Mpl, alpha_M, cs2;
	int pvecback_size;
	int pvec_stable_params_size; // size of output vector, controlled by input parameter return_format
	/** end experimental bit*/
	pbpaw = parameters_and_workspace;
	pba =  pbpaw->pba;
	pvec_stable_params_smg = pbpaw->pvec_stable_params_smg;
	/** experimental bit*/
	pvecback = pbpaw->pvecback;
	pvecback_B = pbpaw->pvecback_B;
	pvecback_size = pba->bg_size;
	pvec_stable_params_size = pba->num_stable_params;
	/** end experimental bit*/
	int last_index;
	double D_kin;

	// double Btilde = y[pba->index_bibw_B_tilde_smg];
	// double dBtilde = y[pba->index_bibw_dB_tilde_smg];

	// pba->stable_params_derived_smg[((pba->bt_bw_size-1) - index_loga)*pba->num_stable_params_derived + pba->index_derived_braiding_smg] = 2. * (1. - dBtilde/Btilde);
	// double alpha_B = pba->stable_params_derived_smg[((pba->bt_bw_size-1) - index_loga)*pba->num_stable_params_derived + pba->index_derived_braiding_smg];

	// To be used with 1st order ODE
	double alpha_B = y[pba->index_bibw_B_tilde_smg];
	pba->stable_params_derived_smg[((pba->bt_bw_size-1) - index_loga)*pba->num_stable_params_derived + pba->index_derived_braiding_smg] = alpha_B;

	class_call(array_interpolate_spline(
                                        pba->stable_params_lna_smg,
                                        pba->stable_params_size_smg,
                                        pba->stable_params_smg,
                                        pba->ddstable_params_smg,
                                        pba->num_stable_params,
                                        loga,
                                        &last_index,
                                        pvec_stable_params_smg,
                                        pba->num_stable_params,
                                        error_message),
               error_message,
               error_message);

	D_kin = pvec_stable_params_smg[pba->index_stable_Dkin_smg];
	pba->stable_params_derived_smg[((pba->bt_bw_size-1) - index_loga)*pba->num_stable_params_derived + pba->index_derived_kineticity_smg] = D_kin - 1.5 * alpha_B*alpha_B;

	/** Experimental bit to store alpha_B', works only with 1st order ODE for now */
	// Assign pvecback_B elements based on interpolation background table above
	if (pba->has_dcdm == _TRUE_ || pba->has_dr == _TRUE_) {
		class_call(array_interpolate_spline(
                                        pba->loga_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dloga2_table,
                                        pba->bg_size,
                                        loga,
                                        &last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
		/* dcdm */
		if (pba->has_dcdm == _TRUE_) {
			pvecback_B[pba->index_bi_rho_dcdm] = pvecback[pba->index_bg_rho_dcdm];
		}
		/* dr */
		if (pba->has_dr == _TRUE_) {
			pvecback_B[pba->index_bi_rho_dr] = pvecback[pba->index_bg_rho_dr];
		}
	}

	// call directly background_functions to speed up entire code and increase accuracy
	a = exp(loga);
	class_call(background_functions(pba, a, pvecback_B, normal_info, pvecback),
             pba->error_message,
             error_message);

	if (pba->expansion_model_smg == wext || pba->expansion_model_smg == rho_de) {
		class_call(interpolate_rho_smg_p_smg(pba, log(a), pba->loga_final_bw_integration, pvecback),
                 pba->error_message,
                 pba->error_message
      	);
		// Update quantities depending on rho_smg and p_smg
		H = sqrt(pvecback[pba->index_bg_rho_tot_wo_smg] + pvecback[pba->index_bg_rho_smg] - pba->K/a/a);
		dH = (-1.5*a*(pvecback[pba->index_bg_rho_tot_wo_smg] + pvecback[pba->index_bg_rho_smg]
													+ pvecback[pba->index_bg_p_tot_wo_smg] + pvecback[pba->index_bg_p_smg]) + pba->K/a)/a/H; // dH/dloga = 1/aH * dH/dtau
	}
	else {
		H = pvecback[pba->index_bg_H];
		dH = pvecback[pba->index_bg_H_prime]/a/H; // dH/dloga = 1/aH * dH/dtau
	}
	
	rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg]; // all matter exluding scalar field
	P_tot = pvecback[pba->index_bg_p_tot_wo_smg]; // all matter excluding scalar field

	Delta_Mpl = pvec_stable_params_smg[pba->index_stable_Delta_Mpl_smg];
	cs2 = pvec_stable_params_smg[pba->index_stable_cs2_smg];
	alpha_M = pvec_stable_params_smg[pba->index_stable_Mpl_running_smg];

	pba->stable_params_derived_smg[((pba->bt_bw_size-1) - index_loga)*pba->num_stable_params_derived + pba->index_derived_braiding_prime_smg] = (alpha_B - 2.) * (0.5 * alpha_B + alpha_M - dH/H) + D_kin * cs2 + 3. * (rho_tot + P_tot)/(H * H * (1.+Delta_Mpl));
	/** End experimental bit*/

	return _SUCCESS_;

}

/**
 * At some step during the backward integration of the background equation for rho_smg,
 * this function evaluates rho_smg at some specific times, and stores these vaules in stable_rho_smg.
 *
 * This is one of the few functions in the code which is passed to the generic_integrator() routine.
 * Since generic_integrator() should work with functions passed from various modules, the format of the arguments
 * is a bit special:
 * - fixed parameters and workspaces are passed through a generic pointer.
 *   generic_integrator() doesn't know the content of this pointer.
 * - the error management is a bit special: errors are not written as usual to pba->error_message, but to a generic
 *   error_message passed in the list of arguments.
 *
 * @param loga                     Input: current value of log(a)
 * @param y                        Input: current vector of backward integrated quantities (with index_bibw)
 * @param dy                       Input: current derivative of y w.r.t log(a)
 * @param index_loga               Input: index of the log(a) value within loga_bw_table
 * @param parameters_and_workspace Input/output: fixed parameters (e.g. indices), workspace, background structure used to derive relevant quantities
 * @param error_message            Output: error message
 */

int background_sources_bw_rho_smg(
                       double loga,
                       double * y,
                       double * dy,
                       int index_loga,
                       void * parameters_and_workspace,
                       ErrorMsg error_message
                       ) {

	struct background_parameters_and_workspace * pbpaw;
	struct background * pba;

	pbpaw = parameters_and_workspace;
	pba =  pbpaw->pba;

	pba->stable_rho_smg[(pba->bt_bw_rho_smg_size-1) - index_loga] = y[0];

	return _SUCCESS_;

}

/**
* Stability tests for smg.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input: vector of background quantities
* @param pvecback_integration Input: vector of background quantities to be integrated
* @return the error status
*/
int stability_tests_smg(
        								struct background *pba,
        								double * pvecback,
												double * pvecback_integration
											  ) {

	/* Horndeski stability tests
	* only if not overriden
	* and model is tuned!
	*/
	if ((pba->parameters_tuned_smg == _TRUE_) &&
	    (pba->skip_stability_tests_smg == _FALSE_)){

	  class_test(pba->min_D_smg <= -fabs(pba->D_safe_smg),
	      pba->error_message,
	      "Ghost instability for scalar field perturbations with minimum D=%g at a=%e\n", pba->min_D_smg, pba->a_min_D_smg);
	  class_test(pba->min_cs2_smg < -fabs(pba->cs2_safe_smg),
	      pba->error_message,
	      "Gradient instability for scalar field perturbations with minimum c_s^2=%g at a=%e\n", pba->min_cs2_smg, pba->a_min_cs2_smg);
	  class_test(pba->min_M2_smg < -fabs(pba->M2_safe_smg),
	      pba->error_message,
	      "Ghost instability for metric tensor perturbations with minimum M*^2=%g at a=%e\n", pba->min_M2_smg, pba->a_min_M2_smg);
	  class_test(pba->min_ct2_smg < -fabs(pba->ct2_safe_smg),
	      pba->error_message,
	      "Gradient instability for metric tensor perturbations with minimum c_t^2=%g at a=%e\n",pba->min_ct2_smg,pba->a_min_ct2_smg);

	 }

	return _SUCCESS_;
}


/**
* Numerical derivatives of the alphas.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input: vector of background quantities
* @param pvecback_derivs      Output: vector of derivatives
* @param i                    Input: counter for the time step
* @return the error status
*/
int derivatives_alphas_smg(
													 struct background *pba,
													 double * pvecback,
													 double * pvecback_derivs,
													 int i
													 ) {

	/* -> write in the table (overwrite the alpha time derivatives, which were set to nan in background_functions)
	 * direction of copy: add the corresponding indices to the coordinates
	 * thing to be copied: the address (&) to the element of pvecback corresponding to the index we want
	 * size: just a single double number
	 * -> repeat for all necessary quantities
	 */

	/* necessary for calling array_interpolate(), but never used */
	int last_index;

	// TODO_EB: note that the derivative is now calculated w.r.t. loga, while our _prime are w.r.t. tau
	double factor = pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H];
	double d_over_dtau;

	/*NOTE: here we compute the derivatives of quantities coputed in background_gravity_functions_smg during the integration.
	 * for quantities that depend on these derivatives (e.g. the gamma_i functions determining the effective mass)
	 * there is an additional loop at the end of background_solve
	 */

	// Kineticity'
	d_over_dtau = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
	copy_to_background_table_smg(pba, i, pba->index_bg_kineticity_prime_smg, d_over_dtau);

   //Braiding'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_braiding_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_braiding_prime_smg, d_over_dtau);

   //Planck mass run rate'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_mpl_running_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_mpl_running_prime_smg, d_over_dtau);

   //Tensor excess'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_tensor_excess_prime_smg, d_over_dtau);

  //Beyond horndeski'
	d_over_dtau = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
	copy_to_background_table_smg(pba, i, pba->index_bg_beyond_horndeski_prime_smg, d_over_dtau);

	if (pba->field_evolution_smg == _TRUE_) {
		//Braiding_over_phi'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_braiding_over_phi_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_braiding_over_phi_prime_smg, d_over_dtau);

		//Beyond_horndeski_over_phi'
		d_over_dtau = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_over_phi_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_beyond_horndeski_over_phi_prime_smg, d_over_dtau);
	}

   //H''
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_H_prime];
	 copy_to_background_table_smg(pba, i, pba->index_bg_H_prime_prime, d_over_dtau);

   // p_tot_wo_smg'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_p_tot_wo_prime_smg, d_over_dtau);

   // p_smg'
	 d_over_dtau = factor*pvecback_derivs[pba->index_bg_p_smg];
	 copy_to_background_table_smg(pba, i, pba->index_bg_p_prime_smg, d_over_dtau);

	// Planck's mass running
	// Only need to compute it if neither self consistent field evolution nor evolving M_pl in terms of alpha_M
	// check equation 3.3 of Bellini & Sawicki 2014

	if(pba->gravity_model_smg != stable_params && pba->field_evolution_smg == _FALSE_ && pba->M_pl_evolution_smg == _FALSE_){
		double alpha_M = pvecback_derivs[pba->index_bg_delta_M2_smg]/pvecback[pba->index_bg_M2_smg];
		copy_to_background_table_smg(pba, i, pba->index_bg_mpl_running_smg, alpha_M);
	}

	if(pba->background_verbose > 15 && fabs(1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]/factor)>1e-8)
printf("a = %g, (delta H')/H' = %g \n", pvecback[pba->index_bg_a], 1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]/factor);

	return _SUCCESS_;
}


/**
* Copy to the background table _smg quantities.
*
* @param pba                  Input: pointer to background structure
* @param row                  Input: table row
* @param column               Input: table column
* @param value                Input: value to copy
* @return the error status
*/
int copy_to_background_table_smg(
																 struct background *pba,
                                 int row,
                                 int column,
                                 double value
															   ) {

	/* needed for growing table */
	void * memcopy_result;

	memcopy_result = memcpy(pba->background_table + row*pba->bg_size + column,
	                        &value, 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + row*pba->bg_size + column,
	           pba->error_message, "cannot copy data back to pba->background_table");

  return _SUCCESS_;
}
