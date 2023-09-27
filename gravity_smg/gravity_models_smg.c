#include "gravity_models_smg.h"


/**
* Fix gravity_model properties.
* This includes two cases:
* - in scalar tensor theories this is the only place
*   were you specify the theory
* - in parametrizations this regulates the evolution
*   of the perturbations, while the background is fixed
*   in "_expansion_properties_"
*
* @param pfc                   Input: pointer to local structure
* @param pba                   Input/output: pointer to background structure
* @param string1               Input: name of the gravity model
* @param has_tuning_index_smg  Input: index of parameters_smg which should be tuned
* @param has_dxdy_guess_smg    Input: parameter variation range
* @param errmsg                Input: Error message
* @return the error status
*/
int gravity_models_gravity_properties_smg(
                                          struct file_content * pfc,
                                          struct background *pba,
                                          char string1[_ARGUMENT_LENGTH_MAX_],
                                          int has_tuning_index_smg,
                                          int has_dxdy_guess_smg,
                                          ErrorMsg errmsg
                                          ) {

  int flag1, flag2, flag3;
  int flag4, flag5, flag6, flag7;
  int int1;
  double * pointer1; 
  double * pointer2; 
  double * pointer3;
  char string2[_ARGUMENT_LENGTH_MAX_];
  char string3[_ARGUMENT_LENGTH_MAX_];
  /** Loop over the different models
  * flag1 checks whether external input file has been provided
  * flag2 keeps track of whether model has been identified
  */

  flag1=_FALSE_;
  flag2=_FALSE_;
  flag4=_FALSE_;
  flag5=_FALSE_;
  flag6=_FALSE_;
  flag7=_FALSE_;

  /** read the model and loop over models to set several flags and variables
  * field_evolution_smg: for self-consistent scalar tensor theories, need to evolve the background equations
  * M_pl_evolution_smg: for some parameterizations, need to integrate M_pl from alpha_M
  * Primary and secondary parameters: The tuning is alway in terms of a value in parameters_smg, therefore
  *  -> real models: "parameters_smg" to pba->parameters_smg
  *  -> parameterizations: "parameters_smg" to pba->parameters_2_smg
  *                        "expansion_smg" to pba->parameters_smg
  * NOTE: can change class_read_list_of_doubles_or_default <-> class_read_list_of_doubles
  * to make it mandatory or allow for default values
  */

  if (strcmp(string1,"propto_omega") == 0) {
     pba->gravity_model_smg = propto_omega;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"propto_scale") == 0) {
     pba->gravity_model_smg = propto_scale;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"constant_alphas") == 0) {
     pba->gravity_model_smg = constant_alphas;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 5;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"eft_alphas_power_law") == 0) {
     pba->gravity_model_smg = eft_alphas_power_law;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 8;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"eft_gammas_power_law") == 0) {
     pba->gravity_model_smg = eft_gammas_power_law;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 8;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

  if (strcmp(string1,"eft_gammas_exponential") == 0) {
     pba->gravity_model_smg = eft_gammas_exponential;
     pba->field_evolution_smg = _FALSE_;
     pba->M_pl_evolution_smg = _TRUE_;
     flag2=_TRUE_;
     pba->parameters_2_size_smg = 8;
     class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
   }

// MC parametrization based on external input file
  if (strcmp(string1,"external_alphas") == 0) {
	pba->gravity_model_smg = external_alphas;
	pba->field_evolution_smg = _FALSE_;
	pba->M_pl_evolution_smg = _TRUE_;
	flag2=_TRUE_;
	pba->parameters_2_size_smg = 1;
	class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);

  class_call(parser_read_string(pfc,
                                  "smg_file_name",
                                  &string2,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);

    if (flag1 == _TRUE_) {
        pba->has_smg_file = _TRUE_;
        class_read_string("smg_file_name",pba->smg_file_name);
    }else{
        class_stop(errmsg,"external_alphas parameterization requested. Please specify full path to file containing alpha functions");
    }

  /* for reading alpha_X functions */
  FILE * input_file;
  int row,status;
  double tmp1,tmp2,tmp3,tmp4,tmp5; // as many as the number of columns in input file
  int index_alpha;

  pba->ext_alphas_size_smg = 0;

  if (pba->has_smg_file == _TRUE_) {

    input_file = fopen(pba->smg_file_name,"r");
    class_test(input_file == NULL,
               errmsg,
               "Could not open file %s",pba->smg_file_name);

    /* Find size of table */
    for (row=0,status=5; status==5; row++){
      status = fscanf(input_file,"%lf %lf %lf %lf %lf",&tmp1,&tmp2,&tmp3,&tmp4,&tmp5);
    }
    rewind(input_file);
    pba->ext_alphas_size_smg = row-1;

    /** - define indices and allocate arrays for interpolation table*/

    class_alloc(pba->ext_alphas_lna_smg,pba->ext_alphas_size_smg*sizeof(double),errmsg);

    index_alpha = 0;
    class_define_index(pba->index_bg_ext_kineticity_smg,_TRUE_,index_alpha,1);
    class_define_index(pba->index_bg_ext_braiding_smg,_TRUE_,index_alpha,1);
    class_define_index(pba->index_bg_ext_mpl_running_smg,_TRUE_,index_alpha,1);
    class_define_index(pba->index_bg_ext_tensor_excess_smg,_TRUE_,index_alpha,1);
    pba->ext_num_alphas = index_alpha;
    class_alloc(pba->ext_alphas_smg,pba->ext_alphas_size_smg*pba->ext_num_alphas*sizeof(double),errmsg);
    class_alloc(pba->ext_ddalphas_smg,pba->ext_alphas_size_smg*pba->ext_num_alphas*sizeof(double),errmsg);

    /* - fill a table of ln(a) and alphas*/
    for (row=0; row<pba->ext_alphas_size_smg; row++){
      status = fscanf(input_file,"%lf %lf %lf %lf %lf",
                      &pba->ext_alphas_lna_smg[row],
                      &pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_kineticity_smg],
                      &pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_braiding_smg],
                      &pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_mpl_running_smg],
                      &pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_tensor_excess_smg]);
      // printf("%d: (lna, alpha_K, alpha_B, alpha_M, alpha_T) = (%e,%e,%e,%e,%e)\n",row,
      //         pba->ext_alphas_lna_smg[row],
      //         pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_kineticity_smg],
      //         pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_braiding_smg],
      //         pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_mpl_running_smg],
      //         pba->ext_alphas_smg[row*pba->ext_num_alphas + pba->index_bg_ext_tensor_excess_smg]);
    }
    fclose(input_file);


    /** - spline the table for later interpolation */

    class_call(array_spline_table_lines(
                                        pba->ext_alphas_lna_smg,
                                        pba->ext_alphas_size_smg,
                                        pba->ext_alphas_smg,
                                        pba->ext_num_alphas,
                                        pba->ext_ddalphas_smg,
                                        _SPLINE_NATURAL_,
                                        errmsg),
              errmsg,errmsg);
  }

  }// MC end of external_alphas

  // MC stable parametrization with external input file
  if (strcmp(string1,"stable_params") == 0) {
    pba->gravity_model_smg = stable_params;
    pba->field_evolution_smg = _FALSE_;
    pba->M_pl_evolution_smg = _FALSE_;
    flag2=_TRUE_;
    pba->parameters_2_size_smg = 1;
    // read alpha_B(z=0) to set initial conditions
    class_read_list_of_doubles("parameters_smg",pba->parameters_2_smg,pba->parameters_2_size_smg);
    // read z_gr_smg
    class_read_double("z_gr_smg",pba->z_gr_smg);

    class_call(parser_read_string(pfc,
                                    "smg_file_name",
                                    &string2,
                                    &flag1,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "lna_smg",
                                    &string2,
                                    &flag4,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "Delta_M2",
                                    &string2,
                                    &flag5,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "D_kin",
                                    &string2,
                                    &flag6,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "cs2",
                                    &string2,
                                    &flag7,
                                    errmsg),
                errmsg,
                errmsg);

    if (flag1 == _TRUE_ && flag4 == _FALSE_ && flag5 == _FALSE_ && flag6 == _FALSE_ && flag7 == _FALSE_) {
        pba->has_smg_file = _TRUE_;
        class_read_string("smg_file_name",pba->smg_file_name);
    } else if (flag1 == _FALSE_ && flag4 == _TRUE_ && flag5 == _TRUE_ && flag6 == _TRUE_ && flag7 == _TRUE_) {
        pba->has_smg_file = _FALSE_;
    } else {
        class_stop(errmsg,"stable_params parameterization requested. Either specify full path to file containing Delta(M*^2), D_kin and cs^2 functions and leave the parameters lna, DeltaM2, Dkin and cs2 unspecified, or leave filename unspecified and provide arrays for stable parametrization in input file.");
    }

    /* for reading Delta_M_pl^2, D_kin and cs^2 functions */
    FILE * input_file;
    int row,status;
    double tmp1,tmp2,tmp3,tmp4; // as many as the number of columns in input file
    int int1;
    int index_stable_params, index_stable_params_derived, index_stable_params_aux;

    pba->stable_params_size_smg = 0;

    index_stable_params = 0;
    class_define_index(pba->index_stable_Delta_Mpl_smg,_TRUE_,index_stable_params,1);
    class_define_index(pba->index_stable_Dkin_smg,_TRUE_,index_stable_params,1);
    class_define_index(pba->index_stable_cs2_smg,_TRUE_,index_stable_params,1);
    class_define_index(pba->index_stable_Mpl_running_smg,_TRUE_,index_stable_params,1);
    pba->num_stable_params = index_stable_params;
    // parameters derived from ODE integration
    index_stable_params_derived = 0;
    class_define_index(pba->index_derived_braiding_smg,_TRUE_,index_stable_params_derived,1);
    class_define_index(pba->index_derived_kineticity_smg,_TRUE_,index_stable_params_derived,1);
    pba->num_stable_params_derived = index_stable_params_derived;
    // auxiliary array to compute alpha_M from input Delta_M_pl^2
    index_stable_params_aux = 0;
    class_define_index(pba->index_aux_Delta_Mpl_smg,_TRUE_,index_stable_params_aux,1);
    class_define_index(pba->index_aux_dMpl_smg,_TRUE_,index_stable_params_aux,1);
    class_define_index(pba->index_aux_ddMpl_smg,_TRUE_,index_stable_params_aux,1);
    pba->num_stable_params_aux = index_stable_params_aux;

    if (pba->has_smg_file == _TRUE_) { // read from file

      input_file = fopen(pba->smg_file_name,"r");
      class_test(input_file == NULL,
                errmsg,
                "Could not open file %s",pba->smg_file_name);

      /* Find size of table */
      for (row=0,status=4; status==4; row++){
        status = fscanf(input_file,"%lf %lf %lf %lf",&tmp1,&tmp2,&tmp3,&tmp4);
      }
      rewind(input_file);
      pba->stable_params_size_smg = row-1;

      /** - define indices and allocate arrays for interpolation table*/
      class_alloc(pba->stable_params_lna_smg,pba->stable_params_size_smg*sizeof(double),errmsg);

    } else { // provide directly arrays for lna_smg, Delta(M*^2), D_kin and cs2

      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"lna_smg",&pba->stable_params_size_smg,&pba->stable_params_lna_smg,&flag1,errmsg),
             errmsg,
             errmsg);
      class_call(parser_read_list_of_doubles(pfc,"Delta_M2",&int1,&pointer1,&flag1,errmsg),
             errmsg,
             errmsg);
    // printf("gravity_models_gravity_properties_smg(): pba->stable_params_size_smg=%d \t int1=%d\n",pba->stable_params_size_smg,int1);
    // for (row=0; row<int1; row++){
    //   printf("lna_smg[%d] = %e \t Delta_M2[%d] = %e \n",
    //           row,pba->stable_params_lna_smg[row],
    //           row,pointer1[row]);
    // }
    // class_stop(errmsg,"stop here for now.\n");
      class_test(int1 != pba->stable_params_size_smg,errmsg,"Size of Delta_M2 array must match that of lna_smg.\n");
      class_call(parser_read_list_of_doubles(pfc,"D_kin",&int1,&pointer2,&flag1,errmsg),
             errmsg,
             errmsg);
      class_test(int1 != pba->stable_params_size_smg,errmsg,"Size of D_kin array must match that of lna_smg.\n");
      class_call(parser_read_list_of_doubles(pfc,"cs2",&int1,&pointer3,&flag1,errmsg),
             errmsg,
             errmsg);
      class_test(int1 != pba->stable_params_size_smg,errmsg,"Size of cs2 array must match that of lna_smg.\n");

    }

    /** - allocate various vectors */
    class_alloc(pba->stable_params_smg,pba->stable_params_size_smg*pba->num_stable_params*sizeof(double),errmsg);
    class_alloc(pba->ddstable_params_smg,pba->stable_params_size_smg*pba->num_stable_params*sizeof(double),errmsg);
    class_alloc(pba->stable_params_derived_smg,pba->stable_params_size_smg*pba->num_stable_params_derived*sizeof(double),errmsg);
    class_alloc(pba->ddstable_params_derived_smg,pba->stable_params_size_smg*pba->num_stable_params_derived*sizeof(double),errmsg);
    class_alloc(pba->stable_params_aux_smg,pba->stable_params_size_smg*pba->num_stable_params_aux*sizeof(double),errmsg);

    if (pba->has_smg_file == _TRUE_) {
      /* - fill a table of ln(a) and stable parameters*/
      for (row=0; row<pba->stable_params_size_smg; row++){
        status = fscanf(input_file,"%lf %lf %lf %lf",
                        &pba->stable_params_lna_smg[row],
                        &pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Delta_Mpl_smg],
                        &pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Dkin_smg],
                        &pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_cs2_smg]);
        // Initialize alpha_M and update value below as dMpl/Mpl
        pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Mpl_running_smg] = 0.;
      }
      fclose(input_file);
    } else {
      for (row=0; row<pba->stable_params_size_smg; row++){
        /* - fill a table of stable parameters*/
        pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Delta_Mpl_smg] = pointer1[row];
        pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Dkin_smg] = pointer2[row];
        pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_cs2_smg] = pointer3[row];
        // Initialize alpha_M and update value below as dMpl/Mpl
        pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Mpl_running_smg] = 0.;
      }
      free(pointer1);
      free(pointer2);
      free(pointer3);
    }

    // printf("stable_params_size_smg = %d\n",pba->stable_params_size_smg);
    // for (row=0; row<pba->stable_params_size_smg; row++){
    //   printf("lna_smg[%d] = %e \t Delta_M2[%d] = %e \t D_kin[%d] = %e \t cs2[%d] = %e \n",
    //           row,pba->stable_params_lna_smg[row],
    //           row,pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Delta_Mpl_smg],
    //           row,pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Dkin_smg],
    //           row,pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_cs2_smg]);
    // }

    // class_stop(errmsg,"stop here for now.\n");

    // Assign value to a_file_gr_smg
    pba->a_file_gr_smg = exp(pba->stable_params_lna_smg[0]);
    // Check that largest provided scale factor is >= 1. for stable numerical derivatives
    double a_max = 1.;
    class_test((exp(pba->stable_params_lna_smg[pba->stable_params_size_smg-1]) < a_max),
          errmsg,
          "maximum scale factor in input file is %f. For stable parametrization it must be >= %f for accurate numerical derivatives of smg parameters at late times", exp(pba->stable_params_lna_smg[pba->stable_params_size_smg-1]), a_max);

    /** - spline stable input parameters for later interpolation */
    class_call(array_spline_table_lines(pba->stable_params_lna_smg,
                                        pba->stable_params_size_smg,
                                        pba->stable_params_smg,
                                        pba->num_stable_params,
                                        pba->ddstable_params_smg,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
              pba->error_message,
              pba->error_message);

    /** - copy Delta_Mpl and ddMpl to auxiliary array */
    for (row=0; row<pba->stable_params_size_smg; row++){
      copy_to_aux_array_smg(pba,row,pba->index_aux_Delta_Mpl_smg,pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Delta_Mpl_smg]);
      copy_to_aux_array_smg(pba,row,pba->index_aux_ddMpl_smg,pba->ddstable_params_smg[row*pba->num_stable_params + pba->index_stable_Delta_Mpl_smg]);
    }
    /* Differentiate splined M_pl^2 to obtain alpha_m over the range of lna used by the ODE for alpha_b. */
    // First step gives dM_pl^2/dlna
    class_call(array_derive_spline_table_line_to_line(pba->stable_params_lna_smg,
                                                      pba->stable_params_size_smg,
                                                      pba->stable_params_aux_smg,
                                                      pba->num_stable_params_aux,
                                                      pba->index_aux_Delta_Mpl_smg,
                                                      pba->index_aux_ddMpl_smg,
                                                      pba->index_aux_dMpl_smg,
                                                      pba->error_message),
              pba->error_message,
              pba->error_message);

    // Update value of alpha_M in pba->stable_params_smg
    double Mpl, dMpl; // M_pl^2, dM_pl^2/dlna
    for (row=0; row<pba->stable_params_size_smg; row++){
                    Mpl = pba->stable_params_aux_smg[row*pba->num_stable_params_aux + pba->index_aux_Delta_Mpl_smg] + 1.;
                    dMpl = pba->stable_params_aux_smg[row*pba->num_stable_params_aux + pba->index_aux_dMpl_smg];
                    pba->stable_params_smg[row*pba->num_stable_params + pba->index_stable_Mpl_running_smg] = dMpl/Mpl;
    }

    /** - re-spline stable input parameters for later interpolation */
    class_call(array_spline_table_lines(pba->stable_params_lna_smg,
                                        pba->stable_params_size_smg,
                                        pba->stable_params_smg,
                                        pba->num_stable_params,
                                        pba->ddstable_params_smg,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
              pba->error_message,
              pba->error_message);

  }// MC end of stable_params

  if (strncmp("quintessence", string1, strlen("quintessence")) == 0) {
     // Check if gravity_model has quintessence as prefix.
     // Add here all variables common to quintessence.
     pba->is_quintessence_smg = _TRUE_;
     class_read_double("quintessence_w_safe_smg", pba->quintessence_w_safe_smg);
   }

  if (strcmp(string1,"quintessence_monomial") == 0) {
     pba->gravity_model_smg = quintessence_monomial;
     pba->field_evolution_smg = _TRUE_;
     pba->is_quintessence_smg = _TRUE_;
     flag2=_TRUE_;

     pba->parameters_size_smg = 4;
     class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

     /* Guess for the parameter variation range.
        *
        * For the initial parameter one can use:
      *
      * 	rho_smg = 1/2*a_ini^-2*phi_prime_ini^2 + V0*3*H0^2/h^2*phi_ini^N
      *
      * However, for the range of variation it is better to use
      *
      * 	Omega = rho_smg/(rho_smg + rho_m)
      *
      * => dOmega/dx_i = rho_m/(rho_smg+rho_m)^2 drho_smg/dx_i
      * => tuning_dxdy_guess_smg = (dOmega/dx_i)^{-1}
      * where we use rho_m ~ H_0^2
        *
        * drho_smg/dV0 = 10^-7*phi_ini^N
      */

     double N = pba->parameters_smg[0];
     double V0 = pba->parameters_smg[1];
     double phi_prime_ini_smg = pba->parameters_smg[2];
     double phi_ini_smg =  pba->parameters_smg[3];

     double P_ini = pow(phi_ini_smg, N);  // V=cte*P(phi)

     double phi_end_guess = fmax(phi_ini_smg,2); //guess the final value of the field

     // class_test( ((abs(N)<1) || (abs(N)>7)), errmsg, "Exponent out of range. N must be a interger in (1,7)-range" );

     if (has_tuning_index_smg == _FALSE_)
       pba->tuning_index_smg = 1; //use V0 for default tuning

     if (has_dxdy_guess_smg == _FALSE_) {

       if(pba->tuning_index_smg == 1) {
         // if(phi_ini_smg != 0){
         V0 = pba->Omega0_smg/pow(phi_end_guess,N);
         pba->tuning_dxdy_guess_smg = 1./pow(phi_end_guess,N);
         pba->parameters_smg[1] = V0;
         // }
         // else{
         //   V0 = pba->Omega0_smg/pow(1.e-40,N);
         //   pba->tuning_dxdy_guess_smg = 1./pow(1.e-40,N);
         //   pba->parameters_smg[1] = V0;
         //
         // }
       }

       if(pba->tuning_index_smg == 3){
          phi_ini_smg = pow(pba->Omega0_smg/V0, 1./N);
          pba->parameters_smg[3] = phi_ini_smg;
          pba->tuning_dxdy_guess_smg = phi_ini_smg/(pba->Omega0_smg)/N;
       }
     }
     //end of no has_dxdy_guess_smg
   }
  //end of quintessence_monomial

  if (strcmp(string1,"quintessence_tracker") == 0) {
     pba->gravity_model_smg = quintessence_tracker;
     pba->field_evolution_smg = _TRUE_;
     pba->is_quintessence_smg = _TRUE_;
     flag2=_TRUE_;

     pba->parameters_size_smg = 6;
     class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

     double K_ini = pba->parameters_smg[0];
     double P_ini =  pba->parameters_smg[1];
     double V0 = pba->parameters_smg[2];
     double n = pba->parameters_smg[3];
     double m = pba->parameters_smg[4];
     double lambda = pba->parameters_smg[5];

     /* Guess for the parameter variation range.
     *
     * For the initial parameter one can use:
     *  V = H0^2/h^2* V0 * phi^-n exp(lambda*phi^m)
     *
     *  minimum at phi0 = (n/(lambda*m))^(1/m)
     *  -> choose V0 ~ V(phi0)~Omega_smg H_0^2
     *
     * Initial conditions: see background.c
     *
     *  choose phi_ini so V = P_ini*sqrt(rho_rad)
     *  choose phi_prime_ini so K = K_ini*sqrt(rho_rad)
     *
     */

     double phi_0 = pow(n/lambda/m,1./m); /* minimum of the potential */
     double v_0_guess = (pow(phi_0,-n) * exp(lambda*pow(phi_0,m))); /*V/V0 at the minimum*/

     if (has_tuning_index_smg == _FALSE_)
       pba->tuning_index_smg = 2; //use V0 for default tuning

     if (has_dxdy_guess_smg == _FALSE_){
       if(pba->tuning_index_smg == 2){

         V0 = 3* pba->h * (pba->Omega0_smg)/v_0_guess;
         pba->tuning_dxdy_guess_smg = 3. * pba->h/ (v_0_guess); //*(1-pba->Omega0_smg) -> removed, lead to instability!
         pba->parameters_smg[2] = V0;
       }
     }
     //end of no has_dxdy_guess_smg
   }
  //end of tracker

  if (strcmp(string1,"alpha_attractor_canonical") == 0) {
     pba->gravity_model_smg = alpha_attractor_canonical;
     pba->field_evolution_smg = _TRUE_;
     pba->is_quintessence_smg = _TRUE_;
     flag2=_TRUE_;

     pba->parameters_size_smg = 6;
     class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

     class_call(parser_read_string(pfc,"log_10_param_alpha",&string2,&flag2,errmsg),
                errmsg,
                errmsg);

     if(flag2 == _TRUE_ && ((strstr(string2,"y") != NULL) || (strstr(string2,"Y") != NULL))){
       pba->parameters_smg[2] = pow(10, pba->parameters_smg[2]);
     }

     class_call(parser_read_string(pfc,"use_phi_no_f",&string2,&flag2,errmsg),
                errmsg,
                errmsg);

     if(flag2 == _TRUE_ && ((strstr(string2,"y") != NULL) || (strstr(string2,"Y") != NULL))){
       pba->parameters_smg[1] =  pba->parameters_smg[1]/sqrt(pba->parameters_smg[2]);
     }

     /* Guess for the parameter variation range. Copied from galileons.
        *
        * For the initial parameter one can use:
      *
      * However, for the range of variation it is better to use
      *
      * 	Omega = rho_smg/(rho_smg + rho_m)
      *
      * => dOmega/dx_i = rho_m/(rho_smg+rho_m)^2 drho_smg/dx_i
      * => tuning_dxdy_guess_smg = (dOmega/dx_i)^{-1}
      * where we use rho_m ~ H_0^2
        *
      */

     //Note: f = phi/sqrt(alpha)

     if (pba->Omega_smg_debug == 0.) {
       double phi_prime_ini = pba->parameters_smg[0];
       double f_ini = pba->parameters_smg[1];
       double alpha = pba->parameters_smg[2];
       double c = pba->parameters_smg[3];
       double p = pba->parameters_smg[4];
       double n = pba->parameters_smg[5];
       double x = tanh(f_ini/(sqrt(6)));
       double v = alpha* pow(x,p)/pow(1+x, 2*n); // v = V/c^2
       double rho_c = pow(pba->H0, 2);


       if (has_tuning_index_smg == _FALSE_)
         pba->tuning_index_smg = 3; //use V0 for default tuning

       if (has_dxdy_guess_smg == _FALSE_){
         if(pba->tuning_index_smg == 3){
             c = sqrt(pba->Omega0_smg * rho_c / v);
             pba->tuning_dxdy_guess_smg = pow((1 + 3*pba->Omega0_smg) * pba->H0,2)/(2*c*v);
             pba->parameters_smg[3] = c;
         }
       }//end of no has_dxdy_guess_smg
     }
     //end Omega_smg_debug
   }
  //end of  alpha_attractor_canonical

  if (strcmp(string1,"galileon") == 0) {
     pba->gravity_model_smg = galileon;
     pba->field_evolution_smg = _TRUE_;
     pba->parameters_size_smg = 7;
     flag2=_TRUE_;

     /* Galileon dynamics pulls towards the shift-symmetry attractor n = 0 with
      *
      * n/H0 = xi(c2 - 6c3 xi + 18c4 xi^2 + 5c5 xi^4)
      *
      * and xi = \dot\phi H /H0^2
      *
      * If attractor_ic_smg => n=0 is set in background_initial_conditions
     */


     /* Guess for the parameter variation range. For the initial parameter one can use
      *
      * 	rho_smg*H^2/H_0^4 = c2 xi^2/6 - 2 c3 xi^3 + 15/2 c4 xi^4 + 7/3 c5 xi^5
      *
      * (Barreira+ '14 2.22 for z\neq 0), which equals Omega_smg at z=0 only if tuned.
      *
      * There are three submodels (taken on tracker):
      * 	1) rogue mode: user sets all (if Omega_smg_debug or NO attractor_ic_smg)
      * 	2) cubic attractor: c3, xi set for attractor
      * 	3) quartic/quintic attractor: c3, c4 set xi for attractor
      */


     // read submodel: remember flag2 used for test over gravity models
     class_call(parser_read_string(pfc,"gravity_submodel",&string2,&flag3,errmsg),
            errmsg,
            errmsg);

     /*1) base galileon, user specifies everything!  */
     if (flag3==_FALSE_) {

       class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

     }
     else {
       //a submodel is given

       /*  temporary allocation, will be rewritten latter
        *  order is xi, phi, c2, c3, c4, c5  */
       class_alloc(pba->parameters_smg, sizeof(double*)*7,pba->error_message);
       double * input_params_gal; //dummy allocation vector

       double xi, c2, c3, c4, c5;
       double phi0 = 0, c1 = 0;

       /*2) cubic Galileon in the attractor
         * Omega = -c2^3/(6^3 c3^2) and xi = c2/(6c3)
         */
       if (strcmp(string2,"cubic") == 0) {

         c2 = -1.;
         c3 = -sqrt(-pow(c2/6.,3)/pba->Omega0_smg);
         xi = c2/6./c3;
         c4 = 0;
         c5 = 0;

       }//end of cubic
        /* 3) quartic Galileon on the attractor
         */
       else if (strcmp(string2,"quartic") == 0) {

         class_read_list_of_doubles("parameters_smg",input_params_gal,1);
         xi = input_params_gal[0];
         c2 = -1;
         c3 = (4.*pba->Omega0_smg + c2*pow(xi,2))/(2.*pow(xi,3));
         c4 = (6.*pba->Omega0_smg + c2*pow(xi,2))/(9.*pow(xi,4));
         c5 = 0;

       }/* 4) quartic Galileon on the attractor
         */
       else if (strcmp(string2,"quintic") == 0) {//Quintic case

         class_read_list_of_doubles("parameters_smg",input_params_gal,2);
         xi = input_params_gal[0];
         c2 = -1;
         c3 = input_params_gal[1];
         c4 = -(10*pba->Omega0_smg + 3*c2*pow(xi,2) - 8*c3*pow(xi,3))/(9.*pow(xi,4));
         c5 = (4*pba->Omega0_smg + pow(xi,2)*(c2 - 2*c3*xi))/pow(xi,5);

       }//end of quintic
       else {
             class_test(flag3 == _TRUE_,
        errmsg,
        "Galileon: you specified a gravity_submodel that could not be identified. \n Options are: cubic, quartic, quintic");
       };

       /* Set parameters for submodels */
       pba->parameters_smg[0] = xi;
       pba->parameters_smg[1] = c1;
       pba->parameters_smg[2] = c2;
       pba->parameters_smg[3] = c3;
       pba->parameters_smg[4] = c4;
       pba->parameters_smg[5] = c5;
       pba->parameters_smg[6] = phi0;

     }
     //end of submodels

     /* default tuning index is 3 */
     if (has_tuning_index_smg == _FALSE_){
       pba->tuning_index_smg = 3; //use c3 for default tuning
       //Use the tracker condition for the cubic to define xi_0, in case xi is used as an IC.
         pba->tuning_dxdy_guess_smg = 2./pow(pba->parameters_smg[2]/6./pba->parameters_smg[3],3);
       //pba->tuning_dxdy_guess_smg = 2./pow(pba->parameters_smg[0],3); // d(c3)/d(Omega_smg) = 2/xi^3 and xi = c2/6./c3;
     }
     class_test(has_dxdy_guess_smg == _TRUE_ && has_tuning_index_smg == _FALSE_,
                errmsg,
                "Galileon: you gave dxdy_guess_smg but no tuning_index_smg. You need to give both if you want to tune the model yourself");

     class_call(parser_read_string(pfc,"attractor_ic_smg",&string3,&flag3,errmsg),
          errmsg,
          errmsg);

     if(flag3 == _TRUE_ && ((strstr(string3,"y") != NULL) || (strstr(string3,"Y") != NULL))){
       pba->attractor_ic_smg = _TRUE_;
     }
     else{
       pba->attractor_ic_smg = _FALSE_;
     }
   }
  //end of Galileon

  if (strcmp(string1,"brans dicke") == 0 || strcmp(string1,"Brans Dicke") == 0 || strcmp(string1,"brans_dicke") == 0) {
    pba->gravity_model_smg = brans_dicke;
    pba->field_evolution_smg = _TRUE_;
    flag2=_TRUE_;

    pba->parameters_size_smg = 4;
    class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
    pba->parameters_smg[0] = 2*pba->Omega0_smg;
    pba->tuning_dxdy_guess_smg = 0.5;
    pba->tuning_index_2_smg = 2;

    /** Read the desired Planck mass and check that the necessary information is provided.
    *  if needed re-assign shooting parameter for the Planck mass
    */
    class_call(parser_read_string(pfc, "M_pl_tuning_smg", &string3, &flag3, errmsg),
         errmsg,
         errmsg);

    if (flag3 == _TRUE_){
      if((strstr(string3,"y") != NULL) || (strstr(string3,"Y") != NULL)){

        pba->M_pl_tuning_smg = _TRUE_;

        class_read_double("M_pl_today_smg",pba->M_pl_today_smg);

        class_call(parser_read_string(pfc,"normalize_G_NR",
         		&string3,
         		&flag3,
         		errmsg),
         	errmsg,
         	errmsg);

        if (flag3 == _TRUE_){
           if((strstr(string3,"y") != NULL) || (strstr(string3,"Y") != NULL)){
             double omega_BD = pba->parameters_smg[1];
             pba->M_pl_today_smg = (4.+2.*omega_BD)/(3.+2.*omega_BD);
          }
        }

        class_read_double("param_shoot_M_pl_smg",pba->parameters_smg[pba->tuning_index_2_smg]);
           // printf("updating param = %e to tune M_pl \n",pba->parameters_smg[pba->tuning_index_2_smg]);
       }
     }
   }

  if (strcmp(string1,"nkgb") == 0 || strcmp(string1,"n-kgb") == 0 || strcmp(string1,"N-KGB") == 0 || strcmp(string1,"nKGB") == 0) {
     // This is self-accelerating KGB with K=-X and G(X)=1/n g^(2n-1)/2 * X^n
     pba->gravity_model_smg = nkgb;
     pba->field_evolution_smg = _TRUE_;
     if (has_tuning_index_smg == _FALSE_ && pba->Omega_smg_debug == 0){
       pba->tuning_index_smg = 0; //use g for default tuning
     }
     class_test(has_dxdy_guess_smg == _TRUE_ && has_tuning_index_smg == _FALSE_,
                errmsg,
                "nKGB: you gave dxdy_guess_smg but no tuning_index_smg. You need to give both if you want to tune the model yourself");
     if(has_dxdy_guess_smg == _FALSE_){
       pba->tuning_dxdy_guess_smg = -0.5;
     }
     flag2=_TRUE_;

     pba->parameters_size_smg = 3; // g, n, xi0 == rho_DE_0(shift charge)/rho_DE_0(total)
     class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
     class_test(pba->parameters_smg[1]<=0.5,errmsg,"In n-KGB G(X)=X^n n>1/2 for acceleration. Note that limit n->1/2 is singular and models become badly behaved for n<0.7");
     class_test(pba->parameters_smg[2]>=1.,errmsg,"In n-KGB, Rshift0<1 for positive energy density today.");
     class_test(pba->parameters_smg[2]<0.,errmsg,"In n-KGB, Rshift0>=0, or ICs for background can't be set.");
   }
  //end of nKGB

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify gravity_theory value, check that it is one of 'propto_omega', 'propto_scale', 'constant_alphas', 'eft_alphas_power_law', 'eft_gammas_power_law', 'eft_gammas_exponential', 'external_alphas', 'stable_params', 'brans_dicke', 'galileon', 'nKGB', 'quintessence_monomial', 'quintessence_tracker', 'alpha_attractor_canonical' ...");

  return _SUCCESS_;
}

/**
* Fix expansion history properties.
*
* @param pfc                   Input: pointer to local structure
* @param pba                   Input/output: pointer to background structure
* @param string1               Input: name of the gravity model
* @param errmsg                Input: Error message
* @return the error status
*/
int gravity_models_expansion_properties_smg(
                                            struct file_content * pfc,
                                            struct background *pba,
                                            char string1[_ARGUMENT_LENGTH_MAX_],
                                            ErrorMsg errmsg
                                            ) {

  /* Define local variables */
  int flag1 = _FALSE_;
  int flag2 = _FALSE_;
  int flag3 = _FALSE_;
  int flag4 = _FALSE_;
  char string2[_ARGUMENT_LENGTH_MAX_];
  int n;

  //possible expansion histories. Can make tests, etc...
  if (strcmp(string1,"lcdm") == 0) {
    pba->expansion_model_smg = lcdm;
    flag2=_TRUE_;
    pba->parameters_size_smg = 1;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);
  }
  //accept different names

  if (strcmp(string1,"wowa") == 0 || strcmp(string1,"w0wa") == 0 || strcmp(string1,"cpl") == 0 ) {
    pba->expansion_model_smg = wowa;
    flag2=_TRUE_;
    pba->parameters_size_smg = 3;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);
  }

  if (strcmp(string1,"wowa_w") == 0 || strcmp(string1,"w0wa_w") == 0 || strcmp(string1,"cpl_w") == 0 ) {
    pba->expansion_model_smg = wowa_w;
    flag2=_TRUE_;
    pba->parameters_size_smg = 3;
    pba->rho_evolution_smg=_TRUE_;
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);
  }

  if (strcmp(string1,"wede") == 0) {
    //ILSWEDE
    pba->expansion_model_smg = wede;
    flag2=_TRUE_;
    pba->parameters_size_smg = 3;
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);
    // 	//optimize the guessing BUG: eventually leads to problem in the MCMC, perhaps the guess is too good?
    // 	if(pba->tuning_index_smg == 0){
    // 	  pba->parameters_smg[0] = pba->Omega0_smg;
    // 	}
  }

  if (strcmp(string1,"wext") == 0) {
    // wext only works with stable parametrisation
    class_test(pba->gravity_model_smg != stable_params,
               pba->error_message,
               "wext background parametrisation only works for Horndeski gravity and gravity_model = stable_params")

    pba->expansion_model_smg = wext;
    flag2=_TRUE_;
    pba->parameters_size_smg = 1;
    pba->rho_evolution_smg=_FALSE_; // set here to FALSE because rho_smg will be integrated backward in time and before integration of other background quantities
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);

    class_call(parser_read_string(pfc,
                                  "expansion_file_name",
                                  &string2,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);

    class_call(parser_read_string(pfc,
                                    "lna_de",
                                    &string2,
                                    &flag3,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "de_evo",
                                    &string2,
                                    &flag4,
                                    errmsg),
                errmsg,
                errmsg);

    if (flag1 == _TRUE_ && flag3 == _FALSE_ && flag4 == _FALSE_) {
        pba->has_expansion_file = _TRUE_;
        class_read_string("expansion_file_name",pba->expansion_file_name);
    } else if (flag1 == _FALSE_ && flag3 == _TRUE_ && flag4 == _TRUE_) {
        pba->has_expansion_file = _FALSE_;
    } else {
        class_stop(errmsg,"wext expansion model requested. Either specify full path to file containing w function and leave the parameters lna_de and de_evo unspecified, or leave filename unspecified and provide arrays for DE background evolution in input file.");
    }

    /* for reading w from file */
    FILE * input_file;
    int row,status;
    int int1;
    double tmp1,tmp2; // as many as the number of columns in input file

    pba->stable_wext_size_smg = 0;

    if (pba->has_expansion_file == _TRUE_) { // read from file

      input_file = fopen(pba->expansion_file_name,"r");
      class_test(input_file == NULL,
                errmsg,
                "Could not open file %s",pba->expansion_file_name);

      /* Find size of table */
      for (row=0,status=2; status==2; row++){
        status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      pba->stable_wext_size_smg = row-1;

      /** - allocate various vectors */
      class_alloc(pba->stable_wext_lna_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
      class_alloc(pba->stable_wext_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
      /* - fill a table of ln(a) and stable parameters*/
      for (row=0; row<pba->stable_wext_size_smg; row++){
        status = fscanf(input_file,"%lf %lf",
                        &pba->stable_wext_lna_smg[row],
                        &pba->stable_wext_smg[row]);
      }
      fclose(input_file);      

    } else { // provide directly arrays for lna_de and DE eos
      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"lna_de",&pba->stable_wext_size_smg,&pba->stable_wext_lna_smg,&flag1,errmsg),
             errmsg,
             errmsg);
      class_call(parser_read_list_of_doubles(pfc,"de_evo",&int1,&pba->stable_wext_smg,&flag1,errmsg),
             errmsg,
             errmsg);
      class_test(int1 != pba->stable_wext_size_smg,errmsg,"Size of de_evo array must match that of lna_de.\n");
    }

    class_alloc(pba->ddstable_wext_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
    class_alloc(pba->stable_rho_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
    class_alloc(pba->ddstable_rho_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);

    // Assign value to a_file_lcdm_smg
    pba->a_file_lcdm_smg = exp(pba->stable_wext_lna_smg[0]);
    // Check that largest provided scale factor is >= 1. for stable numerical derivatives
    double a_max = 1.;
    class_test((exp(pba->stable_wext_lna_smg[pba->stable_wext_size_smg-1]) < a_max),
          errmsg,
          "maximum scale factor in input file is %f. For stable parametrization it must be >= %f for accurate numerical derivatives of smg parameters at late times", exp(pba->stable_wext_lna_smg[pba->stable_wext_size_smg-1]), a_max);

      /** - spline stable input parameters for later interpolation */
    class_call(array_spline_table_lines(pba->stable_wext_lna_smg,
                                        pba->stable_wext_size_smg,
                                        pba->stable_wext_smg,
                                        1,
                                        pba->ddstable_wext_smg,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
              pba->error_message,
              pba->error_message);

    // printf("stable_wext_size_smg = %d\n",pba->stable_wext_size_smg);
    // for (row=0; row<pba->stable_wext_size_smg; row++){
    //   printf("lna_de[%d] = %e \t wext[%d] = %e\n",
    //           row,pba->stable_wext_lna_smg[row],
    //           row,pba->stable_wext_smg[row]);
    // }

    // class_stop(errmsg,"stop here for now.\n");

  } // MC end of wext

  if (strcmp(string1,"rho_de") == 0) {
    // rho_de only works with stable parametrisation
    class_test(pba->gravity_model_smg != stable_params,
               pba->error_message,
               "rho_de background parametrisation only works for Horndeski gravity and gravity_model = stable_params")

    pba->expansion_model_smg = rho_de;
    flag2=_TRUE_;
    pba->parameters_size_smg = 1;
    pba->rho_evolution_smg=_FALSE_;
    class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);

    class_call(parser_read_string(pfc,
                                  "expansion_file_name",
                                  &string2,
                                  &flag1,
                                  errmsg),
               errmsg,
               errmsg);

    class_call(parser_read_string(pfc,
                                    "lna_de",
                                    &string2,
                                    &flag3,
                                    errmsg),
                errmsg,
                errmsg);
    class_call(parser_read_string(pfc,
                                    "de_evo",
                                    &string2,
                                    &flag4,
                                    errmsg),
                errmsg,
                errmsg);

    if (flag1 == _TRUE_ && flag3 == _FALSE_ && flag4 == _FALSE_) {
        pba->has_expansion_file = _TRUE_;
        class_read_string("expansion_file_name",pba->expansion_file_name);
    } else if (flag1 == _FALSE_ && flag3 == _TRUE_ && flag4 == _TRUE_) {
        pba->has_expansion_file = _FALSE_;
    } else {
        class_stop(errmsg,"rho_de expansion model requested. Either specify full path to file containing normalised rho_de(lna)/rho_de(0) function and leave the parameters lna_de and de_evo unspecified, or leave filename unspecified and provide arrays for DE background evolution in input file.");
    }

    /* for reading w from file */
    FILE * input_file;
    int row,status;
    int int1;
    double tmp1,tmp2; // as many as the number of columns in input file

    pba->stable_wext_size_smg = 0; // for simplicity use same variables and vectors of wext parametrization whenever it makes sense

    if (pba->has_expansion_file == _TRUE_) { // read from file

      input_file = fopen(pba->expansion_file_name,"r");
      class_test(input_file == NULL,
                errmsg,
                "Could not open file %s",pba->expansion_file_name);

      /* Find size of table */
      for (row=0,status=2; status==2; row++){
        status = fscanf(input_file,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(input_file);
      pba->stable_wext_size_smg = row-1;

      /** - allocate various vectors */
      class_alloc(pba->stable_wext_lna_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
      class_alloc(pba->stable_rho_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);
      // class_alloc(pba->ddstable_rho_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);

      /* - fill in vectors for ln(a) and DE density normilised at z=0, i.e. rho_de(lna)/rho_de(0)*/
      for (row=0; row<pba->stable_wext_size_smg; row++){
        status = fscanf(input_file,"%lf %lf",
                        &pba->stable_wext_lna_smg[row],
                        &pba->stable_rho_smg[row]);
      }
      fclose(input_file);      

    } else { // provide directly arrays for lna_de and rho_de
      /* Read */
      class_call(parser_read_list_of_doubles(pfc,"lna_de",&pba->stable_wext_size_smg,&pba->stable_wext_lna_smg,&flag1,errmsg),
             errmsg,
             errmsg);
      class_call(parser_read_list_of_doubles(pfc,"de_evo",&int1,&pba->stable_rho_smg,&flag1,errmsg),
             errmsg,
             errmsg);
      class_test(int1 != pba->stable_wext_size_smg,errmsg,"Size of de_evo array must match that of lna_de.\n");
    }
    /** - allocate vector for second derivatives */
    class_alloc(pba->ddstable_rho_smg,pba->stable_wext_size_smg*sizeof(double),errmsg);

    // Assign value to a_file_lcdm_smg
    pba->a_file_lcdm_smg = exp(pba->stable_wext_lna_smg[0]);
    // Check that largest provided scale factor is >= 1. for stable numerical derivatives
    double a_max = 1.;
    class_test((exp(pba->stable_wext_lna_smg[pba->stable_wext_size_smg-1]) < a_max),
          errmsg,
          "maximum scale factor in input file is %f. For stable parametrization it must be >= %f for accurate numerical derivatives of smg parameters at late times", exp(pba->stable_wext_lna_smg[pba->stable_wext_size_smg-1]), a_max);

      /** - spline stable input parameters for later interpolation */
    class_call(array_spline_table_lines(pba->stable_wext_lna_smg,
                                        pba->stable_wext_size_smg,
                                        pba->stable_rho_smg,
                                        1,
                                        pba->ddstable_rho_smg,
                                        _SPLINE_EST_DERIV_,
                                        pba->error_message),
              pba->error_message,
              pba->error_message);

    // printf("stable_wext_size_smg = %d\n",pba->stable_wext_size_smg);
    // for (row=0; row<pba->stable_wext_size_smg; row++){
    //   printf("lna_de[%d] = %e \t rho_de[%d] = %e\n",
    //           row,pba->stable_wext_lna_smg[row],
    //           row,pba->stable_rho_smg[row]);
    // }

    // class_stop(errmsg,"stop here for now.\n");

  } // MC end of rho_de

  class_test(flag2==_FALSE_,
             errmsg,
             "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede, wext, rho_de ...");

  return _SUCCESS_;
}


/**
* Overwrite the default values of the Gs depending on the model.
* This is the only place where the shape of the Gs is specified.
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback_B    Input: vector containing all {B} quantities
* @param pgf           Input: pointer to G_functions_and_derivs structure
* @return the error status
*/
int gravity_models_get_Gs_smg(
                              struct background *pba,
                              double a,
                              double * pvecback_B,
                              struct G_functions_and_derivs *pgf
                              ) {

  double phi = pvecback_B[pba->index_bi_phi_smg];
  double phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];
  double X = 0.5*pow(phi_prime/a,2);

  /* Overwrite the Horneski functions and their partial derivatives depending on the model */
  if (pba->gravity_model_smg == quintessence_monomial) {

    double N = pba->parameters_smg[0];
    double V0 = pba->parameters_smg[1];

    pgf->G2 = X - V0*pow(pba->H0/pba->h,2.)*pow(phi,N); // V written as in arXiv:1406.2301 in CLASS units
    pgf->G2_X = 1.;
    pgf->G2_phi = -N*V0*pow(pba->H0/pba->h,2.)*pow(phi,N-1.);
  }

  else if (pba->gravity_model_smg == quintessence_tracker) {
      // V written as in arXiv:1406.2301 in CLASS units

    double V0 = pba->parameters_smg[2];
    double n = pba->parameters_smg[3];
    double m = pba->parameters_smg[4];
    double lambda = pba->parameters_smg[5];
    double V, V_phi;

    V = pow(pba->H0/pba->h,2)* V0 * pow(phi, -n) * exp(lambda*pow(phi, m));
    V_phi = (lambda * m * pow(phi, m) - n) /phi * V;

    pgf->G2 = X - V;
    pgf->G2_X = 1.;
    pgf->G2_phi = -V_phi;
  }

  else if (pba->gravity_model_smg == alpha_attractor_canonical) {
      // V written as in eq. 12 of arXiv:1505.00815 in CLASS units with canonical field.

    double alpha = pba->parameters_smg[2];
    double c2 = pow(pba->parameters_smg[3], 2);
    double p = pba->parameters_smg[4];
    double n = pba->parameters_smg[5];
    double V, V_phi;
    double x = tanh(phi/(sqrt(6*alpha)));
    double y = sinh(sqrt(2/(3*alpha))*phi);
    V = alpha* c2 * pow(x,p)/pow(1+x, 2*n);

    V_phi = sqrt(2*alpha/3)*c2* pow(x,p)* (p +(p-2*n)*x)/(y * pow(1+x,2*n +1));

    pgf->G2 = X - V;
    pgf->G2_X = 1.;
    pgf->G2_phi = -V_phi;

   class_test((phi < 0. ) && ( phi_prime > 0 )
           && (pba->parameters_tuned_smg == _TRUE_)
           && (pba->skip_stability_tests_smg == _FALSE_),
           pba->error_message,
           "The model has started oscillating with first minimum at a = %e. Since <w> = 0 yields matter, it cannot make the universe expand accelerately.", a);
  }

  else if(pba->gravity_model_smg == galileon){
    //TODO: change name with ccc_galileon

    double M3 = pow(pba->H0,2);

    double Lambda1 = pba->parameters_smg[1]*M3;
    double Lambda2 = pba->parameters_smg[2];
    double Lambda3 = pba->parameters_smg[3]/M3;
    double Lambda4 = pba->parameters_smg[4]*pow(M3,-2);
    double Lambda5 = pba->parameters_smg[5]*pow(M3,-3);

    pgf->G2 = Lambda2*X - 0.5*Lambda1*phi;
    pgf->G2_X = Lambda2;
    pgf->G2_phi = - 0.5*Lambda1;
    /*  pgf->G_3 = -2*Lambda3*X */
    pgf->G3_X = -2.*Lambda3;
    /* pgf->G_4 = 1/2 + Lambda4 X^2 */
    pgf->DG4 = Lambda4*pow(X,2);
    pgf->G4 = 1/2. + pgf->DG4;
    pgf->G4_X = 2.*Lambda4*X;
    pgf->G4_XX = 2.*Lambda4;
    /* pgf->G_5 = Lambda5*pow(X,2) */
    pgf->G5_X = 2.*Lambda5*X;
    pgf->G5_XX = 2.*Lambda5;
  }

  else if(pba->gravity_model_smg == brans_dicke){

    /* Brans-Dicke can't cause acceleration:
     * - V is a constant potential, basically a cosmological constant
     * - omega is the Brans-Dicke parameter in charge of the fifth force
     */
    double V = 3.*pba->parameters_smg[0]*pow(pba->H0,2);
    double omega = pba->parameters_smg[1];

    pgf->G2 = -V + omega*X/phi;
    pgf->G2_X = omega/phi;
    pgf->G2_Xphi = -omega/pow(phi,2);
    pgf->G2_phi = -omega*X/pow(phi,2);

    pgf->DG4 = (phi-1.)/2.;
    pgf->G4 = phi/2.;
    pgf->G4_phi = 1./2.;

  }

  else if(pba->gravity_model_smg == nkgb){

    /* Action is

      -X + 1/n * g^[(2n-1)/2] Lambda (X/Lambda^4)^n box(phi)

      g was picked like this so that it approx. remains g*Omega_smg0 ~ O(1) for all n
      Note that for n=1/4 the Lambda mass scales cancels out, so we set it to 1.
    */

    double g = pba->parameters_smg[0];
    double npow = pba->parameters_smg[1];
    double ngpow = copysign(1.,g)*pow(fabs(g),(2.*npow-1.)/2.)/npow;
    double H0=pba->H0;

    pgf->G2    = -X;
    pgf->G2_X  = -1.;

    // pgf->G3 = 1/n g^[(2n-1)/2] Lambda (X/Lambda^4)^n

    pgf->G3_X = npow*ngpow*pow(X,npow-1)/pow(H0,2*npow);
    pgf->G3_XX = npow*(npow-1.)*ngpow*pow(X,npow-2)/pow(H0,2*npow);

  }

  return _SUCCESS_;
}


/**
* Get the background parametrizations, either parametrizing
* rho_smg and p_smg or w_smg.
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_B    Input: vector containing all {B} quantities
* @return the error status
*/
int gravity_models_get_back_par_smg(
                                    struct background *pba,
                                    double a,
                                    double * pvecback,
                                    double * pvecback_B
                                    ) {


  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

  if (pba->expansion_model_smg == lcdm){
    
    double Omega_const_smg = pba->parameters_smg[0];

    pvecback[pba->index_bg_rho_smg] = Omega_const_smg*pow(pba->H0,2);
    pvecback[pba->index_bg_p_smg] = -Omega_const_smg*pow(pba->H0,2);
  }

  else if (pba->expansion_model_smg == wowa){

    double Omega_const_smg = pba->parameters_smg[0];
    double w0 = pba->parameters_smg[1];
    double wa = pba->parameters_smg[2];

    pvecback[pba->index_bg_rho_smg] = Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1. + w0 + wa)) * exp(3.*wa*(a-1.));
    pvecback[pba->index_bg_p_smg] = (w0+(1-a)*wa) * Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1.+w0+wa)) * exp(3.*wa*(a-1.));
  }

  else if (pba->expansion_model_smg == wowa_w){

    //double Omega_const_smg = pba->parameters_smg[0];
    double w0 = pba->parameters_smg[1];
    double wa = pba->parameters_smg[2];

    pvecback[pba->index_bg_w_smg] = w0+(1-a)*wa;

    //DT: All commented out parts are now before th definition of pvecback_integration[pba->index_bi_rho_smg] (L2784)

    //if (pba->initial_conditions_set_smg == _FALSE_) {
    // Here we provide wi wf from w= (1-a)*wi+a*wf.
    // This is useful to set the initial conditions  for the energy density.
    // The value inferred here is just a guess, since then the shooting will modify it.
    //  double wi = w0+wa;
    //  double wf = w0;

    //  pvecback[pba->index_bg_rho_smg] = Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1. + wi)) * exp(3.*(wi-wf)*(a-1.));
    //}
    //else {
    pvecback[pba->index_bg_rho_smg] = pvecback_B[pba->index_bi_rho_smg];
    //}
    pvecback[pba->index_bg_p_smg] = pvecback[pba->index_bg_w_smg] * pvecback[pba->index_bg_rho_smg];

  }

  else if (pba->expansion_model_smg == wede){

    //Doran-Robbers model astro-ph/0601544
    //as implemented in Pettorino et al. 1301.5279
    //TODO: check these expressions, they probably assume the standard evolution/friedmann eqs, etc...
    //TODO: rewrite the expressions integrating the equation of state

    double Om0 = pba->parameters_smg[0];
    double w0 = pba->parameters_smg[1];
    double Ome = pba->parameters_smg[2] + 1e-10;

    //NOTE: I've regularized the expression adding a tiny Omega_e
    double Om = ((Om0 - Ome*(1.-pow(a,-3.*w0)))/(Om0 + (1.-Om0)*pow(a,3*w0)) + Ome*(1.-pow(a,-3*w0)));
    double dOm_da = (3*pow(a,-1 - 3*w0)*(-1 + Om0)*(-2*pow(a,3*w0)*(-1 + Om0)*Ome + Om0*Ome + pow(a,6*w0)*(Om0 - 2*Ome + Om0*Ome))*w0)/pow(-(pow(a,3*w0)*(-1 + Om0)) + Om0,2); //from Mathematica
    //I took a_eq = a*rho_r/rho_m, with rho_r = 3*p_tot_wo_smg
    double a_eq = 3.*a*p_tot/(pvecback[pba->index_bg_rho_b]+pvecback[pba->index_bg_rho_cdm]); //tested!
    double w = a_eq/(3.*(a+a_eq)) -a/(3.*(1-Om)*Om)*dOm_da;

    pvecback[pba->index_bg_rho_smg] = rho_tot*Om/(1.-Om);
    //pow(pba->H0,2)/pow(a,3)*Om*(Om-1.)/(Om0-1.)*(1.+a_eq/a)/(1.+a_eq); //this eq is from Pettorino et al, not working
    pvecback[pba->index_bg_p_smg] = w*pvecback[pba->index_bg_rho_smg];

//       if (a>0.9)
// 	printf("a = %e, w = %f, Om_de = %e, rho_de/rho_t = %e \n",a,w,Om,
// 	       pvecback[pba->index_bg_rho_smg]/(pvecback[pba->index_bg_rho_smg]+rho_tot));
  }

  else if (pba->expansion_model_smg == wext) {
    // Set here to zero and update their values later on depending on GR->MG transition redshift
    pvecback[pba->index_bg_rho_smg] = 0.;
    pvecback[pba->index_bg_p_smg] = 0.;
  }

  else if (pba->expansion_model_smg == rho_de) {
    // Set here to zero and update their values later on depending on GR->MG transition redshift
    pvecback[pba->index_bg_rho_smg] = 0.;
    pvecback[pba->index_bg_p_smg] = 0.;
  }

  return _SUCCESS_;
}


/**
* Get the alphas parametrizations.
*
* @param pba           Input: pointer to background structure
* @param a             Input: scale factor
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_B    Input: vector containing all {B} quantities
* @return the error status
*/
int gravity_models_get_alphas_par_smg(
                                      struct background *pba,
                                      double a,
                                      double * pvecback,
                                      double * pvecback_B
                                      ) {


  double * ext_alphas;
  double lna = log(a);
  int last_index=0;

  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg]+pvecback[pba->index_bg_rho_smg];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg]+pvecback[pba->index_bg_p_smg];
  double delta_M_pl = pvecback_B[pba->index_bi_delta_M_pl_smg];
  double Omega_smg = pvecback[pba->index_bg_rho_smg]/rho_tot;

  if (pba->gravity_model_smg == propto_omega) {

    double c_k = pba->parameters_2_smg[0];
    double c_b = pba->parameters_2_smg[1];
    double c_m = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];

    pvecback[pba->index_bg_kineticity_smg] = c_k*Omega_smg;
    pvecback[pba->index_bg_braiding_smg] = c_b*Omega_smg;
    pvecback[pba->index_bg_tensor_excess_smg] = c_t*Omega_smg;
    pvecback[pba->index_bg_mpl_running_smg] = c_m*Omega_smg;
    pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;
  }

  else if (pba->gravity_model_smg == propto_scale) {

    double c_k = pba->parameters_2_smg[0];
    double c_b = pba->parameters_2_smg[1];
    double c_m = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];

    pvecback[pba->index_bg_kineticity_smg] = c_k*a;
    pvecback[pba->index_bg_braiding_smg] = c_b*a;
    pvecback[pba->index_bg_tensor_excess_smg] = c_t*a;
    pvecback[pba->index_bg_mpl_running_smg] = c_m*a;
    pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;
  }

  else if (pba->gravity_model_smg == constant_alphas) {

    double c_k = pba->parameters_2_smg[0];
    double c_b = pba->parameters_2_smg[1];
    double c_m = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];

    pvecback[pba->index_bg_kineticity_smg] = c_k;
    pvecback[pba->index_bg_braiding_smg] = c_b;
    pvecback[pba->index_bg_tensor_excess_smg] = c_t;
    pvecback[pba->index_bg_mpl_running_smg] = c_m;
    pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;
  }

  else if (pba->gravity_model_smg == eft_alphas_power_law) {

    double M_0 = pba->parameters_2_smg[0];
    double c_k = pba->parameters_2_smg[1];
    double c_b = pba->parameters_2_smg[2];
    double c_t = pba->parameters_2_smg[3];
    double M_0_exp = pba->parameters_2_smg[4];
    double c_k_exp = pba->parameters_2_smg[5];
    double c_b_exp = pba->parameters_2_smg[6];
    double c_t_exp = pba->parameters_2_smg[7];

    pvecback[pba->index_bg_kineticity_smg] = c_k*pow(a, c_k_exp);
    pvecback[pba->index_bg_braiding_smg] = c_b*pow(a, c_b_exp);
    pvecback[pba->index_bg_tensor_excess_smg] = c_t*pow(a, c_t_exp);
    pvecback[pba->index_bg_mpl_running_smg] = M_0*M_0_exp*pow(a, M_0_exp)/(1. + M_0*pow(a, M_0_exp));
    pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;
  }

  else if ((pba->gravity_model_smg == eft_gammas_power_law) || (pba->gravity_model_smg == eft_gammas_exponential)) {

    double Omega=0., g1=0., g2=0., g3=0., Omega_p=0., Omega_pp=0., g3_p=0.;
    double Omega_0=0., g1_0=0., g2_0=0., g3_0=0.;
    double Omega_exp=0., g1_exp=0., g2_exp=0., g3_exp=0.;

    Omega_0 = pba->parameters_2_smg[0];
    g1_0 = pba->parameters_2_smg[1];
    g2_0 = pba->parameters_2_smg[2];
    g3_0 = pba->parameters_2_smg[3];
    Omega_exp = pba->parameters_2_smg[4];
    g1_exp = pba->parameters_2_smg[5];
    g2_exp = pba->parameters_2_smg[6];
    g3_exp = pba->parameters_2_smg[7];

    if (pba->gravity_model_smg == eft_gammas_power_law) {
      Omega = Omega_0*pow(a,Omega_exp);
      Omega_p = Omega_0*Omega_exp*pow(a,Omega_exp-1.); // Derivative w.r.t. the scale factor
      Omega_pp = Omega_0*Omega_exp*(Omega_exp-1.)*pow(a,Omega_exp-2.); // Derivative w.r.t. the scale factor
      g1 = g1_0*pow(a,g1_exp);
      g2 = g2_0*pow(a,g2_exp);
      g3 = g3_0*pow(a,g3_exp);
      g3_p = g3_0*g3_exp*pow(a,g3_exp-1.); // Derivative w.r.t. the scale factor
    }
    else { //(pba->gravity_model_smg == eft_gammas_exponential)
      Omega = exp(Omega_0*pow(a,Omega_exp))-1.;
      Omega_p = Omega_0*Omega_exp*pow(a,Omega_exp-1.)*exp(Omega_0*pow(a,Omega_exp)); // Derivative w.r.t. the scale factor
      Omega_pp = Omega_0*Omega_exp*pow(a,Omega_exp-2.)*exp(Omega_0*pow(a,Omega_exp))*(Omega_exp-1.+Omega_0*Omega_exp*pow(a,Omega_exp)); // Derivative w.r.t. the scale factor
      g1 = exp(g1_0*pow(a,g1_exp))-1.;
      g2 = exp(g2_0*pow(a,g2_exp))-1.;
      g3 = exp(g3_0*pow(a,g3_exp))-1.;
      g3_p = g3_0*g3_exp*pow(a,g3_exp-1.)*exp(g3_0*pow(a,g3_exp)); // Derivative w.r.t. the scale factor
    }


    double c_over_H2 = (-pow(a,2.)*Omega_pp + 3.*(Omega + a*Omega_p/2.)*(rho_tot + p_tot)/rho_tot + 3.*(pvecback[pba->index_bg_rho_smg]+pvecback[pba->index_bg_p_smg])/rho_tot)/2.;

    pvecback[pba->index_bg_kineticity_smg] = 2.*(2.*g1*pow(pba->H0,2.)/rho_tot + c_over_H2)/(1. + Omega + g3);
    pvecback[pba->index_bg_braiding_smg] = -(g2*pba->H0/sqrt(rho_tot) + a*Omega_p)/(1. + Omega + g3);
    pvecback[pba->index_bg_tensor_excess_smg] = -g3/(1. + Omega + g3);
    pvecback[pba->index_bg_mpl_running_smg] = a*(Omega_p + g3_p)/(1. + Omega + g3);
    pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
    pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;

  }
  else if(pba->gravity_model_smg == external_alphas){
      pvecback[pba->index_bg_delta_M2_smg] = delta_M_pl; //M2-1
      pvecback[pba->index_bg_M2_smg] = 1.+delta_M_pl;

      if(lna < pba->ext_alphas_lna_smg[0]){
        // MC set all alphas to zero for scale factors smaller than those in the external input file
        // pvecback[pba->index_bg_kineticity_smg] = 1.e-16;
        // pvecback[pba->index_bg_braiding_smg] = 0.;
        // pvecback[pba->index_bg_tensor_excess_smg] = 0.;
        // pvecback[pba->index_bg_mpl_running_smg] = 1.e-16;
        // MC set all alphas to first element in array for scale factors smaller than those in the external input file
        pvecback[pba->index_bg_kineticity_smg] = pba->ext_alphas_smg[pba->index_bg_ext_kineticity_smg];
        pvecback[pba->index_bg_braiding_smg] = pba->ext_alphas_smg[pba->index_bg_ext_braiding_smg];
        pvecback[pba->index_bg_tensor_excess_smg] = pba->ext_alphas_smg[pba->index_bg_ext_tensor_excess_smg];
        pvecback[pba->index_bg_mpl_running_smg] = pba->ext_alphas_smg[pba->index_bg_ext_mpl_running_smg];

        // printf("lna=%e \t alpha_K=%e \t alpha_B=%e \t alpha_M=%e \t alpha_T=%e\n",
        //       lna,
        //       pvecback[pba->index_bg_kineticity_smg],
        //       pvecback[pba->index_bg_braiding_smg],
        //       pvecback[pba->index_bg_mpl_running_smg],
        //       pvecback[pba->index_bg_tensor_excess_smg]);
      }
      else{
        // MC interpolate/extrapolate alphas from values in input table
        class_alloc(ext_alphas,pba->ext_num_alphas*sizeof(double),pba->error_message);

        class_call(array_interpolate_extrapolate_spline(
                                            pba->ext_alphas_lna_smg,
                                            pba->ext_alphas_size_smg,
                                            pba->ext_alphas_smg,
                                            pba->ext_ddalphas_smg,
                                            pba->ext_num_alphas,
                                            lna,
                                            &last_index,
                                            ext_alphas,
                                            pba->ext_num_alphas,
                                            pba->error_message),
                  pba->error_message,
                  pba->error_message);

        pvecback[pba->index_bg_kineticity_smg] = ext_alphas[pba->index_bg_ext_kineticity_smg];
        pvecback[pba->index_bg_braiding_smg] = ext_alphas[pba->index_bg_ext_braiding_smg];
        pvecback[pba->index_bg_tensor_excess_smg] = ext_alphas[pba->index_bg_ext_tensor_excess_smg];
        pvecback[pba->index_bg_mpl_running_smg] = ext_alphas[pba->index_bg_ext_mpl_running_smg];

        // printf("lna=%e \t alpha_K=%e \t alpha_B=%e \t alpha_M=%e \t alpha_T=%e\n",
        //       lna,
        //       pvecback[pba->index_bg_kineticity_smg],
        //       pvecback[pba->index_bg_braiding_smg],
        //       pvecback[pba->index_bg_mpl_running_smg],
        //       pvecback[pba->index_bg_tensor_excess_smg]);

        free(ext_alphas);
      }
  }


  return _SUCCESS_;
}


/**
* Set the right initial conditions for each gravity model.
*
* @param pba                  Input: pointer to background structure
* @param a                    Input: scale factor
* @param pvecback             Output: vector of background quantities (assumed to be already allocated)
* @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
* @param ptr_rho_rad          Input: pointer to the density of radiation
* @return the error status
*/
int gravity_models_initial_conditions_smg(
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

	switch (pba->gravity_model_smg) {

	  case quintessence_monomial:
	  	pvecback_integration[pba->index_bi_phi_smg] = pba->parameters_smg[3];
	  	pvecback_integration[pba->index_bi_phi_prime_smg] = pba->parameters_smg[2]*pba->H0;
			break;

	  case quintessence_tracker:

	    /* Tracker quintessence at early times
	    *  V = H0^2/h^2* V0 * phi^-n exp(lambda*phi^m)
	    *
	    *  log(V/V0) = lambda phi^m - n log (phi)
	    *
	    *  choose phi_ini so V = P_ini*sqrt(rho_rad)
	    *  choose phi_prime_ini so K = K_ini*sqrt(rho_rad)
	    *
	    * initial guess: phi_ini = (log(rho_rad/H0^2)/lambda)^(1/m))
	    * then adjust to the right potential using Newton's method
	    *
	    */

      p1 = pba->parameters_smg[3];// n
      p2 = pba->parameters_smg[4]; //m
      p3 = pba->parameters_smg[5]; //lambda
      // phi_scale = gsl_sf_lambert_W0(10);

      //initial guess
      phi_scale = pow(log(rho_rad/pba->parameters_smg[2]/pow(pba->H0/pba->h,2)/p3), 1./pba->parameters_smg[4]);

      //log of the potential
      V_scale =100;//start off at a high value

			//do a newton root finding
			while (i < 100){

			    V_scale = p3*pow(phi_scale,p2) - p1*log(phi_scale) - pba->parameters_smg[1]*log(rho_rad/pba->parameters_smg[2]/pow(pba->H0/pba->h,2));

			    if (fabs(V_scale) < 1e-5)
			        break;

			    phi_scale -= (V_scale)/((p3*p2*pow(phi_scale,p2)-p1)/phi_scale);
			    i++;
			}
			printf("V_scale %e, i %i \n",V_scale,i);

      pvecback_integration[pba->index_bi_phi_smg] = phi_scale;
      pvecback_integration[pba->index_bi_phi_prime_smg] = -a*sqrt(pba->parameters_smg[0]*2*rho_rad);

			break;

	  case alpha_attractor_canonical:
	  	// Note: we are using as input parameters f = phi/sqrt(alpha)
	  	pvecback_integration[pba->index_bi_phi_smg] = pba->parameters_smg[1]*sqrt(pba->parameters_smg[2]);
	  	pvecback_integration[pba->index_bi_phi_prime_smg] = pba->parameters_smg[0]*pba->H0;
			break;


		/* Attractor IC: H\dot\phi = constant = H_0^2 \xi
		 * phi_prime = xi*a*H_0^2/H (\dot\phi = phi_prime/a)
		 * assumes H = sqrt(rho_rad) initially
		 *
		 * if attractor ICs change xi to fit some attractor value n = 0 with
		 *
		 * n/H0 = xi(c2 - 6c3 xi + 18c4 xi^2 + 5c5 xi^4)
		 *
		 * This is done at the level of background_initial_conditions
		 * (this function does not seem to get access to the parameter being varied!)
		 *
		 * We pick the nonzero xi closest to the user's given value!
		 */

		case galileon:
			if(pba->attractor_ic_smg == _TRUE_){

			  double xi = pba->parameters_smg[0];
			  double c2 = pba->parameters_smg[2];
			  double c3 = pba->parameters_smg[3];
			  double c4 = pba->parameters_smg[4];
			  double c5 = pba->parameters_smg[5];

			  double x[3];
			  double Omega_smg;
			  int complex;
			  int i;
			  double mindiff = 1e100;

			  double xi_start = xi;

			  //Don't use poly_3_split, is bad unless c3 tiny!
			  rf_solve_poly_3(5.*c5,18.*c4,-6.*c3,1.*c2,x,&complex);

			  if (pba->background_verbose > 2)
			  {
			    printf(" galileon attractor \n");
			    printf(" c5 = %.2e, c4 = %.2e, c3 = %.3e, c2 = %.3e \n ",c5,c4,c3,c2);
			    printf(" Solutions (complex = %i): \n",complex);
			  }

			  for (i=0; i<3;i+=1){
			    Omega_smg = pow(x[i],2)*(c2/6. -2.*c3*x[i] +15./2.*c4*pow(x[i],2) + 7./3.*c5*pow(x[i],3));
			    if (x[i]!=0 && fabs(x[i]-xi)<mindiff){
			      xi_start = x[i];
			      mindiff = fabs(x[i]-xi);
			    }
			    if (pba->background_verbose > 2)
			      printf("  xi_%i = %e, Omega_* = %.2e \n",i, x[i],Omega_smg);
			  }
			  if (pba->background_verbose > 2)
			    printf(" Dealer's pick: xi = %e (wish = %e) \n",xi_start,xi);
			  pvecback_integration[pba->index_bi_phi_prime_smg] = a*xi_start*pow(pba->H0,2)/sqrt(rho_rad);
			}
			else {
				/* non attractor ICs */
			  pvecback_integration[pba->index_bi_phi_prime_smg] = a*pba->parameters_smg[0]*pow(pba->H0,2)/sqrt(rho_rad);
			}
			//phi is irrelevant
			pvecback_integration[pba->index_bi_phi_smg] = pba->parameters_smg[6];
			break;

    /* BD IC: note that the field value is basically the planck mass,
     * so its initial value should be around 1
     * the derivative we take to be zero, but this can be extended
     */
    case brans_dicke:
			pvecback_integration[pba->index_bi_phi_smg] = pba->parameters_smg[2];
			pvecback_integration[pba->index_bi_phi_prime_smg] = pba->parameters_smg[3];
			break;

	  case nkgb:
		  {
				/* Action is

		      -X + 1/n * g^[(2n-1)/2] Lambda (X/Lambda^4)^n box(phi)

		      with Lambda^(4n-1)=MPl^(2n-1)*H0^2n

		      g was picked like this so that it approx. remains g*Omega_smg0 ~ O(1) for all n

		      Since the energy density in KGB must be >0, then -inf<xi0<1 and the combination xicomb>0
		      The alpha descrition breaks down for phidot=0, so we assume this is never crossed

		      This all implies that phidot on the attractor has the same sign as g, and therefore
		      phi dot always has the same sign as g.

		      Rshift0 = (phidot0*J0)/rho_smg0, i.e. fraction of the DE energy-density today in the shift-charge as opposed to the vacuum part

		    */

		    double g = pba->parameters_smg[0];
		    double n = pba->parameters_smg[1];
		    double Rshift0 = pba->parameters_smg[2];

	      double H = sqrt(rho_rad); //TODO: for low n -> need to solve tracker + Constraint simultaneously. Increasing H (e.g. double H = 10*sqrt(rho_rad);) works for n~0.65.
	      double H0 = pba->H0;

	      double signg = copysign(1.,g);
	      g=fabs(g);
	      double Rshiftcomb = (2.-Rshift0)/(2.*(1.-Rshift0));

	      double phidot0=0., phidot_attr_init= 0., charge_init=0., phidot_init=0.;

	      phidot0 = signg * sqrt(2./g)*H0 * pow(Rshiftcomb/(3*sqrt(2)),1./(2.*n-1.));             // value of phidot today, if xi0=0, then this is attractor
	      phidot_attr_init = signg * sqrt(2./g)*H0 * pow(H0/(3.*sqrt(2.)*H),1./(2.*n-1.));    // value of phidot on attractor at initial time
	      charge_init  = phidot0*Rshift0/(2*(1-Rshift0))*pow(a,-3);    // implied value of required shift charge initially

	      if(fabs(charge_init/phidot_attr_init)<1.){
	        /* test if initial shift charge is large c.f. the on-attractor phidot. If no, we are nearly on attractor
	        at initial time anyway, so just correct the attractor solution slightly
	        if yes, then we are off-attractor and then we have approximate analytic solution in limit of
	        n_init >> phidot. For the range n_init ~ phidot just use the solution for large shift charge.
	        by the late universe, this will be an irrelevant error. */

	        phidot_init = phidot_attr_init + charge_init/(2.*n-1.);
	      }
	      else{
	        phidot_init = signg * pow( fabs(charge_init) * pow(fabs(phidot_attr_init),2.*n-1.),1./(2.*n));
	      }

	      pvecback_integration[pba->index_bi_phi_smg] = 0.0; //shift-symmetric, i.e. this is irrelevant
	      pvecback_integration[pba->index_bi_phi_prime_smg] = a*phidot_init ;
			}
	    break;

	  case propto_omega:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[4]-1.;
			break;

	  case propto_scale:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[4]-1.;
			break;

	  case constant_alphas:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[4]-1.;
			break;

	  case eft_alphas_power_law:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[0]*pow(a, pba->parameters_2_smg[4]);
			break;

	  case eft_gammas_power_law:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[0]*pow(a,pba->parameters_2_smg[4]) + pba->parameters_2_smg[3]*pow(a,pba->parameters_2_smg[7]);
			break;

	  case eft_gammas_exponential:
			pvecback_integration[pba->index_bi_delta_M_pl_smg] = exp(pba->parameters_2_smg[0]*pow(a,pba->parameters_2_smg[4])) + exp(pba->parameters_2_smg[3]*pow(a,pba->parameters_2_smg[7])) -2.;
			break;

    case external_alphas:
	    pvecback_integration[pba->index_bi_delta_M_pl_smg] = pba->parameters_2_smg[0]-1.;
	    break;

    case stable_params:
      /* sets IC at a_initial for integration of alpha_M. It's only needed to initially fill in array and it will be adjusted
      after backward integration*/
	    // pvecback_integration[pba->index_bi_delta_M_pl_smg] = 0.;
      /* For each backward integrated variable fix the initial condition at a_final. */
      pvecback_bw_integration[pba->index_bibw_B_tilde_smg] = 1.; // Arbitrary constant. We set it to 1
	    pvecback_bw_integration[pba->index_bibw_dB_tilde_smg] = 1. - 0.5 * pba->parameters_2_smg[0]; // 1 - alpha_B0/2
	    break;
	  }

  /* expansion_smg when w is parametrized and rho obtained with integration */
	if (pba->rho_evolution_smg == _TRUE_){

		//DT: moved here from the definition of wowa_w model, to avoid pvecback[pba->index_bg_rho_smg] mix up
    // TODO_EB: rethink this
		if (pba->expansion_model_smg == wowa_w) {
	    double Omega_const_smg = pba->parameters_smg[0];
	    double w0 = pba->parameters_smg[1];
	    double wa = pba->parameters_smg[2];

	    double w_smg = w0+(1.-a)*wa;
	    double wi = w0+wa;
	    double wf = w0;
	    double rho_smg_today = Omega_const_smg * pow(pba->H0,2);
	    double integral_smg = 3.*((1.+ wi)*log(1./a) + (wi-wf)*(a-1.));
	    pvecback_integration[pba->index_bi_rho_smg] = rho_smg_today * exp(integral_smg);
	  }
		else {
	  	pvecback_integration[pba->index_bi_rho_smg] = pvecback[pba->index_bg_rho_smg];
		}
	}

	return _SUCCESS_;
}


/**
* Send _smg information about specific models to the standard output.
*
* @param pba                  Input: pointer to background structure
* @return the error status
*/
int gravity_models_print_stdout_smg(
				  													struct background *pba
				  													){

  switch (pba->gravity_model_smg) {

    case quintessence_monomial:
      printf("Modified gravity: quintessence_monomial with parameters: \n");
      printf("-> N = %g, V0 = %g, V0* = %g, phi_prime_ini = %g, phi_ini = %g \n",
	      pba->parameters_smg[0], pba->parameters_smg[1], pba->parameters_smg[1]*pow(pba->H0/pba->h,2),
	      pba->parameters_smg[2], pba->parameters_smg[3]);
    break;

    case quintessence_tracker:
      printf("Modified gravity: quintessence_tracker with parameters: \n");
      printf("-> K_ini = %g, P_ini = %g, V0 = %g, V0* = %g, n = %g, m = %g, lambda=%g  \n",
	    pba->parameters_smg[0], pba->parameters_smg[1], pba->parameters_smg[2]*pow(pba->H0/pba->h,2),
	    pba->parameters_smg[2], pba->parameters_smg[3], pba->parameters_smg[4], pba->parameters_smg[5]);
    break;

    case alpha_attractor_canonical:
      printf("Modified gravity: alpha_attractor_canonical with parameters: \n");
      printf("-> f = phi/sqrt(alpha) \n");
      printf("-> phi_prime_ini = %g, f_ini = %g, alpha = %g, c = %g, p = %g, n = %g \n",
	      pba->parameters_smg[0], pba->parameters_smg[1], pba->parameters_smg[2],
	      pba->parameters_smg[3], pba->parameters_smg[4], pba->parameters_smg[5]);
    break;

    case galileon:
      printf("Modified gravity: covariant Galileon with parameters: \n");
      printf(" -> c_1 = %g, c_2 = %g, c_3 = %g \n    c_4 = %g, c_5 = %g, xi_ini = %g (xi_end = %g) \n",
	    pba->parameters_smg[1],pba->parameters_smg[2],pba->parameters_smg[3],pba->parameters_smg[4],pba->parameters_smg[5],pba->parameters_smg[0], pba->xi_0_smg);
    break;

    case brans_dicke:
      printf("Modified gravity: Brans Dicke with parameters: \n");
      printf(" -> Lambda = %g, omega_BD = %g, \n    phi_ini = %g (phi_0 = %g), phi_prime_ini = %g \n",
	    pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2],pba->phi_0_smg,pba->parameters_smg[3]);
    break;

    case nkgb:
     printf("Modified gravity: Kinetic Gravity Braiding with K=-X and G=1/n g^(2n-1)/2 * X^n with parameters: \n");
     printf(" -> g = %g, n = %g, phi_ini = 0.0, smg density fraction from shift charge term = %g. \n",
	    pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2]);
    break;

    case propto_omega:
      printf("Modified gravity: propto_omega with parameters: \n");
      printf(" -> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	     pba->parameters_2_smg[4]);
    break;

    case propto_scale:
      printf("Modified gravity: propto_scale with parameters: \n");
      printf(" -> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	     pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	     pba->parameters_2_smg[4]);
    break;

    case constant_alphas:
      printf("Modified gravity: constant_alphas with parameters: \n");
      printf(" -> c_K = %g, c_B = %g, c_M = %g, c_T = %g, M_*^2_init = %g \n",
	      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	      pba->parameters_2_smg[4]);
    break;

    case eft_alphas_power_law:
      printf("Modified gravity: eft_alphas_power_law with parameters: \n");
      printf(" -> M_*^2_0 = %g, c_K = %g, c_B = %g, c_T = %g, M_*^2_exp = %g, c_K_exp = %g, c_B_exp = %g, c_T_exp = %g\n",
	      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],
	      pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
    break;

    case eft_gammas_power_law:
      printf("Modified gravity: eft_gammas_power_law with parameters: \n");
      printf(" -> Omega_0 = %g, gamma_1 = %g, gamma_2 = %g, gamma_3 = %g, Omega_0_exp = %g, gamma_1_exp = %g, gamma_2_exp = %g, gamma_3_exp = %g \n",
	      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
    break;

    case eft_gammas_exponential:
      printf("Modified gravity: eft_gammas_exponential with parameters: \n");
      printf(" -> Omega_0 = %g, gamma_1 = %g, gamma_2 = %g, gamma_3 = %g, Omega_0_exp = %g, gamma_1_exp = %g, gamma_2_exp = %g, gamma_3_exp = %g \n",
	      pba->parameters_2_smg[0],pba->parameters_2_smg[1],pba->parameters_2_smg[2],pba->parameters_2_smg[3],pba->parameters_2_smg[4],pba->parameters_2_smg[5],pba->parameters_2_smg[6],pba->parameters_2_smg[7]);
    break;

    case external_alphas:
      printf("Modified gravity: external_alphas from file: \n");
      printf("%s\n",
	      pba->smg_file_name);
     break;

    case stable_params:
      printf("Modified gravity: stable_params from file: \n");
      printf("%s\n",
	      pba->smg_file_name);
     break;

    default:
      printf("Modified gravity: output not implemented in print_stdout_gravity_parameters_smg() \n");
  }

  if(pba->field_evolution_smg==_FALSE_) {
    switch (pba->expansion_model_smg) {

      case lcdm:
        printf("Parameterized model with LCDM expansion \n");
        printf("-> Omega_smg = %f \n",pba->parameters_smg[0]);
      break;

      case wowa:
        printf("Parameterized model with CPL expansion \n");
        printf("-> Omega_smg = %f, w0 = %f, wa = %e \n",
	       pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2]);
      break;

      case wowa_w:
        printf("Parameterized model with CPL expansion \n");
        printf("-> Omega_smg = %f, w0 = %f, wa = %e \n",
	       pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2]);
      break;

      case wede:    //ILSWEDE
        printf("Parameterized model with variable EoS + Early DE \n");
        printf("-> Omega_smg = %f, w = %f, Omega_e = %f \n",pba->parameters_smg[0],pba->parameters_smg[1],pba->parameters_smg[2]);
      break;

      case wext:
        printf("Parameterized model with external smg equation of state from file: \n");
        printf("%s\n",pba->expansion_file_name);
      break;

      case rho_de:
        printf("Parameterized model with external smg normalised energy density from file: \n");
        printf("%s\n",pba->expansion_file_name);
      break;

      default:
        printf("Parameterized model: expansion history output not implemented in print_stdout_gravity_parameters_smg() \n");

    }
  }

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
int copy_to_aux_array_smg(
																 struct background *pba,
                                 int row,
                                 int column,
                                 double value
															   ) {

	/* needed for growing table */
	void * memcopy_result;

	memcopy_result = memcpy(pba->stable_params_aux_smg + row*pba->num_stable_params_aux + column,
	                        &value, 1*sizeof(double));
	class_test(memcopy_result != pba->stable_params_aux_smg + row*pba->num_stable_params_aux + column,
	           pba->error_message, "cannot copy data to pba->stable_params_aux_smg");

  return _SUCCESS_;

  }
