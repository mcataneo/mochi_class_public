#include "input_smg.h"

#include "background.h"
#include "perturbations.h"

int input_warnings_smg(
  struct perturbations * ppt,
  int input_verbose
) {

  /* Here we put a warning as we want to encourage hi_class users to get
  h_prime from the trace of the Einstein ij equation than from the Einstein 00
  equation. This is because the Einstein 00 equation has a gauge dependent
  singularity that can be removed using the trace of the Einstein ij equation.
  */
  if (input_verbose > 0) {
    if (ppt->get_h_from_trace == _FALSE_) {
      printf("\n");
      printf("WARNING: you set get_h_from_trace to False.\n");
      printf("While this is still accepted in hi_class, it can cause gauge dependent\n");
      printf("singularities if your model crosses alphaB=2. For this reason in\n");
      printf("future versions of the code this option will be removed and the\n");
      printf("Einstein 00 equation will be used only to set the initial conditions\n");
      printf("for h_prime and as a test to check that it is satisfied during the evolution.\n");
      printf("On the other hand this is a safe choice if you want very large k modes\n");
      printf("(typically k>10 Mpc^-1), where the constraint and the dynamical equations\n");
      printf("disagree by a non negligible amount in some of the models studied.\n");
      printf("\n");
    }
    else if (ppt->get_h_from_trace == _TRUE_) {
      printf("\n");
      printf("WARNING: you set get_h_from_trace to True.\n");
      printf("While this will be the default option in future versions of hi_class it might\n");
      printf("be safer to set get_h_from_trace to False if you want very large k modes\n");
      printf("(typically k>10 Mpc^-1). In this regime the constraint and the dynamical \n");
      printf("equations disagree by a non negligible amount some of the models studied.\n");
      printf("\n");
    }
  }

  return _SUCCESS_;
}

int input_read_parameters_smg(
  struct file_content * pfc,
  struct precision * ppr,
  struct background * pba,
  struct perturbations * ppt,
  ErrorMsg errmsg
) {

  int flag1, flag2, flag3;
  double param1;
  int entries_read;
  int int1, n;
  char string1[_ARGUMENT_LENGTH_MAX_];
  char string2[_ARGUMENT_LENGTH_MAX_];

  pba->has_smg = _TRUE_;  //default is _FALSE_

  /** Main flag for the quasi-static approximation scheme */

  class_call(parser_read_string(pfc,"method_qs_smg",&string1,&flag1,errmsg),
             errmsg,
             errmsg);

  if (flag1 == _TRUE_) {

    if ((strstr(string1,"automatic") != NULL) || (strstr(string1,"a") != NULL) || (strstr(string1,"A") != NULL)) {
      ppt->method_qs_smg = automatic;
    }

    if ((strstr(string1,"fully_dynamic") != NULL) || (strstr(string1,"fd") != NULL) || (strstr(string1,"FD") != NULL)) {
      ppt->method_qs_smg = fully_dynamic;
    }

    if ((strstr(string1,"quasi_static") != NULL) || (strstr(string1,"qs") != NULL) || (strstr(string1,"QS") != NULL)) {
      ppt->method_qs_smg = quasi_static;
    }

    if ((strstr(string1,"fully_dynamic_debug") != NULL) || (strstr(string1,"fdd") != NULL) || (strstr(string1,"FDD") != NULL)) {
      ppt->method_qs_smg = fully_dynamic_debug;
    }

    if ((strstr(string1,"quasi_static_debug") != NULL) || (strstr(string1,"qsd") != NULL) || (strstr(string1,"QSD") != NULL)) {
      ppt->method_qs_smg = quasi_static_debug;
    }
  }

  class_call(parser_read_string(pfc, "use_pert_var_deltaphi_smg", &string1, &flag1, errmsg),
    errmsg,
    errmsg);

  if (flag1 == _TRUE_){
    if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
      ppt->use_pert_var_deltaphi_smg = _TRUE_;
    }
    else{
      ppt->use_pert_var_deltaphi_smg = _FALSE_;
    }
  }

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

  class_call(parser_read_string(pfc,"gravity_model",&string1,&flag1,errmsg),
       errmsg,
       errmsg);


  if (flag1 == _FALSE_) {
    printf(" gravity_model not read, default will be used \n");
  }
  else {
  /** Read tuning parameter and guess for the parameter variation range
   * These can be adjusted latter on a model basis
   */
  int has_tuning_index_smg, has_dxdy_guess_smg;

  class_read_int("tuning_index_smg",pba->tuning_index_smg);
  has_tuning_index_smg = flag1;

  class_read_double("tuning_dxdy_guess_smg",pba->tuning_dxdy_guess_smg);
  has_dxdy_guess_smg = flag1;
  if (has_dxdy_guess_smg == _FALSE_)
    pba->tuning_dxdy_guess_smg = 1;

  /** Loop over the different models
   * flag2 keeps track of whether model has been identified
   */
    flag2=_FALSE_;

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

  if (strncmp("quintessence", string1, strlen("quintessence")) == 0){
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

if (has_dxdy_guess_smg == _FALSE_){

  if(pba->tuning_index_smg == 1){
//           if(phi_ini_smg != 0){
          V0 = pba->Omega0_smg/pow(phi_end_guess,N);
          pba->tuning_dxdy_guess_smg = 1./pow(phi_end_guess,N);
          pba->parameters_smg[1] = V0;
//           }
//           else{
//             V0 = pba->Omega0_smg/pow(1.e-40,N);
//             pba->tuning_dxdy_guess_smg = 1./pow(1.e-40,N);
//             pba->parameters_smg[1] = V0;
//
//           }
  }

  if(pba->tuning_index_smg == 3){
     phi_ini_smg = pow(pba->Omega0_smg/V0, 1./N);
     pba->parameters_smg[3] = phi_ini_smg;
     pba->tuning_dxdy_guess_smg = phi_ini_smg/(pba->Omega0_smg)/N;
  }
}//end of no has_dxdy_guess_smg
    }//end of quintessence_monomial


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
    }//end of no has_dxdy_guess_smg
    } //end of tracker


    if (strcmp(string1,"alpha_attractor_canonical") == 0) {
pba->gravity_model_smg = alpha_attractor_canonical;
pba->field_evolution_smg = _TRUE_;
  pba->is_quintessence_smg = _TRUE_;
flag2=_TRUE_;

pba->parameters_size_smg = 6;
class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

  class_call(parser_read_string(pfc,"log_10_param_alpha",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

    if(flag1 == _TRUE_ && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))){
      pba->parameters_smg[2] = pow(10, pba->parameters_smg[2]);
    }

  class_call(parser_read_string(pfc,"use_phi_no_f",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

    if(flag1 == _TRUE_ && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))){
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
    }//end Omega_smg_debug

    } //endif  alpha_attractor_canonical



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
if (flag3==_FALSE_){

  class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);

}
else {//a submodel is given

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

}//end of submodels


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

  class_call(parser_read_string(pfc,"attractor_ic_smg",&string1,&flag1,errmsg),
       errmsg,
       errmsg);

    if(flag1 == _TRUE_ && ((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL))){
        pba->attractor_ic_smg = _TRUE_;
    }
    else{
        pba->attractor_ic_smg = _FALSE_;
    }

    }//end of Galileon
    if (strcmp(string1,"brans dicke") == 0 || strcmp(string1,"Brans Dicke") == 0 || strcmp(string1,"brans_dicke") == 0) {
pba->gravity_model_smg = brans_dicke;
pba->field_evolution_smg = _TRUE_;
flag2=_TRUE_;

pba->parameters_size_smg = 4;
class_read_list_of_doubles("parameters_smg",pba->parameters_smg,pba->parameters_size_smg);
pba->parameters_smg[0] = 2*pba->Omega0_smg;
pba->tuning_dxdy_guess_smg = 0.5;
pba->tuning_index_2_smg = 2;
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

    class_test(flag2==_FALSE_,
   errmsg,
   "could not identify gravity_theory value, check that it is one of 'propto_omega', 'propto_scale', 'constant_alphas', 'eft_alphas_power_law', 'eft_gammas_power_law', 'eft_gammas_exponential', 'brans_dicke', 'galileon', 'nKGB', 'quintessence_monomial', 'quintessence_tracker', 'alpha_attractor_canonical' ...");

  }// end of loop over models

  if(pba->field_evolution_smg == _TRUE_){

    //TODO: include generic stuff for covariant theories

  }
  else { //if no self-consistent evolution, need a parameterization for Omega_smg

    class_test(ppt->use_pert_var_deltaphi_smg==_TRUE_,
      errmsg,
      "It is not consistent to evolve delta_phi_smg and choose parametrized models.");

    class_call(parser_read_string(pfc,"expansion_model",&string1,&flag1,errmsg),
   errmsg,
   errmsg);
    if (flag1 == _FALSE_)
printf("No expansion model specified, will take default one \n");

    flag2 = _FALSE_;

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


    if (strcmp(string1,"wede") == 0) {    //ILSWEDE
      pba->expansion_model_smg = wede;
      flag2=_TRUE_;
      pba->parameters_size_smg = 3;
      class_read_list_of_doubles_or_default("expansion_smg",pba->parameters_smg,0.0,pba->parameters_size_smg);
      // 	//optimize the guessing BUG: eventually leads to problem in the MCMC, perhaps the guess is too good?
      // 	if(pba->tuning_index_smg == 0){
      // 	  pba->parameters_smg[0] = pba->Omega0_smg;
      // 	}
    }
    class_test(flag2==_FALSE_,
   errmsg,
   "could not identify expansion_model value, check that it is either lcdm, wowa, wowa_w, wede ...");

  }

  /** Other generic specifications:
   * - whether stability tests are skipped (skip_stability_tests_smg) or softened (cs2_safe_smg)
   * - thresholds for approximations in the cubic Friedmann equation
   * - add a value to have better behaved perturbations
   * - approximations in the perturbations
   */

  class_read_double("cs2_safe_smg",pba->cs2_safe_smg);
  class_read_double("D_safe_smg",pba->D_safe_smg);
  class_read_double("ct2_safe_smg",pba->ct2_safe_smg);
  class_read_double("M2_safe_smg",pba->M2_safe_smg);

  class_read_double("pert_ic_tolerance_smg",ppr->pert_ic_tolerance_smg);
  class_read_double("pert_ic_ini_z_ref_smg",ppr->pert_ic_ini_z_ref_smg);
  class_read_double("pert_ic_regulator_smg",ppr->pert_ic_regulator_smg);
  class_read_double("pert_qs_ic_tolerance_test_smg",ppr->pert_qs_ic_tolerance_test_smg);

  class_read_double("a_min_stability_test_smg",pba->a_min_stability_test_smg);

  class_read_double("kineticity_safe_smg",pba->kineticity_safe_smg); // minimum value of the kineticity (to avoid trouble)
  class_read_double("min_a_pert_smg",ppr->min_a_pert_smg);

  class_call(parser_read_string(pfc,
			  "skip_stability_tests_smg",
			  &string1,
			  &flag1,
			  errmsg),
	errmsg,
	errmsg);

  if (flag1 == _TRUE_){
    if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
pba->skip_stability_tests_smg = _TRUE_;
    }
    else{
pba->skip_stability_tests_smg = _FALSE_;
    }
  }

  //IC for perturbations
  class_call(parser_read_string(pfc,
			  "pert_initial_conditions_smg",
			  &string1,
			  &flag1,
			  errmsg),
	errmsg,
	errmsg);

  if (strcmp(string1,"single_clock") == 0) {
    ppt->pert_initial_conditions_smg = single_clock;
    }
  else if (strcmp(string1,"gravitating_attr") == 0) {
    ppt->pert_initial_conditions_smg = gravitating_attr;
    }
  else if (strcmp(string1,"zero") == 0) {
    ppt->pert_initial_conditions_smg = zero;
    }
  else if (strcmp(string1,"kin_only") == 0) {
    ppt->pert_initial_conditions_smg = kin_only;
  }
  else if (strcmp(string1,"ext_field_attr") == 0 ){//this is the default
    ppt->pert_initial_conditions_smg = ext_field_attr;
  }

//     else {
//       if (ppt->perturbations_verbose > 1)
// 	printf(" Initial conditions for Modified gravity perturbations not specified, using default \n");
//     }

  /** re-assign shooting parameter (for no-tuning debug mode) */
  if (pba->Omega_smg_debug == 0)
    class_read_double("shooting_parameter_smg",pba->parameters_smg[pba->tuning_index_smg]);

  // test that the tuning is correct
  class_test(pba->tuning_index_smg >= pba->parameters_size_smg,
       errmsg,
       "Tuning index tuning_index_smg = %d is larger than the number of entries %d in parameters_smg. Check your .ini file.",
       pba->tuning_index_smg,pba->parameters_size_smg);

   /** Read the desired Planck mass and check that the necessary information is provided.
    *  if needed re-assign shooting parameter for the Planck mass
    */
   flag1==_FALSE_;
   class_read_double("M_pl_today_smg",pba->M_pl_today_smg);
   if (flag1==_TRUE_){

     class_test(pba->gravity_model_smg!=brans_dicke,
 	 errmsg,
 	 "You asked to tune M_pl(today) to %e but currently this is only allowed for Brans-Dicke\n",
 	 pba->M_pl_today_smg);

     class_call(parser_read_string(pfc,"normalize_G_NR",
 			  &string1,
 			  &flag1,
 			  errmsg),
 	errmsg,
 	errmsg);

   if (flag1 == _TRUE_){
     if((strstr(string1,"y") != NULL) || (strstr(string1,"Y") != NULL)){
         double omega_BD = pba->parameters_smg[1];
         pba->M_pl_today_smg = (4.+2.*omega_BD)/(3.+2.*omega_BD);
     }
   }

     class_read_double("param_shoot_M_pl_smg",pba->parameters_smg[pba->tuning_index_2_smg]);
 //       printf("updating param = %e to tune M_pl \n",pba->parameters_smg[pba->tuning_index_2_smg]);
   }

  //how much info on background.dat?
  class_read_double("output_background_smg",pba->output_background_smg);

  return _SUCCESS_;
}

int input_readjust_precision(
  struct precision * ppr
) {

  /** readjust some precision parameters for modified gravity */

  //otherwise problems with ISW effect
  if (ppr->perturbations_sampling_stepsize > 0.05)
    ppr->perturbations_sampling_stepsize=0.05;

  return _SUCCESS_;
}

int input_default_params_smg(
  struct background * pba,
  struct perturbations * ppt
) {

  /** - background structure */

  pba->gravity_model_smg = propto_omega; /* gravitational model */
  pba->expansion_model_smg = lcdm; /*expansion model (only for parameterizations*/
  pba->Omega0_smg = 0.; /* Scalar field defaults */
  pba->M_pl_today_smg = 1.; //*Planck mass today*/
  pba->M_pl_tuning_smg = _FALSE_; //* Tune Planck mass?*/
  pba->Omega_smg_debug = 0;
  pba->field_evolution_smg = _FALSE_; /* does the model require solving the background equations? */
  pba->M_pl_evolution_smg = _FALSE_; /* does the model require integrating M_pl from alpha_M? */
  pba->skip_stability_tests_smg = _FALSE_; /*if you want to skip the stability tests for the perturbations */
  pba->a_min_stability_test_smg = 0; /** < skip stability tests for a < a_min */

  pba->hubble_evolution = _TRUE_; /** dynamical evolution of Friedmann eq. */
  pba->hubble_friction = 3.; /** friction coefficient in H' equation: H' = ... + H_friction*(H^2 - rho_crit) [NOT ONLY IN SMG!] */
  pba->rho_evolution_smg = _FALSE_; /*does the model require to evolve the background energy density? (only for parameterizations)*/

  pba->kineticity_safe_smg = 0; /* value added to the kineticity, useful to cure perturbations at early time in some models */
  pba->cs2_safe_smg = 0; /* threshold to consider the sound speed of scalars negative in the stability check */
  pba->D_safe_smg = 0; /* threshold to consider the kinetic term of scalars negative in the stability check */
  pba->ct2_safe_smg = 0; /* threshold to consider the sound speed of tensors negative in the stability check */
  pba->M2_safe_smg = 0; /* threshold to consider the kinetic term of tensors (M2) negative in the stability check */


  /*set stability quantities to nonzero values*/
  pba->min_M2_smg = 1e10;
  pba->min_ct2_smg = 1e10;
  pba->min_D_smg = 1e10;
  pba->min_cs2_smg = 1e10;

  pba->min_bra_smg = 4.;
  pba->max_bra_smg = 0.;
  pba->quintessence_w_safe_smg = 0;

  pba->parameters_smg = NULL;
  pba->parameters_size_smg = 0;
  pba->tuning_index_smg = 0;
  pba->tuning_dxdy_guess_smg = 1;

  pba->output_background_smg = 1; /**< amount of information printed onbackground.dat output */

  pba->has_smg= _FALSE_;
  pba->parameters_tuned_smg = _FALSE_;
  pba->shooting_failed = _FALSE_;
  pba->is_quintessence_smg = _FALSE_;
  pba->attractor_ic_smg = _TRUE_;  /* only read for those models in which it is implemented */
  pba->initial_conditions_set_smg = _FALSE_;


  /** - perturbations structure */

  ppt->get_h_from_trace=_FALSE_; /* Get h' from Einstein trace rather than 00 (not only _smg!!) */

  ppt->method_qs_smg=fully_dynamic;
  ppt->pert_initial_conditions_smg = ext_field_attr; /* default IC for perturbations in the scalar */

  ppt->use_pert_var_deltaphi_smg=_FALSE_;

  return _SUCCESS_;
}
