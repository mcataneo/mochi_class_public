#include "gravity_models_smg.h"


int gravity_models_properties_smg() {

  return _SUCCESS_;
}


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


// input_read_parameters_smg
// background_gravity_functions_smg
// background_initial_conditions_smg
// print_stdout_gravity_parameters_smg
