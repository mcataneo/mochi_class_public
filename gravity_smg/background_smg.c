#include "background_smg.h"

/** This function fills the modified gravity part of background_functions.
 * First all the Horndeski functions G_i(X,phi) are computed.
 * A loop is used to allow different implementations without erasing previous ones
 * Note that in CLASS units a canonical field has G2 = X ~ [Mpc]^-2
 * This implies that phi_prime ~ [Mpc]^-1
 * and phi is dimensionless (we can think of phi as given in units of the Planck mass
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * Numerical derivatives may further serve as a consistency check.
 * TODO: Add a background_write_alpha_primes
 */

int background_gravity_functions_smg(
				 struct background *pba,
				 double * pvecback_B,
				 short return_format,
				 double * pvecback
				 ){

      // scalar field + curvature not yet implemented
  class_test(pba->K !=0 ,
	     pba->error_message,
	     "has_smg with curvature K = %e not yet implemented",pba->K);

  if (pba->field_evolution_smg == _TRUE_) {

    /* declare variables and set defaults to zero */
    double a, phi, phi_prime, H;
    double x,f,df;
    int n, n_max=100;
    double X,rho_tot,p_tot;
    double E0,E1,E2,E3;
    double P0,P1,P2;
    double R0,R1,R2;
    double G2=0;
    double G2_X=0, G2_XX=0, G2_XXX=0;
    double G2_phi=0, G2_Xphi=0, G2_XXphi=0;
    double G2_phiphi=0, G2_Xphiphi=0;
    double G3_X=0, G3_XX=0, G3_XXX=0;
    double G3_phi=0, G3_Xphi=0, G3_XXphi=0;
    double G3_phiphi=0, G3_Xphiphi=0;
    double G3_phiphiphi=0;
    double G4 = 1./2., DG4=0;
    double G4_X=0, G4_XX=0, G4_XXX=0, G4_XXXX=0;
    double G4_phi=0, G4_Xphi=0, G4_XXphi=0, G4_XXXphi=0;
    double G4_phiphi=0, G4_Xphiphi=0, G4_XXphiphi=0;
    double G4_phiphiphi=0, G4_Xphiphiphi=0;
    double G5_X=0, G5_XX=0, G5_XXX=0, G5_XXXX=0;
    double G5_phi=0, G5_Xphi=0, G5_XXphi=0, G5_XXXphi=0;
    double G5_phiphi=0, G5_Xphiphi=0, G5_XXphiphi=0;
    double G5_phiphiphi=0, G5_Xphiphiphi=0;
    double F4=0;
    double F4_X=0, F4_XX=0, F4_XXX=0;
    double F4_phi=0, F4_Xphi=0, F4_XXphi=0;
    double F4_phiphi=0, F4_Xphiphi=0;
    double F5=0;
    double F5_X=0, F5_XX=0, F5_XXX=0;
    double F5_phi=0, F5_Xphi=0, F5_XXphi=0;
    double F5_phiphi=0, F5_Xphiphi=0;

    a = pvecback_B[pba->index_bi_a];
    if (pba->hubble_evolution == _TRUE_)
      H = pvecback_B[pba->index_bi_H];

    phi = pvecback_B[pba->index_bi_phi_smg];
    phi_prime = pvecback_B[pba->index_bi_phi_prime_smg];

    X = 0.5*pow(phi_prime/a,2);
    rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
    p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

    /* Overwrite the Horneski functions and their partial derivatives depending on the model */
    if (pba->gravity_model_smg == quintessence_monomial) {

      double N = pba->parameters_smg[0];
      double V0 = pba->parameters_smg[1];

      G2 = X - V0*pow(pba->H0/pba->h,2.)*pow(phi,N); // V written as in arXiv:1406.2301 in CLASS units
      G2_X = 1.;
      G2_phi = -N*V0*pow(pba->H0/pba->h,2.)*pow(phi,N-1.);
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

      G2 = X - V;
      G2_X = 1.;
      G2_phi = -V_phi;
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

      G2 = X - V;
      G2_X = 1.;
      G2_phi = -V_phi;

     class_test((phi < 0. ) && ( phi_prime > 0 )
             && (pba->parameters_tuned_smg == _TRUE_)
             && (pba->skip_stability_tests_smg == _FALSE_),
             pba->error_message,
             "The model has started oscillating with first minimum at a = %e. Since <w> = 0 yields matter, it cannot make the universe expand accelerately.", a);
    }

    else if(pba->gravity_model_smg == galileon){//TODO: change name with ccc_galileon

      double M3 = pow(pba->H0,2);

      double Lambda1 = pba->parameters_smg[1]*M3;
      double Lambda2 = pba->parameters_smg[2];
      double Lambda3 = pba->parameters_smg[3]/M3;
      double Lambda4 = pba->parameters_smg[4]*pow(M3,-2);
      double Lambda5 = pba->parameters_smg[5]*pow(M3,-3);

      G2 = Lambda2*X - 0.5*Lambda1*phi;
      G2_X = Lambda2;
      G2_phi = - 0.5*Lambda1;
      /*  G_3 = -2*Lambda3*X */
      G3_X = -2.*Lambda3;
      /* G_4 = 1/2 + Lambda4 X^2 */
      DG4 = Lambda4*pow(X,2);
      G4 = 1/2. + DG4;
      G4_X = 2.*Lambda4*X;
      G4_XX = 2.*Lambda4;
      /* G_5 = Lambda5*pow(X,2) */
      G5_X = 2.*Lambda5*X;
      G5_XX = 2.*Lambda5;
    }
    else if(pba->gravity_model_smg == brans_dicke){

      /* Brans-Dicke can't cause acceleration:
       * - V is a constant potential, basically a cosmological constant
       * - omega is the Brans-Dicke parameter in charge of the fifth force
       */
      double V = 3.*pba->parameters_smg[0]*pow(pba->H0,2);
      double omega = pba->parameters_smg[1];

      G2 = -V + omega*X/phi;
      G2_X = omega/phi;
      G2_Xphi = -omega/pow(phi,2);
      G2_phi = -omega*X/pow(phi,2);

      DG4 = (phi-1.)/2.;
      G4 = phi/2.;
      G4_phi = 1./2.;

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

      G2    = -X;
      G2_X  = -1.;

      // G3 = 1/n g^[(2n-1)/2] Lambda (X/Lambda^4)^n

      G3_X = npow*ngpow*pow(X,npow-1)/pow(H0,2*npow);
      G3_XX = npow*(npow-1.)*ngpow*pow(X,npow-2)/pow(H0,2*npow);

    }



    //TODO: Write the Bellini-Sawicki functions and other information to pvecback

    pvecback[pba->index_bg_phi_smg] = phi; // value of the scalar field phi
    pvecback[pba->index_bg_phi_prime_smg] = phi_prime; // value of the scalar field phi derivative wrt conformal time

    /** - Modified time-time Friedmann equation
     * E0 + E1 H + E3 H^3 = E2 H^2
     * NOTE: H is NOT the curly H!!
     * NOTE: added rho_smg, p_smg separately
     */

    E0 =
    1./3.*(
      + 3.*rho_tot - G2 + 2.*X*(G2_X - G3_phi)
    );

    E1 =
    - 2.*(
      + G4_phi - X*(G3_X - 2.*G4_Xphi)
    )*phi_prime/a;

    E2 =
    2.*(
      1./2. + DG4 - X*(4.*G4_X - 3.*G5_phi)
      - 2.*pow(X,2)*(2.*G4_XX - G5_Xphi + 10.*F4) - 8.*pow(X,3)*F4_X
    );

    E3 =
    2./3.*X*(
      + 5.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
    )*phi_prime/a;

    /* Rewrite if ICs not set or no evolution */
    if (pba->initial_conditions_set_smg == _FALSE_ || pba->hubble_evolution == _FALSE_){
      class_test_except(E3*pow(E0,1./2.) > 1e-10 && pba->initial_conditions_set_smg == _FALSE_,
		  pba->error_message,
      free(pvecback);free(pvecback_B),
           " E3=%e is large in Friedmann constraint when setting ICs ",  E3);
      /* Use Newton's method */
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
      if (pba->background_verbose > 5 && pba->initial_conditions_set_smg == _FALSE_ )
	printf(" Initial H = %e, sqrt(rho) = %e, ratio = %e, n=%i \n", H, sqrt(rho_tot),sqrt(rho_tot)/H,n);
    }

    pvecback[pba->index_bg_E0_smg] = E0;
    pvecback[pba->index_bg_E1_smg] = E1;
    pvecback[pba->index_bg_E2_smg] = E2;
    pvecback[pba->index_bg_E3_smg] = E3;

    /* Rewritten by the constraint if ICs have not been set */
    pvecback[pba->index_bg_H] = H;

    class_test_except(isnan(pvecback[pba->index_bg_H]),
	       pba->error_message,
         free(pvecback);free(pvecback_B),
               " H=%e is not a number at a = %e. phi = %e, phi_prime = %e, E0=%e, E1=%e, E2=%e, E3=%e ",
	       pvecback[pba->index_bg_H],a,phi,phi_prime,E0,E1,E2,E3);


    /** TODO: add references
      * - Modified space-space Friedmann equation and scalar field equation
      * P0 phi'' + P1 H' + P2 = 0
      * R0 phi'' + R1 H' + R2 = 0
      * They will be mixed in general: write the coefficients and solve for phi'', H'
      */

    P0 =
    - 2./3.*(
      + G4_phi - X*(G3_X - 2.*G4_Xphi)
      - 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*H*phi_prime/a
      - X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,2)
    )/a;

    P1 =
    - 2./3.*(
      + 1. + 2.*DG4 - 2.*X*(2.*G4_X - G5_phi)
      - 8.*F4*pow(X,2) - 2.*X*(G5_X - 12.*X*F5 )*H*phi_prime/a
    );

    P2 =
    - 1./3.*a*(
      + 3.*(rho_tot + p_tot)
      + 2.*X*(G2_X - 2.*G3_phi + 2.*G4_phiphi)
      - 4.*(
        + G4_phi - X*(2.*G3_X - 6.*G4_Xphi + G5_phiphi) + 4.*pow(X,2)*F4_phi
      )*H*phi_prime/a
      + 4.*X*(
        + 5.*(G4_X - G5_phi) + 2.*X*(5.*G4_XX - 3.*G5_Xphi + 20.*F4)
        + 4.*pow(X,2)*(5.*F4_X + 3.*F5_phi)
      )*pow(H,2)
      + 4.*X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,3)*phi_prime/a
    );

    R0 =
    1./3.*(
      + (
        + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
      )/a/H
      + 6.*(
        + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
      )*pow(a,-2)*phi_prime
      + 2.*(
        + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
        + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
      )*pow(a,-2)*pow(H,2)*phi_prime
      + 6.*(
        + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
        + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
      )*H/a
    );

    R1 =
    2.*(
      - (G4_phi - X*(G3_X - 2.*G4_Xphi))/H
      + 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*phi_prime/a
      + X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )*H
    );

    R2 =
    1./3.*(
      - (G2_phi - 2.*X*(G2_Xphi - G3_phiphi))*a/H
      + 2.*(
        + G2_X - 2.*G3_phi - X*(G2_XX - 4.*G3_Xphi + 6.*G4_Xphiphi)
      )*phi_prime
      - 2.*(
        + 6.*G4_phi - 3.*X*(G3_X - G5_phiphi)
        + 6.*pow(X,2)*(G3_XX - 4.*G4_XXphi + G5_Xphiphi - 6.*F4_phi)
        - 24.*pow(X,3)*F4_Xphi
      )*a*H
      + 4.*(
        + 3.*(G4_X - G5_phi) - X*(3.*G4_XX - 4.*G5_Xphi)
        - 2.*pow(X,2)*(3.*G4_XXX - 2.*G5_XXphi + 18.*F4_X + 12.*F5_phi)
        - 12.*pow(X,3)*(F4_XX + F5_Xphi)
      )*pow(H,2)*phi_prime
      + 2.*X*(
        + 3.*G5_X - 4.*X*(2.*G5_XX - 15.*F5)
        - 4.*pow(X,2)*(G5_XXX - 48.*F5_X) + 48.*pow(X,3)*F5_XX
      )*a*pow(H,3)
    );

    pvecback[pba->index_bg_P0_smg] = P0;
    pvecback[pba->index_bg_P1_smg] = P1;
    pvecback[pba->index_bg_P2_smg] = P2;
    pvecback[pba->index_bg_R0_smg] = R0;
    pvecback[pba->index_bg_R1_smg] = R1;
    pvecback[pba->index_bg_R2_smg] = R2;

    class_test_except((P1*R0 - P0*R1) == 0 ,
	       pba->error_message,
         free(pvecback);free(pvecback_B),
               "scalar field mixing with metric has degenerate denominator at a = %e, phi = %e, phi_prime = %e \n with P1 = %e, R0 =%e, P0=%e, R1=%e, \n H=%e, E0=%e, E1=%e, E2=%e, E3=%e",
	       a,phi,phi_prime, P1, R0, P0, R1,
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


    /* Computing density, pressure (also rewritten as shift and current) */

    /* Energy density of the field */
    pvecback[pba->index_bg_rho_smg] =
      - (G2 - 2.*X*(G2_X - G3_phi))/3.
      - 2.*(G4_phi - X*(G3_X - 2.*G4_Xphi))*H*phi_prime/a
      - 2.*(
        + DG4 - X*(4.*G4_X - 3.*G5_phi) - 2.*pow(X,2)*(2.*G4_XX - G5_Xphi + 10.*F4)
        - 8.*pow(X,3)*F4_X
      )*pow(H,2)
      + 2./3.*X*(
        + 5.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,3)*phi_prime/a;

    /* Pressure of the field */
    pvecback[pba->index_bg_p_smg] =
    (
      + G2 - 2.*X*(G3_phi - 2.*G4_phiphi)
      + 2.*(
        + 3.*DG4 - X*(2.*G4_X + G5_phi)
        + 2.*pow(X,2)*(4.*G4_XX - 3.*G5_Xphi + 10.*F4)
        + 8.*pow(X,3)*(2.*F4_X + 3.*F5_phi)
      )*pow(H,2)
      + 2.*(
        + G4_phi + X*(G3_X - 6.*G4_Xphi + 2.*G5_phiphi) - 8.*pow(X,2)*F4_phi
      )*H*phi_prime/a
      + 2.*X*(
        + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,3)*phi_prime/a
      + 4.*(
        + DG4 - X*(2.*G4_X - G5_phi) - 4.*pow(X,2)*F4
        - X*(G5_X - 12.*X*F5)*H*phi_prime/a
      )*pvecback[pba->index_bg_H_prime]/a
      + 2.*(
        + G4_phi - X*(G3_X - 2.*G4_Xphi)
        - X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*pow(H,2)
        - 2.*(
          + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/3.;

    /* Current of the field */
    pvecback[pba->index_bg_current_smg] =
      + (G2_X - G3_phi)*phi_prime/a + 6.*X*(G3_X - 2.*G4_Xphi)*H
      + 3.*(
        + 2.*G4_X - G5_phi + 2.*X*(2.*G4_XX - G5_Xphi + 8.*F4)
        + 8.*pow(X,2)*F4_X
      )*pow(H,2)*phi_prime/a
      + 2.*X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,3);

      /* Shift of the field */
      pvecback[pba->index_bg_shift_smg] =
        + G2_phi + 2.*G3_phi*H*phi_prime/a
        + 12.*(G4_phi + 2.*pow(X,2)*F4_phi)*pow(H,2)
        + 2.*(
          + 3.*G5_phi - 2.*X*G5_Xphi - 12.*pow(X,2)*F5_phi
        )*pow(H,3)*phi_prime/a
        + 6.*(G4_phi + G5_phi*H*phi_prime/a)*pvecback[pba->index_bg_H_prime]/a
        + (
          + G3_phi + 6.*G4_Xphi*H*phi_prime/a
          + 3.*(G5_phi + 2.*X*G5_Xphi)*pow(H,2)
        )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg];


    /* Computing alphas at the end (alpha_T, alpha_M depend on phi'') */

    /* Planck mass */
    pvecback[pba->index_bg_delta_M2_smg] =
      + 2.*(DG4 - X*(2.*G4_X - G5_phi) - 4.*F4*pow(X,2))
      - 2.*X*(G5_X - 12.*X*F5)*H*phi_prime/a;

    pvecback[pba->index_bg_M2_smg] = 1. + pvecback[pba->index_bg_delta_M2_smg];

    /* alpha_K kineticity */
    pvecback[pba->index_bg_kineticity_over_phiphi_smg] =
    (
      + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
      + 6.*(
        + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
      )*H*phi_prime/a
      + 2.*(
        + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
        + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
      )*pow(H,3)*phi_prime/a
      + 6.*(
        + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
        + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
      )*pow(H,2)
    )/pvecback[pba->index_bg_M2_smg];

    pvecback[pba->index_bg_kineticity_smg] = 2.*pow(H,-2)*X*pvecback[pba->index_bg_kineticity_over_phiphi_smg];

    /* alpha_B braiding */
    pvecback[pba->index_bg_braiding_over_phi_smg] =
    2.*(
      - G4_phi + X*(G3_X - 2.*G4_Xphi)
      + X*(
        + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,2)
      + 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
      )*H*phi_prime/a
    )/pvecback[pba->index_bg_M2_smg];

    pvecback[pba->index_bg_braiding_smg] = phi_prime/a/H*pvecback[pba->index_bg_braiding_over_phi_smg];

    /* alpha_T: tensor speed excess */
    pvecback[pba->index_bg_tensor_excess_smg] =
    2.*X*(
      + 2*(G4_X - G5_phi + 2.*X*F4)
      + 2*(G5_X - 6.*X*F5)*H*phi_prime/a
      - G5_X*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* alpha_M: Planck mass running */
    pvecback[pba->index_bg_mpl_running_smg] =
    2.*(
      + 2.*X*(
        + G4_X - G5_phi + 2.*X*(G4_XX - G5_Xphi + 4.*F4)
        + 4.*pow(X,2)*(F4_X + 3.*F5_phi)
      )
      + (G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi)*phi_prime/a/H
      + X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*H*phi_prime/a
      - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]*pow(a,-2)*phi_prime/H
      + (
        - 3.*X*G5_X
        + 2.*pow(X,2)*(30.*F5 - G5_XX)
        + 24.*pow(X,3)*F5_X
        - (
          + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
        )*phi_prime/a/H
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]

    )/pvecback[pba->index_bg_M2_smg];

    /* alpha_H: Beyond Horndeski */
    pvecback[pba->index_bg_beyond_horndeski_over_phi_smg] =
    4.*X*(
      F4*phi_prime/a - 6.*X*H*F5
    )*H/pvecback[pba->index_bg_M2_smg];

    pvecback[pba->index_bg_beyond_horndeski_smg] = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg]*phi_prime/a/H;



    /* Computing B functions (intermediate step for perturbations) */

    /* B0_smg */
    pvecback[pba->index_bg_B0_smg] =
    (
      + G4_phi - X*(3.*G3_X - 10.*G4_Xphi + 2.*G5_phiphi) + 8.*pow(X,2)*F4_phi
      - 3.*(
        + G4_X - G5_phi + 2.*X*(G4_XX - 2./3.*G5_Xphi + 4.*F4)
        + 4.*pow(X,2)*(F4_X + F5_phi)
      )*H*phi_prime/a
      - X*(3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X)*pow(H,2)
      - (G2_X/2. - G3_phi + G4_phiphi)*phi_prime/a/H
    )/pvecback[pba->index_bg_M2_smg];

    /* B1_smg */
    pvecback[pba->index_bg_B1_smg] =
    6.*(
      + (
        + G4_phi + X*(3.*G3_X - 16.*G4_Xphi + 6.*G5_phiphi)
        + 2.*pow(X,2)*(G3_XX - 20.*F4_phi - 6.*G4_XXphi + 2.*G5_Xphiphi)
        - 16.*pow(X,3)*F4_Xphi
      )
      + X*(
        + 3.*G5_X + 12.*X*(G5_XX - 15.*F5)
        + 4.*pow(X,2)*(G5_XXX - 60.*F5_X) - 48.*pow(X,3)*F5_XX
      )*pow(H,2)
      + (
        + (G2_X - 2.*G3_phi + 4.*G4_phiphi)/2. - X*(G3_Xphi - 2.*G4_Xphiphi)
      )*phi_prime/a/H
      + (
        + G4_X - G5_phi + X*(40.*F4 + 14.*G4_XX + (-13.)*G5_Xphi)
        + 2.*pow(X,2)*(34.*F4_X + 4.*G4_XXX + 36.*F5_phi + (-3.)*G5_XXphi)
        + 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
      )*H*phi_prime/a
      - 2.*(
        + X*(
          + 3.*G5_X + 2.*X*(G5_XX - 30.*F5) - 24.*pow(X,2)*F5_X
        )
        + (
          + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
        )*phi_prime/H/a
      )/a*pvecback[pba->index_bg_H_prime]
      + (
        - 2.*(G4_X - G5_phi)
        - 2.*X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
        - 4.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X)
        - 16.*pow(X,3)*F4_XX
        - (G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi))*phi_prime/a/H
        - (
          + 3.*G5_X - X*(120.*F5 - 7.*G5_XX)
          - 2.*pow(X,2)*(66.*F5_X - G5_XXX) - 24.*pow(X,3)*F5_XX
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B2_smg */
    pvecback[pba->index_bg_B2_smg] =
    6.*(
      + (G2_phi/2. - X*(G3_phiphi - 2.*G4_phiphiphi))*pow(H,-2)
      + 3.*G4_phi - X*(2.*G4_Xphi + G5_phiphi)
      + 2.*pow(X,2)*(4.*G4_XXphi - 3.*G5_Xphiphi + 10.*F4_phi)
      + 8.*pow(X,3)*(2.*F4_Xphi + 3.*F5_phiphi)
      + X*(
        + G5_Xphi + 2.*X*(G5_XXphi - 18.*F5_phi) - 24.*pow(X,2)*F5_Xphi
      )*H*phi_prime/a
      + (
        + G4_phiphi + X*(G3_Xphi - 6.*G4_Xphiphi + 2.*G5_phiphiphi)
        - 8.*pow(X,2)*F4_phiphi
      )*phi_prime/a/H
      - 2.*(
        + X*(G5_Xphi - 12.*X*F5_phi)*phi_prime/a
        - (G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi)/H
      )*pvecback[pba->index_bg_H_prime]/a/H
      + (
        - X*(3.*G5_Xphi + 2.*X*(G5_XXphi - 30.*F5_phi) - 24.*pow(X,2)*F5_Xphi)
        + (G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi))*pow(H,-2)
        - 2.*(
          + G4_Xphi - G5_phiphi + X*(2.*G4_XXphi - G5_Xphiphi + 8.*F4_phi)
          + 4.*pow(X,2)*F4_Xphi
        )*phi_prime/a/H
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B3_smg */
    pvecback[pba->index_bg_B3_smg] =
    4.*(
      + G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi
      + X*(
        + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,2)
      + 2.*X*(
        + G4_XX - G5_Xphi + 2.*F4 + 2.*X*(F4_X + 3.*F5_phi)
      )*H*phi_prime/a
      - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]/a
      - (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 6.*F4) + 4.*pow(X,2)*F4_X
        + (G5_X + X*(G5_XX - 24.*F5) - 12.*pow(X,2)*F5_X)*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B4_smg */
    pvecback[pba->index_bg_B4_smg] =
    2.*(
      + G4_phi - X*(2.*G4_Xphi - G5_phiphi) - 4.*pow(X,2)*F4_phi
      + 2.*X*(
        + 2.*F4 + G4_XX - G5_Xphi + 2.*X*(F4_X + 3.*F5_phi)
      )*H*phi_prime/a
      + X*(
        + G5_X + 2.*X*(G5_XX - 18.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,2)
      - X*(G5_X - 12.*X*F5)*pvecback[pba->index_bg_H_prime]/a
      - (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 6.*F4) + 4.*pow(X,2)*F4_X
        + (
          + G5_X - X*(24.*F5 - G5_XX) - 12.*pow(X,2)*F5_X
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B5_smg */
    pvecback[pba->index_bg_B5_smg] =
    (
      - 3.*G4_phi + X*(3.*G3_X - 4.*G4_Xphi - 2.*G5_phiphi)
      - 2.*pow(X,2)*(G3_XX - 6.*G4_XXphi + 2.*G5_Xphiphi - 12.*F4_phi) + 16.*pow(X,3)*F4_Xphi
      + (
        + G2_X - 2.*G3_phi + 2.*X*(G3_Xphi - 2.*G4_Xphiphi)
      )*phi_prime/a/H/2.
      + (
        + 5.*(G4_X - G5_phi) - X*(2.*G4_XX - 5.*G5_Xphi - 8.*F4)
        - 2.*pow(X,2)*( + 4.*G4_XXX - 3.*G5_XXphi + 22.*F4_X + 24.*F5_phi)
        - 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
      )*H*phi_prime/a
      + X*(
        + 3.*G5_X - 4.*X*(2.*G5_XX - 15.*F5)
        - 4.*pow(X,2)*(G5_XXX - 48.*F5_X) + 48.*pow(X,3)*F5_XX
      )*pow(H,2)
      + 2.*(
        + X*(3.*G5_X - 2.*X*(30.*F5 - G5_XX) - 24.*pow(X,2)*F5_X)
        + (
          + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 8.*F4) + 4.*pow(X,2)*F4_X
        )*phi_prime/a/H
      )*pvecback[pba->index_bg_H_prime]/a
      + (
        + 2.*(
          + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
          + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
        )
        + (
          + 3.*G5_X + X*(7.*G5_XX - 120.*F5)
          + 2.*pow(X,2)*(G5_XXX - 66.*F5_X) - 24.*pow(X,3)*F5_XX
        )*H*phi_prime/a
        + (
          + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
        )*phi_prime/a/H
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B6_smg */
    pvecback[pba->index_bg_B6_smg] =
    4.*(
      + G4_phi - X*(2.*G4_Xphi - G5_phiphi)
      + 2.*X*(G4_XX - G5_Xphi)*H*phi_prime/a
      + X*(G5_X + 2.*X*G5_XX)*pow(H,2)
      - X*G5_X*pvecback[pba->index_bg_H_prime]/a
      - (
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi)
        + (G5_X + X*G5_XX)*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B7_smg */
    pvecback[pba->index_bg_B7_smg] =
    - 2.*(
      + G2_X - 2.*G3_phi
      - X*(G2_XX - 8.*G3_Xphi + 18.*G4_Xphiphi)
      - 2.*pow(X,2)*(G2_XXX - 4.*G3_XXphi + 6.*G4_XXphiphi)
      + (
        + G2_Xphi - 2.*G3_phiphi + 2.*X*(G2_XXphi - G3_Xphiphi)
      )*phi_prime/a/H/2.
      + 6.*(
        + G4_X - G5_phi - X*(G4_XX - 2.*G5_Xphi)
        - 4./3.*pow(X,2)*(9.*G4_XXX - 7.*G5_XXphi + 45.*F4_X + 30.*F5_phi)
        - 4./3.*pow(X,3)*(3.*G4_XXXX - 2.*G5_XXXphi + 39.*F4_XX + 33.*F5_Xphi)
        - 8.*pow(X,4)*(F4_XXX + F5_XXphi)
      )*pow(H,2)
      + (
        + 3.*G5_X - X*(13.*G5_XX - 120.*F5) - 4.*pow(X,2)*(5.*G5_XXX - 159.*F5_X)
        - 4.*pow(X,3)*(G5_XXXX - 96.*F5_XX) + 48.*pow(X,4)*F5_XXX
      )*pow(H,3)*phi_prime/a
      + 3.*(
        + G3_X - 2.*G4_Xphi - G5_phiphi
        - X*(3.*G3_XX - 16.*G4_XXphi + 5.*G5_Xphiphi - 24.*F4_phi)
        - 2.*pow(X,2)*(G3_XXX - 4.*G4_XXXphi + G5_XXphiphi - 18.*F4_Xphi)
        + 8.*pow(X,3)*F4_XXphi
      )*H*phi_prime/a
      + 6.*(
        + G4_X - G5_phi + X*(8.*G4_XX - 5.*G5_Xphi + 24.*F4)
        + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 18.*F4_X) + 8.*pow(X,3)*F4_XX
        - (
          - 3.*G5_X - X*(7.*G5_XX - 120.*F5) - 2.*pow(X,2)*(G5_XXX - 66.*F5_X)
          + 24.*pow(X,3)*F5_XX
        )*H*phi_prime/a/2.
        + (
          + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
        )*phi_prime/a/H/2.
      )*pvecback[pba->index_bg_H_prime]/a
      + 3.*(
        + (G3_X - 3.*G4_Xphi)
        + X*(5.*G3_XX - 12.*G4_XXphi)
        + 2.*pow(X,2)*(G3_XXX - 2.*G4_XXXphi)
        + (
          + G5_X + 3.*X*(3.*G5_XX - 40.*F5) + 4.*pow(X,2)*(2.*G5_XXX - 75.*F5_X)
          + 4./3.*pow(X,3)*(G5_XXXX - 108.*F5_XX) - 16.*pow(X,4)*F5_XXX
        )*pow(H,2)
        + (
          + G2_XX/2. - 2./3.*G3_Xphi + X*(G2_XXX - G3_XXphi)/3.
        )*phi_prime/a/H
        + (
          + 3.*(3.*G4_XX - 2.*G5_Xphi + 8.*F4)
          + X*(16.*G4_XXX - 9.*G5_XXphi + 96.*F4_X)
          + 2.*pow(X,2)*(2.*G4_XXXX - G5_XXXphi + 30.*F4_XX)
          + 8.*pow(X,3)*F4_XXX
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B8_smg */
    pvecback[pba->index_bg_B8_smg] =
    - 2.*(
      + G2_X/2. - G3_phi + X*(G3_Xphi - 2.*G4_Xphiphi)
      + (
        + 3.*(G4_X - G5_phi) - X*(2.*G4_XX - 3.*G5_Xphi + 8.*F4)
        - 2.*pow(X,2)*(4.*G4_XXX - 3.*G5_XXphi + 26.*F4_X + 24.*F5_phi)
        - 8.*pow(X,3)*(2.*F4_XX + 3.*F5_Xphi)
      )*pow(H,2)
      + (
        + G3_X - 3.*G4_Xphi - X*(G3_XX - 6.*G4_XXphi + 2.*G5_Xphiphi - 12.*F4_phi)
        + 8.*pow(X,2)*F4_Xphi
      )*H*phi_prime/a
      + (
        + G5_X - 3.*X*(G5_XX - 20.*F5) - 2.*pow(X,2)*(G5_XXX - 54.*F5_X)
        + 24.*pow(X,3)*F5_XX
      )*pow(H,3)*phi_prime/a
      + 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 10.*F4) + 4.*pow(X,2)*F4_X
        + (
          + G5_X + X*(G5_XX - 36.*F5) - 12.*pow(X,2)*F5_X
        )*H*phi_prime/a
      )*pvecback[pba->index_bg_H_prime]/a
      + (
        + G3_X - 3.*G4_Xphi + X*(G3_XX - 2.*G4_XXphi)
        + (
          + G5_X + 5.*X*(G5_XX - 24.*F5) + 2.*pow(X,2)*(G5_XXX - 66.*F5_X)
          - 24.*pow(X,3)*F5_XX
        )*pow(H,2)
        + 2.*(
          + 3.*G4_XX - 2.*G5_Xphi + 12.*F4 + X*(2.*G4_XXX - G5_XXphi + 18.*F4_X)
          + 4.*pow(X,2)*F4_XX
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B9_smg */
    pvecback[pba->index_bg_B9_smg] =
    2.*(
      + 6.*G4_phiphi
      - 3.*X*(G3_Xphi - G5_phiphiphi)
      + 6.*pow(X,2)*(G3_XXphi - 4.*G4_XXphiphi + G5_Xphiphiphi - 6.*F4_phiphi)
      - 24.*pow(X,3)*F4_Xphiphi
      + (
        + G2_phiphi - (G2_Xphiphi - G3_phiphiphi)*2.*X
      )*pow(H,-2)/2.
      - (
        + G2_Xphi - 2.*G3_phiphi - X*(G2_XXphi - 4.*G3_Xphiphi + 6.*G4_Xphiphiphi)
      )*phi_prime/a/H
      - 2.*(
        + 3.*(G4_Xphi - G5_phiphi) - X*(3.*G4_XXphi - 4.*G5_Xphiphi)
        - 2.*pow(X,2)*(3.*G4_XXXphi - 2.*G5_XXphiphi + 18.*F4_Xphi + 12.*F5_phiphi)
        - 12.*pow(X,3)*(F4_XXphi + F5_Xphiphi)
      )*H*phi_prime/a
      - X*(
        + 3.*G5_Xphi - 4.*X*(2.*G5_XXphi - 15.*F5_phi)
        - 4.*pow(X,2)*(G5_XXXphi - 48.*F5_Xphi) + 48.*pow(X,3)*F5_XXphi
      )*pow(H,2)
      - 3.*(
        + X*(3.*G5_Xphi + 2.*X*(G5_XXphi - 30.*F5_phi) - 24.*pow(X,2)*F5_Xphi)
        - (G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi))*pow(H,-2)
        + 2.*(
          + G4_Xphi - G5_phiphi + X*(8.*F4_phi + 2.*G4_XXphi - G5_Xphiphi)
          + 4.*pow(X,2)*F4_Xphi
        )*phi_prime/a/H
      )*pvecback[pba->index_bg_H_prime]/a
      - (
        + 3.*(G4_Xphi - G5_phiphi)
        + 3.*X*(8.*G4_XXphi - 5.*G5_Xphiphi + 24.*F4_phi)
        + 6.*pow(X,2)*(2.*G4_XXXphi - G5_XXphiphi + 18.*F4_Xphi)
        + 24.*pow(X,3)*F4_XXphi
        + (
          + G2_Xphi - 2.*G3_phiphi + 2.*X*(G2_XXphi - G3_Xphiphi)
        )*pow(H,-2)/2.
        + 3.*(
          + G3_Xphi - 3.*G4_Xphiphi + X*(G3_XXphi - 2.*G4_XXphiphi)
        )*phi_prime/a/H
        + (
          + 3.*G5_Xphi + X*(7.*G5_XXphi - 120.*F5_phi)
          - 2.*pow(X,2)*(66.*F5_Xphi - G5_XXXphi) - 24.*pow(X,3)*F5_XXphi
        )*H*phi_prime/a
      )*pow(a,-2)*pvecback[pba->index_bg_phi_prime_prime_smg]
    )/pvecback[pba->index_bg_M2_smg];

    /* B10_smg */
    pvecback[pba->index_bg_B10_smg] =
    (
      + 3.*(
        + G4_phi - X*(3.*G3_X - 8.*G4_Xphi) - 2.*pow(X,2)*(G3_XX - 2.*G4_XXphi)
      )
      - (
        + G2_X - 2.*G3_phi + 2.*X*(G2_XX - G3_Xphi)
      )*phi_prime/a/H/2.

      - X*(
        + 15.*G5_X + 20.*X*(G5_XX - 21.*F5) + 4.*pow(X,2)*(G5_XXX - 84.*F5_X)
        - 48.*pow(X,3)*F5_XX
      )*pow(H,2)
      - 3.*(
        + 3.*(G4_X - G5_phi) + X*(12.*G4_XX - 7.*G5_Xphi + 40.*F4)
        + 2.*pow(X,2)*(2.*G4_XXX - G5_XXphi + 22.*F4_X) + 8.*pow(X,3)*F4_XX
      )*H*phi_prime/a
    )/pvecback[pba->index_bg_M2_smg];

    /* B11_smg */
    pvecback[pba->index_bg_B11_smg] =
    (
      + G4_phi
      - X*(G3_X - 2.*G4_Xphi)
      - X*(
        + 3.*G5_X + 2.*X*(G5_XX - 42.*F5) - 24.*pow(X,2)*F5_X
      )*pow(H,2)
      - 2.*(
        + G4_X - G5_phi + X*(2.*G4_XX - G5_Xphi + 10.*F4) + 4.*pow(X,2)*F4_X
      )*H*phi_prime/a
    )/pvecback[pba->index_bg_M2_smg];

    /* B12_smg */
    pvecback[pba->index_bg_B12_smg] =
    (
      + 3.*G4_phi
      - 3.*X*(4.*G4_Xphi - 3.*G5_phiphi)
      - 6.*pow(X,2)*(10.*F4_phi + 2.*G4_XXphi - G5_Xphiphi)
      - 24.*pow(X,3)*F4_Xphi
      + (
        + G2_phi/2. - X*(G2_Xphi - G3_phiphi)
      )*pow(H,-2)
      - X*(
        + 5.*G5_Xphi + 2.*X*(G5_XXphi - 42.*F5_phi) - 24.*pow(X,2)*F5_Xphi
      )*H*phi_prime/a
      + 3.*(
        + G4_phiphi - X*(G3_Xphi - 2.*G4_Xphiphi)
      )*phi_prime/a/H
    )/pvecback[pba->index_bg_M2_smg];

  }// end of if pba->field_evolution_smg
  else{

    double a, delta_M_pl;
    double rho_tot, p_tot;
    double Omega_smg;

    a = pvecback_B[pba->index_bi_a];
    delta_M_pl = pvecback_B[pba->index_bi_delta_M_pl_smg];

    rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
    p_tot = pvecback[pba->index_bg_p_tot_wo_smg];

    //initialize the values to the defaults
    pvecback[pba->index_bg_kineticity_smg] = 0;
    pvecback[pba->index_bg_braiding_smg] = 0.;
    pvecback[pba->index_bg_tensor_excess_smg] = 0.;
    pvecback[pba->index_bg_beyond_horndeski_smg] = 0.;
    pvecback[pba->index_bg_M2_smg] = 1.;
    pvecback[pba->index_bg_delta_M2_smg] = 0.;
    pvecback[pba->index_bg_mpl_running_smg] = 0.;

     if (pba->expansion_model_smg == lcdm){

      double Omega_const_smg = pba->parameters_smg[0];

      pvecback[pba->index_bg_rho_smg] = Omega_const_smg*pow(pba->H0,2);
      pvecback[pba->index_bg_p_smg] = -Omega_const_smg*pow(pba->H0,2);
    }

    if (pba->expansion_model_smg == wowa){

      double Omega_const_smg = pba->parameters_smg[0];
      double w0 = pba->parameters_smg[1];
      double wa = pba->parameters_smg[2];

      pvecback[pba->index_bg_rho_smg] = Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1. + w0 + wa)) * exp(3.*wa*(a-1.));
      pvecback[pba->index_bg_p_smg] = (w0+(1-a)*wa) * Omega_const_smg * pow(pba->H0,2)/pow(a,3.*(1.+w0+wa)) * exp(3.*wa*(a-1.));
    }
    if (pba->expansion_model_smg == wowa_w){

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

    }//ILSWEDE
    if (pba->expansion_model_smg == wede){

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

    rho_tot += pvecback[pba->index_bg_rho_smg];
    p_tot += pvecback[pba->index_bg_p_smg];


    Omega_smg = pvecback[pba->index_bg_rho_smg]/rho_tot; //used for some parameterizations


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

    }//end of look pover models


    pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);
    /** - compute derivative of H with respect to conformal time */
    pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

    //add friction term
    if (pba->hubble_evolution == _TRUE_ && pba->initial_conditions_set_smg == _TRUE_){
      pvecback[pba->index_bg_H] = pvecback_B[pba->index_bi_H];
      /** - compute derivative of H with respect to conformal time */
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(pvecback_B[pba->index_bi_H]*pvecback_B[pba->index_bi_H] - rho_tot - pba->K/a/a);
    }

    // Compute time derivative of rho_smg
    if (pba->rho_evolution_smg == _TRUE_){
        pvecback[pba->index_bg_rho_prime_smg] = -3.*a*pvecback[pba->index_bg_H]*(1.+pvecback[pba->index_bg_w_smg])*pvecback[pba->index_bg_rho_smg];
    }

}//end of parameterized mode

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
  pvecback[pba->index_bg_p_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_prime_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] = 0.;
  pvecback[pba->index_bg_kinetic_D_over_phiphi_prime_smg] = 0.;
  pvecback[pba->index_bg_cs2num_smg] = 0.;
  pvecback[pba->index_bg_cs2num_prime_smg] = 0.;
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
  if ( (pba->skip_stability_tests_smg == _FALSE_) && (pba->parameters_tuned_smg == _TRUE_) && (pba->Omega_smg_debug == 0) ){
      double a = pvecback_B[pba->index_bi_a];
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

  return _SUCCESS_;

}
