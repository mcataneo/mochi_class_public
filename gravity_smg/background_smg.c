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
				 double a,
				 double * pvecback_B,
				 short return_format,
				 double * pvecback,
         double * ptr_rho_tot,
				 double * ptr_p_tot,
				 double * ptr_rho_de
				 ){

	 /* Scalar field */
   /** if (has_smg) do all the mess to compute the Hubble rate, etc... */

	 /** - save rho_tot and p_tot without smg (it is important that smg is the last thing to be computed) */
   pvecback[pba->index_bg_rho_tot_wo_smg] = *ptr_rho_tot;
   pvecback[pba->index_bg_p_tot_wo_smg] = *ptr_p_tot;
   //NOTE: add the field energy and pressure after the debug has been added

  // scalar field + curvature not yet implemented
  class_test(pba->K !=0 ,
	     pba->error_message,
	     "has_smg with curvature K = %e not yet implemented",pba->K);

  if (pba->field_evolution_smg == _TRUE_) {

    /* declare variables and set defaults to zero */
    double phi, phi_prime, H;
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

    if (pba->hubble_evolution == _TRUE_)
      H = exp(pvecback_B[pba->index_bi_logH]);

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

    double delta_M_pl;
    double rho_tot, p_tot;
    double Omega_smg;

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
      pvecback[pba->index_bg_H] = exp(pvecback_B[pba->index_bi_logH]);
      /** - compute derivative of H with respect to conformal time */
      pvecback[pba->index_bg_H_prime] += - a*pba->hubble_friction*(pvecback[pba->index_bg_H]*pvecback[pba->index_bg_H] - rho_tot - pba->K/a/a);
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
	if (pba->rho_evolution_smg == _FALSE_) {
		pvecback[pba->index_bg_w_smg] = pvecback[pba->index_bg_p_smg] / pvecback[pba->index_bg_rho_smg];
	}

  return _SUCCESS_;

}


/* - indices for scalar field (modified gravity) */
int hi_class_define_indices_bg(
				 struct background *pba,
				 int * index_bg
			 ) {

	class_define_index(pba->index_bg_phi_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_phi_prime_prime_smg,pba->field_evolution_smg,*index_bg,1);

	class_define_index(pba->index_bg_rho_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_rho_prime_smg,pba->rho_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_w_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_current_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_shift_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_delta_M2_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kineticity_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_tensor_excess_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_mpl_running_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_kineticity_over_phiphi_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_braiding_over_phi_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_beyond_horndeski_over_phi_prime_smg,_TRUE_,*index_bg,1);

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
	class_define_index(pba->index_bg_kinetic_D_over_phiphi_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_kinetic_D_over_phiphi_prime_smg,_TRUE_,*index_bg,1);
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

	class_define_index(pba->index_bg_E0_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_E1_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_E2_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_E3_smg,pba->field_evolution_smg,*index_bg,1);

	class_define_index(pba->index_bg_P0_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_P1_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_P2_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_R0_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_R1_smg,pba->field_evolution_smg,*index_bg,1);
	class_define_index(pba->index_bg_R2_smg,pba->field_evolution_smg,*index_bg,1);

	class_define_index(pba->index_bg_rho_tot_wo_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_H_prime_prime,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_tot_wo_prime_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_p_prime_smg,_TRUE_,*index_bg,1);

	class_define_index(pba->index_bg_G_eff_smg,_TRUE_,*index_bg,1);
	class_define_index(pba->index_bg_slip_eff_smg,_TRUE_,*index_bg,1);

  return _SUCCESS_;
}


/* - indices for scalar field (modified gravity) */
int hi_class_define_indices_bi(
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


/* -> write in the table (overwrite the alpha time derivatives, which were set to nan in background_functions)
 * direction of copy: add the corresponding indices to the coordinates
 * thing to be copied: the address (&) to the element of pvecback corresponding to the index we want
 * size: just a single double number
 * -> repeat for all necessary quantities
 */
int background_derivs_alphas_smg(
				struct background *pba,
				double * pvecback,
				double * pvecback_derivs,
				int i
			) {

	/* needed for growing table */
  void * memcopy_result;
	/* necessary for calling array_interpolate(), but never used */
	int last_index;

	//Need to update pvecback
	class_call(background_at_tau(pba,
			 pba->tau_table[i],
			 long_info,
			 inter_normal,
			 &last_index, //should be no problem to use the same one as for the derivatives
			 pvecback),
	pba->error_message,
	pba->error_message);

	/*NOTE: here we compute the derivatives of quantities coputed in background_gravity_functions_smg during the integration.
	 * for quantities that depend on these derivatives (e.g. the gamma_i functions determining the effective mass)
	 * there is an additional loop at the end of background_solve
	 */

	// Kineticity'
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kineticity_prime_smg,
				&pvecback_derivs[pba->index_bg_kineticity_smg],
				1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kineticity_prime_smg,
					 pba->error_message,
					 "cannot copy data back to pba->background_table");

   //Braiding'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_braiding_prime_smg,
		      &pvecback_derivs[pba->index_bg_braiding_smg],
		      1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_braiding_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

   //Planck mass run rate'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_prime_smg,
		      &pvecback_derivs[pba->index_bg_mpl_running_smg],
		      1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

   //Tensor excess'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_tensor_excess_prime_smg,
		      &pvecback_derivs[pba->index_bg_tensor_excess_smg],
		      1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_tensor_excess_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

  //Beyond horndeski'
  memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_beyond_horndeski_prime_smg,
		      &pvecback_derivs[pba->index_bg_beyond_horndeski_smg],
		      1*sizeof(double));
  class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_beyond_horndeski_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

   //Braiding_over_phi'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_braiding_over_phi_prime_smg,
		      &pvecback_derivs[pba->index_bg_braiding_over_phi_smg],
		      1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_braiding_over_phi_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

  //Beyond_horndeski_over_phi'
  memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_beyond_horndeski_over_phi_prime_smg,
		      &pvecback_derivs[pba->index_bg_beyond_horndeski_over_phi_smg],
		      1*sizeof(double));
  class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_beyond_horndeski_over_phi_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

   //H''
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_H_prime_prime,
      		&pvecback_derivs[pba->index_bg_H_prime],
      		1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_H_prime_prime,
           pba->error_message,
           "cannot copy data back to pba->background_table");

   // p_tot_wo_smg'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_p_tot_wo_prime_smg,
         &pvecback_derivs[pba->index_bg_p_tot_wo_smg],
         1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_p_tot_wo_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

   // p_smg'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_p_prime_smg,
         &pvecback_derivs[pba->index_bg_p_smg],
         1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_p_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

	// Planck's mass running
	// Only need to compute it if neither self consistent field evolution nor evolving M_pl in terms of alpha_M
	// check equation 3.3 of Bellini & Sawicki 2014

	if (pba->field_evolution_smg == _FALSE_ && pba->M_pl_evolution_smg == _FALSE_){

		double alpha_M = pvecback_derivs[pba->index_bg_delta_M2_smg]/pvecback[pba->index_bg_M2_smg]/pvecback[pba->index_bg_a]/pvecback[pba->index_bg_H];

		memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_smg,
			&alpha_M, //write using the address
			1*sizeof(double));
		class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_mpl_running_smg,
		 pba->error_message,
		 "cannot copy data back to pba->background_table");
	}

	if(pba->background_verbose > 15 && fabs(1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H])>1e-8)
printf("a = %g, (delta H')/H' = %g \n", pvecback[pba->index_bg_a], 1. - pvecback[pba->index_bg_H_prime]/pvecback_derivs[pba->index_bg_H]);

	return _SUCCESS_;
}

int background_gravity_functions_A_C_smg(
        struct background *pba,
        double * pvecback,
        double * pvecback_derivs,
        int i
			) {

	/* needed for growing table */
  void * memcopy_result;

	double a = pvecback[pba->index_bg_a];
	double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
	double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
	double p_smg = pvecback[pba->index_bg_p_smg];
	double rho_smg = pvecback[pba->index_bg_rho_smg];

	//TODO: clean up this!!
	// speed of sound, wanted for output

	double H = pvecback[pba->index_bg_H];
	double H_p = pvecback[pba->index_bg_H_prime];

	double M2 = pvecback[pba->index_bg_M2_smg];
	double DelM2 = pvecback[pba->index_bg_delta_M2_smg];
	double kin = pvecback[pba->index_bg_kineticity_smg];
	double bra = pvecback[pba->index_bg_braiding_smg];
	double run = pvecback[pba->index_bg_mpl_running_smg];
	double ten = pvecback[pba->index_bg_tensor_excess_smg];
	double beh = pvecback[pba->index_bg_beyond_horndeski_smg];
	double dM2 = pvecback[pba->index_bg_delta_M2_smg];

	double kin_ss = pvecback[pba->index_bg_kineticity_over_phiphi_smg];
	double bra_s = pvecback[pba->index_bg_braiding_over_phi_smg];
	double beh_s = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg];

	double B0 = pvecback[pba->index_bg_B0_smg];
	double B1 = pvecback[pba->index_bg_B1_smg];
	double B2 = pvecback[pba->index_bg_B2_smg];
	double B3 = pvecback[pba->index_bg_B3_smg];
	double B4 = pvecback[pba->index_bg_B4_smg];
	double B5 = pvecback[pba->index_bg_B5_smg];
	double B6 = pvecback[pba->index_bg_B6_smg];
	double B7 = pvecback[pba->index_bg_B7_smg];
	double B8 = pvecback[pba->index_bg_B8_smg];
	double B9 = pvecback[pba->index_bg_B9_smg];
	double B10 = pvecback[pba->index_bg_B10_smg];
	double B11 = pvecback[pba->index_bg_B11_smg];
	double B12 = pvecback[pba->index_bg_B12_smg];

	//need to update the time derivatives of the interesting functions

	double kin_p = pvecback_derivs[pba->index_bg_kineticity_smg];
	double bra_p = pvecback_derivs[pba->index_bg_braiding_smg];
	double run_p = pvecback_derivs[pba->index_bg_mpl_running_smg];
	double ten_p = pvecback_derivs[pba->index_bg_tensor_excess_smg];
	double beh_p = pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
	double p_tot_p = pvecback_derivs[pba->index_bg_p_tot_wo_smg];
	double p_smg_p = pvecback_derivs[pba->index_bg_p_smg];


  // kinetic term D
  pvecback[pba->index_bg_kinetic_D_smg] = kin + 3./2.*pow(bra,2);
  memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_smg,
                          &pvecback[pba->index_bg_kinetic_D_smg], 1*sizeof(double));
  class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_smg,
             pba->error_message, "cannot copy data back to pba->background_table");

	// A0
	pvecback[pba->index_bg_A0_smg] =
	1./2.*(
	 + bra - 3.*(rho_smg + p_smg + (rho_tot + p_tot)*DelM2/M2)*pow(H,-2)
	);
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A0_smg,
	                      &pvecback[pba->index_bg_A0_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A0_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A1
	pvecback[pba->index_bg_A1_smg] =
	+ (1. + ten)*kin
	- 3.*(beh*(1. + run) + run - ten + beh_p/a/H)*bra;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A1_smg,
	                      &pvecback[pba->index_bg_A1_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A1_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A2
	pvecback[pba->index_bg_A2_smg] =
	- (kin + 3./2.*pow(bra,2))*(2. + run)
	- 9./4.*bra*(
	 + (2. - bra)*(rho_smg+p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
	)*pow(H,-2)
	- 3./2.*bra*bra_p/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A2_smg,
	                      &pvecback[pba->index_bg_A2_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A2_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A3
	pvecback[pba->index_bg_A3_smg] = bra*beh;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A3_smg,
	                      &pvecback[pba->index_bg_A3_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A3_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A4
	pvecback[pba->index_bg_A4_smg] =
	9./2.*kin*(
	 + (2. - bra)*(rho_smg+p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
	)*pow(H,-2)
	+ 3.*(bra*kin_p - kin*bra_p)/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A4_smg,
	                      &pvecback[pba->index_bg_A4_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A4_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A5
	pvecback[pba->index_bg_A5_smg] = - beh*kin;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A5_smg,
	                      &pvecback[pba->index_bg_A5_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A5_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A6
	pvecback[pba->index_bg_A6_smg] =
	 + 9./4.*(
	   + (2.*kin + 9.*bra)*(2.*DelM2/M2 - bra)
	   + 4.*(kin + 3./2.*pow(bra,2))*run
	 )
	 + 9.*(kin + 9./2.*bra)*rho_smg*pow(H,-2)/M2
	 + 9./2.*(
	   + (kin + 9.*bra)*(2.*DelM2/M2 - bra)
	   + 2.*(kin + 3./2.*pow(bra,2))*run
	 )*pow(H,-2)*p_tot
	 + 81./4.*bra*(
	   + 2.*rho_smg*(p_tot + p_smg)/M2
	   - 2.*(1./M2 - 2. + bra)*p_tot*p_smg
	   + (2. - bra)*pow(p_smg,2)
	   + (2.*DelM2 - bra*M2)*pow(p_tot,2)/M2
	 )*pow(H,-4)
	 + 9./2.*(
	   - 9.*bra*(1./M2 - 2. + bra)
	   + kin*(2. - bra)
	   + 2.*(kin + 3./2.*pow(bra,2))*run
	 )*pow(H,-2)*p_smg
	 + 3.*(
	   + bra*kin_p
	   - (kin - 9./2.*bra - 9./2.*bra*pow(H,-2)*(p_tot + p_smg))*bra_p
	 )/a/H
	 + 9.*(
	   + (kin*DelM2/M2 + 3./2.*pow(bra,2))*p_tot_p
	   + (kin + 3./2.*pow(bra,2))*p_smg_p
	 )*pow(H,-3)/a;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A6_smg,
	                      &pvecback[pba->index_bg_A6_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A6_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A7
	pvecback[pba->index_bg_A7_smg] =
	- 2.*kin*beh
	+ 3.*bra*(bra + 2.*beh)*(1. + run)
	+ 2.*(kin + 3.*bra)*(run - ten)
	+ 9./2.*bra*(
	 + (2. - bra - 2.*beh)*(rho_smg + p_smg)
	 + (2.*DelM2/M2 - bra - 2.*beh)*(rho_tot + p_tot)
	)*pow(H,-2)
	+ 3.*bra*(bra_p + 2.*beh_p)/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A7_smg,
	                      &pvecback[pba->index_bg_A7_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A7_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A8
	pvecback[pba->index_bg_A8_smg] = run - ten - beh;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A8_smg,
	                      &pvecback[pba->index_bg_A8_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A8_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A9
	pvecback[pba->index_bg_A9_smg] =
	+ 3./4.*(
	 + (2. - bra)*(rho_smg + p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2)
	+ 1./2.*bra_p/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A9_smg,
	                      &pvecback[pba->index_bg_A9_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A9_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A10
	pvecback[pba->index_bg_A10_smg] =
	 bra + 2.*run - (2. - bra)*ten + 2.*(1. + run)*beh + 2.*beh_p/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A10_smg,
	                      &pvecback[pba->index_bg_A10_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A10_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A11
	pvecback[pba->index_bg_A11_smg] =
	- (kin + 3./2.*pow(bra,2))*(4. + run)
	+ 3./4.*(
	 + (4.*kin + 6.*bra + 3.*pow(bra,2))*(rho_smg + p_smg)
	 + (4.*kin + 6.*bra*DelM2/M2 + 3.*pow(bra,2))*(rho_tot + p_tot)
	)*pow(H,-2)
	- (kin_p + 3./2.*bra*bra_p)/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A11_smg,
	                      &pvecback[pba->index_bg_A11_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A11_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A12
	pvecback[pba->index_bg_A12_smg] =
	+ kin/2. - 3.*bra*(3./2. - bra) - run*(kin + 3./2.*pow(bra,2))
	- 9./4.*(
	 + (6.*DelM2/M2 + bra*(2./M2 - 7. + 2.*bra))*pow(H,-2)*rho_tot
	 + (6.*DelM2/M2  - 2.*kin + bra*(2./M2 - 5. - 2.*bra))*pow(H,-2)*p_tot
	 + (6. - 5.*bra - 2.*pow(bra,2) - 2.*kin)*pow(H,-2)*p_smg
	 - 6.*(1./M2 - 2. + bra)*pow(H,-4)*p_tot*p_smg
	 + 3.*(2. - bra)*pow(H,-4)*pow(p_smg,2)
	 + 3.*(2.*DelM2/M2 - bra)*(
	   + rho_tot*p_tot + pow(p_tot,2) + rho_tot*p_smg
	 )*pow(H,-4)
	 + (2. - bra)*(
	   + 3. - 2.*bra + 3.*pow(H,-2)*(p_tot + p_smg)
	 )*pow(H,-2)*rho_smg
	 + 2.*bra*pow(H,-3)*p_tot_p/a/M2
	)
	- (
	 + kin_p
	 + 3./2.*(3. + bra + 3.*pow(H,-2)*(p_tot + p_smg))*bra_p
	)/a/H;

	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A12_smg,
	                      &pvecback[pba->index_bg_A12_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A12_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A13
	pvecback[pba->index_bg_A13_smg] =
	- bra - 2.*run + (2. - bra)*ten - (2. + bra + 2.*run)*beh
	- 3./2.*(
	 + (2. - bra - 2.*beh)*(rho_smg + p_smg)*pow(H,-2)
	 + (2.*DelM2/M2 - bra - 2.*beh)*(rho_tot + p_tot)*pow(H,-2)
	)
	- (bra_p + 2.*beh_p)/a/H;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A13_smg,
	                      &pvecback[pba->index_bg_A13_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A13_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A14
	pvecback[pba->index_bg_A14_smg] = - (kin + 3.*bra)/2.;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A14_smg,
	                      &pvecback[pba->index_bg_A14_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A14_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A15
	pvecback[pba->index_bg_A15_smg] = - 1./2.*bra - beh;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A15_smg,
	                      &pvecback[pba->index_bg_A15_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A15_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// A16
	pvecback[pba->index_bg_A16_smg] =
	- 1./2.*(kin + 3.*bra)
	+ 9./4.*(
	 + (2. - bra)*(rho_smg + p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2);
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A16_smg,
	                      &pvecback[pba->index_bg_A16_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A16_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// kinetic term D over phiphi
	pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] =
	+ kin_ss + 3./2.*pow(bra_s,2);
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_over_phiphi_smg,
	                       &pvecback[pba->index_bg_kinetic_D_over_phiphi_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_over_phiphi_smg,
	          pba->error_message, "cannot copy data back to pba->background_table");

	// C0
	pvecback[pba->index_bg_C0_smg] = B0;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C0_smg,
	                      &pvecback[pba->index_bg_C0_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C0_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C1
	pvecback[pba->index_bg_C1_smg] = kin_ss*(1. + ten) - 3./2.*bra_s*B6;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C1_smg,
	                      &pvecback[pba->index_bg_C1_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C1_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C2
	pvecback[pba->index_bg_C2_smg] = - kin_ss*(2. + run) - 3.*bra_s*B5;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C2_smg,
	                      &pvecback[pba->index_bg_C2_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C2_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C3
	pvecback[pba->index_bg_C3_smg] = bra_s*beh_s;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C3_smg,
	                      &pvecback[pba->index_bg_C3_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C3_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C4
	pvecback[pba->index_bg_C4_smg] = kin_ss*B1 - 3.*bra_s*B7;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C4_smg,
	                      &pvecback[pba->index_bg_C4_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C4_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C5
	pvecback[pba->index_bg_C5_smg] = - kin_ss*beh_s;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C5_smg,
	                      &pvecback[pba->index_bg_C5_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C5_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C6
	pvecback[pba->index_bg_C6_smg] = kin_ss*B2 - 3.*bra_s*B9;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C6_smg,
	                      &pvecback[pba->index_bg_C6_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C6_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C7
	pvecback[pba->index_bg_C7_smg] = kin_ss*B3 - 3.*bra_s*B8;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C7_smg,
	                      &pvecback[pba->index_bg_C7_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C7_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C8
	pvecback[pba->index_bg_C8_smg] = B4;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C8_smg,
	                      &pvecback[pba->index_bg_C8_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C8_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C9
	pvecback[pba->index_bg_C9_smg] = - bra_s - bra_s*run/2. + B5;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C9_smg,
	                      &pvecback[pba->index_bg_C9_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C9_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C10
	pvecback[pba->index_bg_C10_smg] = bra_s*(1. + ten) + B6;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C10_smg,
	                      &pvecback[pba->index_bg_C10_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C10_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C11
	pvecback[pba->index_bg_C11_smg] = bra_s*B1/2. + B7;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C11_smg,
	                      &pvecback[pba->index_bg_C11_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C11_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C12
	pvecback[pba->index_bg_C12_smg] = bra_s*B2/2. + B9;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C12_smg,
	                      &pvecback[pba->index_bg_C12_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C12_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C13
	pvecback[pba->index_bg_C13_smg] = bra_s*B3/2. + B8;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C13_smg,
	                      &pvecback[pba->index_bg_C13_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C13_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C14
	pvecback[pba->index_bg_C14_smg] = B10;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C14_smg,
	                      &pvecback[pba->index_bg_C14_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C14_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C15
	pvecback[pba->index_bg_C15_smg] = B11;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C15_smg,
	                      &pvecback[pba->index_bg_C15_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C15_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	// C16
	pvecback[pba->index_bg_C16_smg] = B12;
	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C16_smg,
	                      &pvecback[pba->index_bg_C16_smg], 1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C16_smg,
	         pba->error_message, "cannot copy data back to pba->background_table");

	pvecback[pba->index_bg_lambda_1_smg] = (run + (-1.)*ten)*(-3.)*bra + (1. + ten)*kin;

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_1_smg,
	     &pvecback[pba->index_bg_lambda_1_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_1_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_2_smg] = (- 2.*dM2 + bra*M2)*(rho_tot + p_tot)*(-3.)/2.*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/2.*pow(H,-2) + pow(H,-1)*bra_p*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_smg,
	     &pvecback[pba->index_bg_lambda_2_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_3_smg] = (2. + run)*(-1.)/2.*pvecback[pba->index_bg_kinetic_D_smg] + (-3.)/4.*bra*pvecback[pba->index_bg_lambda_2_smg];

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_3_smg,
	     &pvecback[pba->index_bg_lambda_3_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_3_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_4_smg] = kin*pvecback[pba->index_bg_lambda_2_smg] + (2.*kin*bra_p + (-1.)*bra*kin_p)*(-1.)*pow(H,-1)*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_4_smg,
	     &pvecback[pba->index_bg_lambda_4_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_4_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_5_smg] = (bra + 2.*run + (-2.)*ten + bra*ten)*3./2.*bra + (run + (-1.)*ten)*pvecback[pba->index_bg_kinetic_D_smg] + 3./2.*bra*pvecback[pba->index_bg_lambda_2_smg];

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_5_smg,
	     &pvecback[pba->index_bg_lambda_5_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_5_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_6_smg] = 3./2.*(((9./2.*bra + kin)*dM2*pow(M2,-1) + (-9.)/4.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*pow(rho_tot,2) + ((9.*bra + kin)*dM2*pow(M2,-1) + (-9.)/2.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_tot*p_tot + 9./2.*bra*(dM2 - M2*bra/2.)*pow(M2,-1)*pow(p_tot,2) + (kin*dM2*pow(M2,-1) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*(rho_tot + p_tot)*rho_smg + ((kin - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_smg + ((9.*bra + kin)*(2. - bra)/2. + pvecback[pba->index_bg_kinetic_D_smg]*run - 9./2.*bra*pow(M2,-1))*rho_tot + 9.*bra*(1. - bra/2. - pow(M2,-1)/2.)*p_tot)*(rho_smg + p_smg) + 9./2.*bra*(1. - bra/2.)*pow(rho_smg + p_smg,2))*pow(H,-4) + (((9.*bra*(rho_tot + p_tot) - 2.*kin*(rho_tot + rho_smg)) + (rho_smg + p_smg)*9.*bra)*bra_p/2. + (rho_tot + rho_smg)*bra*kin_p + (2.*dM2*kin + 3.*pow(bra,2)*M2)*3./2.*pow(M2,-1)*p_tot_p + 3.*pvecback[pba->index_bg_kinetic_D_smg]*p_smg_p)*pow(H,-3)*pow(a,-1)/2.;

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_6_smg,
	     &pvecback[pba->index_bg_lambda_6_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_6_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_7_smg] = ((-2.) + bra)*(4. + run)*(-1.)/8.*pvecback[pba->index_bg_kinetic_D_smg] + ((-2.)*(2. + dM2) + bra*M2)*(rho_tot + p_tot)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg]*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (pvecback[pba->index_bg_kinetic_D_smg]*bra_p + ((-2.) + bra)*((-3.)*bra*bra_p + (-1.)*kin_p))*1./8.*pow(H,-1)*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_7_smg,
	     &pvecback[pba->index_bg_lambda_7_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_7_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_8_smg] = ((-2.) + bra)*(4. + run)*1./8.*pvecback[pba->index_bg_kinetic_D_smg] + 3./8.*(rho_tot + p_tot)*(((-9.)*bra + (-2.)*pvecback[pba->index_bg_kinetic_D_smg]*(3. + 2.*dM2 - bra*M2))*(-1.)/2. + (-rho_tot*dM2 - (p_smg + rho_smg*M2))*9.*pow(H,-2)*pow(M2,-1))*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/8.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (-2.*dM2 + bra*M2)*(rho_tot + p_tot)*(p_tot + p_smg)*27./16.*pow(H,-4)*pow(M2,-2) + ((-9.)*(rho_tot + p_tot) + (-6.)*bra*pow(H,2)*M2 + 3.*pow(bra,2)*pow(H,2)*M2 + (-1.)*pow(H,2)*pvecback[pba->index_bg_kinetic_D_smg]*M2)*1./8.*pow(H,-3)*pow(M2,-1)*bra_p*pow(a,-1) + ((-2.) + bra)*1./8.*pow(H,-1)*kin_p*pow(a,-1) + ((-2.) + bra)*9./16.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_smg,
	     &pvecback[pba->index_bg_lambda_8_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_9_smg] = ((-2.) + 3.*bra)*pvecback[pba->index_bg_kinetic_D_smg] + 2.*pvecback[pba->index_bg_lambda_3_smg] + (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_2_smg])*(((-3.) + 2.*bra)*(-3.)/2. + (p_tot + p_smg)*9./2.*pow(H,-2)) + (3.*bra*bra_p + kin_p)*(-1.)*pow(H,-1)*pow(a,-1) + (-9.)/2.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_9_smg,
	     &pvecback[pba->index_bg_lambda_9_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_9_smg,
	         pba->error_message,
	         "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_lambda_10_smg] = (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_3_smg])*(-2.) + (3.*bra*dM2 + kin*M2)*(rho_tot + p_tot)*3.*pow(H,-2)*pow(M2,-1) + (3.*bra + kin)*(rho_smg + p_smg)*3.*pow(H,-2) + (-1.)*pow(H,-1)*kin_p*pow(a,-1);

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_10_smg,
	     &pvecback[pba->index_bg_lambda_10_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_10_smg,
	         pba->error_message,
	         "cannot copy data back to pba->background_table");


	 pvecback[pba->index_bg_lambda_11_smg] = bra + 2.*run - (2.-bra)*ten;

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_11_smg,
	     &pvecback[pba->index_bg_lambda_11_smg],
	     1*sizeof(double));
	 class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_11_smg,
	           pba->error_message,
	           "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_cs2num_smg] = ((-2.) + bra)*((-1.)*bra + (-2.)*run + 2.*ten + (-1.)*bra*ten)*1./2. + pvecback[pba->index_bg_lambda_2_smg];

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_smg,
	     &pvecback[pba->index_bg_cs2num_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	pvecback[pba->index_bg_cs2_smg] = pvecback[pba->index_bg_cs2num_smg]/pvecback[pba->index_bg_kinetic_D_smg];

	memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_cs2_smg,
	     &pvecback[pba->index_bg_cs2_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_cs2_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	// TODO_EB: revisit G_eff and slip_eff for beyond horndeski
	double beta_1 = (run + (-1.)*ten)*2. + (1. + ten)*bra;
	double beta_2 = 2.*beta_1 + (2. + (-2.)*M2 + bra*M2)*(rho_tot + p_tot)*(-3.)*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)*pow(H,-2) + 2.*pow(H,-1)*bra_p*pow(a,-1);

	if (bra*beta_1 == 0.) {
		pvecback[pba->index_bg_G_eff_smg] = 1./M2;
	}
	else {
		pvecback[pba->index_bg_G_eff_smg] = (1. - bra*beta_1*pow(bra*beta_1 - beta_2,-1))/M2;
	}

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_G_eff_smg,
	     &pvecback[pba->index_bg_G_eff_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_G_eff_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


				 if (2.*(run - ten)*beta_1 + ten*beta_2 == 0.) {
					 pvecback[pba->index_bg_slip_eff_smg] = 1.;
			 	}
			 	else {
					pvecback[pba->index_bg_slip_eff_smg] = 1. - (2.*(run - ten)*beta_1 + ten*beta_2)*pow((run - ten)*2.*beta_1 + (1. + ten)*beta_2,-1);
			 	}

	     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_slip_eff_smg,
	     &pvecback[pba->index_bg_slip_eff_smg],
	     1*sizeof(double));
	class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_slip_eff_smg,
	       pba->error_message,
	       "cannot copy data back to pba->background_table");


	/* Here we update the minimum values of the stability quantities
	* test will be performed based on the lowest values
	*/
	//TODO_EB: move this to an appropriate place
	if (a > pba->a_min_stability_test_smg && pba->parameters_tuned_smg == _TRUE_){
		if (pvecback[pba->index_bg_kinetic_D_smg] < pba->min_D_smg){
			pba->min_D_smg = pvecback[pba->index_bg_kinetic_D_smg];
		}

		if (pvecback[pba->index_bg_cs2_smg] <= pba->min_cs2_smg){
			pba->min_cs2_smg = pvecback[pba->index_bg_cs2_smg];
		}

		if (pvecback[pba->index_bg_M2_smg] < pba->min_M2_smg){
			pba->min_M2_smg = pvecback[pba->index_bg_M2_smg];
		}

		if (pvecback[pba->index_bg_tensor_excess_smg] + 1. < pba->min_ct2_smg){
			pba->min_ct2_smg = 1. + pvecback[pba->index_bg_tensor_excess_smg];
		}

		if (pvecback[pba->index_bg_braiding_smg] < pba->min_bra_smg){
			pba->min_bra_smg = pvecback[pba->index_bg_braiding_smg];
		}
		if (pvecback[pba->index_bg_braiding_smg] > pba->max_bra_smg){
			pba->max_bra_smg = pvecback[pba->index_bg_braiding_smg];
		}
	}

	return _SUCCESS_;
}

int background_stability_tests_smg(
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

	  class_test_except(pba->min_D_smg <= -fabs(pba->D_safe_smg),
	      pba->error_message,
	      free(pvecback);free(pvecback_integration);background_free(pba),
	      "Ghost instability for scalar field perturbations with minimum D=%g \n",pba->min_D_smg);
	  class_test_except(pba->min_cs2_smg < -fabs(pba->cs2_safe_smg),
	      pba->error_message,
	      free(pvecback);free(pvecback_integration);background_free(pba),
	      "Gradient instability for scalar field perturbations with minimum c_s^2=%g \n",pba->min_cs2_smg);
	  class_test_except(pba->min_M2_smg < -fabs(pba->M2_safe_smg),
	      pba->error_message,
	      free(pvecback);free(pvecback_integration);background_free(pba),
	      "Ghost instability for metric tensor perturbations with minimum M*^2=%g \n",pba->min_M2_smg);
	  class_test_except(pba->min_ct2_smg < -fabs(pba->ct2_safe_smg),
	      pba->error_message,
	      free(pvecback);free(pvecback_integration);background_free(pba),
	      "Gradient instability for metric tensor perturbations with minimum c_t^2=%g \n",pba->min_ct2_smg);

	 }

	return _SUCCESS_;
}

int background_hi_class_second_loop(
        struct background *pba,
        double * pvecback
			) {

	/* necessary for calling array_interpolate(), but never used */
	int last_index;
	int i;

	/** - second loop over lines, overwrite derivatives that can't be analytically computed from background_functions
	 * Fill the derivatives of the Bellini-Sawicki functions in pvecback
	 * This is done just by overwriting the pvecback entries corresponding to the relevant indice
	 */

	double * pvecback_derivs;
	class_alloc(pvecback_derivs,pba->bg_size*sizeof(double),pba->error_message);

	for (i=0; i < pba->bt_size; i++) {

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

	  /* - indices for scalar field (modified gravity) */
	  class_call(background_derivs_alphas_smg(pba, pvecback, pvecback_derivs, i),
	    pba->error_message,
	    pba->error_message
	  );

	  class_call(background_gravity_functions_A_C_smg(pba,pvecback,pvecback_derivs,i),
	    pba->error_message,
	    pba->error_message
	  );

	}

	free(pvecback_derivs);

	return _SUCCESS_;
}

int background_hi_class_third_loop(
        struct background *pba,
        double * pvecback,
				double * pvecback_integration
			) {

	/* needed for growing table */
  void * memcopy_result;
	/* necessary for calling array_interpolate(), but never used */
	int last_index;
	int i;

	double * pvecback_derivs;
	class_alloc(pvecback_derivs,pba->bg_size*sizeof(double),pba->error_message);
	 /* Yet another (third!) loop to make sure the background table makes sense
	 */
	for (i=0; i < pba->bt_size; i++) {

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
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_prime_smg,
		      &pvecback_derivs[pba->index_bg_cs2num_smg],
		      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_cs2num_prime_smg,
             pba->error_message,
             "cannot copy data back to pba->background_table");

     //D'
     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_prime_smg,
			      &pvecback_derivs[pba->index_bg_kinetic_D_smg],
			      1*sizeof(double));
     class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_prime_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");

     //D_over_phiphi'
     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_over_phiphi_prime_smg,
			      &pvecback_derivs[pba->index_bg_kinetic_D_over_phiphi_smg],
			      1*sizeof(double));
     class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_kinetic_D_over_phiphi_prime_smg,
              pba->error_message,
              "cannot copy data back to pba->background_table");

     //A9'
     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A9_prime_smg,
     	      &pvecback_derivs[pba->index_bg_A9_smg],
     	      1*sizeof(double));
     class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A9_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

    //A10'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A10_prime_smg,
    	      &pvecback_derivs[pba->index_bg_A10_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A10_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

   //A12'
   memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A12_prime_smg,
   	      &pvecback_derivs[pba->index_bg_A12_smg],
   	      1*sizeof(double));
   class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A12_prime_smg,
          pba->error_message,
          "cannot copy data back to pba->background_table");

    //A13'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_A13_prime_smg,
    	      &pvecback_derivs[pba->index_bg_A13_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_A13_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

     //C9'
     memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C9_prime_smg,
     	      &pvecback_derivs[pba->index_bg_C9_smg],
     	      1*sizeof(double));
     class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C9_prime_smg,
            pba->error_message,
            "cannot copy data back to pba->background_table");

    //C10'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C10_prime_smg,
    	      &pvecback_derivs[pba->index_bg_C10_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C10_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

    //C12'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C12_prime_smg,
    	      &pvecback_derivs[pba->index_bg_C12_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C12_prime_smg,
          pba->error_message,
          "cannot copy data back to pba->background_table");

    //C13'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_C13_prime_smg,
    	      &pvecback_derivs[pba->index_bg_C13_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_C13_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

    //lambda_2'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_prime_smg,
    	      &pvecback_derivs[pba->index_bg_lambda_2_smg],
    	      1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_2_prime_smg,
           pba->error_message,
           "cannot copy data back to pba->background_table");

    //lambda_8'
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_prime_smg,
          &pvecback_derivs[pba->index_bg_lambda_8_smg],
          1*sizeof(double));
    class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_8_prime_smg,
          pba->error_message,
         "cannot copy data back to pba->background_table");

         //lambda_9'
         memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_9_prime_smg,
         	      &pvecback_derivs[pba->index_bg_lambda_9_smg],
         	      1*sizeof(double));
         class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_9_prime_smg,
                pba->error_message,
                "cannot copy data back to pba->background_table");

        //lambda_11'
        memcopy_result = memcpy(pba->background_table + i*pba->bg_size + pba->index_bg_lambda_11_prime_smg,
        	      &pvecback_derivs[pba->index_bg_lambda_11_smg],
        	      1*sizeof(double));
        class_test(memcopy_result != pba->background_table + i*pba->bg_size + pba->index_bg_lambda_11_prime_smg,
               pba->error_message,
               "cannot copy data back to pba->background_table");


	   class_call(background_at_tau(pba,
				   pba->tau_table[i],
				   long_info,
				   inter_normal,
				   &last_index, //should be no problem to use the same one as for the derivatives
				   pvecback),
		 pba->error_message,
		 pba->error_message);

	   // check if any of the values becomes nan
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

	return _SUCCESS_;
}

int background_initial_conditions_smg(
        struct background *pba,
				double a,
				double * pvecback,
        double * pvecback_integration,
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
	//         phi_scale = gsl_sf_lambert_W0(10);

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
	//         printf("V_scale %e, i %i \n",V_scale,i);

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
	else
	  {/* non attractor ICs */
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
	  }

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

	    if (pba->rho_evolution_smg == _TRUE_){

	//DT: moved here from the definition of wowa_w model, to avoid pvecback[pba->index_bg_rho_smg] mix up
	if (pba->expansion_model_smg == wowa_w)
	  {
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
	      else
	        {
	    pvecback_integration[pba->index_bi_rho_smg] = pvecback[pba->index_bg_rho_smg];
	  }
	    }

	return _SUCCESS_;
}

int background_store_columntitles_smg(
				struct background *pba,
				char titles[_MAXTITLESTRINGLENGTH_]
			) {

	class_store_columntitle(titles,"(.)rho_smg",_TRUE_);
  class_store_columntitle(titles,"(.)p_smg",_TRUE_);

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
  }

	return _SUCCESS_;
}

int background_store_doubles_smg(
				struct background *pba,
        double * pvecback,
				double * dataptr,
				int * ptr_storeidx
			) {

	int storeidx = *ptr_storeidx;

	class_store_double(dataptr,pvecback[pba->index_bg_rho_smg],_TRUE_,storeidx);
  class_store_double(dataptr,pvecback[pba->index_bg_p_smg],_TRUE_,storeidx);

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
		class_store_double(dataptr,pvecback[pba->index_bg_B0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B3_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B4_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B5_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B6_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B7_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B8_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B9_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B10_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B11_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_B12_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C0_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C1_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C2_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C3_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C4_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C5_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C6_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C7_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C8_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C9_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C10_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C11_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C12_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C13_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C15_smg],pba->field_evolution_smg,storeidx);
		class_store_double(dataptr,pvecback[pba->index_bg_C16_smg],pba->field_evolution_smg,storeidx);
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
	}

	*ptr_storeidx = storeidx;

	return _SUCCESS_;
}

int background_gravity_parameters(
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

   default:
       printf("Modified gravity: output not implemented in background_gravity_parameters() \n");


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

      default:
       printf("Parameterized model: expansion hisotry output not implemented in background_gravity_parameters() \n");

    }

  }

  return _SUCCESS_;

}

int background_free_smg(
			  struct background *pba
			) {

	if (pba->parameters_smg != NULL)
	  free(pba->parameters_smg);
	//dealocate parameters_2_smg only for parameterizations
	if (pba->field_evolution_smg == _FALSE_ && pba->parameters_2_smg != NULL)
	  free(pba->parameters_2_smg);

	return _SUCCESS_;
}

int background_print_smg(
			  struct background *pba,
				double * pvecback,
				double * pvecback_integration
			) {

	printf(" -> Omega_smg = %f, wanted %f ",pvecback[pba->index_bg_rho_smg]/pvecback[pba->index_bg_rho_crit], pba->Omega0_smg);
	if(pba->has_lambda == _TRUE_)
		printf(", Omega_Lambda = %f", pba->Omega0_lambda);
	printf("\n");
	if (pba->background_verbose > 3) {
		printf("Minimal stability values: cs2 = %g, ct2 = %g, D = %g, M2 = %g \n",pba->min_cs2_smg,pba->min_ct2_smg,pba->min_D_smg,pba->min_M2_smg);
	}

	if (pba->field_evolution_smg == _TRUE_){
	    pba->xi_0_smg = pvecback_integration[pba->index_bi_phi_prime_smg]*pvecback[pba->index_bg_H]/pow(pba->H0,2);
	    pba->phi_0_smg = pvecback_integration[pba->index_bi_phi_smg];
	}

	background_gravity_parameters(pba);

	return _SUCCESS_;
}

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
