/** @file gravity_functions_smg.c Documented gravity_functions_smg module
 *
 * Emilio Bellini, Ignacy Sawicki, Miguel Zumalacarregui, TODO_EB: date here xx.xx.xxxx
 *
 * This module contains all the complicated expressions
 * used in hi_class. In particular, they are casted in
 * different functions, depending on their type:
 * - background as a functions of the Gs
 * - alphas as a functions of the Gs
 * - As as a functions of the alphas
 * - Bs as a functions of the Gs
 * - Cs as a functions of the As
 *
 * The purpose of this module is twofold:
 * - isolate long expressions in one single file
 * - have a file that can be optimized for specific theories.
 *   Indeed, the problem of numerical errors can be alleviated
 *   when analytic cancellations greatly simplify long expressions
 */

#include "gravity_functions_smg.h"
// background_gravity_functions_smg


int gravity_functions_back_from_Gs_smg(

) {

  return _SUCCESS_;
}

int gravity_functions_alphas_from_Gs_smg(

) {

  return _SUCCESS_;
}


/**
* Get gravity functions As from alphas.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input/Output: vector of background quantities
* @param pvecback_derivs      Input: vector of derivatives
* @return the error status
*/
int gravity_functions_As_from_alphas_smg(
                                         struct background *pba,
                                         double * pvecback,
                                         double * pvecback_derivs
                                         ) {

  // basic background quantities
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double H_p = pvecback[pba->index_bg_H_prime];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_smg = pvecback[pba->index_bg_p_smg];
  double rho_smg = pvecback[pba->index_bg_rho_smg];
  // factor to convert lna derivatives to tau derivatives
  double factor = a*H;

  // alphas
  double M2 = pvecback[pba->index_bg_M2_smg];
  double DelM2 = pvecback[pba->index_bg_delta_M2_smg];
  double kin = pvecback[pba->index_bg_kineticity_smg];
  double bra = pvecback[pba->index_bg_braiding_smg];
  double run = pvecback[pba->index_bg_mpl_running_smg];
  double ten = pvecback[pba->index_bg_tensor_excess_smg];
  double beh = pvecback[pba->index_bg_beyond_horndeski_smg];
  double dM2 = pvecback[pba->index_bg_delta_M2_smg];

  // need to update the time derivatives of the interesting functions
  double kin_p = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
  double bra_p = factor*pvecback_derivs[pba->index_bg_braiding_smg];
  double run_p = factor*pvecback_derivs[pba->index_bg_mpl_running_smg];
  double ten_p = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
  double beh_p = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
  double p_tot_p = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
  double p_smg_p = factor*pvecback_derivs[pba->index_bg_p_smg];

  // kinetic term D
  pvecback[pba->index_bg_kinetic_D_smg] = kin + 3./2.*pow(bra,2);

  // A0
	pvecback[pba->index_bg_A0_smg] =
	1./2.*(
    + bra - 3.*(rho_smg + p_smg + (rho_tot + p_tot)*DelM2/M2)*pow(H,-2)
	);

  // A1
	pvecback[pba->index_bg_A1_smg] =
	+ (1. + ten)*kin
	- 3.*(beh*(1. + run) + run - ten + beh_p/a/H)*bra;

  // A2
	pvecback[pba->index_bg_A2_smg] =
	- (kin + 3./2.*pow(bra,2))*(2. + run)
	- 9./4.*bra*(
	  + (2. - bra)*(rho_smg+p_smg)
	  + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
	)*pow(H,-2)
	- 3./2.*bra*bra_p/a/H;

  // A3
	pvecback[pba->index_bg_A3_smg] = bra*beh;

  // A4
  pvecback[pba->index_bg_A4_smg] =
  9./2.*kin*(
    + (2. - bra)*(rho_smg+p_smg)
    + (2.*DelM2/M2 - bra)*(rho_tot+p_tot)
  )*pow(H,-2)
  + 3.*(bra*kin_p - kin*bra_p)/a/H;

  // A5
	pvecback[pba->index_bg_A5_smg] = - beh*kin;

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

  // A8
	pvecback[pba->index_bg_A8_smg] = run - ten - beh;

  // A9
	pvecback[pba->index_bg_A9_smg] =
	+ 3./4.*(
	  + (2. - bra)*(rho_smg + p_smg)
	  + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2)
	+ 1./2.*bra_p/a/H;

  // A10
	pvecback[pba->index_bg_A10_smg] =
	bra + 2.*run - (2. - bra)*ten + 2.*(1. + run)*beh + 2.*beh_p/a/H;

  // A11
	pvecback[pba->index_bg_A11_smg] =
	- (kin + 3./2.*pow(bra,2))*(4. + run)
	+ 3./4.*(
	  + (4.*kin + 6.*bra + 3.*pow(bra,2))*(rho_smg + p_smg)
	  + (4.*kin + 6.*bra*DelM2/M2 + 3.*pow(bra,2))*(rho_tot + p_tot)
	)*pow(H,-2)
	- (kin_p + 3./2.*bra*bra_p)/a/H;

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

  // A13
	pvecback[pba->index_bg_A13_smg] =
	- bra - 2.*run + (2. - bra)*ten - (2. + bra + 2.*run)*beh
	- 3./2.*(
	 	+ (2. - bra - 2.*beh)*(rho_smg + p_smg)*pow(H,-2)
	 	+ (2.*DelM2/M2 - bra - 2.*beh)*(rho_tot + p_tot)*pow(H,-2)
	)
	- (bra_p + 2.*beh_p)/a/H;

  // A14
	pvecback[pba->index_bg_A14_smg] = - (kin + 3.*bra)/2.;

  // A15
	pvecback[pba->index_bg_A15_smg] = - 1./2.*bra - beh;

  // A16
	pvecback[pba->index_bg_A16_smg] =
	- 1./2.*(kin + 3.*bra)
	+ 9./4.*(
	 + (2. - bra)*(rho_smg + p_smg)
	 + (2.*DelM2/M2 - bra)*(rho_tot + p_tot)
	)*pow(H,-2);


  // TODO_EB: remove below this when IC are recalculated

  pvecback[pba->index_bg_lambda_1_smg] = (run + (-1.)*ten)*(-3.)*bra + (1. + ten)*kin;

  pvecback[pba->index_bg_lambda_2_smg] = (- 2.*dM2 + bra*M2)*(rho_tot + p_tot)*(-3.)/2.*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/2.*pow(H,-2) + pow(H,-1)*bra_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_3_smg] = (2. + run)*(-1.)/2.*pvecback[pba->index_bg_kinetic_D_smg] + (-3.)/4.*bra*pvecback[pba->index_bg_lambda_2_smg];

	pvecback[pba->index_bg_lambda_4_smg] = kin*pvecback[pba->index_bg_lambda_2_smg] + (2.*kin*bra_p + (-1.)*bra*kin_p)*(-1.)*pow(H,-1)*pow(a,-1);

	pvecback[pba->index_bg_lambda_5_smg] = (bra + 2.*run + (-2.)*ten + bra*ten)*3./2.*bra + (run + (-1.)*ten)*pvecback[pba->index_bg_kinetic_D_smg] + 3./2.*bra*pvecback[pba->index_bg_lambda_2_smg];

	pvecback[pba->index_bg_lambda_6_smg] = 3./2.*(((9./2.*bra + kin)*dM2*pow(M2,-1) + (-9.)/4.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*pow(rho_tot,2) + ((9.*bra + kin)*dM2*pow(M2,-1) + (-9.)/2.*pow(bra,2) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_tot*p_tot + 9./2.*bra*(dM2 - M2*bra/2.)*pow(M2,-1)*pow(p_tot,2) + (kin*dM2*pow(M2,-1) - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*(rho_tot + p_tot)*rho_smg + ((kin - bra*kin/2. + pvecback[pba->index_bg_kinetic_D_smg]*run)*rho_smg + ((9.*bra + kin)*(2. - bra)/2. + pvecback[pba->index_bg_kinetic_D_smg]*run - 9./2.*bra*pow(M2,-1))*rho_tot + 9.*bra*(1. - bra/2. - pow(M2,-1)/2.)*p_tot)*(rho_smg + p_smg) + 9./2.*bra*(1. - bra/2.)*pow(rho_smg + p_smg,2))*pow(H,-4) + (((9.*bra*(rho_tot + p_tot) - 2.*kin*(rho_tot + rho_smg)) + (rho_smg + p_smg)*9.*bra)*bra_p/2. + (rho_tot + rho_smg)*bra*kin_p + (2.*dM2*kin + 3.*pow(bra,2)*M2)*3./2.*pow(M2,-1)*p_tot_p + 3.*pvecback[pba->index_bg_kinetic_D_smg]*p_smg_p)*pow(H,-3)*pow(a,-1)/2.;

	pvecback[pba->index_bg_lambda_7_smg] = ((-2.) + bra)*(4. + run)*(-1.)/8.*pvecback[pba->index_bg_kinetic_D_smg] + ((-2.)*(2. + dM2) + bra*M2)*(rho_tot + p_tot)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg]*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*3./16.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (pvecback[pba->index_bg_kinetic_D_smg]*bra_p + ((-2.) + bra)*((-3.)*bra*bra_p + (-1.)*kin_p))*1./8.*pow(H,-1)*pow(a,-1);

	pvecback[pba->index_bg_lambda_8_smg] = ((-2.) + bra)*(4. + run)*1./8.*pvecback[pba->index_bg_kinetic_D_smg] + 3./8.*(rho_tot + p_tot)*(((-9.)*bra + (-2.)*pvecback[pba->index_bg_kinetic_D_smg]*(3. + 2.*dM2 - bra*M2))*(-1.)/2. + (-rho_tot*dM2 - (p_smg + rho_smg*M2))*9.*pow(H,-2)*pow(M2,-1))*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)/8.*pow(H,-2)*pvecback[pba->index_bg_kinetic_D_smg] + (-2.*dM2 + bra*M2)*(rho_tot + p_tot)*(p_tot + p_smg)*27./16.*pow(H,-4)*pow(M2,-2) + ((-9.)*(rho_tot + p_tot) + (-6.)*bra*pow(H,2)*M2 + 3.*pow(bra,2)*pow(H,2)*M2 + (-1.)*pow(H,2)*pvecback[pba->index_bg_kinetic_D_smg]*M2)*1./8.*pow(H,-3)*pow(M2,-1)*bra_p*pow(a,-1) + ((-2.) + bra)*1./8.*pow(H,-1)*kin_p*pow(a,-1) + ((-2.) + bra)*9./16.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_9_smg] = ((-2.) + 3.*bra)*pvecback[pba->index_bg_kinetic_D_smg] + 2.*pvecback[pba->index_bg_lambda_3_smg] + (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_2_smg])*(((-3.) + 2.*bra)*(-3.)/2. + (p_tot + p_smg)*9./2.*pow(H,-2)) + (3.*bra*bra_p + kin_p)*(-1.)*pow(H,-1)*pow(a,-1) + (-9.)/2.*bra*pow(H,-3)*pow(M2,-1)*p_tot_p*pow(a,-1);

	pvecback[pba->index_bg_lambda_10_smg] = (pvecback[pba->index_bg_kinetic_D_smg] + (-1.)*pvecback[pba->index_bg_lambda_3_smg])*(-2.) + (3.*bra*dM2 + kin*M2)*(rho_tot + p_tot)*3.*pow(H,-2)*pow(M2,-1) + (3.*bra + kin)*(rho_smg + p_smg)*3.*pow(H,-2) + (-1.)*pow(H,-1)*kin_p*pow(a,-1);

  pvecback[pba->index_bg_lambda_11_smg] = bra + 2.*run - (2.-bra)*ten;


	pvecback[pba->index_bg_cs2num_smg] = ((-2.) + bra)*((-1.)*bra + (-2.)*run + 2.*ten + (-1.)*bra*ten)*1./2. + pvecback[pba->index_bg_lambda_2_smg];

  // TODO_EB: rewrite cs2, Geff and slip for beyond Horndeski (calculate them in the hi_class.nb Mathematica notebook)
  // TODO_EB: check if there is a better alternative to regularizing these quantities

	if (pvecback[pba->index_bg_cs2num_smg] == pvecback[pba->index_bg_kinetic_D_smg]) {
		pvecback[pba->index_bg_cs2_smg] = 1.;
	}
	else {
		pvecback[pba->index_bg_cs2_smg] = pvecback[pba->index_bg_cs2num_smg]/pvecback[pba->index_bg_kinetic_D_smg];
	}

	double beta_1 = (run + (-1.)*ten)*2. + (1. + ten)*bra;
	double beta_2 = 2.*beta_1 + (2. + (-2.)*M2 + bra*M2)*(rho_tot + p_tot)*(-3.)*pow(H,-2)*pow(M2,-1) + ((-2.) + bra)*(rho_smg + p_smg)*(-3.)*pow(H,-2) + 2.*pow(H,-1)*bra_p*pow(a,-1);

	if (bra*beta_1 == 0.) {
		pvecback[pba->index_bg_G_eff_smg] = 1./M2;
	}
	else {
		pvecback[pba->index_bg_G_eff_smg] = (1. - bra*beta_1*pow(bra*beta_1 - beta_2,-1))/M2;
	}

  if (2.*(run - ten)*beta_1 + ten*beta_2 == 0.) {
		pvecback[pba->index_bg_slip_eff_smg] = 1.;
	}
	else {
		pvecback[pba->index_bg_slip_eff_smg] = 1. - (2.*(run - ten)*beta_1 + ten*beta_2)*pow((run - ten)*2.*beta_1 + (1. + ten)*beta_2,-1);
	}

  return _SUCCESS_;
}





int gravity_functions_Bs_from_Gs_smg(

) {

  return _SUCCESS_;
}





/**
* Get gravity functions Cs from Bs.
*
* @param pba                  Input: pointer to background structure
* @param pvecback             Input/Output: vector of background quantities
* @return the error status
*/
int gravity_functions_Cs_from_Bs_smg(
                                     struct background *pba,
                                     double * pvecback,
                                     double * pvecback_derivs
                                     ) {

  // basic background quantities
  double a = pvecback[pba->index_bg_a];
  double H = pvecback[pba->index_bg_H];
  double H_p = pvecback[pba->index_bg_H_prime];
  double p_tot = pvecback[pba->index_bg_p_tot_wo_smg];
  double rho_tot = pvecback[pba->index_bg_rho_tot_wo_smg];
  double p_smg = pvecback[pba->index_bg_p_smg];
  double rho_smg = pvecback[pba->index_bg_rho_smg];
  // factor to convert lna derivatives to tau derivatives
  double factor = a*H;

  // alphas
  double M2 = pvecback[pba->index_bg_M2_smg];
  double DelM2 = pvecback[pba->index_bg_delta_M2_smg];
  double kin = pvecback[pba->index_bg_kineticity_smg];
  double bra = pvecback[pba->index_bg_braiding_smg];
  double run = pvecback[pba->index_bg_mpl_running_smg];
  double ten = pvecback[pba->index_bg_tensor_excess_smg];
  double beh = pvecback[pba->index_bg_beyond_horndeski_smg];
  double dM2 = pvecback[pba->index_bg_delta_M2_smg];

  // need to update the time derivatives of the interesting functions
  double kin_p = factor*pvecback_derivs[pba->index_bg_kineticity_smg];
  double bra_p = factor*pvecback_derivs[pba->index_bg_braiding_smg];
  double run_p = factor*pvecback_derivs[pba->index_bg_mpl_running_smg];
  double ten_p = factor*pvecback_derivs[pba->index_bg_tensor_excess_smg];
  double beh_p = factor*pvecback_derivs[pba->index_bg_beyond_horndeski_smg];
  double p_tot_p = factor*pvecback_derivs[pba->index_bg_p_tot_wo_smg];
  double p_smg_p = factor*pvecback_derivs[pba->index_bg_p_smg];

  // alphas over scalar field
  double kin_ss = pvecback[pba->index_bg_kineticity_over_phiphi_smg];
  double bra_s = pvecback[pba->index_bg_braiding_over_phi_smg];
  double beh_s = pvecback[pba->index_bg_beyond_horndeski_over_phi_smg];

  // Bs
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

  // kinetic term D over phiphi
  pvecback[pba->index_bg_kinetic_D_over_phiphi_smg] =
  + kin_ss + 3./2.*pow(bra_s,2);

  // C0
  pvecback[pba->index_bg_C0_smg] = B0;

  // C1
  pvecback[pba->index_bg_C1_smg] = kin_ss*(1. + ten) - 3./2.*bra_s*B6;

  // C2
  pvecback[pba->index_bg_C2_smg] = - kin_ss*(2. + run) - 3.*bra_s*B5;

  // C3
  pvecback[pba->index_bg_C3_smg] = bra_s*beh_s;

  // C4
  pvecback[pba->index_bg_C4_smg] = kin_ss*B1 - 3.*bra_s*B7;

  // C5
  pvecback[pba->index_bg_C5_smg] = - kin_ss*beh_s;

  // C6
  pvecback[pba->index_bg_C6_smg] = kin_ss*B2 - 3.*bra_s*B9;

  // C7
  pvecback[pba->index_bg_C7_smg] = kin_ss*B3 - 3.*bra_s*B8;

  // C8
  pvecback[pba->index_bg_C8_smg] = B4;

  // C9
  pvecback[pba->index_bg_C9_smg] = - bra_s - bra_s*run/2. + B5;

  // C10
  pvecback[pba->index_bg_C10_smg] = bra_s*(1. + ten) + B6;

  // C11
  pvecback[pba->index_bg_C11_smg] = bra_s*B1/2. + B7;

  // C12
  pvecback[pba->index_bg_C12_smg] = bra_s*B2/2. + B9;

  // C13
  pvecback[pba->index_bg_C13_smg] = bra_s*B3/2. + B8;

  // C14
  pvecback[pba->index_bg_C14_smg] = B10;

  // C15
  pvecback[pba->index_bg_C15_smg] = B11;

  // C16
  pvecback[pba->index_bg_C16_smg] = B12;

  return _SUCCESS_;
}
