
class_precision_parameter(tol_einstein00_reldev,double,1.e-3) /**< tolerance to deviations w.r.t. the Einstein 00 equation. Useful if get_h_from_trace,int==_TRUE_. (not only _smg!!) */

class_precision_parameter(einstein00_friction,double,1.) /**< friction term muliplying the Einstein 00 equation to correct for h''. (not only _smg!!) */

/*
 * hi_class (_smg) parameters
 * */

/* This is the time at which we test for the initial value of the QS approximation.
It has to be at least a_ini_over_a_today_default */
class_precision_parameter(a_ini_test_qs_smg,double,1.e-14)

class_precision_parameter(n_min_qs_smg,int,1e2) /**< minimum number of steps used to sample the quantities in the quasi-static approximation (qs_smg) */
class_precision_parameter(n_max_qs_smg,int,1e4) /**< maximum number of steps used to sample the quantities in the quasi-static approximation (qs_smg) */
class_precision_parameter(z_fd_qs_smg,double,10.) /**< minimum redshift after which the user requires the full-dynamic evolution */
class_precision_parameter(trigger_mass_qs_smg,double,1.e3) /**< if the mass is above this trigger the quasi-static approximation is switched on */
class_precision_parameter(trigger_rad_qs_smg,double,1.e3) /**< if the radiation component is still important w.r.t.\ the scalar field the quasi-static approximation can not be used */
class_precision_parameter(eps_s_qs_smg,double,0.01) /**< when the system enters the quasi-static evolution this parameter measures how much the oscillation are decaying with time */


class_precision_parameter(min_a_pert_smg,double,1.) /**< minimum value of scale factor to start integration (important to test some ede models */
class_precision_parameter(pert_ic_tolerance_smg,double,2.e-2) /**< tolerance to deviations from n=2 for IC h~tau^n. Negative values override test */
class_precision_parameter(pert_ic_ini_z_ref_smg,double,1.e10) /**<Reference z to carry out test for conservation of curvature before pert evolution*/
class_precision_parameter(pert_ic_regulator_smg,double,1.e-15)  /* minumum size of denominator in IC expressions: regulate to prevent infinities. Negative => off */
class_precision_parameter(pert_qs_ic_tolerance_test_smg,double,1.) /* maximal fractional contribution to (0i) equation of SMG terms in QS initial condition */
