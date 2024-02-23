/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

/* this main runs only the background part */

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermodynamics th;           /* for thermodynamics */
  struct perturbations pt;         /* for source functions */
  struct transfer tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct harmonic hr;          /* for output spectra */
  struct fourier fo;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct distortions sd;      /* for spectral distortions */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */

  if (input_init(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&hr,&fo,&le,&sd,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  /****** here you can output the evolution of any background
	  quanitity you are interested in ******/

  int index_tau;
  FILE * output_file;
  // output_file = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/output/designer_fR_small_fr_stable_params_full_grid.dat","w");
	// output_file = fopen("/Users/matteoc/Documents/Projects/Pseudo_Emulator_v2/test_hiclass/designer_fR/designer_fR_stable_params_full_grid_hiclass_v2.dat","w");
  // output_file = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/brans_dicke_stable_params_hiclass.dat","w");
  // output_file = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/cubic_galileon_mnu_0p4_background.dat","w");
  output_file = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/nkgb_background.dat","w");
  // output_file = fopen("/Users/matteoc/Documents/Projects/Pseudo_Emulator_v2/test_hiclass/designer_fR/propto_omega_stable_params_full_grid_hiclass.dat","w");
  // FILE * output_file_w;
  // output_file_w = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/w_de_stable_nkgb_mnu_0p4_hiclass.dat","w");
  // output_file_w = fopen("/Users/matteoc/Documents/Projects/class_env/hi_class_pub_devel/rho_de_stable_brans_dicke_hiclass.dat","w");
  
  for (index_tau=0; index_tau<ba.bt_size; index_tau++) {

    // fprintf(stdout,
	  //   "tau=%e z=%e a=%e H=%e\n",
	  //   ba.tau_table[index_tau],
	  //   ba.z_table[index_tau],
	  //   ba.background_table[index_tau*ba.bg_size+ba.index_bg_a],
	  //   ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]);

      // fprintf(stdout,
	    // "%e \t %e \t %e \t %e \t %e\n",
	    // log(ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]),
      // ba.background_table[index_tau*ba.bg_size+ba.index_bg_kineticity_smg],
	    // ba.background_table[index_tau*ba.bg_size+ba.index_bg_braiding_smg],
      // ba.background_table[index_tau*ba.bg_size+ba.index_bg_mpl_running_smg],
      // ba.background_table[index_tau*ba.bg_size+ba.index_bg_tensor_excess_smg]);

		// fprintf(output_file,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
		// log(ba.background_table[index_tau*ba.bg_size + ba.index_bg_a]),
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_mpl_running_smg],
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_braiding_smg],
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_kineticity_smg],
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_delta_M2_smg],
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_cs2_smg],
		// ba.background_table[index_tau*ba.bg_size + ba.index_bg_kinetic_D_smg],
    // ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]*ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]*ba.background_table[index_tau*ba.bg_size + ba.index_bg_mpl_running_prime_smg],
		// ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]*ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]*ba.background_table[index_tau*ba.bg_size + ba.index_bg_braiding_prime_smg],
		// ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]*ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]*ba.background_table[index_tau*ba.bg_size + ba.index_bg_kineticity_prime_smg],
		// ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]*ba.background_table[index_tau*ba.bg_size+ba.index_bg_H]*ba.background_table[index_tau*ba.bg_size + ba.index_bg_kinetic_D_prime_smg],
    // ba.background_table[index_tau*ba.bg_size + ba.index_bg_cs2num_prime_smg]);

    double a = ba.background_table[index_tau*ba.bg_size+ba.index_bg_a];
    
    if (a > exp(-5)) {
      // fprintf(output_file,"%.15e %.15e %.15e %.15e\n",
      // log(ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]),
      // ba.background_table[index_tau*ba.bg_size + ba.index_bg_delta_M2_smg],
      // ba.background_table[index_tau*ba.bg_size + ba.index_bg_kinetic_D_smg],
      // // ba.background_table[index_tau*ba.bg_size + ba.index_bg_cs2_smg]
      // 1.
      // );

      // fprintf(output_file_w,"%.15e %.15e\n",
      // log(ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]),
      // // ba.background_table[index_tau*ba.bg_size + ba.index_bg_p_smg]/ba.background_table[index_tau*ba.bg_size + ba.index_bg_rho_smg]
      // ba.background_table[index_tau*ba.bg_size + ba.index_bg_rho_smg]/ba.background_table[(ba.bt_size-1)*ba.bg_size + ba.index_bg_rho_smg]
      // );

      fprintf(output_file,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
      log(ba.background_table[index_tau*ba.bg_size + ba.index_bg_a]),
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_H],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_rho_tot],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_p_tot],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_rho_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_p_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_M2_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_braiding_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_mpl_running_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_cs2_smg],
      ba.background_table[index_tau*ba.bg_size + ba.index_bg_kinetic_D_smg]
      );

      // fprintf(output_file,"%.15e %.15e\n",
      // log(ba.background_table[index_tau*ba.bg_size+ba.index_bg_a]),
      // ba.background_table[index_tau*ba.bg_size + ba.index_bg_braiding_smg]
      // );
    }

    if (a==1.) {
      printf("alpha_B(z=0) = %.15e\n",ba.background_table[index_tau*ba.bg_size + ba.index_bg_braiding_smg]);
    }
	}
	fclose(output_file);
  // fclose(output_file_w);
	// class_stop(pba->error_message,"Stop here, for now.");

  /****** all calculations done, now free the structures ******/

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
