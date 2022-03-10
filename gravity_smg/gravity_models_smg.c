#include "gravity_models_smg.h"


int gravity_models_properties_smg() {

  return _SUCCESS_;
}


int gravity_models_get_Gs_smg(
                              struct background *pba,
                              struct G_functions_and_derivs *g_fun
                              ) {

  // g_fun->G4 = 0.;
  // g_fun->F4 = 0.;
  printf("%f, %f\n", g_fun->G4, g_fun->F4);

  return _SUCCESS_;
}


// input_read_parameters_smg
// background_gravity_functions_smg
// background_initial_conditions_smg
// print_stdout_gravity_parameters_smg
