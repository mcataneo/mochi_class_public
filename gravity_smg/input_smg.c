#include "input_smg.h"

int input_warnings_smg(
  int get_h_from_trace,
  int input_verbose
) {

  /* Here we put a warning as we want to encourage hi_class users to get
  h_prime from the trace of the Einstein ij equation than from the Einstein 00
  equation. This is because the Einstein 00 equation has a gauge dependent
  singularity that can be removed using the trace of the Einstein ij equation.
  */
  if (input_verbose > 0) {
    if (get_h_from_trace == _FALSE_) {
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
    else if (get_h_from_trace == _TRUE_) {
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
