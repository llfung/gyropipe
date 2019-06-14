#ifndef INTERP1IN_H
#define INTERP1IN_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#define ONE_N_X 205

/* Function Declarations */
extern void interp1in_(const double x[ONE_N_X], const double y[ONE_N_X],
  const double xi[], double yi[], const int *insiz);

#endif

/* End of code generation (chosen_interpolation.h) */
