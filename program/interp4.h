/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * chosen_interpolation3.h
 *
 * Code generation for function 'chosen_interpolation3'
 *
 */

#ifndef INTERP4_H
#define INTERP4_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#define N_X 97
#define N_Y 13
#define N_Z 13
#define N_W 13

#define N_XYZW 213109

#define N_XP 99
#define N_YP 15
#define N_ZP 15
#define N_WP 15

#define N_XYP 1485
#define N_XYZP 22275
#define N_XYZWP 334125


/* Function Declarations */
extern void interp4_libgen_(const double x[N_X],  const double y[N_Y],    const double z[N_Z],  const double w[N_W], const double v[N_XYZW],
             double VV[N_XYZWP]);
extern void interp4_interp_(const double x[N_X],  const double y[N_Y],    const double z[N_Z],  const double w[N_W], const double VV[N_XYZWP],
             const double xq[], const double yq[],   const double zq[], const double wq[],double vq[], const int *insiz);

#endif

/* End of code generation (chosen_interpolation3.h) */
