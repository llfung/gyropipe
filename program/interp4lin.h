/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * chosen_interpolation4.h
 *
 * Code generation for function 'chosen_interpolation4'
 *
 */

#ifndef INTERP4LIN_H
#define INTERP4LIN_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>

#define N_X 97
#define N_Y 13
#define N_Z 13
#define N_W 13

#define N_XYZW 213109

#define N_XY 1261
#define N_XYZ 16393


/* Function Declarations */
extern void interp4lin_(const double x[N_X], const double
  y[N_Y], const double z[N_Z], const double w[N_W],const double V[N_XYZW], const double
  xq[], const double yq[], const double zq[], const double wq[],
  double Vq[], const int *insiz);

#endif

/* End of code generation (chosen_interpolation4.h) */
