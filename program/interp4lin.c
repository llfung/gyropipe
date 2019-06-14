/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * chosen_interpolation4.c
 *
 * Code generation for function 'chosen_interpolation4'
 *
 */

/* Include files */

#include "interp4lin.h"


/* Function Definitions */
void interp4lin_(const double x[N_X], const double
  y[N_Y], const double z[N_Z], const double w[N_W],const double V[N_XYZW], const double
  xq[], const double yq[], const double zq[],
  const double wq[],
  double Vq[], const int *insiz)
{
  double xmin[4];
  double xmax[4];
  int k;
  int high_i;
  double Xq[4];
  int low_i;
  int low_ip1;
  int b_low_i;
  int mid_i;
  int c_low_i;
  int d_low_i;
  double vbox[16];
  static const short vbbidx[16] = { 1, 2, N_X+1, N_X+2, N_XY+1, N_XY+2, N_XY+N_X+1, N_XY+N_X+2, N_XYZ+1, N_XYZ+2,
    N_XYZ+N_X+1, N_XYZ+N_X+2, N_XYZ+N_XY+1, N_XYZ+N_XY+2, N_XYZ+N_XY+N_X+1, N_XYZ+N_XY+N_X+2 };

  double r;
  xmin[0] = x[0];
  xmax[0] = x[N_X-1];
  xmin[1] = y[0];
  xmax[1] = y[N_Y-1];
  xmin[2] = z[0];
  xmax[2] = z[N_Z-1];
  xmin[3] = w[0];
  xmax[3] = w[N_W-1];


  for (k = 0; k < *insiz; k++) {
    Xq[0] = xq[k];
    Xq[1] = yq[k];
    Xq[2] = zq[k];
    Xq[3] = wq[k];

      low_i = 1;
      low_ip1 = 2;
      high_i = N_X;
      while (high_i > low_ip1) {
        mid_i = (low_i + high_i) >> 1;
        if (Xq[0] >= x[mid_i - 1]) {
          low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

      b_low_i = 1;
      low_ip1 = 2;
      high_i = N_Y;
      while (high_i > low_ip1) {
        mid_i = (b_low_i + high_i) >> 1;
        if (Xq[1] >= y[mid_i - 1]) {
          b_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

      c_low_i = 1;
      low_ip1 = 2;
      high_i = N_Z;
      while (high_i > low_ip1) {
        mid_i = (c_low_i + high_i) >> 1;
        if (Xq[2] >= z[mid_i - 1]) {
          c_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

      d_low_i = 1;
      low_ip1 = 2;
      high_i = N_W;
      while (high_i > low_ip1) {
        mid_i = (d_low_i + high_i) >> 1;
        if (Xq[3] >= w[mid_i - 1]) {
          d_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

      high_i = ((low_i + (b_low_i - 1) * N_X) + (c_low_i - 1) * N_XY) + (d_low_i -
        1) * N_XYZ;
      for (low_ip1 = 0; low_ip1 < 16; low_ip1++) {
        vbox[low_ip1] = V[(vbbidx[low_ip1] + high_i) - 2];
      }

      if (xq[k] == x[low_i - 1]) {
        for (high_i = 0; high_i < 8; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 2];
        }
      } else if (xq[k] == x[low_i]) {
        for (high_i = 0; high_i < 8; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 1];
        }
      } else {
        r = (xq[k] - x[low_i - 1]) / (x[low_i] -
          x[low_i - 1]);
        for (high_i = 0; high_i < 8; high_i++) {
          low_ip1 = (high_i + 1) << 1;
          vbox[high_i] = (1.0 - r) * vbox[low_ip1 - 2] + r * vbox[low_ip1 - 1];
        }
      }

      if (yq[k] == y[b_low_i - 1]) {
        for (high_i = 0; high_i < 4; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 2];
        }
      } else if (yq[k] == y[b_low_i]) {
        for (high_i = 0; high_i < 4; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 1];
        }
      } else {
        r = (yq[k] - y[b_low_i - 1]) / (y[b_low_i] -
          y[b_low_i - 1]);
        for (high_i = 0; high_i < 4; high_i++) {
          vbox[high_i] = (1.0 - r) * vbox[((high_i + 1) << 1) - 2] + r * vbox
            [((high_i + 1) << 1) - 1];
        }
      }

      if (zq[k] == z[c_low_i - 1]) {
        for (high_i = 0; high_i < 2; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 2];
        }
      } else if (zq[k] == z[c_low_i]) {
        for (high_i = 0; high_i < 2; high_i++) {
          vbox[high_i] = vbox[((high_i + 1) << 1) - 1];
        }
      } else {
        r = (zq[k] - z[c_low_i - 1]) / (z[c_low_i] -
          z[c_low_i - 1]);
        for (high_i = 0; high_i < 2; high_i++) {
          vbox[high_i] = (1.0 - r) * vbox[((high_i + 1) << 1) - 2] + r * vbox
            [((high_i + 1) << 1) - 1];
        }
      }

      if (!(wq[k] == w[d_low_i - 1])) {
        if (wq[k] == w[d_low_i]) {
          for (high_i = 0; high_i < 1; high_i++) {
            vbox[0] = vbox[1];
          }
        } else {
          r = (wq[k] - w[d_low_i - 1]) / (w[d_low_i] -
            w[d_low_i - 1]);
          for (high_i = 0; high_i < 1; high_i++) {
            vbox[high_i] = (1.0 - r) * vbox[0] + r * vbox[1];
          }
        }
      }

      Vq[k] = vbox[0];

  }
}




/* End of code generation (chosen_interpolation4.c) */
