/* Include files */
#include <string.h>
#include <math.h>
#include "interp4.h"

void interp4_libgen_(const double x[N_X],  const double y[N_Y],    const double z[N_Z],  const double w[N_W], const double v[N_XYZW],double VV[N_XYZWP])
{
  int ix;
  int iy;
  int iz;
  int iw;
  int VV_tmp;
  unsigned mem1;
  unsigned mem2;


  mem1=N_XYZWP;
  mem2=N_X;

  for (iw= 0; iw<N_W; iw++ ) {
	  for (iz = 0; iz < N_Z; iz++) {
		for (iy = 0; iy < N_Y; iy++) {
		  memcpy(&VV[(iw * N_XYZP + iz * N_XYP + iy * N_XP) + N_XYZP+N_XYP+N_XP+1], &v[iw * N_Z*N_Y*N_X + iz * N_Y*N_X + iy
				 * N_X], mem2 * sizeof(double));
		}
	  }
  }
  // Boundary cases
  for (iw = 0; iw<N_WP; iw++ ) {
	  for (iz = 0; iz < N_ZP; iz++) {
		for (iy = 0; iy < N_YP; iy++) {
		  VV_tmp = N_XP * iy + N_XYP * iz + N_XYZP * iw;
		  VV[VV_tmp] = (3.0 * VV[1 + VV_tmp] - 3.0 * VV[2 + VV_tmp]) + VV[3 + VV_tmp];
		  VV[N_X+1 + VV_tmp] = (3.0 * VV[N_X + VV_tmp] - 3.0 * VV[N_X-1 + VV_tmp]) + VV[N_X-2 +
			VV_tmp];
		}
	  }
  }
  for (iw = 0; iw<N_WP; iw++ ) {
	  for (iz = 0; iz < N_ZP; iz++) {
		for (ix = 0; ix < N_XP; ix++) {
		  VV_tmp = ix + N_XYP * iz + N_XYZP * iw;
		  VV[VV_tmp] = (3.0 * VV[N_XP + VV_tmp] - 3.0 * VV[2*N_XP + VV_tmp]) + VV[3*N_XP  + VV_tmp];
		  VV[(N_Y+1)*N_XP + VV_tmp] = (3.0 * VV[N_XP*N_Y + VV_tmp] - 3.0 * VV[N_XP*(N_Y-1) + VV_tmp]) + VV[N_XP*(N_Y-2) + VV_tmp];
		}
	  }
  }
  for (iw = 0; iw<N_WP; iw++ ) {
	  for (iy = 0; iy < N_YP; iy++) {
		for (ix = 0; ix < N_XP; ix++) {
		  VV_tmp = ix + N_XP * iy + N_XYZP * iw ;
		  VV[VV_tmp] = (3.0 * VV[N_XYP + VV_tmp] - 3.0 * VV[N_XYP*2 + VV_tmp]) + VV[N_XYP*3 + VV_tmp];
		  VV[N_XYP*(N_Z+1) + VV_tmp] = (3.0 * VV[N_XYP*N_Z + VV_tmp] - 3.0 * VV[N_XYP*(N_Z-1) + VV_tmp]) + VV[N_XYP*(N_Z-2) + VV_tmp];
		}
	  }
  }
  for (iz = 0; iz<N_ZP; iz++ ) {
	  for (iy = 0; iy < N_YP; iy++) {
		for (ix = 0; ix < N_XP; ix++) {
		  VV_tmp = ix + N_XP * iy + N_XYP * iz ;
		  VV[VV_tmp] = (3.0 * VV[N_XYZP + VV_tmp] - 3.0 * VV[N_XYZP*2 + VV_tmp]) + VV[N_XYZP*3 + VV_tmp];
		  VV[N_XYZP*(N_W+1) + VV_tmp] = (3.0 * VV[N_XYZP*N_W + VV_tmp] - 3.0 * VV[N_XYZP*(N_W-1) + VV_tmp]) + VV[N_XYZP*(N_W-2) + VV_tmp];
		}
	  }
}
}
void interp4_interp_(const double x[N_X],  const double y[N_Y],    const double z[N_Z],  const double w[N_W], const double VV[N_XYZWP],
             const double xq[], const double yq[],   const double zq[], const double wq[],double vq[], const int *insiz)
{
  double dx;
  double dy;
  double dz;
  double dw;
  int VV_tmp;
  int k;
  int b_low_i;
  int a_low_i;
  int c_low_i;
  int d_low_i;
  int low_ip1;
  int mid_i;
  int high_i;
  double g;
  double t;
  double s;
  double h;
  double hh;
  double tt_tmp;
  int vik_tmp_tmp;
  int vik_tmp;
  double b_vik_tmp;
  double c_vik_tmp;
  double d_vik_tmp;
  double e_vik_tmp;
  double b_tt_tmp;
  double c_tt_tmp;
  int b_vik_tmp_tmp;
  int c_vik_tmp_tmp;
  double vik;
  double vik2;

  dx=(x[N_X-1]-x[0])/(N_X-1);
  dy=(y[N_Y-1]-y[0])/(N_Y-1);
  dz=(z[N_Z-1]-z[0])/(N_Z-1);
  dw=(w[N_W-1]-w[0])/(N_W-1);

  for (k = 0; k < *insiz; k++) {
	  // Finding a_low_i (high_i discarded)
	  a_low_i = 1;
      low_ip1 = 2;
      high_i = N_X;
      while (high_i > low_ip1) {
        mid_i = (a_low_i + high_i) >> 1;
        if (xq[k] >= x[mid_i - 1]) {
          a_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

	  // Finding b_low_i (high_i discarded)
      b_low_i = 1;
      low_ip1 = 2;
      high_i = N_Y;
      while (high_i > low_ip1) {
        mid_i = (b_low_i + high_i) >> 1;
        if (yq[k] >= y[mid_i - 1]) {
          b_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

	  // Finding c_low_i (high_i discarded)
      c_low_i = 1;
      low_ip1 = 2;
      high_i = N_Z;
      while (high_i > low_ip1) {
        mid_i = (c_low_i + high_i) >> 1;
        if (zq[k] >= z[mid_i - 1]) {
          c_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

	  // Finding d_low_i (high_i discarded)
      d_low_i = 1;
      low_ip1 = 2;
      high_i = N_W;
      while (high_i > low_ip1) {
        mid_i = (d_low_i + high_i) >> 1;
        if (wq[k] >= w[mid_i - 1]) {
          d_low_i = mid_i;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

	  s = (xq[k] - x[a_low_i - 1]) / dx;
      t = (yq[k] - y[b_low_i - 1]) / dy;
      h = (zq[k] - z[c_low_i - 1]) / dz;
	  g = (wq[k] - w[d_low_i - 1]) / dw;

      hh = ((2.0 - h) * h - 1.0) * h;
      tt_tmp = ((2.0 - t) * t - 1.0) * t;

      mid_i = N_XYP * (c_low_i - 1)+ N_XYZP * (d_low_i-1);;
      vik_tmp_tmp = a_low_i + N_XP * (b_low_i - 1);
      vik_tmp = vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i - 1) + N_XYP * (c_low_i - 1) + NXYZP * (d_low_i - 1);
      b_vik_tmp = ((2.0 - s) * s - 1.0) * s;
      c_vik_tmp = (3.0 * s - 5.0) * s * s + 2.0;
      d_vik_tmp = ((4.0 - 3.0 * s) * s + 1.0) * s;
      e_vik_tmp = (s - 1.0) * s * s;

      b_tt_tmp = (3.0 * t - 5.0) * t * t + 2.0;
      b_vik_tmp_tmp = a_low_i + N_XP * b_low_i;
      low_ip1 = b_vik_tmp_tmp + mid_i; // a_low_i + N_XP * b_low_i + N_XYP * (c_low_i - 1) + NXYZP * (d_low_i - 1);
      c_tt_tmp = ((4.0 - 3.0 * t) * t + 1.0) * t;
      c_vik_tmp_tmp = a_low_i + N_XP * (b_low_i + 1);
      high_i = c_vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i+1) + N_XYP * (c_low_i - 1) + NXYZP * (d_low_i - 1);

      vik  = VV[vik_tmp - 1] * tt_tmp   * b_vik_tmp * hh + VV[vik_tmp]     * tt_tmp   * c_vik_tmp * hh +
	         VV[vik_tmp + 1] * tt_tmp   * d_vik_tmp * hh + VV[vik_tmp + 2] * tt_tmp   * e_vik_tmp * hh ;
      vik += VV[low_ip1 - 1] * b_tt_tmp * b_vik_tmp * hh + VV[low_ip1]     * b_tt_tmp * c_vik_tmp * hh +
	         VV[low_ip1 + 1] * b_tt_tmp * d_vik_tmp * hh + VV[low_ip1 + 2] * b_tt_tmp * e_vik_tmp * hh ;

	  vik += VV[ high_i - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[high_i] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[high_i + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[high_i + 2] * c_tt_tmp * e_vik_tmp * hh;
      s = (t - 1.0) * t * t;
      low_ip1 = a_low_i + N_XP * (b_low_i + 2);
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (3.0 * h - 5.0) * h * h + 2.0;
	  mid_i = N_XYP * c_low_i + N_XYZP * (d_low_i-1);;
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = ((4.0 - 3.0 * h) * h + 1.0) * h;
      mid_i = N_XYP * (c_low_i + 1) + N_XYZP * (d_low_i-1);;
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (h - 1.0) * h * h;
      mid_i = N_XYP * (c_low_i + 2) + N_XYZP * (d_low_i-1);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

	  vik2 = vik*((2.0 - g) * g - 1.0) * g;



      hh = ((2.0 - h) * h - 1.0) * h;
      mid_i = N_XYP * (c_low_i - 1) + N_XYZP * d_low_i;
      vik_tmp = vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i - 1) + N_XYP * (c_low_i - 1) + N_XYZP * d_low_i;
      low_ip1 = b_vik_tmp_tmp + mid_i; // a_low_i + N_XP * b_low_i + N_XYP * (c_low_i - 1) + N_XYZP * d_low_i;
      high_i = c_vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i+1) + N_XYP * (c_low_i - 1) + N_XYZP * d_low_i;

      vik = VV[vik_tmp - 1] * tt_tmp   * b_vik_tmp * hh + VV[vik_tmp]     * tt_tmp   * c_vik_tmp * hh +
	         VV[vik_tmp + 1] * tt_tmp   * d_vik_tmp * hh + VV[vik_tmp + 2] * tt_tmp   * e_vik_tmp * hh ;
      vik += VV[low_ip1 - 1] * b_tt_tmp * b_vik_tmp * hh + VV[low_ip1]     * b_tt_tmp * c_vik_tmp * hh +
	         VV[low_ip1 + 1] * b_tt_tmp * d_vik_tmp * hh + VV[low_ip1 + 2] * b_tt_tmp * e_vik_tmp * hh ;

	  vik += VV[ high_i - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[high_i] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[high_i + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[high_i + 2] * c_tt_tmp * e_vik_tmp * hh;

      low_ip1 = a_low_i + N_XP * (b_low_i + 2);
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (3.0 * h - 5.0) * h * h + 2.0;
	  mid_i = N_XYP * c_low_i + N_XYZP * d_low_i;
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = ((4.0 - 3.0 * h) * h + 1.0) * h;
      mid_i = N_XYP * (c_low_i + 1) + N_XYZP * d_low_i;
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (h - 1.0) * h * h;
      mid_i = N_XYP * (c_low_i + 2) + N_XYZP * d_low_i;
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

	  vik2 += vik * ((3.0 * g - 5.0) * g * g + 2.0);




	  hh = ((2.0 - h) * h - 1.0) * h;
      mid_i = N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+1);
      vik_tmp = vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i - 1) + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+1);
      low_ip1 = b_vik_tmp_tmp + mid_i; // a_low_i + N_XP * b_low_i + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+1);
      high_i = c_vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i+1) + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+1);

      vik = VV[vik_tmp - 1] * tt_tmp   * b_vik_tmp * hh + VV[vik_tmp]     * tt_tmp   * c_vik_tmp * hh +
	         VV[vik_tmp + 1] * tt_tmp   * d_vik_tmp * hh + VV[vik_tmp + 2] * tt_tmp   * e_vik_tmp * hh ;
      vik += VV[low_ip1 - 1] * b_tt_tmp * b_vik_tmp * hh + VV[low_ip1]     * b_tt_tmp * c_vik_tmp * hh +
	         VV[low_ip1 + 1] * b_tt_tmp * d_vik_tmp * hh + VV[low_ip1 + 2] * b_tt_tmp * e_vik_tmp * hh ;

	  vik += VV[ high_i - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[high_i] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[high_i + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[high_i + 2] * c_tt_tmp * e_vik_tmp * hh;

      low_ip1 = a_low_i + N_XP * (b_low_i + 2);
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (3.0 * h - 5.0) * h * h + 2.0;
	  mid_i = N_XYP * c_low_i + N_XYZP * (d_low_i+1);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = ((4.0 - 3.0 * h) * h + 1.0) * h;
      mid_i = N_XYP * (c_low_i + 1) + N_XYZP * (d_low_i+1);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (h - 1.0) * h * h;
      mid_i = N_XYP * (c_low_i + 2) + N_XYZP * (d_low_i+1);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

	  vik2 += vik * ((4.0 - 3.0 * g) * g + 1.0) * g;



	  hh = ((2.0 - h) * h - 1.0) * h;
      mid_i = N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+2);
      vik_tmp = vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i - 1) + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+2);
      low_ip1 = b_vik_tmp_tmp + mid_i; // a_low_i + N_XP * b_low_i + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+2);
      high_i = c_vik_tmp_tmp + mid_i; // a_low_i + N_XP * (b_low_i+1) + N_XYP * (c_low_i - 1) + N_XYZP * (d_low_i+2);

      vik = VV[vik_tmp - 1] * tt_tmp   * b_vik_tmp * hh + VV[vik_tmp]     * tt_tmp   * c_vik_tmp * hh +
	         VV[vik_tmp + 1] * tt_tmp   * d_vik_tmp * hh + VV[vik_tmp + 2] * tt_tmp   * e_vik_tmp * hh ;
      vik += VV[low_ip1 - 1] * b_tt_tmp * b_vik_tmp * hh + VV[low_ip1]     * b_tt_tmp * c_vik_tmp * hh +
	         VV[low_ip1 + 1] * b_tt_tmp * d_vik_tmp * hh + VV[low_ip1 + 2] * b_tt_tmp * e_vik_tmp * hh ;

	  vik += VV[ high_i - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[high_i] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[high_i + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[high_i + 2] * c_tt_tmp * e_vik_tmp * hh;

      low_ip1 = a_low_i + N_XP * (b_low_i + 2);
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (3.0 * h - 5.0) * h * h + 2.0;
	  mid_i = N_XYP * c_low_i + N_XYZP * (d_low_i+2);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = ((4.0 - 3.0 * h) * h + 1.0) * h;
      mid_i = N_XYP * (c_low_i + 1) + N_XYZP * (d_low_i+2);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

      hh = (h - 1.0) * h * h;
      mid_i = N_XYP * (c_low_i + 2) + N_XYZP * (d_low_i+2);
      vik_tmp = vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * tt_tmp * e_vik_tmp * hh;
      vik_tmp = b_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * b_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * b_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * b_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * b_tt_tmp * e_vik_tmp * hh;
      vik_tmp = c_vik_tmp_tmp + mid_i;
      vik += VV[vik_tmp - 1] * c_tt_tmp * b_vik_tmp * hh;
      vik += VV[vik_tmp] * c_tt_tmp * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * c_tt_tmp * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * c_tt_tmp * e_vik_tmp * hh;
      vik_tmp = low_ip1 + mid_i;
      vik += VV[vik_tmp - 1] * s * b_vik_tmp * hh;
      vik += VV[vik_tmp] * s * c_vik_tmp * hh;
      vik += VV[vik_tmp + 1] * s * d_vik_tmp * hh;
      vik += VV[vik_tmp + 2] * s * e_vik_tmp * hh;

	  vik2 += vik * (g - 1.0) * g * g;

      vq[k] = vik2 / 16.0;
  }
}
