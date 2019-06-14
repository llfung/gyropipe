#include <string.h>
#include <math.h>
#include "interp1in.h"

/* Function Definitions */
void interp1in_(const double x[ONE_N_X], const double y[ONE_N_X], const double xi[], double
  yi[], const int *insiz)
{
  double minx;
  double secx;
  double penx;
  double maxx;
  double h;
  int k;
  int b_k;
  double s;
  double sd2;
  double ssd2;
  double s3m4;
  double c0;
  double c1;
  minx = x[0];
  secx = x[1];
  penx = x[ONE_N_X-2];
  maxx = x[ONE_N_X-1];
  h = (x[ONE_N_X-1] - x[0]) / (double)(ONE_N_X-1);

	for (k = 0; k < *insiz; k++) {
	  if ((xi[k] >= minx) && (xi[k] <= maxx)) {
		if (xi[k] < secx) {
		  b_k = 1;
		  s = (xi[k] - minx) / h;
		} else if (xi[k] >= penx) {
		  b_k = ONE_N_X-1;
		  s = (xi[k] - penx) / h;
		} else {
		  b_k = (int)(1.0 + floor((xi[k] - minx) / h));
		  s = (xi[k] - x[b_k - 1]) / h;
		}

		sd2 = s / 2.0;
		ssd2 = s * sd2;
		s3m4 = 3.0 * s - 4.0;
		c0 = -s * (s * (sd2 - 1.0) + 0.5);
		c1 = ssd2 * (s3m4 - 1.0) + 1.0;
		s3m4 = -sd2 * (s * s3m4 - 1.0);
		sd2 = ssd2 * (s - 1.0);
		if (b_k < 2) {
		  yi[k] = ((c0 * ((3.0 * y[b_k - 1] - 3.0 * y[b_k]) + y[b_k + 1]) + c1 *
					y[b_k - 1]) + s3m4 * y[b_k]) + sd2 * y[b_k + 1];
		} else if (b_k < (ONE_N_X-1)) {
		  yi[k] = ((c0 * y[b_k - 2] + c1 * y[b_k - 1]) + s3m4 * y[b_k]) + sd2 *
			y[b_k + 1];
		} else {
		  yi[k] = ((c0 * y[b_k - 2] + c1 * y[b_k - 1]) + s3m4 * y[b_k]) + sd2 *
			((3.0 * y[b_k] - 3.0 * y[b_k - 1]) + y[b_k - 2]);
		}
	  }
	}
}
