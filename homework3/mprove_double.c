
#define NRANSI
#include "nrutil.h"
#include "myheader.h"
#include <stdlib.h>

void mprove_double(double **a, double **alud, int n, int indx[], double b[], double x[])
{
	void lubksb_double(double **a, int n, int *indx, double b[]);
	int j,i;
	double sdp;
	double *r;

	r=(double*)malloc(sizeof(double)*(n+1));
	for (i=1;i<=n;i++) {
		sdp = -b[i];
		for (j=1;j<=n;j++) sdp += a[i][j]*x[j];
		r[i]=sdp;
	}
	lubksb_double(alud,n,indx,r);
	for (i=1;i<=n;i++) x[i] -= r[i];
	free(r);
}
#undef NRANSI
