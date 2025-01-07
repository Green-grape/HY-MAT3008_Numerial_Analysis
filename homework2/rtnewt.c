
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define JMAX 20

float rtnewt_t(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc, int* iter_cnt)
{
	void nrerror(char error_text[]);
	int j;
	float df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=JMAX;j++) {
		*iter_cnt = j;
		(*funcd)(rtn,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) return rtn;
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}
#undef JMAX
