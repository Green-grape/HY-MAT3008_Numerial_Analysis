#ifndef CUSTOM_HEADER
#define CUSTOM_HEADER

void ludcmp_double(double **a, int n, int *indx, double *d);
void lubksb_double(double **a, int n, int *indx, double b[]);
void mprove_double(double **a, double **alud, int n, int indx[], double b[], double x[]);

#endif CUSTOM_HEADER