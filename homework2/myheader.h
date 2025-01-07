#ifndef CUSTOM_HEADER
#define CUSTOM_HEADER

float rtbis_t(float (*func)(float), float x1, float x2, float xacc, int*);
float rtsec_t(float (*func)(float), float x1, float x2, float xacc, int*);
float rtflsp_t(float (*func)(float), float x1, float x2, float xacc, int*);
float rtnewt_t(void (*funcd)(float, float *, float *), float x1, float x2, float xacc, int*);
float rtsafe_t(void (*funcd)(float, float *, float *), float x1, float x2, float xacc, int*);
float rtmuller(float(*func)(float), float x1, float x2, float xacc, int*);
float legendre0(int n, float x);
float legendre1(int n, float x);

#endif