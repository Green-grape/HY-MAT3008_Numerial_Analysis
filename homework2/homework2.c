#include "nr.h"
#include "myheader.h"
#include <stdlib.h>

// 일반적인 함수에서 답을 찾는 함수
float* find_root(float (*target_func)(float), float(*solve_func)(float(*)(float), float, float, float, int*),float start, float end, float eps, int* root_cnt, int *iter_cnt){
    int nb=100, iter=100, nb_length=nb+1;
    float xb1[nb_length], xb2[nb_length];
    zbrak(target_func, start, end, iter, xb1, xb2, &nb);
    float* result = (float*)malloc(sizeof(float)*(nb));
    for(int i=1; i<=nb; i++){
        result[i-1] = solve_func(target_func, xb1[i], xb2[i], eps, iter_cnt);
    }
    *root_cnt = nb;
    return result;
}

// Newton's method 종류의 함수에서 답을 찾는 함수
float* find_root_nt(float (*target_func)(float), void(*target_func_with_diff)(float, float*, float*), float(*solve_func)(void (*)(float, float *, float *), float, float, float, int*),float start, float end, float eps, int* root_cnt, int *iter_cnt){
    int nb=100, iter=100, nb_length=nb+1;
    float xb1[nb_length], xb2[nb_length];
    zbrak(target_func, start, end, iter, xb1, xb2, &nb);
    float* result = (float*)malloc(sizeof(float)*(nb));
    for(int i=1; i<=nb; i++){
        result[i-1] = solve_func(target_func_with_diff, xb1[i], xb2[i], eps, iter_cnt);
    }
    *root_cnt = nb;
    return result;
}

void besse(float x, float *y, float *dy){
    *y = bessj0(x);
    *dy = -bessj1(x);
}

float target_legendre0(float x){
    return legendre0(10, x);
}

void legendre(float x, float *y, float *dy){
    *y = legendre0(10, x);
    *dy = legendre1(10, x);
}

int main(void){
    // find possible return values for the following function
    float start=1.0, end=10.0, eps=1e-6;
    int iter_cnt=0;
    float(*solve_funcs[])(float(*)(float), float, float, float, int*)={rtbis_t, rtsec_t, rtflsp_t, rtmuller};
    char* solve_funcs_name[]={"rtbis_t", "rtsec_t", "rtflsp_t", "rtmuller"};
    float(*solve_funcs_nt[])(void (*)(float, float *, float *), float, float, float, int*)={rtsafe_t, rtnewt_t};
    char* solve_funcs_nt_name[]={"rtsafe_t", "rtnewt_t"};
    printf("Solve Bessel Function\n");
    for(int i=0;i<4;i++){
        int nb=-1;
        float* result = find_root(bessj0, solve_funcs[i], start, end, eps, &nb, &iter_cnt);
        for(int j=0; j<nb; j++){
            printf("%s Root %d: %f with iter_cnt:%d\n", solve_funcs_name[i],j, result[j], iter_cnt);
        }
        iter_cnt=0;
    }
    for(int i=0;i<2;i++){
        int nb=-1;
        float* result = find_root_nt(bessj0, besse, solve_funcs_nt[i], start, end, eps, &nb, &iter_cnt);
        for(int j=0; j<nb; j++){
            printf("%s Root %d: %f with iter_cnt:%d\n", solve_funcs_nt_name[i],j, result[j], iter_cnt);
        }
        iter_cnt=0;
    }
    printf("Solve Legendre Polynomial\n");
    int nb=-1;
    start=-1.0, end=1.0;
    float* result = find_root_nt(target_legendre0, legendre, rtsafe_t,start,end, eps, &nb, &iter_cnt);
    for(int j=0; j<nb; j++){
        printf("rtsafe_t Root %d: %f with iter_cnt:%d\n",j, result[j], iter_cnt);
    }
    return 0;
}