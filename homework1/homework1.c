#include "nr.h"
#include "myheader.h"

double get_double_eps(){
    double one = 1.0;
    double eps = 1.0;
    double two = 2.0;
    while (one + eps != one){
        eps /= two;
    }
    return eps*two;
}

float get_eps(){
    float one = 1.0;
    float eps = 1.0;
    float two = 2.0;
    while (one + eps != one){
        eps /=two;
    }
    return eps*two;
}

int main(){
    int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
    float eps, epsneg, xmin, xmax;
    double eps_double, epsneg_double, xmin_double, xmax_double;

    machar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax);
    machar_double(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &eps_double, &epsneg_double, &xmin_double, &xmax_double);
    double custom_double_eps= get_double_eps();
    float custom_eps = get_eps();
    printf("Accuracy of floating point number: %0.20f\n", eps);
    printf("Accuracy of double floating point number: %0.20f\n", eps_double);
    printf("Accuracy of custom floating point number: %0.20f\n", custom_eps);
    printf("Accuracy of custom double floating point number: %0.20f\n", custom_double_eps);
    return 0;
}

