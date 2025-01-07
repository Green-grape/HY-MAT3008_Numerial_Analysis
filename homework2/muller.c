#include <math.h>
#include <stdio.h>

// muller method for root finding
float rtmuller(float(*func)(float), float x1, float x2, float xacc, int* iter_cnt){
    float x3, f1, f2, f3,x4;
    float temp=x2;
    x2=(x1+x2)/2;
    x3=temp;
    f1=func(x1);
    f2=func(x2);
    f3=func(x3);
    for(int j=1;j<=100;j++){
        *iter_cnt = j;
        float c = func(x3);
        float b=((x1-x3)*(x1-x3)*(f2-f3)-(x2-x3)*(x2-x3)*(f1-f3))/((x1-x3)*(x2-x3)*(x1-x2));
        float a=((x2-x3)*(f1-f3)-(x1-x3)*(f2-f3))/((x1-x3)*(x2-x3)*(x1-x2));
        int sign_b = b>0?1:-1;
        float x4=x3-2*c/(b+sign_b*sqrt(b*b-4*a*c));
        if(fabs(x4-x3)<xacc || c==0.0) return x4;
        x1=x2;
        x2=x3;
        x3=x4;
    }
    return 0.0;
}