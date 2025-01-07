// Legendre Polynomial

float legendre0(int n, float x){
    if(n==0) return 1.0;
    if(n==1) return x;
    return ((2.0*n-1.0)*x*legendre0(n-1, x) - (n-1)*legendre0(n-2, x))/(float)n;
}

// Legendre 다항식의 일차 미분 함수
float legendre1(int n, float x) {
    if (n == 0) return 0.0;
    return (n / (1 - x * x)) * (x * legendre0(n, x) - legendre0(n - 1, x));
}