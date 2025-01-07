#include <stdio.h>
#include <stdlib.h>
#include "nr.h"
#include "myheader.h"

void read_linear_equation(float ***A, float **b, int* n, FILE* fp);
void solve_linear_equation_gaussj(float** A, float* b, int n, float* x);
void solve_linear_equation_lu(float** A, float* b, int n, float* x);
void solve_linear_equation_svd(float** A, float* b, int n, float* x);
void solve_linear_equation_mprove(float** A, float* b, int n, float* x);
void solve_linear_equation_lu_double(double** A, double* b, int n, double* x);
void solve_linear_equation_mprove_double(double** A, double* b, int n, double* x);
void read_linear_equation_double(double*** A, double** b, int* n, FILE* fp);
double get_determinant_double(double** A, int n);
float get_determinant(float** A, int n);
float** get_inverse(float** A, int n);

int main(){
    const char* filenames[]={"./linear_equation/lineq1.dat", "./linear_equation/lineq2.dat", "./linear_equation/lineq3.dat"};
    const char* solve_names[]={"Gauss-Jordan", "LU", "SVD", "M-Prove"};
    const char* solve_names_double[]={"Double LU", "Double M-Prove"};
    void(*solve_funcs[])(float** A, float* b, int n, float* x)={solve_linear_equation_gaussj, solve_linear_equation_lu, solve_linear_equation_svd, solve_linear_equation_mprove};
    void(*solve_funcs_double[])(double** A, double* b, int n, double* x)={solve_linear_equation_lu_double, solve_linear_equation_mprove_double};
    FILE* fp;
    int file_cnt=3;
    for(int i=0; i<file_cnt; i++){
        fp=fopen(filenames[i], "r");
        if(fp==NULL){
            printf("File open error\n");
            return 1;
        }
        float **A, *b;
        int n;
        read_linear_equation(&A, &b, &n, fp);
        printf("Linear Equation %d\n", i+1);
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                printf("%f ", A[j][k]);
            }
            printf("= %f\n", b[j]);
        }
        float *x;
        x=(float*)malloc(sizeof(float)*n);
        for(int j=0; j<4; j++){
            solve_funcs[j](A, b, n, x);
            printf("Solution with %s\n", solve_names[j]);
            for(int k=0; k<n; k++){
                printf("%.20f ", x[k]);
            }
            printf("\n");
        }
        rewind(fp);
        double **A_double, *b_double;
        read_linear_equation_double(&A_double, &b_double, &n, fp);
        printf("Linear Equation Double %d\n", i+1);
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                printf("%lf ", A_double[j][k]);
            }
            printf("= %lf\n", b_double[j]);
        }
        double *x_double;
        x_double=(double*)malloc(sizeof(double)*n);
        for(int j=0; j<2; j++){
            solve_funcs_double[j](A_double, b_double, n, x_double);
            printf("Solution with %s\n", solve_names_double[j]);
            for(int k=0; k<n; k++){
                printf("%.20f ", x_double[k]);
            }
            printf("\n");
        }
        free(x_double);
        free(x);
        for(int j=0; j<n; j++){
            free(A[j]);
            free(A_double[j]);
        }
        free(A);
        free(A_double);
        free(b);
        free(b_double);
        fclose(fp);
    }
    for(int i=0; i<file_cnt; i++){
        fp=fopen(filenames[i], "r");
        if(fp==NULL){
            printf("File open error\n");
            return 1;
        }
        float **A, *b;
        int n;
        read_linear_equation(&A, &b, &n, fp);
        printf("Linear Equation %d\n", i+1);
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                printf("%f ", A[j][k]);
            }
            printf("= %f\n", b[j]);
        }
        float **A_inv;
        printf("Determinant: %f\n", get_determinant(A, n));
        A_inv=get_inverse(A, n);
        if(A_inv==NULL){
            printf("Can't find inverse matrix\n");
            continue;
        }
        printf("Inverse Matrix\n");
        for(int j=0; j<n; j++){
            for(int k=0; k<n; k++){
                printf("%f ", A_inv[j][k]);
            }
            printf("\n");
        }
        fclose(fp);
        for(int j=0; j<n; j++){
            free(A[j]);
            free(A_inv[j]);
        }
        free(A_inv);
        free(A);
        free(b);
    }
    return 0;
}

void read_linear_equation(float ***A, float **b, int* n, FILE* fp){
    fscanf(fp, "%d", n);
    (*A)=(float**)malloc(sizeof(float*)*(*n));
    int m;
    fscanf(fp, "%d", &m);
    for(int i=0; i<(*n); i++){
        (*A)[i]=(float*)malloc(sizeof(float)*(*n));
        for(int j=0; j<(*n); j++){
            fscanf(fp, "%f", &(*A)[i][j]);
        }
    }
    (*b)=(float*)malloc(sizeof(float)*(*n));
    for(int i=0; i<(*n); i++){
        fscanf(fp, "%f", &(*b)[i]);
    }
}

void read_linear_equation_double(double ***A, double **b, int* n, FILE* fp){
    fscanf(fp, "%d", n);
    (*A)=(double**)malloc(sizeof(double*)*(*n));
    int m;
    fscanf(fp, "%d", &m);
    for(int i=0; i<(*n); i++){
        (*A)[i]=(double*)malloc(sizeof(double)*(*n));
        for(int j=0; j<(*n); j++){
            fscanf(fp, "%lf", &(*A)[i][j]);
        }
    }
    (*b)=(double*)malloc(sizeof(double)*(*n));
    for(int i=0; i<(*n); i++){
        fscanf(fp, "%lf", &(*b)[i]);
    }
}


void solve_linear_equation_gaussj(float **A, float *b, int n, float *x){
    if(get_determinant(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
    }
    float **A_f, **b_f;
    A_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    b_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        b_f[i]=(float*)malloc(sizeof(float)*2);
    }
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i][1]=b[i-1];
    }
    gaussj(A_f, n, b_f, 1);
    for(int i=1; i<=n; i++){
        x[i-1]=b_f[i][1];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
    }
    free(A_f);
    free(b_f);
}

void solve_linear_equation_lu(float **A, float *b, int n, float *x){
    if(get_determinant(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
    }
    float **A_f, *b_f, *d;
    int *indx;
    A_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    b_f=(float*)malloc(sizeof(float)*(n+1));
    d=(float*)malloc(sizeof(float)*(n+1));
    indx=(int*)malloc(sizeof(int)*(n+1));
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i]=b[i-1];
    }
    ludcmp(A_f, n, indx, d);
    lubksb(A_f, n, indx, b_f);
    for(int i=1; i<=n; i++){
        x[i-1]=b_f[i];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
    }
    free(A_f);
    free(b_f);
    free(d);
}

void solve_linear_equation_svd(float **A, float *b, int n, float *x){
    if(get_determinant(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
    }
    float **A_f, *b_f, *w, **v, *x_f;
    A_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    b_f=(float*)malloc(sizeof(float)*(n+1));
    w=(float*)malloc(sizeof(float)*(n+1));
    v=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        v[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i]=b[i-1];
    }
    x_f=(float*)malloc(sizeof(float)*(n+1));
    for(int i=1; i<=n; i++){
        x_f[i]=0.0;
    }
    svdcmp(A_f, n, n, w, v);
    svbksb(A_f, w, v, n, n, b_f, x_f);
    for(int i=1; i<=n; i++){
        x[i-1]=x_f[i];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
        free(v[i]);
    }
    free(A_f);
    free(b_f);
    free(w);
    free(v);
    free(x_f);
}

void solve_linear_equation_mprove(float **A, float *b, int n, float *x){
    if(get_determinant(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
    }
    FILE *fp=fopen("mprove.txt", "w");
    float **A_f, **alud, *b_f, *x_f;
    int *indx;
    A_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    b_f=(float*)malloc(sizeof(float)*(n+1));
    alud=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        alud[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    indx=(int*)malloc(sizeof(int)*(n+1));
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i]=b[i-1];
    }
    x_f=(float*)malloc(sizeof(float)*(n+1));
    for(int i=1; i<=n; i++){
        x_f[i]=0.0;
    }
    int iter=20;
    for(int i=0;i<iter;i++){
        ludcmp(A_f, n, indx, b_f);
        for(int j=1; j<=n; j++){
            for(int k=1; k<=n; k++){
                alud[j][k]=A_f[j][k];
                A_f[j][k]=A[j-1][k-1];
            }
        }
        mprove(A_f, alud, n, indx, b_f, x_f);
        for(int j=1;j<=n;j++){
            fprintf(fp, "%.20f ", x_f[j]);
        }
        fprintf(fp, "\n");
    }
    for(int i=1; i<=n; i++){
        x[i-1]=x_f[i];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
        free(alud[i]);
    }
    free(A_f);
    free(b_f);
    free(alud);
    free(indx);
    free(x_f);
    fclose(fp);
}

void solve_linear_equation_lu_double(double **A, double *b, int n, double* x){
    if(get_determinant_double(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
        return;
    }
    double **A_f, *b_f, *d;
    int *indx;
    A_f=(double**)malloc(sizeof(double*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(double*)malloc(sizeof(double)*(n+1));
    }
    b_f=(double*)malloc(sizeof(double)*(n+1));
    d=(double*)malloc(sizeof(double)*(n+1));
    indx=(int*)malloc(sizeof(int)*(n+1));
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i]=b[i-1];
    }
    ludcmp_double(A_f, n, indx, d);
    lubksb_double(A_f, n, indx, b_f);
    for(int i=1; i<=n; i++){
        x[i-1]=b_f[i];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
    }
    free(A_f);
    free(b_f);
    free(d);
}

void solve_linear_equation_mprove_double(double **A, double *b, int n, double* x){
    if(get_determinant_double(A, n)==0){
        printf("This is Singular matrix, So we can have many solutions\n");
        return;
    }
    double **A_f, **alud, *b_f, *x_f;
    int *indx;
    FILE *fp=fopen("mprove_double.txt", "w");
    A_f=(double**)malloc(sizeof(double*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(double*)malloc(sizeof(double)*(n+1));
    }
    b_f=(double*)malloc(sizeof(double)*(n+1));
    alud=(double**)malloc(sizeof(double*)*(n+1));
    for(int i=1; i<=n; i++){
        alud[i]=(double*)malloc(sizeof(double)*(n+1));
    }
    indx=(int*)malloc(sizeof(int)*(n+1));
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
        b_f[i]=b[i-1];
    }
    x_f=(double*)malloc(sizeof(double)*(n+1));
    for(int i=1; i<=n; i++){
        x_f[i]=0.0;
    }
    int iter=20;
    for(int i=0;i<iter;i++){
        ludcmp_double(A_f, n, indx, b_f);
        for(int i=1; i<=n; i++){
            for(int j=1; j<=n; j++){
                alud[i][j]=A_f[i][j];
                A_f[i][j]=A[i-1][j-1];
            }
        }
        mprove_double(A_f, alud, n, indx, b_f, x_f);
        for(int j=1;j<=n;j++){
            fprintf(fp, "%.20f ", x_f[j]);
        }
        fprintf(fp, "\n");
    }
    for(int i=1; i<=n; i++){
        x[i-1]=x_f[i];
    }
    for(int i=1; i<=n; i++){
        free(A_f[i]);
        free(alud[i]);
    }
    free(A_f);
    free(b_f);
    free(alud);
    free(indx);
    free(x_f);
    fclose(fp);
}

float get_determinant(float **A, int n){
    if(n==2){
        return A[0][0]*A[1][1]-A[0][1]*A[1][0];
    }else{
        float det=0.0;
        for(int i=0; i<n; i++){
            float **A_sub;
            A_sub=(float**)malloc(sizeof(float*)*(n-1));
            for(int j=0; j<n-1; j++){
                A_sub[j]=(float*)malloc(sizeof(float)*(n-1));
            }
            for(int j=0; j<n-1; j++){
                for(int k=0; k<n-1; k++){
                    if(k<i){
                        A_sub[j][k]=A[j+1][k];
                    }else{
                        A_sub[j][k]=A[j+1][k+1];
                    }
                }
            }
            det+=A[0][i]*get_determinant(A_sub, n-1)*(i%2==0?1:-1);
            for(int j=0; j<n-1; j++){
                free(A_sub[j]);
            }
            free(A_sub);
        }
        return det;
    }
}


double get_determinant_double(double **A, int n){
    if(n==2){
        return A[0][0]*A[1][1]-A[0][1]*A[1][0];
    }else{
        double det=0.0;
        for(int i=0; i<n; i++){
            double **A_sub;
            A_sub=(double**)malloc(sizeof(double*)*(n-1));
            for(int j=0; j<n-1; j++){
                A_sub[j]=(double*)malloc(sizeof(double)*(n-1));
            }
            for(int j=0; j<n-1; j++){
                for(int k=0; k<n-1; k++){
                    if(k<i){
                        A_sub[j][k]=A[j+1][k];
                    }else{
                        A_sub[j][k]=A[j+1][k+1];
                    }
                }
            }
            det+=A[0][i]*get_determinant_double(A_sub, n-1)*(i%2==0?1:-1);
            for(int j=0; j<n-1; j++){
                free(A_sub[j]);
            }
            free(A_sub);
        }
        return det;
    }
}

float** get_inverse(float **A, int n){
    if(get_determinant(A, n)==0){
        return NULL;
    }
    float **A_f, **A_inv_f, **A_inv;
    A_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_f[i][j]=A[i-1][j-1];
        }
    }
    A_inv_f=(float**)malloc(sizeof(float*)*(n+1));
    for(int i=1; i<=n; i++){
        A_inv_f[i]=(float*)malloc(sizeof(float)*(n+1));
    }
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_inv_f[i][j]=(i==j)?1.0:0.0;
        }
    }
    gaussj(A_f, n, A_inv_f, n);
    A_inv=(float**)malloc(sizeof(float*)*(n));
    for(int i=0; i<n; i++){
        A_inv[i]=(float*)malloc(sizeof(float)*(n));
    }
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            A_inv[i-1][j-1]=A_inv_f[i][j];
        }
    }
    for(int i=1; i<=n; i++){
        free(A_inv_f[i]);
    }
    free(A_inv_f);
    return A_inv;
}