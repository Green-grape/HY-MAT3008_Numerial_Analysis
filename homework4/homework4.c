#include "nr.h"
#include <time.h>

long seed=12345;

float uniform_random(float a, float b){
    return a + (b-a)*ran1(&seed);
}

float normal_random(float mean, float std_dev){
    return mean + std_dev*gasdev(&seed);
}

void save_array(float *array, int N, char *file_name){
    FILE *file = fopen(file_name, "w");
    for(int i=0; i<N; i++){
        fprintf(file, "%f\n", array[i]);
    }
    fclose(file);
}

int main(void){
    float a=-3.0, b=4, mean=0.5, std_dev=1.5;
    int sampling[] = {100,1000,10000, 100000};
    char file_name_uni[50], file_name_norm[50];
    for(int i=0; i<3; i++){
        int N = sampling[i];
        float uniform_samples[N], normal_samples[N];
        for(int j=0; j<N; j++){
            uniform_samples[j] = uniform_random(a, b);
            normal_samples[j] = normal_random(mean, std_dev);
        }
        sprintf(file_name_uni, "uniform_samples_%d.dat", N);
        sprintf(file_name_norm, "normal_samples_%d.dat", N);
        save_array(uniform_samples, N, file_name_uni);
        save_array(normal_samples, N, file_name_norm);
    }
    return 0;
}