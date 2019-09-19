#ifndef HAAR_C_H
#define HAAR_C_H

#include <iostream>
#include <cmath>
#include <omp.h>

void haar_2d_c(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator);
void haar_2d_inverse_c(int m, int n, double u[], int n_level);

#endif
