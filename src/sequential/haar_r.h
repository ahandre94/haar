#ifndef HAAR_R_H
#define HAAR_R_H

#include <iostream>
#include <cmath>
#include <omp.h>

void haar_2d_r(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator);
void haar_2d_inverse_r(int m, int n, double u[], int n_level);

#endif
