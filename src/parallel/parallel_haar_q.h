#ifndef PARALLEL_HAAR_Q_H
#define PARALLEL_HAAR_Q_H

#include <iostream>
#include <cmath>
#include <omp.h>

void p_haar_2d(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator);
void p_haar_2d_inverse(int m, int n, double u[], int n_level);

#endif
