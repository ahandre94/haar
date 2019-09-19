#ifndef HAAR_Q_H
#define HAAR_Q_H

#include <iostream>
#include <cmath>
#include <omp.h>

void haar_2d(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator);
void haar_2d_inverse(int m, int n, double u[], int n_level);

#endif
