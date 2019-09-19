#ifndef PARALLEL_THRESHOLD_H
#define PARALLEL_THRESHOLD_H

#include <iostream>
#include <cmath>
#include <omp.h>

void p_threshold(int m, int n, double u[], int n_level, double threshold);
double p_mean(int m, int n, double u[], int n_level);

#endif
