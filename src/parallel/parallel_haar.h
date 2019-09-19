#ifndef PARALLEL_HAAR_H
#define PARALLEL_HAAR_H

#include <opencv2/core.hpp>
#include <omp.h>
#include "parallel_haar_q.h"
#include "parallel_haar_c.h"
#include "parallel_haar_r.h"

using namespace cv;

void p_haar(int m, int n, Mat image, int n_level, double *time, int iterator);
void p_haar_inverse(int m, int n, Mat image, int n_level);
void p_visualizza_haar(int m, int n, double u[], int n_level);

#endif
