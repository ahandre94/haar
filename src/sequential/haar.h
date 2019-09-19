#ifndef HAAR_H
#define HAAR_H

#include <opencv2/core.hpp>
#include "haar_q.h"
#include "haar_c.h"
#include "haar_r.h"

using namespace cv;

void haar(int m, int n, Mat image, int n_level, double *time, int iterator); 
void haar_inverse(int m, int n, Mat image, int n_level);
void visualizza_haar(int m, int n, double u[], int n_level);

#endif
