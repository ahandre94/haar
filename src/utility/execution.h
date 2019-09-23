#ifndef EXECUTION_H
#define EXECUTION_H

#include <opencv2/core.hpp>
#include <iostream>
#include <string>
#include <omp.h>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "dimension.h"
#include "load_img.h"
#include "../sequential/haar.h"
#include "../sequential/threshold.h"
#include "../parallel/parallel_haar.h"
#include "../parallel/parallel_threshold.h"

using namespace cv;
using namespace std;

void exeggutor(vector<Mat> images, 
			   int n_level,
			   int n_level_stop, 
			   string type,
			   double *time,
			   int iterator,
			   void (*haar)(int, int, Mat, int, double*, int), 
			   void (*visualizza)(int, int, double*, int), 
			   void (*inverse)(int, int, Mat, int, double*, int), 
			   void (*threshold)(int, int, double*, int, double), 
			   double (*mean)(int, int, double*, int),
			   double *total_time,
			   double *total_time_inverse,
			   double *threshold_time);
void execution_parallel_mix_sequential(int n_level_to_do, int n_level_stop, double *time, double *total_time);

#endif
