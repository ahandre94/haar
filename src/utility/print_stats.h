#ifndef PRINT_STATS_H
#define PRINT_STATS_H

#include <opencv2/core.hpp>
#include <iostream>
#include <string>
#include <omp.h>

using namespace cv;
using namespace std;

void print_stats(vector<Mat> images, 
				 int n_level, 
				 int n_level_stop, 
				 double seq, 
				 double prl, 
				 double *total_time_seq, 
				 double *total_time_seq_inverse, 
				 double *total_time_par, 
				 double *total_time_par_inverse,
				 double speed_up,
				 double p,
				 double total_seq,
				 double total_par);

#endif
