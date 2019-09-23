#include "utility/execution.h"
#include "utility/print_stats.h"

#define SHOW_IMAGE(x) {namedWindow(#x, WINDOW_AUTOSIZE); imshow(#x, x); waitKey(0);}

using namespace cv;
using namespace std;

int main(int argc, char **argv)
{
	vector<Mat> images = load_img("4k");
	
	int n_level = 5;
	int n_level_stop = 5;
	int stop = n_level - n_level_stop;

	double start, end, seq, prl;
	double speed_up, p;
	int iterator = 0;
	double *time = (double *)malloc(n_level * sizeof(double));
	double *total_time_seq = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_seq_inverse = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par_first = (double *)malloc(images.size() * n_level_stop * sizeof(double));
	double *total_time_par_second = (double *)malloc(images.size() * stop * sizeof(double));
	double *total_time_par_inverse = (double *)malloc(images.size() * n_level_stop * sizeof(double));
	double threshold_time = 0.0;
	double total_seq = 0.0, total_par = 0.0;
	
	// sequential execution
	start = omp_get_wtime();
	exeggutor(images, n_level, n_level_stop, "seq", time, iterator, haar, visualizza_haar, haar_inverse, threshold, mean, total_time_seq, total_time_seq_inverse, &threshold_time);
	end = omp_get_wtime();
	seq = end - start;
	
	// parallel execution
	start = omp_get_wtime();
	exeggutor(images, n_level_stop, n_level, "par", time, iterator, p_haar, p_visualizza_haar, p_haar_inverse, p_threshold, p_mean, total_time_par_first, total_time_par_inverse, &threshold_time);
	execution_parallel_mix_sequential(stop, n_level_stop, time, total_time_par_second);
	end = omp_get_wtime();
	prl = end - start;

	for (unsigned i = 0; i < images.size(); i++)
	{
		for (int j = 0; j < n_level_stop; j++)
		{
			total_time_par[j + i * n_level] = total_time_par_first[j + i * n_level_stop];
		}
		for (int k = n_level_stop; k < n_level; k++)
		{
			total_time_par[k + i * n_level] = total_time_par_second[k - n_level_stop + i * stop];
		}

		for (int j = 0; j < n_level; j++)
		{
			total_seq += (total_time_seq[j + i * n_level] + total_time_seq_inverse[j + i * n_level]);
			total_par += (total_time_par[j + i * n_level] + total_time_par_inverse[j + i * n_level]);
		}
	}

	speed_up = (1.0 - (prl / seq)) * 100;
	p = (total_seq + threshold_time) / seq;

	print_stats(images, n_level, n_level_stop, seq, prl, total_time_seq, total_time_seq_inverse, total_time_par, total_time_par_inverse, speed_up, p, total_seq, total_par);

	free(time);
	free(total_time_seq);
	free(total_time_seq_inverse);
	free(total_time_par_first);
	free(total_time_par_second);
	free(total_time_par);
	free(total_time_par_inverse);
	
	return 0;
}
