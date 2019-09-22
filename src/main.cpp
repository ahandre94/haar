#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "sequential/haar.h"
#include "sequential/threshold.h"
#include "utility/dimension.h"
#include "utility/load_img.h"
#include "parallel/parallel_haar.h"
#include "parallel/parallel_threshold.h"

#define SHOW_IMAGE(x) {namedWindow(#x, WINDOW_AUTOSIZE); imshow(#x, x); waitKey(0);}

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
			   double *threshold_time)
{
	Mat image, image_double, haar_img, inverse_haar_img;
	Size s;

	int m, n;
	bool change_size;
	double media;
	double start = 0.0;

	for (unsigned i = 0; i < images.size(); i++)
	{
		image = images[i];

		//SHOW_IMAGE(image);

		m = image.rows;
		n = image.cols;
		s = Size(n, m);

		tie(change_size, m, n) = change_dimension(m, n);
		if (change_size)
			resize(image, image, Size(n, m));

		image.convertTo(image_double, CV_64F);
		haar(m, n, image_double, n_level, time, iterator);
		
		for (int j = 0; j < n_level; j++)
		{
			total_time[j + i * n_level] = time[j];
		}

		image_double.copyTo(haar_img);
		if (n_level == n_level_stop)
			visualizza(m, n, (double *)(haar_img.data), n_level);
		else if (n_level != n_level_stop && type == "seq")
			visualizza(m, n, (double *)(haar_img.data), n_level);
		haar_img.convertTo(haar_img, CV_8U);
		if (change_size)
			resize(haar_img, haar_img, s);
		//SHOW_IMAGE(haar_img);
		
		std::ostringstream name;
		if (type == "par")
			name << "output_img/" << type << "_img" << i << ".png";
		else if (type == "seq")
			name << "output_img_seq/" << type << "_img" << i << ".png";
		cv::imwrite(name.str(), haar_img);

		if (type == "seq")
			start = omp_get_wtime();
		media = mean(m, n, (double *)(image_double.data), n_level);
		threshold(m, n, (double *)(image_double.data), n_level, media);
		if (type == "seq")
			*threshold_time += (omp_get_wtime() - start);

		iterator = 0;
		inverse(m, n, image_double, n_level, time, iterator);

		for (int j = 0; j < n_level; j++)
		{
			total_time_inverse[j + i * n_level] = time[j];
		}

		image_double.copyTo(inverse_haar_img);
		inverse_haar_img.convertTo(inverse_haar_img, CV_8U);
		if (change_size)
			resize(inverse_haar_img, inverse_haar_img, s);
		//SHOW_IMAGE(inverse_haar_img);
	}
}

// too much overhead, better not to use this function; set n_level and n_level_stop equal
void execution_parallel_mix_sequential(int n_level_to_do, int n_level_stop, double *time, double *total_time)
{
	vector<Mat> images;

	Mat image;
	int m, n;
	bool change_size;
	
	if (n_level_to_do > 0)
	{
		images = load_img("output_img");
		for (unsigned i = 0; i < images.size(); i++)
		{
			image = images[i];
			m = image.rows;
			n = image.cols;

			tie(change_size, m, n) = change_dimension(m, n);
			if (change_size)
			{			
				resize(image, image, Size(n, m));
				images[i] = image;
			}
		}
		
		#pragma omp parallel for private(m, n, image, change_size) shared(n_level_stop, n_level_to_do, time, total_time, images)
		for (unsigned i = 0; i < images.size(); i++)
		{
			image = images[i];

			m = (int)(image.rows / pow(2, n_level_stop - 1));
			n = (int)(image.cols / pow(2, n_level_stop - 1));
			haar(m, n, image, n_level_to_do, time, n_level_stop);

			for (int j = 0; j < n_level_to_do; j++)
			{
				total_time[j + i * n_level_to_do] = time[j + n_level_stop];
			}
		}
	}
}

void stampa(vector<Mat> images, int n_level, string type, string func, double *total_time, FILE *f)
{
	fprintf(f, "%s\n", type.c_str());
	fprintf(f, "%s\n", func.c_str());
	fprintf(f, "n_level\t|");
	if (func == "HAAR FUNCTION")
	{
		for (int i = 0; i < n_level; i++)
		{
			fprintf(f, "\t%d\t\t|\t", i + 1);
		}
	}
	else
	{
		for (int i = n_level; i > 0; i--)
		{
			fprintf(f, "\t%d\t\t|\t", i);
		}
	}
	fprintf(f, "\n");
	for (unsigned i = 0; i < images.size(); i++)
	{	
		fprintf(f, "%d\t|\t", i);
		for (int j = 0; j < n_level; j++)
		{
			fprintf(f, "%lf\t|\t", total_time[j + i * n_level]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
}

void stampa_time_reduction(vector<Mat> images, int n_level, string type, double *total_time_s, double *total_time_p, FILE *f)
{
	double tmp;

	fprintf(f, "TIME REDUCTION PER LEVEL %s\n", type.c_str());
	fprintf(f, "n_level\t\t|\t");
	if (type == "HAAR")
	{
		for (int i = 0; i < n_level; i++)
		{
			fprintf(f, "\t%d\t\t|\t", i + 1);
		}
	}
	else
	{
		for (int i = n_level; i > 0; i--)
		{
			fprintf(f, "\t%d\t\t|\t", i);
		}
	}
	fprintf(f, "\n");
	for (unsigned i = 0; i < images.size(); i++)
	{
		fprintf(f, "%dx%d\t|\t", images[i].cols, images[i].rows);
		for (int j = 0; j < n_level; j++)
		{
			tmp = total_time_s[j + i * n_level] / total_time_p[j + i * n_level];
			(tmp < 1.0) ? fprintf(f, "[%lf]\t|\t", tmp) : fprintf(f, "%lf\t|\t", tmp);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
}

int main(int argc, char **argv)
{
	vector<Mat> images = load_img("4k");
	
	int n_level = 5;
	int n_level_stop = 5;
	int stop = n_level - n_level_stop;

	FILE *f;
	double start, end, seq, prl;
	double speed_up, p;
	int iterator;
	double *time = (double *)malloc(n_level * sizeof(double));
	double *total_time_seq = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_seq_inverse = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par_first = (double *)malloc(images.size() * n_level_stop * sizeof(double));
	double *total_time_par_second = (double *)malloc(images.size() * stop * sizeof(double));
	double *total_time_par_inverse = (double *)malloc(images.size() * n_level_stop * sizeof(double));
	double threshold_time = 0.0;
	double total_seq = 0.0, total_par = 0.0;

	f = fopen("stats.txt", "a+");
	fprintf(f, "num_core: %d\nn_level: %d\nn_level_stop: %d\n\n", omp_get_num_procs(), n_level, n_level_stop);
	
	// sequential execution
	iterator = 0;
	start = omp_get_wtime();
	exeggutor(images, n_level, n_level_stop, "seq", time, iterator, haar, visualizza_haar, haar_inverse, threshold, mean, total_time_seq, total_time_seq_inverse, &threshold_time);
	end = omp_get_wtime();
	stampa(images, n_level, "SEQUENTIAL WORK", "HAAR FUNCTION", total_time_seq, f);
	stampa(images, n_level, "SEQUENTIAL WORK", "INVERSE HAAR FUNCTION", total_time_seq_inverse, f);
	seq = end - start;
	fprintf(f, "Sequential work took %lf seconds\n\n", seq);
	
	// parallel execution
	iterator = 0;
	start = omp_get_wtime();
	exeggutor(images, n_level_stop, n_level, "par", time, iterator, p_haar, p_visualizza_haar, p_haar_inverse, p_threshold, p_mean, total_time_par_first, total_time_par_inverse, &threshold_time);
	execution_parallel_mix_sequential(stop, n_level_stop, time, total_time_par_second);
	end = omp_get_wtime();

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
	}

	stampa(images, n_level, "PARALLEL WORK", "HAAR FUNCTION", total_time_par, f);
	stampa(images, n_level_stop, "PARALLEL WORK", "INVERSE HAAR FUNCTION", total_time_par_inverse, f);
	prl = end - start;
	fprintf(f, "Parallel work took %lf seconds\n\n", prl);

	speed_up = (1.0 - (prl / seq)) * 100;
	fprintf(f, "Time reduction: %.2lf%%\n", speed_up);

	stampa_time_reduction(images, n_level, "HAAR", total_time_seq, total_time_par, f);
	if (n_level == n_level_stop)
		stampa_time_reduction(images, n_level, "INVERSE HAAR", total_time_seq_inverse, total_time_par_inverse, f);

	for (unsigned i = 0; i < images.size(); i++)
	{
		for (int j = 0; j < n_level; j++)
		{
			total_seq += (total_time_seq[j + i * n_level] + total_time_seq_inverse[j + i * n_level]);
			total_par += (total_time_par[j + i * n_level] + total_time_par_inverse[j + i * n_level]);
		}
	}

	p = (total_seq + threshold_time) / seq;
	fprintf(f, "theoretical speed up: %lf\n", 1 / (1 - p + (p / omp_get_num_procs())));
	fprintf(f, "real speed up: %lf\n", (seq / prl));
	if (n_level == n_level_stop)
		fprintf(f, "haar medium speed up: %lf\n", (total_seq / total_par));
	fprintf(f, "\n\n\n");	

	fclose(f);
	free(time);
	free(total_time_seq);
	free(total_time_seq_inverse);
	free(total_time_par_first);
	free(total_time_par_second);
	free(total_time_par);
	free(total_time_par_inverse);
	
	return 0;
}