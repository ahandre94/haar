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
			   string type,
			   double *time,
			   int iterator,
			   void (*haar)(int, int, Mat, int, double*, int), 
			   void (*visualizza)(int, int, double*, int), 
			   void (*inverse)(int, int, Mat, int), 
			   void (*threshold)(int, int, double*, int, double), 
			   double (*mean)(int, int, double*, int),
			   double *total_time)
{
	Mat image, image_double, haar_img, inverse_haar_img;
	Size s;

	int m, n;
	bool change_size;
	//double media;

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
		if (type == "seq")
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
/*
		media = mean(m, n, (double *)(image_double.data), n_level);
		threshold(m, n, (double *)(image_double.data), n_level, media);

		inverse(m, n, image_double, n_level);

		image_double.copyTo(inverse_haar_img);
		inverse_haar_img.convertTo(inverse_haar_img, CV_8U);
		if (change_size)
			resize(inverse_haar_img, inverse_haar_img, s);
		//SHOW_IMAGE(inverse_haar_img);*/
	}
}

void execution_parallel_mix_sequential(int n_level, int n_level_stop, double *time, double *total_time)
{
	vector<Mat> images;
	int n_level_to_do = n_level - n_level_stop;

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
			m = image.rows;
			n = image.cols;

			m = (int)(m / pow(2, n_level_stop - 1));
			n = (int)(n / pow(2, n_level_stop - 1));
			haar(m, n, image, n_level_to_do, time, n_level_stop);

			for (int j = 0; j < n_level_to_do; j++)
			{
				total_time[j + i * n_level_to_do] = time[j + n_level_stop];
			}
		}
	}
}

void stampa(vector<Mat> images, int n_level, string type, double *total_time, FILE *f)
{
	fprintf(f, "%s\n", type.c_str());
	for (unsigned i = 0; i < images.size(); i++)
	{	
		fprintf(f, "%d\t|\t", i);
		for (int j = 0; j < n_level; j++)
		{
			fprintf(f, "%lf\t|\t", total_time[j + i * n_level]);
		}
		fprintf(f, "\n");
	}
}

int main(int argc, char **argv)
{
	vector<Mat> images = load_img("4k");

	double start, end, seq, prl;
	double speed_up;
	int n_level = 5;
	int n_level_stop = 5;
	int stop = n_level - n_level_stop;
	int iterator;
	double *time = (double *)malloc(n_level * sizeof(double));
	double *total_time_seq = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par = (double *)malloc(images.size() * n_level * sizeof(double));
	double *total_time_par_first = (double *)malloc(images.size() * n_level_stop * sizeof(double));
	double *total_time_par_second = (double *)malloc(images.size() * stop * sizeof(double));
	double tmp, total_seq = 0.0, total_par = 0.0;
	double p;
	FILE *f;

	f = fopen("stats.txt", "a+");
	fprintf(f, "num_core: %d\nn_level: %d\nn_level_stop: %d\n", omp_get_num_procs(), n_level, n_level_stop);
	
	// sequential execution
	iterator = 0;
	start = omp_get_wtime();
	exeggutor(images, n_level, "seq", time, iterator, haar, visualizza_haar, haar_inverse, threshold, mean, total_time_seq);
	end = omp_get_wtime();
	stampa(images, n_level, "SEQUENTIAL WORK", total_time_seq, f);
	seq = end - start;
	fprintf(f, "Sequential work took %lf seconds\n", seq);
	
	// parallel execution
	iterator = 0;
	start = omp_get_wtime();
	exeggutor(images, n_level_stop, "par", time, iterator, p_haar, p_visualizza_haar, p_haar_inverse, p_threshold, p_mean, total_time_par_first);
	execution_parallel_mix_sequential(n_level, n_level_stop, time, total_time_par_second);
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

	stampa(images, n_level, "PARALLEL WORK", total_time_par, f);
	prl = end - start;
	fprintf(f, "Parallel work took %lf seconds\n", prl);

	speed_up = (1.0 - (prl / seq)) * 100;
	fprintf(f, "Time reduction: %.2lf%%\n", speed_up);

	fprintf(f, "TIME REDUCTION PER LEVEL\n");
	for (unsigned i = 0; i < images.size(); i++)
	{
		fprintf(f, "%dx%d\t|\t", images[i].cols, images[i].rows);
		for (int j = 0; j < n_level; j++)
		{
			tmp = total_time_seq[j + i * n_level] / total_time_par[j + i * n_level];
			(tmp < 1.5) ? fprintf(f, "[%lf]\t|\t", tmp) : fprintf(f, "%lf\t|\t", tmp);
			total_seq += total_time_seq[j + i * n_level];
			total_par += total_time_par[j + i * n_level];
		}
		fprintf(f, "\n");
	}

	p = total_seq / seq;
	fprintf(f, "speed up teorico: %lf\n", 1 / (1 - p + (p / omp_get_num_procs())));
	fprintf(f, "speed up reale: %lf\n", (seq / prl));
	fprintf(f, "speed up media: %lf\n", (total_seq / total_par));
	fprintf(f, "\n\n\n");	

	fclose(f);
	free(time);
	free(total_time_seq);
	free(total_time_par_first);
	free(total_time_par_second);
	free(total_time_par);
	
	return 0;
}
