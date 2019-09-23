#include "execution.h"

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

// too much overhead, better not to use this function; set n_level and n_level_stop equal in main function
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

