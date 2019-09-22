#include "print_stats.h"

void indentation(string type, string func, int n_level, FILE *f)
{
	if (type == func)
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
}

void stampa(vector<Mat> images, int n_level, string type, string func, double *total_time, FILE *f)
{
	fprintf(f, "%s\n", type.c_str());
	fprintf(f, "%s\n", func.c_str());
	fprintf(f, "n_level\t|");
	indentation(func, "HAAR FUNCTION", n_level, f);
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
	indentation(type, "HAAR", n_level, f);
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
				   double total_par)
{
	FILE *f;

	f = fopen("stats.txt", "a+");
	fprintf(f, "num_core: %d\nn_level: %d\nn_level_stop: %d\n\n", omp_get_num_procs(), n_level, n_level_stop);
	stampa(images, n_level, "SEQUENTIAL WORK", "HAAR FUNCTION", total_time_seq, f);
	stampa(images, n_level, "SEQUENTIAL WORK", "INVERSE HAAR FUNCTION", total_time_seq_inverse, f);
	fprintf(f, "Sequential work took %lf seconds\n\n", seq);
	stampa(images, n_level, "PARALLEL WORK", "HAAR FUNCTION", total_time_par, f);
	stampa(images, n_level_stop, "PARALLEL WORK", "INVERSE HAAR FUNCTION", total_time_par_inverse, f);
	fprintf(f, "Parallel work took %lf seconds\n\n", prl);
	fprintf(f, "Time reduction: %.2lf%%\n", speed_up);
	stampa_time_reduction(images, n_level, "HAAR", total_time_seq, total_time_par, f);
	if (n_level == n_level_stop)
		stampa_time_reduction(images, n_level, "INVERSE HAAR", total_time_seq_inverse, total_time_par_inverse, f);
	fprintf(f, "theoretical speed up: %lf\n", 1 / (1 - p + (p / omp_get_num_procs())));
	fprintf(f, "real speed up: %lf\n", (seq / prl));
	if (n_level == n_level_stop)
		fprintf(f, "haar medium speed up: %lf\n", (total_seq / total_par));
	fprintf(f, "\n\n\n");
	fclose(f);
}
