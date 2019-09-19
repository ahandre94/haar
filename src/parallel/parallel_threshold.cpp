#include "parallel_threshold.h"

void p_threshold(int m, int n, double u[], int n_level, double threshold)
{
	int i, j;
	int level = (int)pow(2, n_level);
	int l = (int)(m / level);
	int k = (int)(n / level);

	#pragma omp parallel for private(i, j) shared(u, threshold) collapse(2)
	for (i = 0; i < l; i++)
	{
		for (j = k; j < n; j++)
		{
			if (abs(u[j + i * n]) <= threshold)
			{
				u[j + i * n] = 0.0;
			}
		}
	}
	#pragma omp parallel for private(i, j) shared(u, threshold) collapse(2)
	for (i = l; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (abs(u[j + i * n]) <= threshold)
			{
				u[j + i * n] = 0.0;
			}
		}
	}
}

double p_mean(int m, int n, double u[], int n_level)
{
	int i, j;
	int level = (int)pow(2, n_level);
	int l = (int)(m / level);
	int k = (int)(n / level);
	double sum = 0;
	double mean = 0;
	int n_element = (m * n) - (l * k);

	#pragma omp parallel for private(i, j) shared(u) reduction(+:sum) collapse(2)
	for (i = 0; i < l; i++)
	{
		for (j = k; j < n; j++)
		{
			sum += abs(u[j + i * n]);
		}
	}
	#pragma omp parallel for private(i, j) shared(u) reduction(+:sum) collapse(2)
	for (i = l; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			sum += abs(u[j + i * n]);
		}
	}
	mean = sum / n_element;

	return mean;
}
