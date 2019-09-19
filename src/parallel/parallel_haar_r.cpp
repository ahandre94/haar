#include "parallel_haar_r.h"

/**
	Parallel version of 2D-Haar Wavelet transform on rectangular images where 
	the number of columns is less than the number of rows

	@param int m -> number of rows, must be integer power of 2
	@param int n -> number of columns, must be integer power of 2
	@param double u[] -> vector in which the image is read row by row
	@param int n_level -> number of level to perform the haar transform
	@param int n_fix -> must be set as the number of columns and never changed
*/
void p_haar_2d_r(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator)
{
	int i, j;
	int k;
	int flag_giro = (int)(n_fix / n);
	double s = sqrt(2.0);
	double *v, *w;
	double start, end;

	v = (double *)malloc(m * n * sizeof(double));
	w = (double *)malloc(m * n * sizeof(double));

	start = omp_get_wtime();

	#pragma omp parallel for private(i, j) shared(u, v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = u[j + i * n * flag_giro];
			w[j + i * n] = u[j + i * n * flag_giro];
		}
	}

	k = n / 2;

	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			v[j + i * n] = (w[2 * j + i * n] + w[2 * j + 1 + i * n]) / s;
			v[j + i * n + k] = (w[2 * j + i * n] - w[2 * j + 1 + i * n]) / s;
		}
	}
	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			w[j + i * n] = v[j + i * n];
		}
	}

	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m / 2; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = (w[j + 2 * i * n] + w[j + 2 * i * n + n]) / s;
			v[j + i * n + m * k] = (w[j + 2 * i * n] - w[j + 2 * i * n + n]) / s;
		}
	}
	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			u[j + i * n * flag_giro] = v[j + i * n];
		}
	}

	end = omp_get_wtime();
	
	time[iterator] = end - start;
	iterator++;

	free(v);
	free(w);

	if (n_level > 1)
	{
		n_level--;
		p_haar_2d_r(m / 2, n / 2, u, n_level, n_fix, time, iterator);
	}
}

/**
	Parallel version of Inverse 2D-Haar Wavelet transform on rectangular
	images where the number of columns is less than the number of rows
*/
void p_haar_2d_inverse_r(int m, int n, double u[], int n_level)
{
	int i, j;
	int k;
	int flag_giro = (int)pow(2, n_level - 1);
	double s = sqrt(2.0);
	double *v, *w;

	v = (double *)malloc(m * n * sizeof(double));
	w = (double *)malloc(m * n * sizeof(double));

	int m_fix = m;
	int n_fix = n;
	m = m / flag_giro;
	n = n / flag_giro;

	#pragma omp parallel for private(i, j) shared(u, v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = u[j + i * n * flag_giro];
			w[j + i * n] = u[j + i * n * flag_giro];
		}
	}

	k = n / 2;

	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m / 2; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + 2 * i * n] = (w[j + i * n] + w[j + i * n + m * k]) / s;
			v[j + 2 * i * n + n] = (w[j + i * n] - w[j + i * n + m * k]) / s;
		}
	}
	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			w[j + i * n] = v[j + i * n];
		}
	}

	#pragma omp parallel for private(i, j) shared(v, w) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			v[2 * j + i * n] = (w[j + i * n] + w[j + i * n + k]) / s;
			v[2 * j + 1 + i * n] = (w[j + i * n] - w[j + i * n + k]) / s;
		}
	}
	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			u[j + i * n * flag_giro] = v[j + i * n];
		}
	}

	free(v);
	free(w);

	if (n_level > 1)
	{
		n_level--;
		p_haar_2d_inverse_r(m_fix, n_fix, u, n_level);
	}
}
