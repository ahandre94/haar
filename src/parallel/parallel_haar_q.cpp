#include "parallel_haar_q.h"

/**
	Parallel version of 2D-Haar Wavelet transform on square images

	@param int m -> number of rows, must be integer power of 2
	@param int n -> number of columns, must be integer power of 2
	@param double u[] -> vector in which the image is read row by row
						 img = [0 1 2 3			u = [0 1 2 3 4 5 6 7]
								4 5 6 7]
	@param int n_level -> number of level to perform the haar transform
	@param int n_fix -> must be set as the number of columns and never changed
*/
void p_haar_2d(int m, int n, double u[], int n_level, int n_fix, double *time, int iterator)
{
	int i, j;
	int k;
	int flag_giro = (int)(n_fix / n);
	double s = sqrt(2.0);
	double *v;
	double start, end;

	v = (double *)malloc(m * n * sizeof(double));

	start = omp_get_wtime();

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = u[j + i * n * flag_giro];
		}
	}

	k = n / 2;

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			v[j + i * n] = (u[2 * j + i * n * flag_giro] + u[2 * j + 1 + i * n * flag_giro]) / s;
			v[k + j + i * n] = (u[2 * j + i * n * flag_giro] - u[2 * j + 1 + i * n * flag_giro]) / s;
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

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m / 2; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = (u[j + 2 * i * n * flag_giro] + u[j + (2 * i + 1) * n * flag_giro]) / s;
			v[j + (k + i) * n] = (u[j + 2 * i * n * flag_giro] - u[j + (2 * i + 1) * n * flag_giro]) / s;
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

	if (n_level > 1)
	{
		n_level--;
		p_haar_2d(m / 2, n / 2, u, n_level, n_fix, time, iterator);
	}
}


/**
	Parallel version of Inverse 2D-Haar Wavelet transform on square images
*/
void p_haar_2d_inverse(int m, int n, double u[], int n_level)
{
	int i, j;
	int k;
	int flag_giro = (int)pow(2, n_level - 1);
	double s = sqrt(2.0);
	double *v;

	v = (double *)malloc(m * n * sizeof(double));

	int m_fix = m;
	int n_fix = n;
	m = m / flag_giro;
	n = n / flag_giro;

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + i * n] = u[j + i * n * flag_giro];
		}
	}

	k = n / 2;

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m / 2; i++)
	{
		for (j = 0; j < n; j++)
		{
			v[j + 2 * i * n] = (u[j + i * n * flag_giro] + u[j + (k + i) * m * flag_giro]) / s;
			v[j + 2 * i * n + n] = (u[j + i * n * flag_giro] - u[j + (k + i) * m * flag_giro]) / s;
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

	#pragma omp parallel for private(i, j) shared(u, v) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < k; j++)
		{
			v[2 * j + i * n] = (u[j + i * n * flag_giro] + u[k + j + i * n * flag_giro]) / s;
			v[2 * j + 1 + i * n] = (u[j + i * n * flag_giro] - u[k + j + i * n * flag_giro]) / s;
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

	if (n_level > 1)
	{
		n_level--;
		p_haar_2d_inverse(m_fix, n_fix, u, n_level);
	}
}
