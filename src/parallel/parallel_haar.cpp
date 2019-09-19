#include "parallel_haar.h"

void p_haar(int m, int n, Mat image, int n_level, double *time, int iterator)
{
	if (m == n)
	{
		p_haar_2d(m, n, (double *)(image.data), n_level, n, time, iterator);
	}
	if (m < n)
	{
		p_haar_2d_c(m, n, (double *)(image.data), n_level, n, time, iterator);
	}
	if (m > n)
	{
		p_haar_2d_r(m, n, (double *)(image.data), n_level, n, time, iterator);
	}
}

void p_haar_inverse(int m, int n, Mat image, int n_level)
{
	if (m == n)
	{
		p_haar_2d_inverse(m, n, (double *)(image.data), n_level);
	}
	if (m < n)
	{
		p_haar_2d_inverse_c(m, n, (double *)(image.data), n_level);
	}
	if (m > n)
	{
		p_haar_2d_inverse_r(m, n, (double *)(image.data), n_level);
	}
}

void p_visualizza_haar(int m, int n, double u[], int n_level)
{
	int i, j;
	int level = (int)pow(2, n_level);
	int l = (int)(m / level);
	int k = (int)(n / level);

	#pragma omp parallel for private(i, j) shared(u) collapse(2)
	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			u[j + i * n] += (!(i < l && j < k)) ? 127 : -127 * (level - 1);
		}
	}
}
