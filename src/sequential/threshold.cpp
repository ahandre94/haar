#include "threshold.h"

void threshold(int m, int n, double u[], int n_level, double threshold)
{
	int i, j;
	int level = (int)pow(2, n_level);
	int l = (int)(m / level);
	int k = (int)(n / level);

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (!(i < l && j < k))
			{
				if (abs(u[j + i * n]) <= threshold)
				{
					u[j + i * n] = 0;
				}
			}
		}
	}
}

double mean(int m, int n, double u[], int n_level)
{
	int i, j;
	int level = (int)pow(2, n_level);
	int l = (int)(m / level);
	int k = (int)(n / level);
	double sum = 0;
	double mean = 0;
	int n_element = (m * n) - (l * k);

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (!(i < l && j < k))
			{
				sum += abs(u[j + i * n]);
			}
		}
	}
	mean = sum / n_element;

	return mean;
}
