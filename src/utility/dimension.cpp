#include "dimension.h"

/*
	return true if dim is not a integer power of 2
*/
bool check_dimension(int dim)
{
	return !((dim != 0) && ((dim & (dim - 1)) == 0));
}

tuple<bool, int> change_dim(int dim)
{
	bool change = false;
	if (check_dimension(dim))
	{
		int k = 1;
		while (k <= dim)
		{
			k *= 2;
		}
		dim = k;
		change = true;
	}
	return { change, dim };
}

tuple<bool, int, int> change_dimension(int m, int n)
{
	bool change_m, change_n;

	tie(change_m, m) = change_dim(m);
	tie(change_n, n) = change_dim(n);

	return { change_m || change_n, m, n };
}
