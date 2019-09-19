#ifndef DIMENSION_H
#define DIMENSION_H

#include <tuple>

using namespace std;

bool check_dimension(int dim);
tuple<bool, int> change_dim(int dim);
tuple<bool, int, int> change_dimension(int m, int n);

#endif
