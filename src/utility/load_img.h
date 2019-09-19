#ifndef LOAD_IMG_H
#define LOAD_IMG_H

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <iostream>
#include <string>
#include <experimental/filesystem>

using namespace cv;
using namespace std;

vector<Mat> load_img(string path);

#endif
