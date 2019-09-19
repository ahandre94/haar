#include "load_img.h"

namespace fs = std::experimental::filesystem;

vector<Mat> load_img(string path)
{
	vector<Mat> images;

	for (const auto & entry : fs::directory_iterator(path))
	{
		string imageName(entry.path().string());
		Mat img = imread(imageName, IMREAD_GRAYSCALE);

		if (img.empty())
			break;

		images.push_back(img);
	}
	
	return images;
}
