# Haar
Haar Wavelet 2D

A wavelet is a mathematical function useful in digital signal processing and image compression. Wavelets have been used to compress images to a greater extent than is generally possible with other methods. 

The Haar Wavelet is the simplest wavelet and it is efficient to perform both lossless and lossy image compression. It relies on averaging and differencing values in an image matrix to produce a matrix which is sparse or nearly sparse.

The algorithm here presented can perform the Haar transform and the inverse Haar transform on any type of image (the requirement for the Haar is that both the width and the height of the image are integer power of 2, but here it is not strictly necessary thanks to a preventive scaling operation). Two different versions have been implemented: one, which works in sequential, and another one, which works in parallel. 

A further implementation has been made; the idea behind it is the following: as the number of level increases, the size of the image where the Haar function is performed decreases. It is possible to assume that the speedup decreases too. So, the idea is to perform the parallel execution until a certain number of levels (i.e., until the size of the image is big enough) and then to parallelize the images and to perform the sequential version of the Haar Wavelet on them. However, during the execution of this new function there was an excessive amount of overhead, so it is better to avoid using it.

Also a thresholding function has been implemented; this requires a threshold value and so a function which returns the mean of the values on the horizontal, vertical, and diagonal detail coefficients has been implemented too.

Here an example of execution:

Original image
![alt original_lion](https://github.com/ahandre94/haar/blob/master/example/lion.jpg)

4 levels Haar transform
![alt transformed_lion](https://github.com/ahandre94/haar/blob/master/example/4_level_haar.png)

Restored image after thresholding
![alt restored_lion](https://github.com/ahandre94/haar/blob/master/example/restored_lion.png)

# Requirements
OpenCV

OpenMP
