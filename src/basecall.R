// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

// FFTW3
#include "fftw3.h"

//for wavelet
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <algorithm>
#include "wavelet2d.h"
#include <math.h>
#include <iostream>

//for sleep
#include <unistd.h>

#include <Rcpp.h>

using namespace cv;
using namespace Rcpp;

RcppExport SEXP baseCall(SEXP input) {
  BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>

  Rcpp::CharacterVector f(input);
  std::string ff(f[0]);
  Rcpp::Rcout << "Loading image:" << ff << std::endl;
  t0 = clock();
  tobefiltered = imread(ff, -1); //CV_LOAD_IMAGE_GRAYSCALE
  //for debugging  imgSWT('/Users/danielfurth/Documents/mouse_atlas/lighsheet_figure/analysis/Resliceofanalysis.tif', cell.bodies=4, fluorescence.threshold=0, areaSize=0)->mupp2
  t1 = clock();
  double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." << std::endl;


  //check that the image is symmetric. If not then subsample the image by bilinear interpolation.
  bool sumbsample;
  if(rowsWzero!=colsWzero){
    if(rowsWzero>colsWzero){
      Rcpp::Rcout << "Image size "<< colsWzero << "x" << rowsWzero << " is not symmetric, subsampling will be performed." << std::endl;
      Size downsampledsize(static_cast<int>(colsWzero), static_cast<int>(colsWzero));
      resize(tobefiltered, tobefiltered, downsampledsize, INTER_LINEAR);
      sumbsample = true;
    }else{
      Rcpp::Rcout << "Image size "<< colsWzero << "x" << rowsWzero << " is not symmetric, subsampling will be performed." << std::endl;
      Size downsampledsize(static_cast<int>(rowsWzero), static_cast<int>(rowsWzero));
      resize(tobefiltered, tobefiltered, downsampledsize, INTER_LINEAR);
      sumbsample = true;
    }
  }else{
    sumbsample = false;
  }
