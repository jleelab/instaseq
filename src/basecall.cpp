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

//to handle execution time logging
clock_t t0, t1;


/* function to compute the maximum from an image array */
void* maxval(vector<vector<double> > &arr, double &max){
  max = 0;
  for (unsigned int i =0; i < arr.size(); i++) {
    for (unsigned int j =0; j < arr[0].size(); j++) {
      if (max <= arr[i][j]){
        max = arr[i][j];
      }
    }
  }
  return 0;
}

//find closest
int closest(std::vector<int> const& vec, int value) {
  //auto const it = std::upper_bound(vec.begin(), vec.end(), value);

  auto const upper = std::lower_bound(vec.begin(), vec.end(), value);
  auto const lower = upper-1;

  if (upper == vec.end()) { return -1; }

  return *lower;
}


RcppExport SEXP baseCall(SEXP input, SEXP scales, SEXP sigmaR, SEXP outputfilename) {
  BEGIN_RCPP
  Rcpp::RNGScope __rngScope; //this and BEGIN_RCPP and END_RCPP is needed for wrappers such as Rcpp::as<int>

  Mat tobefiltered; // no local
  Mat cellImg;

  int Scale = Rcpp::as<int>(scales);
  int Sigma = Rcpp::as<int>(sigmaR);

  Rcpp::CharacterVector stringi(outputfilename);
  std::string outFilename(stringi[0]);

  Rcpp::CharacterVector f(input);
  std::string ff(f[0]);
  Rcpp::Rcout << "Loading image:" << ff << std::endl;
  t0 = clock();
  std::vector<cv::Mat> multipageImage;
  cv::imreadmulti(ff, multipageImage, cv::IMREAD_ANYDEPTH);
  t1 = clock();
  double time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "LOADED." << " loading took " <<  time_elapsed << " seconds." <<  " Tiles: " << multipageImage.size() << std::endl;

  for (int q=0; q < multipageImage.size(); q++) {
    tobefiltered = multipageImage[q];

    double maxValinput, minValinput;
    minMaxLoc(tobefiltered, &minValinput, &maxValinput);

    //get width and height of original image.
    int rowsWzero = tobefiltered.rows;
    int colsWzero = tobefiltered.cols;

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

  //get width and height of image.
  int rowsW = tobefiltered.rows;
  int colsW = tobefiltered.cols;

  //set the different wavelet resolutions for cases where the image is not sampled by a factor of 2
  int myints[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536};
  std::vector<int> vec (myints, myints + sizeof(myints) / sizeof(int) );

  //resize image to a J size if it is not contained with 2^J sampling.
  if(std::find(vec.begin(), vec.end(), rowsW) != vec.end()) {
    /* resolutions contains rowsW */
  } else {
    /* resolutions does not contain rowsW */
    Rcpp::Rcout << "Image size "<< rowsW << "x" << rowsW << " needs to be adjused to 2^J sampling." << std::endl;
    //save the old row for future reference.
    int oldRowW = rowsW;
    //function closest defined on line 60.
    int newwidth = closest(vec, rowsW);
    Size oldsize = tobefiltered.size();
    Size newsize(static_cast<int>(newwidth), static_cast<int>(newwidth));
    resize(tobefiltered, tobefiltered, newsize, INTER_LINEAR);
    //get width and height of image.
    rowsW = tobefiltered.rows;
    colsW = tobefiltered.cols;
    Rcpp::Rcout << "New size: "<< rowsW << "x" << rowsW << "." << std::endl;
    sumbsample = true;
  }

  //make a 1D vector of the image.
  vector<vector<double> > vec1(rowsW, vector<double>(colsW));

  int k =1;
  for (int i=0; i < rowsW; i++) {
    for (int j =0; j < colsW; j++){
      unsigned char temp;
      temp = ((uchar*) tobefiltered.data + i * tobefiltered.step)[j  * tobefiltered.elemSize() + k ];
      vec1[i][j] = (double) temp;
    }
  }


  //create output vector
  vector<double> output;

  string nm = "db2"; // default "db2"

  //define wavelet scales
  //int J = Rcpp::as<int>(scales);// Rcpp::as<int >(scales); //default is J=6

  //this is the computationally heavy part. Improve swt_2d with OpenCL
  Rcpp::Rcout << "Beginning wavelet transform" << std::endl;
  t0 = clock();
  swt_2d(vec1,6,nm,output);

  t1 = clock();
  time_elapsed = (double)((t1 - t0) / CLOCKS_PER_SEC);
  Rcpp::Rcout << "DONE." << " wavelet took " <<  time_elapsed << " seconds." << std::endl;


  // get the output dimensions
  int rowW,colW;
  dwt_output_dim(vec1, rowW, colW );

  //create OpenCv matrix to store the approximation coefficients in.
  Size imgSize = tobefiltered.size(); // size of output image

  Mat cvImg(imgSize, tobefiltered.depth());

  //define detail coefficients
  vector<vector<double> >  detail(1*rowW,vector<double>(6 * colW));
  //extract detail coefficients from output martix
  Rcpp::Rcout << "Extracting detail coefficients" << std::endl;
  for (int k=0; k < 6; k++) {
    for (int i=0; i < rowW; i++) {
      for(int j=0+ k*colW; j < (k+1)*colW; j++) {
        double temp = (output[(3*k+1)*rowW*colW+ i * colW +j - k*colW]+output[(3*k+2)*rowW*colW+ i * colW +j - k*colW]+output[(3*k+3)*rowW*colW+ i * colW +j - k*colW]);
        detail[i][j]= temp;
      }
    }
  }

  //IplImage *dvImg; // image used for output
  //CvSize imgSz; // size of output image
  cv::Size imgSz; // size of output image

  imgSz.width = 6*colW;
  imgSz.height = 1*rowW;

  //dvImg = cvCreateImage( imgSz, 16, 1 );
  //cv::Mat dvImg( imgSz, CV_16U);

  double max;
  maxval(detail, max);
  /*Mat w(detail[0]),dvImg;
  Rcpp::Rcout << "Detail size:" << detail.size() << std::endl;

  for (int i = 1; i < detail.size(); i++)
    w.push_back(Mat(detail[i]));
  w = w.reshape(detail[0].size());
  w.convertTo(dvImg, CV_16U);
  Rcpp::Rcout << "dvImg size:" << dvImg.size() << std::endl;
  Rcpp::Rcout << "Detail:" << detail[512][512] << std::endl;
   */

  /*
  for (int i = 0; i < imgSz.height; i++ ) {
    for (int j = 0; j < imgSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }
      dvImg.at<ushort>(j,i) = (short) (detail[i][j]);

      // ((ushort*)(dvImg->imageData + dvImg->widthStep*i))[j] = (short) (detail[i][j]);
      // use "max) * 65536.0" instead of max) * maxvue if normalie on single tile.  (short) ( (detail[i][j]/ max) * imageMax);
    }
  }*/

  for (int i = 0; i < imgSz.height; i++ ) {
    for (int j = 0; j < imgSz.width; j++ ){
      if ( detail[i][j] <= 0.0){
        detail[i][j] = 0.0;
      }
    }
  }

  Mat w(detail[0]),dvImg;
  for (int i = 1; i < detail.size(); i++)
    w.push_back(Mat(detail[i]));
  w = w.reshape(0,detail.size());
  w.convertTo(dvImg, CV_16U);



  // }
  //
  //
  // for (int i = 0; i < imgSz.height; i++ ) {
  //   for (int j = 0; j < imgSz.width; j++ ){
  //     if ( detail[i][j] <= 0.0){
  //       detail[i][j] = 0.0;
  //     }
  //
  //     ((ushort*)(dvImg->imageData + dvImg->widthStep*i))[j] =
  //       (short) (detail[i][j]);
  //   }
  // }

  //Mat muppo = cvarrToMat(dvImg);
  cv::Mat muppo = dvImg.clone();

  //take out each separate detail coefficient.
  std::vector< int > x1;
  int x2 = 0;
  for (int j=0; j < 6; j++) {
    x1.push_back(colW*(j+1));
    Mat subImg = muppo(cv::Range(0, rowW), cv::Range(x2,  x1[j]));

    Mat dest(subImg);
    Mat dest2(subImg);


    //check if subsampling was done due to unequal dimesnions, if so correct for it.
    if(sumbsample==true){
      Size upsampledsize(static_cast<int>(colsWzero), static_cast<int>(rowsWzero));
      resize(dest2, dest2, upsampledsize, INTER_LINEAR);
    }

    //check if wavelet scale is the one chosen for segmenting cell bodies
    if(j==Scale){
      cellImg = dest2.clone();
    }
    x2 = x1[j];
  }


  Mat gx, gy;
  Scharr(cellImg, gx, CV_32F, 1, 0); //CV_32F
  Scharr(cellImg, gy, CV_32F, 0, 1);

  /// Generate grad_x and grad_y
  Mat grad_x, grad_y;
  Mat abs_grad_x, abs_grad_y;

  /// Gradient X
  //Scharr( src_gray, grad_x, ddepth, 1, 0, scale, delta, BORDER_DEFAULT );
  Scharr( cellImg, grad_x, CV_32F, 1, 0, 1, 0, BORDER_DEFAULT );
  convertScaleAbs( grad_x, abs_grad_x );

  /// Gradient Y
  //Scharr( src_gray, grad_y, ddepth, 0, 1, scale, delta, BORDER_DEFAULT );
  Scharr( cellImg, grad_y, CV_32F, 0, 1, 1, 0, BORDER_DEFAULT );
  convertScaleAbs( grad_y, abs_grad_y );

  /// Total Gradient (approximate)
  Mat grad;
  addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
  //resize(grad, grad, Size(), 0.25, 0.25);

  //imshow( "Sobel test", grad );


  // Compute the structure tensor, and from it, anisotropy
  Mat gx2, gxy, gy2;
  GaussianBlur(gx.mul(gx), gx2, Size(0, 0), Sigma);
  GaussianBlur(gx.mul(gy), gxy, Size(0, 0), Sigma);
  GaussianBlur(gy.mul(gy), gy2, Size(0, 0), Sigma);
  MatExpr trace = gx2 + gy2;
  MatExpr det = gx2.mul(gy2) - gxy.mul(gxy);
  Mat second_term;
  sqrt(trace.mul(trace) / 4 - det, second_term);
  MatExpr eig1 = trace / 2 + second_term;
  MatExpr eig2 = trace / 2 - second_term;
  //compute the tensors
  //Mat anisotropy = eig1 / eig2;
  Mat energy = trace;

  double maxVal, minVal;

  minMaxLoc(energy, &minVal, &maxVal);
  energy.convertTo(energy, CV_16U, maxValinput/(maxVal - minVal), -minVal * maxValinput/(maxVal - minVal));

  string writetothisfile;
  writetothisfile = "mrd_" + outFilename + std::to_string(q) + ".tif";
  try {
    cv::imwrite(writetothisfile, energy);
    Rcpp::Rcout << writetothisfile << " SAVED" << endl;
  }catch (runtime_error& ex) {
    Rcpp::Rcout << "Cannot save: exception converting image to correct format:\n" << endl;
    return(R_NilValue);
  }


  }

  return(R_NilValue);

//return(R_NilValue);
END_RCPP
}
