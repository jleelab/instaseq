// OpenCV
#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

//for sleep
#include <unistd.h>

#include <Rcpp.h>


using namespace cv;
using namespace Rcpp;
using namespace std; //bad practise remove when you have time.

#include "CThinPlateSpline.h"


RcppExport SEXP forwardWarp(SEXP mx, SEXP my, SEXP transMX, SEXP transMY, SEXP mupp){
  BEGIN_RCPP
  Rcpp::RNGScope __rngScope;


  Rcpp::NumericMatrix mX(mx);
  Rcpp::NumericMatrix mY(my);
  Rcpp::NumericMatrix transmx(transMX);
  Rcpp::NumericMatrix transmy(transMY);
  int nrows = mX.nrow();
  int ncolumns = mX.ncol();

  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncolumns; j++) {
      if(( transmx(i,j)>=0)&(transmx(i,j)<ncolumns)&(transmy(i,j)>=0)&(transmy(i,j)<nrows)){
        mX(transmy(i,j), transmx(i,j) ) = j;
        mY(transmy(i,j), transmx(i,j) ) = i;
      }
    }
  }


  int top;
  int bottom;
  int right;
  int left;

  std::vector<double> vecX;
  std::vector<double> vecY;
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncolumns; j++) {
      if( mX(i,j) == -999 ){
        vecX.clear();
        vecY.clear();
        if( i >= (nrows-1) ){
          top = 0;
        }else{
          top = 1;
        }
        if(mX(i+top, j)!= -999 ){
          vecX.push_back( mX(i+top, j) );
          vecY.push_back( mY(i+top, j) );
        }

        if(i == 0){
          bottom = 0;
        }else{
          bottom = 1;
        }
        if( mX(i-bottom, j)  != -999 ){
          vecX.push_back( mX(i-bottom, j) );
          vecY.push_back( mY(i-bottom, j) );
        }

        if(j>=(ncolumns-1) ){
          right = 0;
        }else{
          right = 1;
        }
        if( mX(i, j+right)  != -999 ){
          vecX.push_back( mX(i, j+right) );
          vecY.push_back( mY(i, j+right) );
        }

        if(j == 0){
          left = 0;
        }else{
          left = 1;
        }
        if( mX(i, j-left)  != -999 ){
          vecX.push_back( mX(i, j-left) );
          vecY.push_back( mY(i, j-left) );
        }

        double averageX = std::accumulate(vecX.begin(), vecX.end(), 0.0) / vecX.size();
        double averageY = std::accumulate(vecY.begin(), vecY.end(), 0.0) / vecY.size();
        mX(i,j) = averageX;
        mY(i,j) = averageY;

      }
    }
  }

  return List::create(
    _["mx"] = mX,
    _["my"] = mY
  );

  END_RCPP
}

RcppExport SEXP ThinPlateRegistration(SEXP srcX, SEXP srcY, SEXP dstX, SEXP dstY){
  BEGIN_RCPP
  Rcpp::RNGScope __rngScope;

  // generate some generic points
  std::vector<cv::Point> iP, iiP;

  NumericVector srX(srcX);
  NumericVector srY(srcY);
  NumericVector dtX(dstX);
  NumericVector dtY(dstY);

  Mat img = Mat::zeros(1760, 1760, CV_32F);

  // push some points into the vector for the source image
  for (unsigned i=0; i < srX.size(); i++) {
    iP.push_back(cv::Point(srX[i],srY[i]));
    iiP.push_back(cv::Point(dtX[i],dtY[i]));
  }
  iP.push_back(cv::Point(img.cols,img.rows));
  iP.push_back(cv::Point(0,0));
  iP.push_back(cv::Point(0,img.rows));
  iP.push_back(cv::Point(img.cols,0));

  iiP.push_back(cv::Point(img.cols,img.rows));
  iiP.push_back(cv::Point(0,0));
  iiP.push_back(cv::Point(0,img.rows));
  iiP.push_back(cv::Point(img.cols,0));

  CThinPlateSpline tps(iP,iiP);

  // warp the image to dst
  Mat dst;

  tps.warpImage(img,dst,0.01,INTER_CUBIC,BACK_WARP);

  Mat mx;
  Mat my;
  tps.getMaps(mx, my);


  vector<float> X;
  vector<float> Y;


  for (int i=0; i < mx.cols; i++) {
    for (int j =0; j < mx.rows; j++){
      float tempX, tempY;
      tempX = mx.at<float>(j,i);
      X.push_back(tempX);
      tempY = my.at<float>(j,i);
      Y.push_back(tempY);

    }
  }



  NumericMatrix Mx(mx.rows, mx.cols, X.begin() );
  NumericMatrix My(my.rows, my.cols, Y.begin() );

  return List::create(
    _["mx"] = Mx,
    _["my"] = My
  );

  END_RCPP

}
