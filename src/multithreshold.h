#ifndef MULTITHRESH_H
#define MULTITHRESH_H

#include "opencv2/opencv.hpp"
#include "opencv2/opencv_modules.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <unistd.h>
#include <Rcpp.h>
/// This useless class does nothing else
/// than to store a number.
class Multithreshold
{
public:
  /// Constructor which sets the number to be stored
  Multithreshold(std::vector<cv::Mat> imgStack, Rcpp::IntegerVector thresholds, Rcpp::IntegerVector allowedareas, Rcpp::IntegerVector channels, int planes);
  /// Returns the stored number
  void runThresh(int start, int stop);
  std::vector<int> thresh;
  std::vector<int> area;
  std::vector<int> ch;
  std::vector<cv::Mat> binaryZstack;
  std::vector<double> zpos;
  std::vector<double> xpos;
  std::vector<double> ypos;
  std::vector<int> rolonychannel;
  std::vector<double> channel01;
  std::vector<double> channel02;
  std::vector<double> channel03;
  std::vector<double> channel04;

  std::vector<double> M0102;
  std::vector<double> M0103;
  std::vector<double> M0104;

  std::vector<double> M0201;
  std::vector<double> M0203;
  std::vector<double> M0204;

  std::vector<double> M0301;
  std::vector<double> M0302;
  std::vector<double> M0304;

  std::vector<double> M0401;
  std::vector<double> M0402;
  std::vector<double> M0403;

  int i4_min ( int i1, int i2 );
  int connectedcomponents3D ( int l, int m, int n, int channel, std::vector<cv::Mat>& src, std::vector<cv::Mat>& dst, std::vector<std::vector<cv::Mat>>& stack, std::vector<double>&  centerofmassZ,  std::vector<double>&  centerofmassX,  std::vector<double>&  centerofmassY, std::vector<int>&  rolonychannel, std::vector<double>& channel01, std::vector<double>& channel02, std::vector<double>& channel03, std::vector<double>& channel04);

private:
  std::vector<cv::Mat> zstack;
  int num_files;
  double sum(std::vector<double> a);
  double mean(std::vector<double> a);
  double sqsum(std::vector<double> a);
  double stdev(std::vector<double> nums);
  double pearsoncoeff(std::vector<double> X, std::vector<double> Y);
  std::vector<double> subtrc(std::vector<double> a, double b);
  std::vector<double> multi(std::vector<double> a, std::vector<double> b);
};
#endif // MULTITHRESH_H
