#include "multithreshold.h"
#include <progress.hpp>
#include <progress_bar.hpp>


#include <cmath>
#include <cstdio>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

Multithreshold::Multithreshold(std::vector<cv::Mat> imgStack, Rcpp::IntegerVector thresholds, Rcpp::IntegerVector allowedareas, Rcpp::IntegerVector channels, int planes){
  zstack = imgStack;

  for (int i=0; i<thresholds.size(); ++i){
    thresh.push_back(thresholds[i]);
    area.push_back(allowedareas[i]);
    ch.push_back(channels[i]);
  }

  num_files = planes;
}


void Multithreshold::runThresh(int startPlane, int stopPlane)
{
  int processPlane = stopPlane - startPlane + 1;
  Progress p(processPlane*4, true);
  cv::Mat src = zstack.at(0).clone();

  cv::Mat dst = cv::Mat::zeros(src.size(), CV_8UC1);

  int rows = src.rows;
  int cols = src.cols;
  std::vector<std::vector<cv::Mat>> srcStack(4, std::vector<cv::Mat>(processPlane, cv::Mat(rows,cols,CV_16U)) );

  Rcpp::Rcout << "Segmenting rolonies" << std::endl;

  for(int plane = startPlane; plane <= stopPlane; ++plane){
    //int planeOfinterest = zstack.size()/8;

    for (int channel = 0; channel < 4; ++channel) {

      cv::Mat tmp;
      src = zstack.at(plane + (num_files/4)*channel ).clone();
      dst = cv::Mat::zeros(src.size(), CV_8UC1);

      std::vector<int> threshOfInterest;

      for(int i=0; i < ch.size(); i++){
        if(ch[i] == channel){
          threshOfInterest.push_back(thresh[i]);
        }
      }


      for(int it=0; it < threshOfInterest.size(); it++) {
        cv::threshold(src, tmp, threshOfInterest[it], 65536, cv::THRESH_BINARY );
        tmp.convertTo(tmp, CV_8UC1);

        std::vector<std::vector<cv::Point> > contours;
        std::vector<cv::Vec4i> hierarchy;
        cv::findContours(tmp, contours, hierarchy, cv::RETR_EXTERNAL, cv::CHAIN_APPROX_SIMPLE );

        cv::Mat mask = cv::Mat::zeros(src.size(), CV_8UC1);
        for (size_t k = 0; k < contours.size(); ++k)
        {
          const std::vector<cv::Point>& contour = contours[k];
          double area0 = cv::contourArea(contour);

          //if area falls within the area limits given by R via "alim" parameter then draw the contour.
          if( (area0 >= area[0]) && (area0 <= area[1]) ){
            cv::Scalar color(255, 255, 255);
            //cv::Scalar blackOutline(0, 0, 0);
            cv::drawContours( mask, contours, k, color, cv::FILLED, 8, hierarchy );
            //cv::drawContours( mask, contours, k, blackOutline, 3, 8, hierarchy );

          }
          //add the contours to the destination Mat.
          cv::bitwise_or(dst, mask, dst);
        }



      }//threshold iterator

      double ScaleFactor = 800/(double)dst.rows;
      cv::Mat show;
      cv::resize(dst, show, cv::Size(), ScaleFactor, ScaleFactor, cv::INTER_CUBIC);
      std::string s = std::to_string(plane);
      cv::putText(show, //target image
                  s, //text
                  cv::Point(10, show.rows / 2), //top-left position
                  cv::FONT_HERSHEY_DUPLEX,
                  1.0,
                  CV_RGB(255, 255, 255), //font color
                  2);
      s = std::to_string(channel);
      cv::putText(show, //target image
                  s, //text
                  cv::Point(10, show.rows / 2 + 20), //top-left position
                  cv::FONT_HERSHEY_DUPLEX,
                  1.0,
                  CV_RGB(255, 255, 255), //font color
                  2);

      cv::imshow("Segment", show);
      cv::waitKey(1);
      srcStack.at(channel).at(plane - startPlane) = zstack.at(plane + (num_files/4)*channel ).clone();
      zstack.at(plane + (num_files/4)*channel ) = dst.clone();
      p.increment();

    }//channel iterator

  }//plane iterator

  //int rows = zstack.at(0).rows;
  //int cols = zstack.at(0).cols;
  std::vector<std::vector<cv::Mat>> binaryZstack(4, std::vector<cv::Mat>(processPlane, cv::Mat(rows,cols,CV_16U)) );
  std::vector<std::vector<cv::Mat>> destination(4, std::vector<cv::Mat>(processPlane, cv::Mat(rows,cols,CV_16U)) );

  int planeCounter;
  for (int channeL = 0; channeL < 4; ++channeL) {

    planeCounter = 0;
  for(int planeL = startPlane; planeL <= stopPlane; ++planeL){
    cv::Mat dsp = zstack.at(planeL + (num_files/4)*channeL ).clone();
    cv::Mat tmp;
    cv::threshold(dsp, tmp, 0, 1, cv::THRESH_BINARY | cv::THRESH_OTSU );
    tmp.convertTo(tmp, CV_16U);


    binaryZstack.at(channeL).at(planeCounter) = tmp;


    cv::Mat dst = cv::Mat::zeros(src.size(), CV_16U);

    destination.at(channeL).at(planeCounter) = dst;

    double ScaleFactor = 800/(double)dsp.rows;
    cv::resize(dsp, dsp, cv::Size(), ScaleFactor, ScaleFactor, cv::INTER_CUBIC);
    cv::imshow("Results", dsp);
    cv::waitKey(80);
    planeCounter++;

  }
  int mupp = connectedcomponents3D ( processPlane, rows, cols, channeL, binaryZstack.at(channeL), destination.at(channeL), srcStack, zpos, ypos, xpos, rolonychannel, channel01, channel02, channel03, channel04);

  }



  //double ScaleFactor = 800/(double)dst.rows;
  //cv::resize(dst, dst, cv::Size(), ScaleFactor, ScaleFactor, cv::INTER_CUBIC);
  //cv::imshow("Segment", dst);
  //cv::waitKey(1);

}


double Multithreshold::sum(std::vector<double> a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += a[i];
  }
  return s;
}

double Multithreshold::mean(std::vector<double> a)
{
  return sum(a) / a.size();
}

double Multithreshold::sqsum(std::vector<double> a)
{
  double s = 0;
  for (int i = 0; i < a.size(); i++)
  {
    s += pow(a[i], 2);
  }
  return s;
}

double Multithreshold::stdev(std::vector<double> nums)
{
  double N = nums.size();
  return pow(sqsum(nums) / N - pow(sum(nums) / N, 2), 0.5);
}

std::vector<double> Multithreshold::subtrc(std::vector<double> a, double b)
{
  std::vector<double> retvect;
  for (int i = 0; i < a.size(); i++)
  {
    retvect.push_back(a[i] - b);
  }
  return retvect;
}

std::vector<double> Multithreshold::multi(std::vector<double> a, std::vector<double> b)
{
  std::vector<double> retvect;
  for (int i = 0; i < a.size() ; i++)
  {
    retvect.push_back(a[i] * b[i]);
  }
  return retvect;
}

double Multithreshold::pearsoncoeff(std::vector<double> X, std::vector<double> Y)
{
  std::vector<double> numerator = multi(subtrc(X, mean(X)), subtrc(Y, mean(Y)));

  std::vector<double> dividend;
  dividend.reserve(X.size());

  double denominator = (X.size()*stdev(X)*stdev(Y));

  std::transform(numerator.begin(), numerator.end(), std::back_inserter( dividend ),
                 std::bind1st(std::divides<double>(), denominator));

  return sum( dividend );
}


//****************************************************************************80

int Multithreshold::i4_min ( int i1, int i2 )

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    I4_MIN returns the minimum of two I4's.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    13 October 1998
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int I1, I2, two integers to be compared.
  //
  //    Output, int I4_MIN, the smaller of I1 and I2.
  //
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}

int Multithreshold::connectedcomponents3D ( int l, int m, int n, int channel, std::vector<cv::Mat>& src, std::vector<cv::Mat>& dst, std::vector<std::vector<cv::Mat>>& stack, std::vector<double>&  centerofmassZ,  std::vector<double>&  centerofmassX,  std::vector<double>&  centerofmassY, std::vector<int>&  rolonychannel, std::vector<double>& channel01, std::vector<double>& channel02, std::vector<double>& channel03, std::vector<double>& channel04) // int a[], int c[]

  //****************************************************************************80
  //
  //  Purpose:
  //
  //    connectedcomponents3D assigns contiguous nonzero pixels to a common component.
  //
  //  Discussion:
  //
  //    On input, the A array contains values of 0 or 1.
  //
  //    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
  //    into connected components.
  //
  //    The pixel A(I,J,K) is "connected" to the pixels:
  //
  //      A(I-1,J,  K  ),  A(I+1,J,  K  ),
  //      A(I,  J-1,K  ),  A(I,  J+1,K  ),
  //      A(I,  J,  K-1),  A(I,  J,  K+1),
  //
  //    so most pixels have 6 neighbors.
  //
  //    On output, COMPONENT_NUM reports the number of components of nonzero
  //    data, and the array C contains the component assignment for
  //    each nonzero pixel, and is 0 for zero pixels.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license.
  //
  //  Modified:
  //
  //    01 February 2012
  //
  //  Author:
  //
  //    John Burkardt
  //
  //  Parameters:
  //
  //    Input, int L, M, N, the order of the array.
  //
  //    Input, int A[L*M*N], the pixel array.
  //
  //    Output, int C[L*M*N], the component array.
  //
  //    Output, int I4BLOCK_COMPONENTS, the number of components
  //    of nonzero data.
  //
{
  int b;
  int c1;
  int component;
  int component_num;
  int i;
  int j;
  int k;
  int north;
  int *p;
  int *q;
  int up;
  int west;


  //
  //  Initialization.  src.at<uchar>(i,j)
  //
  /*for ( k = 0; k < n; k++ )
  {
  for ( j = 0; j < m; j++ )
  {
  for ( i = 0; i < l; i++ )
  {
  c[i+j*l+k*l*m] = 0;
  }
  }
  } */
  component_num = 0;
  //
  //  P is simply used to store the component labels.  The dimension used
  //  here is, of course, usually an absurd overestimate.
  //
  p = new int[l*m*n+1];
  for ( i = 0; i <= l * m * n; i++ )
  {
    p[i] = i;
  }
  //
  //  "Read" the array one pixel at a time.  If a (nonzero) pixel has a north or
  //  west neighbor with a label, the current pixel inherits it.
  //  In case the labels disagree, we need to adjust the P array so we can
  //  later deal with the fact that the two labels need to be merged.
  //
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        //get the array content in each direction N-S, W-E, U-D
        //NORTH-SOUTH
        if ( i == 0 )
        {
          north = 0;
        }
        else
        {
          north = dst.at(i-1).at<ushort>(j,k); //north = c[i-1+j*l+k*l*m];
        }
        //WEST-EAST
        if ( j == 0 )
        {
          west = 0;
        }
        else
        {
          west = dst.at(i).at<ushort>(j-1,k); //west = c[i+(j-1)*l+k*l*m];
        }
        //UP-DOWN
        if ( k == 0 )
        {
          up = 0;
        }
        else
        {
          up = dst.at(i).at<ushort>(j,k-1);//up = c[i+j*l+(k-1)*l*m];
        }
        //check if target pixel is empty  a[i+j*l+k*l*m] != 0
        if ( src.at(i).at<ushort>(j,k) != 0 )
        {
          //
          //  New component?
          //
          if ( north == 0 && west == 0 && up == 0 )
          {
            component_num = component_num + 1;
            dst.at(i).at<ushort>(j,k) = component_num;//  c[i+j*l+k*l*m] = component_num;
          }
          //
          //  One predecessor is labeled.
          //
          else if ( north != 0 && west == 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = north; // c[i+j*l+k*l*m] = north;
          }
          else if ( north == 0 && west != 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = west; //c[i+j*l+k*l*m] = west;
          }
          else if ( north == 0 && west == 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = up; //c[i+j*l+k*l*m] = up;
          }
          //
          //  Two predecessors are labeled.
          //
          else if ( north == 0 && west != 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( west, up );//c[i+j*l+k*l*m] = i4_min ( west, up );
            c1 = i4_min ( p[west], p[up] );
            p[west] = c1;
            p[up] = c1;
          }
          else if ( north != 0 && west == 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, up );//c[i+j*l+k*l*m] = i4_min ( north, up );
            c1 = i4_min ( p[north], p[up] );
            p[north] = c1;
            p[up] = c1;
          }
          else if ( north != 0 && west != 0 && up == 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, west );//c[i+j*l+k*l*m] = i4_min ( north, west );
            c1 = i4_min ( p[north], p[west] );
            p[north] = c1;
            p[west] = c1;
          }
          //
          //  Three predecessors are labeled.
          //
          else if ( north != 0 && west != 0 && up != 0 )
          {
            dst.at(i).at<ushort>(j,k) = i4_min ( north, i4_min ( west, up ) ); //c[i+j*l+k*l*m] = i4_min ( north, i4_min ( west, up ) );
            c1 = i4_min ( p[north], i4_min ( p[west], p[up] ) );
            p[north] = c1;
            p[west] = c1;
            p[up] = c1;
          }
        }
      }
    }
  }
  //
  //  When a component has multiple labels, have the higher labels
  //  point to the lowest one.
  //
  for ( component = component_num; 1 <= component; component-- )
  {
    b = component;
    while ( p[b] != b )
    {
      b = p[b];
    }
    p[component] = b;
  }
  //
  //  Locate the minimum label for each component.
  //  Assign these mininum labels new consecutive indices.
  //
  q = new int[component_num+1];

  for ( j = 0; j <= component_num; j++ )
  {
    q[j] = 0;
  }

  i = 0;
  for ( component = 1; component <= component_num; component++ )
  {
    if ( p[component] == component )
    {
      i = i + 1;
      q[component] = i;
    }
  }

  component_num = i;

  //for statistics
  double X[component_num];
  double Y[component_num];
  double Z[component_num];
  double ch01[component_num];
  double ch02[component_num];
  double ch03[component_num];
  double ch04[component_num];

  double m0102[component_num];
  double m0103[component_num];
  double m0104[component_num];

  double m0201[component_num];
  double m0203[component_num];
  double m0204[component_num];

  double m0301[component_num];
  double m0302[component_num];
  double m0304[component_num];

  double m0401[component_num];
  double m0402[component_num];
  double m0403[component_num];

  uint number_voxels[component_num];
  for(i = 0; i < component_num; i++){
    X[i] = 0;
    Y[i] = 0;
    Z[i] = 0;
    number_voxels[i] = 0;
  }
  //
  //  Replace the labels by consecutive labels.
  //
  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        int index = dst.at(i).at<ushort>(j,k);
        int value = q [ p [  index ] ];
        dst.at(i).at<ushort>(j,k) = value; //q [ p [  index ] ]; // c[i+j*l+k*l*m] = q [ p [ c[i+j*l+k*l*m] ] ];
        Z[value] += (double)i/l;
        X[value] += (double)j/m;
        Y[value] += (double)k/n;
        number_voxels[value]++;
        ch01[value] = 0;
        ch02[value] = 0;
        ch03[value] = 0;
        ch04[value] = 0;

        m0102[value] = 0;
        m0103[value] =  0;
        m0104[value] =  0;

        m0201[value] =  0;
        m0203[value] =  0;
        m0204[value] =  0;

        m0301[value] =  0;
        m0302[value] =  0;
        m0304[value] =  0;

        m0401[value] =  0;
        m0402[value] =  0;
        m0403[value] =  0;

      }
    }
  }

  int maxColor = 4096;

  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      for ( k = 0; k < n; k++ )
      {
        int value = (int)dst.at(i).at<ushort>(j,k);

        ch01[value] += (double)stack.at(0).at(i).at<ushort>(j,k)/maxColor;
        ch02[value] += (double)stack.at(1).at(i).at<ushort>(j,k)/maxColor;
        ch03[value] += (double)stack.at(2).at(i).at<ushort>(j,k)/maxColor;
        ch04[value] += (double)stack.at(3).at(i).at<ushort>(j,k)/maxColor;

        m0102[value] += (int)stack.at(1).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(0).at(i).at<ushort>(j,k)/maxColor : 0;
        m0103[value] += (int)stack.at(2).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(0).at(i).at<ushort>(j,k)/maxColor : 0;
        m0104[value] += (int)stack.at(3).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(0).at(i).at<ushort>(j,k)/maxColor : 0;

        m0201[value] += (int)stack.at(0).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(1).at(i).at<ushort>(j,k)/maxColor : 0;
        m0203[value] += (int)stack.at(2).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(1).at(i).at<ushort>(j,k)/maxColor : 0;
        m0204[value] += (int)stack.at(3).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(1).at(i).at<ushort>(j,k)/maxColor : 0;

        m0301[value] += (int)stack.at(0).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(2).at(i).at<ushort>(j,k)/maxColor : 0;
        m0302[value] += (int)stack.at(1).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(2).at(i).at<ushort>(j,k)/maxColor : 0;
        m0304[value] += (int)stack.at(3).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(2).at(i).at<ushort>(j,k)/maxColor : 0;

        m0401[value] += (int)stack.at(0).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(3).at(i).at<ushort>(j,k)/maxColor : 0;
        m0402[value] += (int)stack.at(1).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(3).at(i).at<ushort>(j,k)/maxColor : 0;
        m0403[value] += (int)stack.at(2).at(i).at<ushort>(j,k) > 30 ? (double)stack.at(3).at(i).at<ushort>(j,k)/maxColor : 0;

      }
    }
  }

  std::cout << (int)stack.at(0).at(28).at<ushort>(772,1070) << std::endl;
  std::cout << (int)stack.at(1).at(28).at<ushort>(772,1070) << std::endl;
  std::cout << (int)stack.at(2).at(28).at<ushort>(772,1070) << std::endl;
  std::cout << (int)stack.at(3).at(28).at<ushort>(772,1070) << std::endl;
  std::cout << ch01[1054] << std::endl;
  std::cout << ch02[1054] << std::endl;
  std::cout << ch03[1054] << std::endl;
  std::cout << ch04[1054] << std::endl;
  std::cout << m0102[1054] << std::endl;
  std::cout << m0103[1054] << std::endl;
  std::cout << m0104[1054] << std::endl;

  for(i = 0; i < component_num; i++){
    centerofmassZ.push_back( ((double)Z[i]/number_voxels[i])  );
    centerofmassX.push_back( ((double)X[i]/number_voxels[i])  );
    centerofmassY.push_back( ((double)Y[i]/number_voxels[i])  );
    channel01.push_back( ((double)ch01[i]/number_voxels[i])  );
    channel02.push_back( ((double)ch02[i]/number_voxels[i])  );
    channel03.push_back( ((double)ch03[i]/number_voxels[i])  );
    channel04.push_back( ((double)ch04[i]/number_voxels[i])  );
    rolonychannel.push_back( channel  );

    M0102.push_back( (double)(m0102[i])/(ch01[i]));
    M0103.push_back( (double)(m0103[i])/(ch01[i]));
    M0104.push_back( (double)(m0104[i])/(ch01[i]));

    M0201.push_back( (double)(m0201[i])/(ch02[i]));
    M0203.push_back( (double)(m0203[i])/(ch02[i]));
    M0204.push_back( (double)(m0204[i])/(ch02[i]));

    M0301.push_back( (double)(m0301[i])/(ch03[i]));
    M0302.push_back( (double)(m0302[i])/(ch03[i]));
    M0304.push_back( (double)(m0304[i])/(ch03[i]));

    M0401.push_back((double)(m0401[i])/(ch04[i]));
    M0402.push_back((double)(m0402[i])/(ch04[i]));
    M0403.push_back((double)(m0403[i])/(ch04[i]));


  }
  std::string filepath;

  for ( i = 0; i < l; i++ )
  {

    filepath = "/Users/danielfurth/Documents/GitHub/PRICKly/prickly/multipoint/out/" + std::to_string(channel) + "_" + std::to_string(i) + "_" + ".tif";

    cv::imwrite(filepath, dst.at(i));

    filepath = "/Users/danielfurth/Documents/GitHub/PRICKly/prickly/multipoint/out/Hello_" + std::to_string(channel) + "_" + std::to_string(i) + "_" + ".tif";
    cv::imwrite(filepath, stack.at(0).at(i));

  }

  delete [] p;
  delete [] q;

  //std::cout << ch01[2] << std::endl;
  //std::cout << channel01[2] << std::endl;

  return component_num;
}

