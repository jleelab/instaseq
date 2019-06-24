// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "multithreshold.h"
// [[Rcpp::export]]

class Stack {
public:
  Stack();
  Stack(Rcpp::StringVector inputFiles);
  //cv::Mat image;
  cv::Mat src;
  cv::Mat dsp;
  Rcpp::NumericVector dim();
  int nrow(), ncol(), nchan();
  std::string depth(), space();
  int num_files;
  int cols;
  int rows;
  int planeOfinterest;
  Rcpp::NumericVector numFiles();
  std::vector<cv::Mat> zpositions;
  bool display(int slice);
  bool render();
  bool updateDisp(int plane, int Max, int Min);
  Rcpp::List click();
  bool threshold();
  Rcpp::List segment(Rcpp::IntegerVector thresh, Rcpp::IntegerVector area, Rcpp::IntegerVector channels, Rcpp::IntegerVector startStop);
  Rcpp::List extract(Rcpp::IntegerVector xpos, Rcpp::IntegerVector ypos, Rcpp::IntegerVector zpos);
private:
  void init();
  void get_quantile(cv::Mat &input, size_t smin, size_t smax,  int &alpha, int &beta );
  std::string imageSpace, imageDepth;
};

Stack::Stack() {

}


Stack::Stack(Rcpp::StringVector inputFiles) {

  num_files = inputFiles.size();
  Progress p(num_files, true);
  Rcpp::Rcout << "Loading " <<  num_files << " stacks into RAM." << std::endl;
  for (int i = 0; i < inputFiles.size(); i++){
    std::string filename(inputFiles[i]);
    //Rcpp::Rcout << "Loading " <<  filename << "\r" << std::endl;

    src = cv::imread(filename, -1);
    if ( src.empty() )      // please, *always check* resource-loading.
    {
      std::cerr << filename << " can't be loaded!" << std::endl;
      continue;
    }
    zpositions.push_back(src);
    cols = src.cols;
    rows = src.rows;

    p.increment();
    // int barWidth = 70;
    // float progress = 0.0;
    //
    // //print progress bar in console
    // Rcpp::Rcout << "  [";
    // int pos = barWidth * progress;
    // for (int j = 0; j < barWidth; ++j) {
    //   if (j < pos) Rcpp::Rcout << "=";
    //   else if (j == pos) Rcpp::Rcout << ">";
    //   else Rcpp::Rcout << " ";
    // }
    // Rcpp::Rcout << "] " << int(progress * 100.0) << "% \r" << std::flush;
    // R_FlushConsole();
    // R_ProcessEvents();
    // progress += (float)1/(num_files-1);

  }
  Rcpp::Rcout << "\n" << std::endl;
  Rcpp::Rcout << "====== LOADING DONE ======" << std::endl;

  /*

  Rcpp::Environment base = Rcpp::Environment::base_env();
  Rcpp::Function pathExpand = base["path.expand"];

  this->image = cv::imread(Rcpp::as<std::string>(pathExpand(inputFile)), cv::IMREAD_UNCHANGED);

  if (!this->image.data) {
  throw std::range_error("Could not open the image file.");
  }
  */
  this->init();
}

void Stack::get_quantile(cv::Mat &input, size_t smin, size_t smax,  int &alpha, int &beta )
{

  std::vector<float> temp_input((float*)input.data, (float*)input.data + input.rows * input.cols);

  std::sort(std::begin(temp_input), std::end(temp_input));

  beta = temp_input[smin];
  alpha =  temp_input[smax];
}

bool Stack::updateDisp(int plane, int Max, int Min) {
  std::vector<cv::Mat> tmpChannels;
  cv::Mat tmp;
  for (int channel = 0; channel < 4; ++channel) {
    tmp = zpositions.at(plane + (num_files/4)*channel ).clone();

    if(tmp.rows>800){
      double ScaleFactor = 800/(double)tmp.rows;
      cv::resize(tmp, tmp, cv::Size(), ScaleFactor, ScaleFactor, cv::INTER_CUBIC);
    }

    tmp.convertTo(tmp, CV_32F, 1.0);
    if(Max == 0){
      this->get_quantile(tmp, tmp.total() * 0.001, tmp.total() * 0.999, Max, Min);
    }
    tmp -= Min;
    tmp.convertTo(tmp, CV_8UC1, 255.0/(Max-0));

    tmpChannels.push_back(tmp);
  }

  std::vector<cv::Mat> redCh;
  redCh.push_back(tmpChannels[1]);
  redCh.push_back(tmpChannels[2]);

  std::vector<cv::Mat> greenCh;
  greenCh.push_back(tmpChannels[0]);
  greenCh.push_back(tmpChannels[3]);

  //dsp = cv::Mat::zeros(tmpChannels[0].rows, tmpChannels[0].cols, CV_8UC3);
  cv::Mat red = cv::Mat::zeros(tmpChannels[0].rows, tmpChannels[0].cols, CV_8UC1);
  cv::Mat blue = cv::Mat::zeros(tmpChannels[0].rows, tmpChannels[0].cols, CV_8UC1);
  cv::Mat green = cv::Mat::zeros(tmpChannels[0].rows, tmpChannels[0].cols, CV_8UC1);

  //cv::merge(redCh, red);
  //cv::merge(greenCh, green);
  cv::add(redCh[0], redCh[1], red);
  cv::add(greenCh[0], greenCh[1], green);
  cv::add(tmpChannels[0], tmpChannels[1], blue);

  std::vector<cv::Mat> colorArray;

  colorArray.push_back(blue);
  colorArray.push_back(green);
  colorArray.push_back(red);

  dsp = cv::Mat::zeros(tmpChannels[0].rows, tmpChannels[0].cols, CV_8UC3);

  cv::merge(colorArray, dsp);
  this->render();
  return true;
}
//cyan RGB = 0 1 1
//magenta RGB = 1 0 1
//mupp[20 + (length(mupp)/length(unique(channels)))*(seq_along(unique(channels))-1)]

void Stack::init() {
  planeOfinterest = zpositions.size()/8;

  this->updateDisp(planeOfinterest, 0, 0);
}

bool Stack::display(int slice) {
  cv::imshow("Display window", zpositions.at(slice));
  cv::waitKey(1);
  return true;
}

Rcpp::NumericVector Stack::numFiles() {
  return Rcpp::NumericVector::create(this->zpositions.size());
}

bool Stack::render() {

  cv::imshow("Display window", dsp);
  cv::waitKey(1);
  return true;
}


// bool Stack::threshold(std::vector<int> intensity, int sizeMax, int sizeMin) {
//
//   for (int channel = 0; channel < 4; ++channel) {
//     cv::Mat tmp = zpositions.at(planeOfinterest + (num_files/4)*channel );
//   }
//
//   cv::imshow("Display window", dsp);
//   cv::waitKey(1);
//   return true;
// }

struct ClickData {
  int x;
  int y;
  int button;
};

void onMouse(int event, int x, int y, int flags, void* data) {
  ClickData* p = (ClickData*) data;

  if  (event == cv::EVENT_LBUTTONDOWN) {
    p->x = x;
    p->y = y;
    p->button = 0;
  } else if (event == cv::EVENT_RBUTTONDOWN) {
    p->x = x;
    p->y = y;
    p->button = 1;
  }
}

Rcpp::List Stack::click() {
  std::vector<ushort> buf;
  std::vector<int> ch;
  std::vector<int> up;
  std::vector<int> selection;

  ClickData d;
  d.button = -1;
  int clickNumber = 0;
  while(d.button != 1){
  clickNumber++;
  d.button = -1;
  cv::setMouseCallback("Display window", onMouse, &d);

  while (d.button == -1) {
    cv::waitKey(10);
  }

  Rcpp::Rcout << "x: " << d.x << std::endl;
  Rcpp::Rcout << "y: " << d.y << std::endl;
  Rcpp::Rcout << "button: " << d.button << std::endl;


  cv::Mat tmp;
  tmp = zpositions.at(planeOfinterest).clone();
  if(tmp.rows>800){
    double ScaleFactor = (double)tmp.rows/800;
    d.x = ScaleFactor*d.x;
    d.y = ScaleFactor*d.y;
  }


  for (int channel = 0; channel < 4; ++channel) {
    cv::Mat tmp = zpositions.at(planeOfinterest + (num_files/4)*channel );

    cv::LineIterator horiz(tmp, cv::Point(d.x-20,d.y), cv::Point(d.x+20,d.y), 8);
    cv::LineIterator verti(tmp, cv::Point(d.x,d.y-20), cv::Point(d.x,d.y+20), 8);

    for(int i = 0; i < horiz.count; i++, ++horiz)
    {
      buf.push_back(tmp.at<ushort>(horiz.pos()));
      ch.push_back(channel);
      up.push_back(0);
      selection.push_back(clickNumber);
    }

    for(int k = 0; k < verti.count; k++, ++verti)
    {
      buf.push_back(tmp.at<ushort>(verti.pos()));
      ch.push_back(channel);
      up.push_back(1);
      selection.push_back(clickNumber);
    }

  }
  }


  return Rcpp::List::create(Rcpp::Named("selection") = selection,
                            Rcpp::Named("values") = buf,
                            Rcpp::Named("channel") = ch,
                            Rcpp::Named("verctical") = up);
}

Rcpp::List Stack::segment(Rcpp::IntegerVector thresh, Rcpp::IntegerVector area, Rcpp::IntegerVector channels, Rcpp::IntegerVector startStop) {
  Multithreshold segmented(zpositions, thresh, area, channels, num_files);
  segmented.runThresh(startStop[0], startStop[1]);

  return Rcpp::List::create(
    Rcpp::_["channel"] = segmented.rolonychannel,
    Rcpp::_["cyan"] = segmented.channel01,
    Rcpp::_["red"] = segmented.channel02,
    Rcpp::_["magenta"] = segmented.channel03,
    Rcpp::_["green"] = segmented.channel04,
    Rcpp::_["x"] = segmented.xpos,
    Rcpp::_["y"] = segmented.ypos,
    Rcpp::_["z"] = segmented.zpos,
    Rcpp::_["m0102"] = segmented.M0102,
    Rcpp::_["m0103"] = segmented.M0103,
    Rcpp::_["m0104"] = segmented.M0104,
    Rcpp::_["m0201"] = segmented.M0201,
    Rcpp::_["m0203"] = segmented.M0203,
    Rcpp::_["m0204"] = segmented.M0204,
    Rcpp::_["m0301"] = segmented.M0301,
    Rcpp::_["m0302"] = segmented.M0302,
    Rcpp::_["m0304"] = segmented.M0304,
    Rcpp::_["m0401"] = segmented.M0401,
    Rcpp::_["m0402"] = segmented.M0402,
    Rcpp::_["m0403"] = segmented.M0403
  );
  //cv::imshow("Display window", tmp);
  //cv::waitKey(1);
}

Rcpp::List Stack::extract(Rcpp::IntegerVector xpos, Rcpp::IntegerVector ypos, Rcpp::IntegerVector zpos)
{
 std::vector<int> value;
 std::vector<int> direction;
 std::vector<int> channelVec;
 std::vector<int> rolony;
 int counter = 1;
 int largest_element = 0;
 int horiz_element = 0;
 int verti_element = 0;
 for(int i = 0; i < xpos.size(); i++){

   for (int channel = 0; channel < 4; ++channel) {

     cv::LineIterator horiz(zpositions.at(zpos[i] + (num_files/4)*channel ), cv::Point(xpos[i]-20, ypos[i]), cv::Point(xpos[i]+20, ypos[i]), 8);
     cv::LineIterator verti(zpositions.at(zpos[i] + (num_files/4)*channel ), cv::Point(xpos[i], ypos[i]-20), cv::Point(xpos[i], ypos[i]+20), 8);

     largest_element = zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(horiz.pos());
     horiz_element = 0;
     for(int j = 0; j < horiz.count; j++, ++horiz)
     {
       if(zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(horiz.pos()) > largest_element){
         largest_element = zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(horiz.pos());
         horiz_element = j;
       }
     }
     horiz_element = horiz_element - 21;

     largest_element = zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(verti.pos());
     verti_element = 0;
     for(int k = 0; k < verti.count; k++, ++verti)
     {
       if(zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(horiz.pos()) > largest_element){
         largest_element = zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(verti.pos());
         verti_element = k;
       }
     }
     verti_element = verti_element - 21;

     cv::LineIterator horiz2(zpositions.at(zpos[i] + (num_files/4)*channel ), cv::Point(xpos[i]+horiz_element-20, ypos[i]), cv::Point(xpos[i]+horiz_element+20, ypos[i]), 8);
     cv::LineIterator verti2(zpositions.at(zpos[i] + (num_files/4)*channel ), cv::Point(xpos[i], ypos[i]+verti_element-20), cv::Point(xpos[i], ypos[i]+verti_element+20), 8);

     for(int j = 0; j < horiz2.count; j++, ++horiz2)
     {
       value.push_back(zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(horiz2.pos()));
       rolony.push_back(counter);
       direction.push_back(1);
       channelVec.push_back(channel);
     }

     for(int k = 0; k < verti2.count; k++, ++verti2)
     {
       value.push_back(zpositions.at(zpos[i] + (num_files/4)*channel ).at<ushort>(verti2.pos()));
       rolony.push_back(counter);
       direction.push_back(0);
       channelVec.push_back(channel);

     }

   }
   counter++;
 }


  return Rcpp::List::create(
    Rcpp::_["id"] = rolony,
    Rcpp::_["value"] = value,
    Rcpp::_["channel"] = channelVec
  );
}


