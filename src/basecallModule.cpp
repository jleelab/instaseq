#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include <progress.hpp>
#include <progress_bar.hpp>
// [[Rcpp::depends(RcppProgress)]]

#include "opencv2/opencv.hpp"
#include "utils.h"


#include "Image.h"
RCPP_EXPOSED_CLASS(Stack);
RCPP_MODULE(class_Stack) {

  class_<Stack>("Stack")

  .constructor()
  .constructor<Rcpp::StringVector>("", &ImageConst1)
  .method("updateDisp", &Stack::updateDisp)
  .method("numFiles", &Stack::numFiles)
  .method("render", &Stack::render)
  .method("display", &Stack::display)
  .method("click", &Stack::click)
  .method("segment", &Stack::segment)
  .method("extract", &Stack::extract);

}

#include "display.h"
RCPP_MODULE(methods_Display) {
  function("_newDisplay", &_newDisplay, List::create(_["window_name"], _["height"],
                                          _["width"]), "");
  function("_destroyDisplay", &_destroyDisplay, List::create(_["window_name"]), "");
           function("_destroyAllDisplays", &_destroyAllDisplays, "", "");

           function("_click", &_click, List::create(_["window_name"]), "");
}
