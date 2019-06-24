void _newDisplay(std::string window_name, int height, int width) {
  cv::Mat blank = cv::Mat::zeros(height, width, CV_8UC3);
  cv::namedWindow(window_name, cv::WINDOW_AUTOSIZE);
  cv::imshow(window_name, blank);
  cv::waitKey(1);
}

void _destroyDisplay(std::string window_name) {
  cv::destroyWindow(window_name);
  cv::waitKey(1);
}

void _destroyAllDisplays() {
  cv::destroyAllWindows();
  cv::waitKey(1);
}

struct sClickData {
  int x;
  int y;
  int button;
};

void sonMouse(int event, int x, int y, int flags, void* data) {
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

Rcpp::DataFrame _click(std::string window_name) {
  ClickData d;
  d.button = -1;
  cv::setMouseCallback(window_name, sonMouse, &d);

  while (d.button == -1) {
    cv::waitKey(10);
  }
  return Rcpp::DataFrame::create(Rcpp::Named("x") = d.x,
                                 Rcpp::Named("y") = d.y,
                                 Rcpp::Named("button") = d.button);
}
