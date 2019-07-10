#ifndef H_SLIDE_MAKER
#define H_SLIDE_MAKER

#include <fstream>
#include <string>
#include <vector>

class SlideMaker{
public:
  SlideMaker(std::string fname, std::string aspect_ratio = "169");
  ~SlideMaker() = default;

  void AddSlide(std::vector<std::string> pnames, unsigned ncols = 0, std::string title = "", 
    std::vector<std::string> col_labels = std::vector<std::string>(),
    std::vector<std::string> row_labels = std::vector<std::string>());
  void AddSlideWithReplace(std::string oldexp, std::string newexp, 
  	std::vector<std::string> pnames, unsigned ncols = 0, std::string title = "",
    std::vector<std::string> col_labels = std::vector<std::string>(),
    std::vector<std::string> row_labels = std::vector<std::string>());
  void Close();

private:
  std::ofstream outfile_;
  std::string filename_;
};

#endif
