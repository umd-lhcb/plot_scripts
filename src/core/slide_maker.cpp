#include "core/slide_maker.hpp"

#include "core/utilities.hpp"

using namespace std;

SlideMaker::SlideMaker(string fname, string aspect_ratio):
  filename_(fname){
  outfile_.open("slides/"+filename_, ofstream::out);
  outfile_<<"\\documentclass[8pt,aspectratio="+aspect_ratio+"]{beamer}\n";
  outfile_<<"\\usepackage{graphicx}\n\n";
  outfile_<<"\\setbeamertemplate{footline}[frame number]{}\n";
  outfile_<<"\\setbeamertemplate{navigation symbols}{}\n";

  // outfile_<<"\\renewcommand{\\arraystretch}{1.2}\n";
  outfile_<<"\\setbeamertemplate{frametitle}[default][center]\n";
  outfile_<<"\\begin{document}\n";
}

void SlideMaker::AddSlide(vector<string> pnames, unsigned ncols, string title, 
                          vector<string> col_labels, vector<string> row_labels){
  //set # of columns and rows
  unsigned nplots = pnames.size();
  unsigned nrows = nplots/ncols;
  unsigned nperrow = nplots/nrows;
  if (nplots%ncols!=0) nrows +=1;

  string width = RoundNumber(0.9/ncols,2).Data();
  string height = RoundNumber(0.85/nrows,2).Data();
  if (title!="") height = RoundNumber(0.9/nrows,2).Data();
  outfile_<<"\\frame{\n";
  if (ncols<=2 || nrows>=3)
    outfile_<<"\\resizebox{\\textheight}{!}{\n";
  else
    outfile_<<"\\resizebox{\\textwidth}{!}{\n";
  outfile_<<"\\begin{tabular}{|c|";
  for (unsigned i(0); i<ncols; i++) outfile_<<"c";
  outfile_<<"|} \n";
  outfile_<<"\\hline \n";
  outfile_<<"\\multicolumn{"<<ncols+1<<"}{|c|}{"+title+"}\\\\ \n";
  outfile_<<"\\hline \n";
  for (unsigned i(0); i<ncols+1; i++) outfile_<<col_labels[i]<<(i==ncols ? " \\\\ \n":" & ");
  outfile_<<"\\hline \n";
  for (size_t i(0); i<nplots; i++){
    if (pnames[i]=="") continue;
    string line = "";
    if (i%unsigned(ncols)==0) line += row_labels[i/nperrow] + "&";
    line += "  \\includegraphics[trim={0.5cm 0.5cm 0.5cm 1cm},clip, height="+height+"\\textheight,width="
                  +width+"\\textwidth,keepaspectratio]{plots/"+pnames[i]+"}";
    if ((i+1)%unsigned(ncols)==0) line +="\\\\ \n \\hline \n";
    else line +="& \n";
    outfile_<<line;
  }
  outfile_<<"\\end{tabular} \n";
  outfile_<<"}}\n";
}

void SlideMaker::AddSlideWithReplace(string oldexp, string newexp, vector<string> pnames, 
                                     unsigned ncols, string title, 
                                     vector<string> col_labels, vector<string> row_labels){
  for (auto &iname: pnames) ReplaceAll(iname, oldexp, newexp);
  AddSlide(pnames, ncols, title, col_labels, row_labels);
}

void SlideMaker::Close(){
  outfile_<<"\\end{document}\n";
  outfile_.close();

  execute("pdflatex -output-directory=slides slides/"+filename_+"  > /dev/null");
  ReplaceAll(filename_, ".tex",".pdf");
  cout<<endl<<"open slides/"+filename_<<endl;
}
