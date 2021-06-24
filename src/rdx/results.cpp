#include "rdx/results.hpp"

using namespace std;

Results::Results(TString iname, float ivalue, vector<float> istatError, vector<float> isystError,
		 int icolor, float icorrel):
  name(iname),
  value(ivalue),
  statError(istatError),
  systError(isystError),
  color(icolor),
  correl(icorrel){
  if(statError.size()==1) statError.push_back(statError[0]);
  if(systError.size()==1) systError.push_back(systError[0]);
  }
  
PadResults::PadResults(TString ititle, float iminX, float imaxX, std::vector<Results> ivresults, 
		       Results ism, Results iaverage):
  title(ititle),
  minX(iminX),
  maxX(imaxX),
  vresults(ivresults),
  sm(ism),
  average(iaverage){

  }

