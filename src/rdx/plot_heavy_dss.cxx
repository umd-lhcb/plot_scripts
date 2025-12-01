#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting

#include "core/utilities.hpp"
#include "core/baby.hpp"
#include "core/process.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/palette.hpp"
#include "core/table.hpp"
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"

using namespace std;
using namespace PlotOptTypes;



int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);

  // Plot types
  PlotOpt log_lumi("txt/plot_styles.txt", "LHCbPaper");
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    //.Stack(StackType::signal_overlay)
    .Overflow(OverflowType::none).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_norm = log_lumi().YAxis(YAxisType::linear).Bottom(BottomType::ratio);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Bottom(BottomType::ratio).Stack(StackType::signal_overlay);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false).Bottom(BottomType::off);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear).Bottom(BottomType::off);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  
  vector<PlotOpt> linplot = {lin_shapes};
  vector<PlotOpt> stackplot = {lin_lumi_norm};
  //vector<PlotOpt> stackplot = {lin_lumi};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "mu_ubdt_ok && (k_p < 200) && (pi_p < 200) && (mu_p < 100) && (iso_p1 < 200) && (iso_p2 < 200) && (iso_p3 < 200) && (nspdhits < 450)";

  string ntuple = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuDst0-12875440-MagUp/D0--25_08_19--mc--12875440--2016--mu--tracker_only.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580)", Process::Type::background, colors("purple"),
                   set<string>({ntuple}), globalCuts + " && truthmatch%100 == 10")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640)", Process::Type::background, colors("red"),
                   set<string>({ntuple}), globalCuts + " && truthmatch%100 == 20")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737)", Process::Type::background, colors("green"),
                   set<string>({ntuple}), globalCuts + " && truthmatch%100 == 30")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978)", Process::Type::background, colors("darkblue"),
                   set<string>({ntuple}), globalCuts + " && truthmatch%100 == 40")); 

     
  PlotMaker pm;

  // 
  pm.Push<Hist1D>(Axis(100,2000,3000, "d_iso1_m", "m(D^{0}#pi) [GeV]"), "wskim_1os", procs, linplot);


  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
