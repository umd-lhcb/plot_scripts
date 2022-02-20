#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

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

  ////////////////////////////////////////// Plot settings //////////////////////////////////////////
  PlotOpt log_lumi_ratio_info("txt/plot_styles.txt", "LHCbPaper");
  log_lumi_ratio_info.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    //.LegendColumns(3)
    .Overflow(OverflowType::both);
  PlotOpt lin_lumi_ratio_info = log_lumi_ratio_info().YAxis(YAxisType::linear);
  
  vector<PlotOpt> linplot = {lin_lumi_ratio_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string gcuts = "d_mass_window_ok && b0_m > 5280";
  string repofolder = "ntuples_rdx/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_step2>("RS USB data", Process::Type::data, colors("data"),
                                                  set<string>({repofolder+"0.9.5-bugfix/Dst_D0-std/Dst_us--22_02_14--std--data--2011--md.root"}), gcuts));
  procs.push_back(Process::MakeShared<Baby_step2>("WS mu USB data", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"0.9.5-bugfix/Dst_D0-std/Dst_ws_Mu_us--22_02_14--std--data--2011--md.root"}), gcuts));


  ////////////////////////////////////////// Making plots //////////////////////////////////////////
  PlotMaker pm;
  // pm.Push<Hist1D>(Axis(4, -0.4, 12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  // pm.Push<Hist1D>(Axis(34, -100, 2650, "el*1000", "E_{#mu}^{*} [MeV]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  // pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");

  pm.Push<Hist1D>(Axis(20, -25, 15, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  pm.Push<Hist1D>(Axis(15, 0, 6000, "el*1000", "E_{#mu}^{*} [MeV]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  pm.Push<Hist1D>(Axis(20, -20, 20, "mm2", "m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");

  pm.Push<Hist1D>(Axis(20, 5280, 10280, "b0_m", "m(D^{*}#mu) [MeV]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  
  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

