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
  //vector<PlotOpt> stackplot = {lin_lumi_norm};
  vector<PlotOpt> stackplot = {lin_lumi};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "(k_p < 200) && (pi_p < 200) && (mu_p < 100) && (iso_p1 < 200) && (iso_p2 < 200) && (iso_p3 < 200) && (nspdhits < 450)";

  string folder = "ntuples/0.9.10-l0_weights/";
  string ntpDDX = folder + "Dst_D0-mc-tracker_only-DDX/Dst--24_05_09--mc--11894610--2016--md--tracker_only.root";
  string ntpDstMu = folder + "Dst_D0-mc-tracker_only-sig_norm/Dst--24_05_10--mc--11574021--2016--md--tracker_only.root";
  string ntpData = "ntuples/0.9.6-2016_production/Dst_D0-std/Dst--23_11_06--std--data--2016--m*.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("data", Process::Type::background, colors("data"),
                   set<string>({ntpData}), globalCuts + "&& is_dd"));
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX", Process::Type::background, colors("green"),set<string>({ntpDDX}), globalCuts)); 
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("D*#mu", Process::Type::background, colors("blue"),set<string>({ntpDstMu}), globalCuts)); 
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("data, #muBDT > 0.25", Process::Type::background, colors("data"),
                   set<string>({ntpData}), globalCuts + "&& is_dd && mu_ubdt_ok",2));
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX, wpid_ubdt", Process::Type::background, colors("green"),
                   set<string>({ntpDDX}), globalCuts, 2)); 
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("D*#mu, wpid_ubdt", Process::Type::background, colors("blue"),
                   set<string>({ntpDstMu}), globalCuts, 2)); 

     
  string basew = "skim_global_ok*wff*wpid_ubdt*wtrg*wtrk*wbr_dd*w_missDDX*wjk";
  vector<NamedFunc> weights({"1", "wdd/wpid_ubdt", "wdd/wpid_ubdt", "1", "wdd", "wdd"});
  PlotMaker pm;

  pm.Push<Hist1D>(Axis(30,0,30, "mu_p", "p_{#mu} [GeV]"), "1", procs, linplot, weights).Tag("uBDT");
  pm.Push<Hist1D>(Axis(30,0,3, "mu_pt", "p_{T,#mu} [GeV]"), "1", procs, linplot, weights).Tag("uBDT");

  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
