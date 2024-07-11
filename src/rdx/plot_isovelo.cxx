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
  string newIso = "((iso_type1>2&&iso_bdt1<0.15) || (iso_type2>2&&iso_bdt2<0.15) || (iso_type3>2&&iso_bdt3<0.15))";
  
  string ntpData = "ntuples/0.9.6-2016_production/Dst_D0-std/Dst--23_11_06--std--data--2016--m*.root";
  string ntpDstMu = "ntuples/0.9.10-l0_weights/Dst_D0-mc-tracker_only-sig_norm/Dst--24_05_10--mc--11574021--2016--md--tracker_only.root";
  string ntpdDDX = "ntuples/0.9.10-l0_weights/Dst_D0-mc-tracker_only-DDX/Dst--24_05_09--mc--11894610--2016--md--tracker_only.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_std>("Data ISO", Process::Type::data, colors("data"),
                                                     set<string>({ntpData}), globalCuts + "&& is_iso"));
  procs.push_back(Process::MakeShared<Baby_run2_std>("Data ISO non-VELO", Process::Type::background, colors("red"),
                                                     set<string>({ntpData}), globalCuts + "&& !is_iso &&"+newIso, 2));
  // procs.push_back(Process::MakeShared<Baby_run2_std>("D*#mu MC ISO", Process::Type::data, colors("blue"),
  //                                                    set<string>({ntpDstMu}), globalCuts));
  // procs.push_back(Process::MakeShared<Baby_run2_std>("D*#mu MC ISO non-VELO", Process::Type::data, colors("purple"),
  //                                                    set<string>({ntpDstMu}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("DDX MC ISO", Process::Type::background, colors("green"),
                                                     set<string>({ntpdDDX}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("DDX MC ISO non-VELO", Process::Type::background, colors("purple"),
                                                     set<string>({ntpdDDX}), globalCuts,2));

     
  vector<NamedFunc> weights({"1", "1", "w*skim_global_ok*(iso_bdt1<0.15)", "w*skim_global_ok*(iso_bdt1>0.15 && ((iso_type1>2&&iso_bdt1<0.15) || (iso_type2>2&&iso_bdt2<0.15) || (iso_type3>2&&iso_bdt3<0.15)))"});
  PlotMaker pm;

  pm.Push<Hist1D>(Axis(30,0,30, "mu_p", "p_{#mu} [GeV]"), "1", procs, linplot, weights);
  pm.Push<Hist1D>(Axis(30,0,3, "mu_pt", "p_{T,#mu} [GeV]"), "1", procs, linplot, weights);
  pm.Push<Hist1D>(Axis(50,0,2, "spi_pt", "#pi_{soft} p_{T} [GeV]"), "1", procs, linplot, weights);
  pm.Push<Hist1D>(Axis(50,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot, weights);
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procs, linplot, weights);
  pm.Push<Hist1D>(Axis(34, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "1", procs, linplot, weights);

  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
