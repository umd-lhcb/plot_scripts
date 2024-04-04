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
    .Overflow(OverflowType::none).Bottom(BottomType::ratio);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Bottom(BottomType::ratio);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear).Bottom(BottomType::ratio);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  
  vector<PlotOpt> linplot = {lin_shapes};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "mu_ubdt_ok && (k_p < 200) && (pi_p < 200) && (mu_p < 100) && (iso_p1 < 200) && (iso_p2 < 200) && (iso_p3 < 200) && (nspdhits < 450) && is_iso";

  string repofolder = "ntuples/";
  string datantp = "ntuples/0.9.6-2016_production/Dst_D0-std/Dst--23_11_06--std--data--2016--mu.root";
  string mcntp = "ntuples/0.9.9-2016_production/Dst_D0-mc-tracker_only-sig_norm/Dst--24_03_27--mc--11574021--2016--mu--tracker_only.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_std>("Data 2016", Process::Type::background, colors("data"),
                                                     set<string>({datantp}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("D*^{+}#mu#nu MC all w", Process::Type::data, colors("blue"),
                                                     set<string>({mcntp}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("D*^{+}#mu#nu MC all w but w_{J/#psi K}", Process::Type::data, colors("green"),
                                                     set<string>({mcntp}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("D*^{+}#mu#nu MC all w but w_{PID}", Process::Type::data, colors("red"),
                                                     set<string>({mcntp}), globalCuts));


     
  NamedFunc logbpt("logbpt", [&](const Baby &b){
    return log(b.b_pt()*1000);
  });
  NamedFunc woverjk("woverjk", [&](const Baby &b){
    if(b.wjk() != 0) return b.wiso()/b.wjk();
    else return 0.;
  });

  string basew = "wskim_iso*skim_global_ok*wff*wtrg*wtrk*wbr_dd*w_missDDX";
  vector<NamedFunc> weights({"1", basew+"*wjk*wpid_ubdt", basew + "*wpid_ubdt", basew + "*wjk"});
  //vector<NamedFunc> weights({"wskim_iso*skim_global_ok*wff*wpid_ubdt*wtrg*wtrk*wbr_dd*w_missDDX*wjk", "1"});
  //vector<NamedFunc> weights({woverjk, "1"});
  PlotMaker pm;

  pm.Push<Hist1D>(Axis(50,0,30, "b_pt", "p_{T}(B^{0}) [GeV]"), "mm2<0.5", procs, linplot, weights).RatioTitle("MC", "Data").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  pm.Push<Hist1D>(Axis(50,1.7,5, "b_eta", "#eta(B^{0})"), "mm2<0.5", procs, linplot, weights).RatioTitle("MC", "Data").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  pm.Push<Hist1D>(Axis(50,0,13, "mu_pt", "p_{T}(#mu^{+}) [GeV]"), "mm2<0.5", procs, linplot, weights).RatioTitle("MC", "Data").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  pm.Push<Hist1D>(Axis(50,1.7,5, "mu_eta", "#eta(#mu^{+})"), "mm2<0.5", procs, linplot, weights).RatioTitle("MC", "Data").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");

  //pm.Push<Hist1D>(Axis(100,5,12.5, logbpt, "log(B^{0} p_{T} [MeV])"), "mm2<0.5", procs, linplot, weights).RatioTitle("MC", "Data").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  //pm.Push<Hist1D>(Axis(100,1.7,6.7, "b_eta", "#eta(B^{0})"), "mm2<0.5", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  //pm.Push<Hist1D>(Axis(100,0,20, "mu_pt", "p_{T}(#mu^{+}) [GeV]"), "mm2<0.5", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  //pm.Push<Hist1D>(Axis(100,0,20, "mu_pt", "p_{T}(#mu^{+}) [GeV]"), "1", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");

  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
