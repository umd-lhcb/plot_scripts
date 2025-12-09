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
  log_lumi.Title(TitleType::info)
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

  string ntpD0 = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuD0-12675011-MagDown/D0--25_12_02--mc--12675011--2016--md--tracker_only.root";
  //string ntpDst = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuDst-11676012-MagDown/Dst--25_12_02--mc--11676012--2016--md--tracker_only.root";
  string ntpDst = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuDst-11676012-MagDown/Dst--25_12_02--mc--11676012--2016--md--tracker_only.root";
  string ntpDst_D0 = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuDst-11676012-MagDown/D0--25_12_02--mc--11676012--2016--md--tracker_only.root";
  string ntpDst0 = "ntuples/0.9.17-all_years/2016/Dstst_heavy/DststHMuDst0-12875440-MagDown/D0--25_12_02--mc--12875440--2016--md--tracker_only.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{0}", Process::Type::background, colors("purple"),
                   set<string>({ntpD0}), globalCuts + " && truthmatch%100 == 10")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{0}", Process::Type::background, colors("red"),
                   set<string>({ntpD0}), globalCuts + " && truthmatch%100 == 20")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{0}", Process::Type::background, colors("green"),
                   set<string>({ntpD0}), globalCuts + " && truthmatch%100 == 30")); 
  procs.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{0}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpD0}), globalCuts + " && truthmatch%100 == 40"));
  
  vector<shared_ptr<Process> > procsDst;
  procsDst.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{*+}", Process::Type::background, colors("purple"),
                   set<string>({ntpDst}), globalCuts + " && truthmatch%100 == 10")); 
  procsDst.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{*+}", Process::Type::background, colors("red"),
                   set<string>({ntpDst}), globalCuts + " && truthmatch%100 == 20")); 
  procsDst.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{*+}", Process::Type::background, colors("green"),
                   set<string>({ntpDst}), globalCuts + " && truthmatch%100 == 30")); 
  procsDst.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{*+}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpDst}), globalCuts + " && truthmatch%100 == 40")); 

  vector<shared_ptr<Process> > procsDst_D0;
  procsDst_D0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{*+}", Process::Type::background, colors("purple"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 10")); 
  procsDst_D0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{*+}", Process::Type::background, colors("red"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 20")); 
  procsDst_D0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{*+}", Process::Type::background, colors("green"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 30")); 
  procsDst_D0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{*+}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 40")); 

  vector<shared_ptr<Process> > procsDst0;
  procsDst0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{*0}", Process::Type::background, colors("purple"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 10")); 
  procsDst0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{*0}", Process::Type::background, colors("red"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 20")); 
  procsDst0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{*0}", Process::Type::background, colors("green"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 30")); 
  procsDst0.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{*0}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 40")); 

  vector<shared_ptr<Process> > procsDsts;
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{*0}", Process::Type::background, colors("purple"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 10")); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{*0}", Process::Type::background, colors("red"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 20")); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{*0}", Process::Type::background, colors("green"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 30")); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{*0}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpDst0}), globalCuts + " && truthmatch%100 == 40")); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2580) #rightarrow D^{*+}", Process::Type::background, colors("purple"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 10",2)); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2640) #rightarrow D^{*+}", Process::Type::background, colors("red"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 20",2)); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2737) #rightarrow D^{*+}", Process::Type::background, colors("green"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 30",2)); 
  procsDsts.push_back(Process::MakeShared<Baby_rdx917>
                  ("D**(2978) #rightarrow D^{*+}", Process::Type::background, colors("darkblue"),
                   set<string>({ntpDst_D0}), globalCuts + " && truthmatch%100 == 40",2)); 

     
  PlotMaker pm;

  // 
  //pm.Push<Hist1D>(Axis(50,2000,3000, "d_iso1_m", "m(D^{0}#pi) [GeV]"), "wskim_1os", procs, linplot);
  pm.Push<Hist1D>(Axis(50,2100,3500, "d_iso1_iso2_m", "m(D^{0}#pi#pi) [GeV]"), "wskim_2os", procs, linplot).Tag("D0");
  pm.Push<Hist1D>(Axis(50,2100,3500, "d_iso1_iso2_m", "m(D^{*+}#pi#pi) [GeV]"), "wskim_2os", procsDst, linplot).Tag("Dst");
  pm.Push<Hist1D>(Axis(50,2100,3500, "d_iso1_iso2_m", "m(D^{0}#pi#pi) [GeV]"), "wskim_2os", procsDst0, linplot).Tag("Dst0");
  pm.Push<Hist1D>(Axis(50,2100,3500, "d_iso1_iso2_m", "m(D^{0}#pi#pi) [GeV]"), "wskim_2os", procsDst_D0, linplot).Tag("Dst_D0");
  pm.Push<Hist1D>(Axis(40,2100,3500, "d_iso1_iso2_m", "m(D^{0}#pi#pi) [GeV]"), "wskim_2os", procsDsts, linplot).Tag("Dsts");


  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
