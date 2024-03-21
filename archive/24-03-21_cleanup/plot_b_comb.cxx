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
    .FileExtensions({"pdf","png"})
    .Overflow(OverflowType::none);
  PlotOpt lin_lumi_ratio_info = log_lumi_ratio_info().YAxis(YAxisType::linear);
  PlotOpt lin_lumishapes_ratio = lin_lumi_ratio_info().YAxis(YAxisType::linear).Stack(StackType::lumi_shapes)
    .Title(TitleType::data);
  PlotOpt lin_shapes_ratio = lin_lumi_ratio_info().YAxis(YAxisType::linear).Stack(StackType::shapes)
    .Title(TitleType::data);
  
  vector<PlotOpt> linplot = {lin_shapes_ratio, lin_lumishapes_ratio};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string usbCuts = "dx_m_ok && b_m > 5400", sigCuts = "dx_m_ok && b_m < 5280";
  string repofolder = "ntuples_rdx/";
  vector<shared_ptr<Process> > procs, procs2, procs3, procs_sig;
  // procs.push_back(Process::MakeShared<Baby_step2>("RS USB", Process::Type::data, colors("data"),
  //                                                 set<string>({repofolder+"0.9.5-bugfix/Dst_D0-std/Dst_us--22_02_14--std--data--2011--md.root"}), usbCuts));
  // procs.push_back(Process::MakeShared<Baby_step2>("WS mu USB", Process::Type::background, colors("run1"),
  //                                               set<string>({repofolder+"0.9.5-bugfix/Dst_D0-std/Dst_ws_Mu_us--22_02_14--std--data--2011--md.root"}), usbCuts));
  procs.push_back(Process::MakeShared<Baby_step2>("WS mu USB #times lin. corr.", Process::Type::data, colors("blue"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));
  procs.push_back(Process::MakeShared<Baby_step2>("WS mu USB #times quad. corr.", Process::Type::data, colors("green"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));
  procs.push_back(Process::MakeShared<Baby_step2>("WS mu USB", Process::Type::data, colors("red"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));
  procs.push_back(Process::MakeShared<Baby_step2>("RS USB", Process::Type::background, colors("data"),
                                                  set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0--*"}), usbCuts));

  procs2.push_back(Process::MakeShared<Baby_step2>("RS USB", Process::Type::background, colors("data"),
                                                  set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0--*"}), usbCuts));
  procs2.push_back(Process::MakeShared<Baby_step2>("WS mu USB", Process::Type::data, colors("red"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));

  procs3.push_back(Process::MakeShared<Baby_step2>("WS mu USB #times lin. corr.", Process::Type::data, colors("blue"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));
  procs3.push_back(Process::MakeShared<Baby_step2>("WS mu USB", Process::Type::data, colors("red"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), usbCuts));
  procs3.push_back(Process::MakeShared<Baby_step2>("RS USB", Process::Type::background, colors("data"),
                                                  set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0--*"}), usbCuts));

  procs_sig.push_back(Process::MakeShared<Baby_step2>("WS mu sig. region #times lin. corr.", Process::Type::data, colors("blue"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), sigCuts));
  procs_sig.push_back(Process::MakeShared<Baby_step2>("WS mu sig. region #times quad. corr.", Process::Type::data, colors("green"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), sigCuts));
  procs_sig.push_back(Process::MakeShared<Baby_step2>("WS mu sig. region", Process::Type::background, colors("red"),
                                                set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/D0_ws_Mu*"}), sigCuts));

  //////// NamedFunc to apply quadratic correction found in rdx-run2-analysis
  NamedFunc quadCorr("quadCorr",[&](const Baby &b){
    return b.b_m()*b.b_m()*(2.04432e-08) + b.b_m()*(-0.000356905) + 2.14111;
  });

  //////// NamedFunc to apply linear correction found in rdx-run2-analysis
  NamedFunc linCorr("linCorr",[&](const Baby &b){
    return b.b_m()*(-5.51244e-05) + 1.05881;
  });

  ////////////////////////////////////////// Making plots //////////////////////////////////////////
  PlotMaker pm;
  // pm.Push<Hist1D>(Axis(4, -0.4, 12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  // pm.Push<Hist1D>(Axis(34, -100, 2650, "el*1000", "E_{#mu}^{*} [MeV]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  // pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("WS mu", "RS");
  
  vector<NamedFunc> weights({linCorr, quadCorr, "1", "1"}), weights3({linCorr, "1", "1"});
  
  pm.Push<Hist1D>(Axis(25, -37, 13, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot,weights).RatioTitle("WS mu", "RS").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(17, 0, 6800, "el*1000", "E_{#mu}^{*} [MeV]"), "1", procs, linplot,weights).RatioTitle("WS mu", "RS").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(20, -20, 20, "mm2", "m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot,weights).RatioTitle("WS mu", "RS").TopRight("Full 2016 data");

  pm.Push<Hist1D>(Axis(20, -0.4, 12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs_sig, linplot,weights).RatioTitle("WS mu", "RS").Tag("Sig").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(34, -100, 2650, "el*1000", "E_{#mu}^{*} [MeV]"), "1", procs_sig, linplot,weights).RatioTitle("WS mu", "RS").Tag("Sig").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m_{miss}^{2} [GeV^{2}]"), "1", procs_sig, linplot,weights).RatioTitle("WS mu", "RS").Tag("Sig").TopRight("Full 2016 data");

  linplot[0] = linplot[0].RatioMinimum(0.8).RatioMaximum(1.39);
  linplot[1] = linplot[1].RatioMinimum(0.8).RatioMaximum(1.19);
  pm.Push<Hist1D>(Axis(23, 5400, 10000, "b_m", "m(D^{0}#mu) [MeV]"), "1", procs2, linplot).RatioTitle("WS mu", "RS").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(23, 5400, 10000, "b_m", "m(D^{0}#mu) [MeV]"), "1", procs, linplot,weights).RatioTitle("WS mu", "RS").Tag("Corrs").TopRight("Full 2016 data");
  pm.Push<Hist1D>(Axis(23, 5400, 10000, "b_m", "m(D^{0}#mu) [MeV]"), "1", procs3, linplot,weights3).RatioTitle("WS mu", "RS").Tag("Corrs3").TopRight("Full 2016 data");

  
  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

