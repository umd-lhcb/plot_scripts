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
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "mu_ubdt_ok && (k_p < 200) && (pi_p < 200) && (mu_p < 100) && (iso_p1 < 200) && (iso_p2 < 200) && (iso_p3 < 200) && (nspdhits < 450) && is_iso";

  string folder = "ntuples/0.9.10-l0_weights/Dst_D0-mc-tracker_only-DDX/";
  string ntpBd = folder + "Dst--24_05_09--mc--12895400--2016--md--tracker_only.root";
  //string ntpBd = folder + "Dst--24_05_09--mc--11894610--2016--md--tracker_only.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body nominal", Process::Type::background, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body: -#sigma_{dalitz}^{lin} (wdal_lm)", Process::Type::background, colors("red"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable")); 
  procs.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body: -#sigma_{dalitz}^{quad} (wdal_qm)", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable")); 

  vector<shared_ptr<Process> > procs2;
  procs2.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procs2.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body nominal", Process::Type::background, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procs2.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body: +#sigma_{dalitz}^{lin} (wdal_lp)", Process::Type::background, colors("red"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable")); 
  procs2.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body: +#sigma_{dalitz}^{quad} (wdal_qp)", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable")); 

  vector<shared_ptr<Process> > procsLin;
  procsLin.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procsLin.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procsLin.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: +#sigma_{dalitz}^{lin} (wdal_lp)", Process::Type::data, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts)); 
  procsLin.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: -#sigma_{dalitz}^{lin} (wdal_lm)", Process::Type::data, colors("red"),
                   set<string>({ntpBd}), globalCuts)); 

  vector<shared_ptr<Process> > procsQuad;
  procsQuad.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procsQuad.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procsQuad.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: +#sigma_{dalitz}^{quad} (wdal_qp)", Process::Type::data, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts)); 
  procsQuad.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: -#sigma_{dalitz}^{quad} (wdal_qm)", Process::Type::data, colors("red"),
                   set<string>({ntpBd}), globalCuts)); 

  vector<shared_ptr<Process> > procs2body;
  procs2body.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procs2body.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procs2body.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: +50% 2-body", Process::Type::data, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts)); 
  procs2body.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: -50% 2-body", Process::Type::data, colors("red"),
                   set<string>({ntpBd}), globalCuts)); 

  vector<shared_ptr<Process> > procsKstar;
  procsKstar.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 3-body", Process::Type::background, colors("purple"),
                   set<string>({ntpBd}), globalCuts + "&& is_dal_variable"));
  procsKstar.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX 2-body", Process::Type::background, colors("yellow"),
                   set<string>({ntpBd}), globalCuts + "&& !is_dal_variable")); 
  procsKstar.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: +#sigma_{K*}", Process::Type::data, colors("darkblue"),
                   set<string>({ntpBd}), globalCuts)); 
  procsKstar.push_back(Process::MakeShared<Baby_run2_std>
                  ("DDX: -#sigma_{K*}", Process::Type::data, colors("red"),
                   set<string>({ntpBd}), globalCuts)); 


     
  string basew = "skim_global_ok*wff*wpid_ubdt*wtrg*wtrk*wbr_dd*w_missDDX*wjk";
  vector<NamedFunc> weights({basew, basew, basew+"*wdal_lm", basew+"*wdal_qm"});
  vector<NamedFunc> weights2({basew, basew, basew+"*wdal_lp", basew+"*wdal_qp"});
  vector<NamedFunc> weightsLin({basew, basew, basew+"*wdal_lp", basew+"*wdal_lm"});
  vector<NamedFunc> weightsQuad({basew, basew, basew+"*wdal_qp", basew+"*wdal_qm"});
  vector<NamedFunc> weights2body({basew, basew, basew+"*(1+0.5*!is_dal_variable)", basew+"*(1-0.5*!is_dal_variable)"});
  vector<NamedFunc> weightsKstar({basew, basew, basew+"*wkst_p", basew+"*wkst_m"});
  PlotMaker pm;

  // Individual contributions
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs, linplot, weights).Tag("linesMinus");
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs2, linplot, weights2).Tag("linesPlus");
  pm.Push<Hist1D>(Axis(50, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procs, linplot, weights).Tag("linesMinus");
  pm.Push<Hist1D>(Axis(50, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procs2, linplot, weights2).Tag("linesPlus");
  pm.Push<Hist1D>(Axis(43, -2.0, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procs, linplot, weights).Tag("linesMinus");
  pm.Push<Hist1D>(Axis(43, -2.0, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procs2, linplot, weights2).Tag("linesPlus");

  // Stacked plots
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procsLin, stackplot, weightsLin).Tag("dalLin");
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procsQuad, stackplot, weightsQuad).Tag("dalQuad");
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procs2body, stackplot, weights2body).Tag("2body");
  pm.Push<Hist1D>(Axis(100,-0.4,12.6, "q2", "q^{2} [GeV^{2}]"), "1", procsKstar, stackplot, weightsKstar).Tag("Kstar");

  pm.Push<Hist1D>(Axis(34, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procsLin, stackplot, weightsLin).Tag("dalLin");
  pm.Push<Hist1D>(Axis(34, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procsQuad, stackplot, weightsQuad).Tag("dalQuad");
  pm.Push<Hist1D>(Axis(34, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procs2body, stackplot, weights2body).Tag("2body");
  pm.Push<Hist1D>(Axis(34, 100,2650, "1000*el", "E^{*}_{#mu} [MeV]"), "q2>10.25", procsKstar, stackplot, weightsKstar).Tag("Kstar");
  
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procsLin, stackplot, weightsLin).Tag("dalLin");
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procsQuad, stackplot, weightsQuad).Tag("dalQuad");
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procs2body, stackplot, weights2body).Tag("2body");
  pm.Push<Hist1D>(Axis(43, -2, 10.9, "mm2", "m^{2}_{miss} [GeV^{2}]"), "1", procsKstar, stackplot, weightsKstar).Tag("Kstar");


  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
