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
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_opt.hpp"

using namespace std;
using namespace PlotOptTypes;


int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  //// User defined colors
  Palette colors("txt/colors.txt", "default");

  //// Typically you would only have one PlotMaker so that you do not run on the same ntuples
  //// several times, but separated them here to make the examples modular and more clear
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// 1D plot mmiss sig/norm/D** ///////////////////////////////////////

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::data)
    .Bottom(BottomType::pull)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::none);
  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log).Title(TitleType::preliminary).Bottom(BottomType::ratio).Overflow(OverflowType::both);
  PlotOpt lin_shapes = lin_lumi().Stack(StackType::shapes).Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes);
  
  vector<PlotOpt> plottypes = {lin_lumi};

   //////// NamedFuncs to select signal, normalization, and D**
  NamedFunc is_dsptau("is_dsptau",[&](const Baby &b){
      return (abs(b.mu_MC_MOTHER_ID())==15 && abs(b.d0_MC_MOTHER_ID())==413 && abs(b.d0_MC_GD_MOTHER_ID())==511);
  });
  NamedFunc is_dspmu("is_dspmu",[&](const Baby &b){
      return (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==413 && abs(b.d0_MC_GD_MOTHER_ID())==511);
  });
  NamedFunc is_dss("is_dss",[&](const Baby &b){
    return ((abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==10411 && abs(b.d0_MC_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==10413 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==20413 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && ((abs(b.d0_MC_MOTHER_ID())==415 && abs(b.d0_MC_GD_MOTHER_ID())==511) |
                                              (abs(b.d0_MC_GD_MOTHER_ID())==415 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511))));
  });
  //////// Signal, normalization, D** processes for mmiss2 plot
  string repofolder = "ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  vector<shared_ptr<Process> > procs_mm;
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("Data (actually MC)",Process::Type::data, colors("data"),
                                                      set<string>({repofolder+run2bare}), "1"));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D^{*+} #tau #nu", Process::Type::signal, colors("green"),
                                                      set<string>({repofolder+run2bare}), is_dsptau));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D^{*+} #mu #nu", Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), is_dspmu));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D^{**} #mu #nu", Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), is_dss));


  PlotMaker pm_mm;
  // Missing mass (in GeV^2, so need to divide by 1e6)
  double iniq2 = -2, endq2 = 10;
  int nplots = 1000;
  for(int ind=0; ind<nplots; ind++){
    double q2cut = iniq2 + ind*(endq2-iniq2);
    string q2cut_s ="FitVar_q2/1000000>" + to_string(q2cut);
    pm_mm.Push<Hist1D>(Axis(75, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"),
                       q2cut_s, procs_mm, plottypes);
  }
  pm_mm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to  


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

