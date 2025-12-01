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
#include "core/hist1d.hpp"
#include "core/plot_opt.hpp"

using namespace std;
using namespace PlotOptTypes;

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  // Defining plot styles
  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::data).Bottom(BottomType::off).Stack(StackType::signal_on_top).Overflow(OverflowType::none);
  vector<PlotOpt> plottypes = {lin_lumi};
  Palette colors("txt/colors.txt", "default");

  // Defining processes (plot components)
  string ntuple = "various/ntuples/Dst--20_03_18--cutflow_mc--cocktail--2016--md.root";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2>("Signal", Process::Type::background, colors("blue"),
                                                 set<string>({ntuple}), "(muplus_MC_MOTHER_ID==15 || muplus_MC_MOTHER_ID==-15) && (D0_MC_MOTHER_ID==413 || D0_MC_MOTHER_ID==-413) && (D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==-511)"));
  procs.push_back(Process::MakeShared<Baby_run2>("Normalization", Process::Type::background, colors("green"),
                                                 set<string>({ntuple}), "(muplus_MC_MOTHER_ID==511 || muplus_MC_MOTHER_ID==-511) && (D0_MC_MOTHER_ID==413 || D0_MC_MOTHER_ID==-413) && (D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==-511)"));
  procs.push_back(Process::MakeShared<Baby_run2>("D**", Process::Type::background, colors("red"),
                                                 set<string>({ntuple}), "(muplus_MC_MOTHER_ID==511 && D0_MC_MOTHER_ID==10411 && D0_MC_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && D0_MC_GD_MOTHER_ID==10413 && D0_MC_GD_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && D0_MC_GD_MOTHER_ID==20413 && D0_MC_GD_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && (D0_MC_MOTHER_ID==415 && D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==415 && D0_MC_GD_GD_MOTHER_ID==511))"));

  // Making plots. Missing mass plot is set to GeV^2 by dividing by 1e6
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(25, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "1", procs, plottypes);
  pm.MakePlots(1); // The "1" is the luminosity scaling

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

