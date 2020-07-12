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

void GetOptions(int argc, char *argv[]);


namespace{
  float lumi = 4.3;
  string example = "search";
}


int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining plot types //////////////////////////////////////////
  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::data)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::none);
  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log).Title(TitleType::simulation).Bottom(BottomType::pull);
  PlotOpt lin_shapes = lin_lumi().Stack(StackType::shapes).Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes);
  
  vector<PlotOpt> plottypes = {lin_lumi, log_lumi, lin_shapes, lin_lumi_shapes};

  Palette colors("txt/colors.txt", "default");


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

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

  //////// Processes
  string repofolder = "ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_bare>("Data (actually MC)",Process::Type::data, colors("data"),
                                                      set<string>({repofolder+run2bare}), "1"));
  procs.push_back(Process::MakeShared<Baby_run2_bare>("Signal", Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), is_dsptau));
  procs.push_back(Process::MakeShared<Baby_run2_bare>("Normalization", Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), is_dspmu));
  procs.push_back(Process::MakeShared<Baby_run2_bare>("D**", Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), is_dss));


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining histograms //////////////////////////////////////////
  

  PlotMaker pm;

  // Missing mass (in GeV^2, so need to divide by 1e6)
  pm.Push<Hist1D>(Axis(75, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>8",
                  procs, plottypes);
  
  pm.min_print_ = true;
  pm.MakePlots(1);

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no data)
      {"example", required_argument, 0, 's'},    // Which example to use: standard, met150, 2015 data
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:l:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'l':
      lumi = atof(optarg);
      break;
    case 's':
      example = optarg;
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
