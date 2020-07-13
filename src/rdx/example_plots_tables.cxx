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
#include "core/hist2d.hpp"
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
    .Bottom(BottomType::pull)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::none);
  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log).Title(TitleType::simulation).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = lin_lumi().Stack(StackType::shapes).Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes);
  
  vector<PlotOpt> plottypes = {lin_lumi, log_lumi, lin_shapes, lin_lumi_shapes};
  //vector<PlotOpt> plottypes = {lin_lumi};

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> scattertype = {style().Stack(StackType::signal_on_top).Title(TitleType::data)};

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

  //////// Signal, normalizatin, D** processes
  string repofolder = "ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  vector<shared_ptr<Process> > procs_mm;
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("Data (actually MC)",Process::Type::data, colors("data"),
                                                      set<string>({repofolder+run2bare}), "1"));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D*^{+} #tau #nu", Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), is_dsptau));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D*^{+} #mu #nu", Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), is_dspmu));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D** #mu #nu", Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), is_dss));

  //////// Processes for scatter plot
  auto data = Process::MakeShared<Baby_run2_bare>("p^{reco}(#mu) > 100 GeV", Process::Type::data, kBlack,
                                                  set<string>({repofolder+run2bare}), "mu_P>100000");
  auto data2 = Process::MakeShared<Baby_run2_bare>("p^{reco}_{T}(#mu) > 8 GeV", Process::Type::data, kBlue,
                                                  set<string>({repofolder+run2bare}), "mu_PT>8000");
  auto bkg = Process::MakeShared<Baby_run2_bare>("MC", Process::Type::background, kBlack,
                                                 set<string>({repofolder+run2bare}), "1");
  data->SetMarkerStyle(20); data->SetMarkerSize(0.4);
  data2->SetMarkerStyle(21);data2->SetMarkerSize(0.4);
  
  vector<shared_ptr<Process> > procs_mu = {data, data2, bkg};


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining histograms //////////////////////////////////////////
  

  PlotMaker pm;

  // Missing mass (in GeV^2, so need to divide by 1e6)
  pm.Push<Hist1D>(Axis(75, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>8",
                  procs_mm, plottypes);

  // Scatter plot of muon
  pm.Push<Hist2D>(Axis(55, -0.6, 0.5, "mu_TRUEP_X/mu_TRUEP_Z", "p_{x}^{true}/p_{z}^{true}(#mu)", {-0.38, 0.38}),
                  Axis(38, -0.38, 0.38, "mu_TRUEP_Y/mu_TRUEP_Z", "p_{y}^{true}/p_{z}^{true}(#mu)", {-0.28, 0.28}),
                  "1", procs_mu, scattertype).TopRight("");
 
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
