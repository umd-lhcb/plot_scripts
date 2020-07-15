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

  time_t begtime, endtime;
  time(&begtime);

  // Plot types
  PlotOpt log_lumi("txt/plot_styles.txt", "LHCbPaper");
  log_lumi.Title(TitleType::simulation)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  
  vector<PlotOpt> lumiplot = {lin_lumi, lin_shapes, log_lumi};
  vector<PlotOpt> shapeplot = {lin_shapes, log_shapes};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string repofolder = "ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  PlotMaker pm;


  /////////////////////////// Cut on muon angle ////////////////////////////////
  NamedFunc mu_pxpz("mu_pxpz",
                          [&](const Baby &b){
                            return (b.mu_PX()/b.mu_PZ());
  });
  NamedFunc mu_pxpz_tru("mu_pxpz_tru",
                          [&](const Baby &b){
                            return (b.mu_TRUEP_X()/b.mu_TRUEP_Z());
  });
  NamedFunc mu_pypz("mu_pypz",
                          [&](const Baby &b){
                            return fabs(b.mu_PY()/b.mu_PZ());
  });
  NamedFunc mu_pypz_tru("mu_pypz_tru",
                          [&](const Baby &b){
                            return fabs(b.mu_TRUEP_Y()/b.mu_TRUEP_Z());
  });

  vector<shared_ptr<Process> > procs_comp_pxpz_mu;
  procs_comp_pxpz_mu.push_back(Process::MakeShared<Baby_run2_bare>("p^{reco}(#mu) > 100 GeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), "mu_P/1000>100"));
  procs_comp_pxpz_mu.push_back(Process::MakeShared<Baby_run2_bare>("10 < p^{reco}(#mu) < 100 GeV",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), "mu_P/1000>10&&mu_P/1000<100"));
  procs_comp_pxpz_mu.push_back(Process::MakeShared<Baby_run2_bare>("5 < p^{reco}(#mu) < 10 GeV",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), "mu_P/1000>5&&mu_P/1000<10"));
  procs_comp_pxpz_mu.push_back(Process::MakeShared<Baby_run2_bare>("p^{reco}(#mu) < 5 GeV",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "mu_P/1000<5"));

  
  
  pm.Push<Hist1D>(Axis(400, -0.2, 0.2, mu_pxpz, "p_{x}^{reco}(#mu)/p_{z}^{reco}(#mu)", {-0.378, 0.378}), "1",
                  procs_comp_pxpz_mu, lumiplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(400, -0.2, 0.2, mu_pxpz_tru, "p_{x}^{true}(#mu)/p_{z}^{true}(#mu)", {-0.378, 0.378}), "1",
                  procs_comp_pxpz_mu, lumiplot).TopRight("13 TeV");


  ////////////////////////////////////// 2D scatter plot mu angles ////////////////////////////////////////

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> scattertype = {style().Stack(StackType::signal_on_top).Title(TitleType::data)};
  //////// Processes for scatter plot
  auto mup_high = Process::MakeShared<Baby_run2_bare>("p^{reco}(#mu) > 100 GeV", Process::Type::data, kBlack,
                                                      set<string>({repofolder+run2bare}), "mu_P>100000");
  auto mupt_high = Process::MakeShared<Baby_run2_bare>("p^{reco}_{T}(#mu) > 8 GeV", Process::Type::data, kBlue,
                                                       set<string>({repofolder+run2bare}), "mu_PT>8000");
  auto all_mu = Process::MakeShared<Baby_run2_bare>("MC", Process::Type::background, kBlack,
                                                    set<string>({repofolder+run2bare}), "1");
  mup_high->SetMarkerStyle(20); mup_high->SetMarkerSize(0.4);
  mupt_high->SetMarkerStyle(21);mupt_high->SetMarkerSize(0.4);
  
  vector<shared_ptr<Process> > procs_mu = {mup_high, mupt_high, all_mu};
  pm.Push<Hist2D>(Axis(55, -0.6, 0.5, "mu_TRUEP_X/mu_TRUEP_Z", "p_{x}^{true}/p_{z}^{true}(#mu)", {-0.38, 0.38}),
                     Axis(45, -0.45, 0.45, "mu_TRUEP_Y/mu_TRUEP_Z", "p_{y}^{true}/p_{z}^{true}(#mu)", {-0.28, 0.28}),
                     "1", procs_mu, scattertype).TopRight("");
  pm.Push<Hist2D>(Axis(80, -100, 100, "mu_TRUEP_X/mu_TRUEP_Z*1000", "p_{x}^{true}/p_{z}^{true}(#mu) [mrad]", {-10, 10}),
                  Axis(80, -100, 100, "mu_TRUEP_Y/mu_TRUEP_Z*1000", "p_{y}^{true}/p_{z}^{true}(#mu) [mrad]", {-10, 10}),
                  "1", procs_mu, scattertype).TopRight("");
  pm.Push<Hist2D>(Axis(80, -100, 100, "mu_PX/mu_PZ*1000", "p_{x}^{reco}/p_{z}^{reco}(#mu) [mrad]", {-10, 10}),
                  Axis(80, -100, 100, "mu_PY/mu_PZ*1000", "p_{y}^{reco}/p_{z}^{reco}(#mu) [mrad]", {-10, 10}),
                  "1", procs_mu, scattertype).TopRight("");


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
