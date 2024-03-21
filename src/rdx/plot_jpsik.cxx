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

  time_t begtime, endtime;
  time(&begtime);

  // Plot types
  PlotOpt log_lumi("txt/plot_styles.txt", "LHCbPaper");
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none).Bottom(BottomType::ratio);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear).Bottom(BottomType::pull);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear).Bottom(BottomType::ratio);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  
  vector<PlotOpt> linplot = {lin_lumi};
  vector<PlotOpt> shapeplot = {lin_shapes};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "b_ownpv_ndof<=250 && ntracks<=450";
  string datafull = "../lhcb-ntuples-gen/ntuples/0.9.6-2016_production/JpsiK-std-step2/JpsiK--22_02_26--std--data--2016-*";
  string datasw = "../lhcb-ntuples-gen/run2-JpsiK/fit/fit_results/JpsiK-22_02_26_23_52-std-fit-2016/fit.root";
  string mcfull = "../lhcb-ntuples-gen/ntuples/0.9.7-rdx_production/JpsiK-mc-step2/*";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC after w", Process::Type::background, colors("green"),
                                                       set<string>({mcfull}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC before w", Process::Type::background, colors("blue"),
                                                       set<string>({mcfull}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_jpsiksw>("J/#psi K data after sWeights", Process::Type::data, colors("data"),
                                                         set<string>({datasw}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K data before sWeights", Process::Type::data, colors("red"),
                                                       set<string>({datafull}), globalCuts));
 
  vector<shared_ptr<Process> > procs2;
  procs2.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC after L0", Process::Type::background, colors("green"),
                                                        set<string>({mcfull}), globalCuts+"&&l0"));
  procs2.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC before L0", Process::Type::background, colors("blue"),
                                                        set<string>({mcfull}), globalCuts));
  procs2.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K data after L0", Process::Type::data, colors("data"),
                                                        set<string>({datafull}), globalCuts+"&&l0"));
  procs2.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K data before L0", Process::Type::data, colors("red"),
                                                        set<string>({datafull}), globalCuts));
 
  vector<shared_ptr<Process> > procs3;
  procs3.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC (#mu p_{T} > 2 GeV)", Process::Type::background, colors("green"),
                                                        set<string>({mcfull}), globalCuts+"&&mu_pt>2000"));
  procs3.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K MC", Process::Type::background, colors("blue"),
                                                        set<string>({mcfull}), globalCuts));
  procs3.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K (#mu p_{T} > 2 GeV)", Process::Type::data, colors("data"),
                                                        set<string>({datafull}), globalCuts+"&&mu_pt>2000"));
  procs3.push_back(Process::MakeShared<Baby_run2_jpsik>("J/#psi K data", Process::Type::data, colors("red"),
                                                        set<string>({datafull}), globalCuts));
 


  vector<NamedFunc> weights({"wpid*wtrk*wjk_occ*wjk_kin", "wpid*wtrk*wjk_occ", "sw_sig", "1"});
  vector<NamedFunc> weights2({"wpid*wtrk*wjk_occ*wjk_kin", "wpid*wtrk*wjk_occ*wjk_kin", "1", "1"});
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(100,1.7,6.7, "b_eta", "#eta(B)"), "1", procs, shapeplot, weights).RatioTitle("Process", "MC after w").Tag("bins100");
  pm.Push<Hist1D>(Axis(9,2,6, "b_eta", "#eta(B)"), "1", procs, shapeplot, weights).RatioTitle("Process", "MC after w").Tag("bins9");
  pm.Push<Hist1D>(Axis(100,1.7,6.7, "b_eta", "#eta(B)"), "1", procs2, shapeplot, weights2).RatioTitle("Process", "MC after L0").Tag("L0");
  pm.Push<Hist1D>(Axis(100,1.7,6.7, "b_eta", "#eta(B)"), "1", procs3, shapeplot, weights2).RatioTitle("Process", "MC after p_{T}").Tag("mupt");

  pm.min_print_ = true;
  pm.MakePlots(1);
  
  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"lumi", required_argument, 0, 'l'},    // Luminosity to normalize MC with (no mc)
      {"example", required_argument, 0, 's'},    // Which example to use: standard, met150, 2015 mc
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
