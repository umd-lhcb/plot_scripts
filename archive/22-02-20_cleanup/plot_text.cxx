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
    .Stack(StackType::data_norm).LegendColumns(3)
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
  
  vector<PlotOpt> linplot = {lin_shapes, lin_lumi};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string repofolder = "ntuples/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2>("Signal", Process::Type::background, colors("dsptau"),
                                                set<string>({repofolder+"pre-0.9.0/Dst-cutflow_mc/Dst--20_03_18--cutflow_mc--cocktail--2016--md.root"}),
                                                 "(muplus_MC_MOTHER_ID==15 || muplus_MC_MOTHER_ID==-15) && (D0_MC_MOTHER_ID==413 || D0_MC_MOTHER_ID==-413) && (D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==-511)"));
  procs.push_back(Process::MakeShared<Baby_run2>("Normalization", Process::Type::background, colors("dspmu"),
                                                set<string>({repofolder+"pre-0.9.0/Dst-cutflow_mc/Dst--20_03_18--cutflow_mc--cocktail--2016--md.root"}),
                                                 "(muplus_MC_MOTHER_ID==511 || muplus_MC_MOTHER_ID==-511) && (D0_MC_MOTHER_ID==413 || D0_MC_MOTHER_ID==-413) && (D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==-511)"));
  procs.push_back(Process::MakeShared<Baby_run2>("D**", Process::Type::background, colors("dss"),
                                                set<string>({repofolder+"pre-0.9.0/Dst-cutflow_mc/Dst--20_03_18--cutflow_mc--cocktail--2016--md.root"}),
                                                 "(muplus_MC_MOTHER_ID==511 && D0_MC_MOTHER_ID==10411 && D0_MC_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && D0_MC_GD_MOTHER_ID==10413 && D0_MC_GD_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && D0_MC_GD_MOTHER_ID==20413 && D0_MC_GD_GD_MOTHER_ID==511) || (muplus_MC_MOTHER_ID==511 && (D0_MC_MOTHER_ID==415 && D0_MC_GD_MOTHER_ID==511 || D0_MC_GD_MOTHER_ID==415 && D0_MC_GD_GD_MOTHER_ID==511))"));

  
  PlotMaker pm;

  pm.Push<Hist1D>(Axis(30, -3, 12,"FitVar_Mmiss2/1000000", " m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot);

  pm.min_print_ = true;
  pm.MakePlots(9);

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
