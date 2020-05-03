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
  
  vector<PlotOpt> linplot = {lin_shapes};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string repofolder = "/Users/manuelf/code/lhcb-ntuples-gen/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run1>("B #rightarrow D^{0}X#mu#nu 2011", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"run1-b2D0MuXB2DMuNuForTauMuLine/ntuples/cutflow-Dst/BCands_Dst-yipeng-cutflow_mc-2011-mag_down.root"}), "1"));
  procs.push_back(Process::MakeShared<Baby_run2>("B #rightarrow D^{0}X#mu#nu 2016", Process::Type::background, colors("run2"),
                                                set<string>({repofolder+"run2-b2D0MuXB2DMuForTauMuLine/ntuples/cutflow-Dst/BCands_Dst-yipeng-cutflow_mc-2016-mag_down.root"}), "1"));

  PlotMaker pm;

  pm.Push<Hist1D>(Axis(100, 0, 10,"D0_PT/1000", " p_{T}(D^{0}) [GeV]"), "1", procs, linplot).RatioTitle("2016 MC", "2011 MC");
  pm.Push<Hist1D>(Axis(100, 0, 10,"(Kplus_PT+piminus0_PT)/1000", "p_{T}(K^{+}) + p_{T}(#pi^{-}) [GeV]",{1.4, 2.5}), "1", procs, linplot).RatioTitle("2016 MC", "2011 MC");

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
