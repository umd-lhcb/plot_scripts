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
  
  vector<PlotOpt> linplot = {lin_shapes, log_shapes};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string repofolder = "../lhcb-ntuples-gen/gen/";
  vector<shared_ptr<Process> > procs;
  // procs.push_back(Process::MakeShared<Baby_phoebe_step2>("Phoebe", Process::Type::background, colors("run2"),
  //                                               set<string>({repofolder+"run1-Dst-step2/Dst--20_07_02--mix--data--2011--md--phoebe-step2.root"}),
  //                                                 "flag_2011 && polarity < 0 && spi_trk_type == 3"));
  procs.push_back(Process::MakeShared<Baby_phoebe_step1>("Phoebe", Process::Type::background, colors("run2"),
                                                set<string>({repofolder+"run1-Dst-step2/Dst--20_09_16--std--data--2011--md--phoebe-step2.root"}), "1"));
  procs.push_back(Process::MakeShared<Baby_step2>("Yipeng", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"run1-Dst_D0-step2/Dst_D0--20_10_12--std--data--2011--md--step2.root"}),"1"));

  ///////// Automatically appending cutflow cuts
  vector<NamedFunc> cuts = {"1", "(b_l0_tis || d0_l0had_tos)", "mu_l0_tis", "mu_pidb > 0",
                            "dst_id_prod > 0", "id_prod > 0", "flag_sel_b0dst"};
  
  vector<string> rownames = {"Data MD 2011", "L0: TIS or D0 TOS", "L0: $\\mu$ TIS",  "$\\mu$ PID",
                             "$D^0\\pi^+$ sign","$D^0\\mu^+$ sign", "BDT$<0.15$, \\texttt{flag\\_b0}"};

  vector<TableRow> table_rows;
  NamedFunc fullcut = "1";
  for(size_t ind = 0; ind < cuts.size(); ind++) {
    string title = (ind==0 ? rownames[ind] : "+ " + rownames[ind]);
    int lines = (ind==0 ? 1 : 0);
    fullcut = fullcut && cuts[ind];
    
    //title = rownames[ind];
    //fullcut = cuts[ind];
    
    table_rows.push_back(TableRow(title,fullcut, 0,lines, "1"));
  }
  
  PlotMaker pm;
  pm.Push<Table>("cutflow", table_rows,procs,0).TotColumn("Ratio"); // Pushing table
  
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
