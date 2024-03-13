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
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string globalCuts = "mu_ubdt_ok && (k_p < 200) && (pi_p < 200) && (mu_p < 100) && (iso_p1 < 200) && (iso_p2 < 200) && (iso_p3 < 200) && (nspdhits < 450) && is_iso";

  string repofolder = "ntuples/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run2_std>("D*^{+}#mu#nu MC 2016", Process::Type::background, colors("darkblue"),
                                                        set<string>({repofolder+"0.9.6-2016_production/Dst_D0-mc-tracker_only-sig_norm/Dst--24_02_12--mc--11574021--2016--mu--tracker_only.root"}), globalCuts));
  procs.push_back(Process::MakeShared<Baby_run2_std>("Data 2016", Process::Type::data, colors("data"),
                                                       set<string>({repofolder+"0.9.6-2016_production/Dst_D0-std/Dst--23_11_06--std--data--2016--mu.root"}), globalCuts));


     
  NamedFunc woverjk("woverjk", [&](const Baby &b){
    if(b.wjk() != 0) return b.wiso()/b.wjk();
    else return 0.;
  });
  //vector<NamedFunc> weights({"wskim_iso*skim_global_ok*wff*wpid_ubdt*wtrg*wtrk*wbr_dd*w_missDDX", "1"});
  vector<NamedFunc> weights({"wskim_iso*skim_global_ok*wff*wpid_ubdt*wtrg*wtrk*wbr_dd*w_missDDX*wjk", "1"});
  //vector<NamedFunc> weights({woverjk, "1"});
  PlotMaker pm;

  pm.Push<Hist1D>(Axis(100,0,20, "mu_pt", "p_{T}(#mu^{+}) [GeV]"), "mm2<0.5", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  pm.Push<Hist1D>(Axis(100,1.7,5, "mu_eta", "#eta(#mu^{+})"), "mm2<0.5", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");
  //pm.Push<Hist1D>(Axis(100,0,20, "mu_pt", "p_{T}(#mu^{+}) [GeV]"), "1", procs, linplot, weights).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");


  // //pm.Push<Table>("cutflow", table_rows,procs,0).TotColumn("Ratio", 614577/1500395.*0.23/0.07*0.080/0.059); // BKK events for the two MC samples
 
  // pm.Push<Hist1D>(Axis(70, -3, 11, "FitVar_q2/1000000", "q^{2} [GeV^{2}]"), "FitVar_Mmiss2/1000000<0.5" && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("dmc");  
  // pm.Push<Hist1D>(Axis(45, -3, 6, "FitVar_q2/1000000", "q^{2} [GeV^{2}]"), "FitVar_q2/1000000<6" && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("q^{2} < 6 GeV^{2}").Tag("dmc");  

  // ///////// FIT VARIABLES
  // pm.Push<Hist1D>(Axis(70, -1, 13, "FitVar_q2/1000000", "q^{2} [GeV^{2}]"), fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("Baseline").Tag("mc");  
  // pm.Push<Hist1D>(Axis(70, 0, 3, "FitVar_El/1000", "E*_{l} [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(80, -2, 10, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("Baseline").Tag("mc");
  // pm.Push<Hist1D>(Axis(80, -2, 10, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(80, -2, 10, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>6" && fullcut1 , procs, linplot).RatioTitle("Data", "MC").SetTitle("q^{2} > 6 GeV^{2}").Tag("mc");
  // pm.Push<Hist1D>(Axis(80, -2, 10, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000<6" && fullcut1 , procs, linplot).RatioTitle("Data", "MC").SetTitle("q^{2} > 6 GeV^{2}").Tag("mc");
  // pm.Push<Hist1D>(Axis(100, -3, 0.5, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_Mmiss2/1000000<0.5" && fullcut1 , procs, linplot).RatioTitle("Data", "MC").SetTitle("m^{2}_{miss} < 0.5 GeV^{2}").Tag("mc");

  // ///////// KINEMATIC VARIABLES (P, PT, ETA)
  // pm.Push<Hist1D>(Axis(80, 1.5, 5.5,mu_eta, "#eta(#mu)"), fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("Baseline").Tag("mc");
  // pm.Push<Hist1D>(Axis(60,0,30, "spi_P/1000", "p(#pi^{-}_{slow}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(70,0,210, "d0_P/1000", "p(D^{0}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,100, "mu_P/1000", "p(#mu^{+}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(50,0,2, "spi_PT/1000", "p_{T}(#pi^{-}_{slow}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,25, "d0_PT/1000", "p_{T}(D^{0}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,10, "mu_PT/1000", "p_{T}(#mu^{+}) [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");

  // ///////// mu, K, pi VARIABLES
  // pm.Push<Hist1D>(Axis(70,0,140, "k_PIDK"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(75,-140,10, "pi_PIDK"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,3, "k_TRACK_CHI2NDOF", "K #chi^{2}_{track}/Ndof"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,3, "pi_TRACK_CHI2NDOF", "#pi #chi^{2}_{track}/Ndof"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,3, "spi_TRACK_CHI2NDOF", "#pi_{slow} #chi^{2}_{track}/Ndof"), selcut && fullcut1, procs,linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,3, "mu_TRACK_CHI2NDOF", "#mu #chi^{2}_{track}/Ndof"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(50,0,0.05, "mu_TRACK_GhostProb"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(70,0,14, "mu_PIDmu"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");

  // ///////// D* VARIABLES
  // pm.Push<Hist1D>(Axis(72,141,150, "dst_MM - d0_MM", "Prefit m(D*)-m(D^{0}) [MeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(60,1835,1895, "d0_MM", "Prefit m(D^{0}) [MeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,100, "dst_FD_ORIVX"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0.9998,1, "d0_DIRA_OWNPV"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0,20000, "d0_FDCHI2_OWNPV"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  
  // ///////// B VARIABLES
  // pm.Push<Hist1D>(Axis(46,-1,0.15, "b0_ISOLATION_BDT"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(100,0.9995,1, "b0_DIRA_OWNPV"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(50,0,25, "b0_ENDVERTEX_CHI2"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(70,2,5.5, "b0_M/1000", "B^{0} mass [GeV]"), selcut && fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle(seltitle).Tag("mc");
  // pm.Push<Hist1D>(Axis(50, 0, 1, "b0_ISOLATION_NNk"), fullcut1, procs, linplot).RatioTitle("Data", "MC").SetTitle("Baseline").Tag("data");  

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
