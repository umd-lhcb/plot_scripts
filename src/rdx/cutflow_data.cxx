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

  string L0 = "muplus_L0Global_TIS && (Y_L0Global_TIS || Dst_2010_minus_L0HadronDecision_TOS)";
  string HLT1_run1 = "Kplus_Hlt1TrackAllL0Decision_TOS || piminus0_Hlt1TrackAllL0Decision_TOS";
  string HLT2_run1 = "D0_Hlt2CharmHadD02HH_D02KPiDecision_TOS";
  string HLT1_run2 = "Kplus_Hlt1Phys_Dec";
  string HLT2_run2 = "D0_Hlt2XcMuXForTauB2XcMuDecision_Dec";
  string trig_run1 = L0 + "&&(" + HLT1_run1 + ")&&"+HLT2_run1;
  string trig_run2 = L0 + "&&(" + HLT1_run2 + ")&&"+HLT2_run2;

  string repofolder = "/Users/manuelf/code/lhcb-ntuples-gen/ntuples/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run1>("Data~2012", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"pre-0.9.0/Dst-cutflow_data/Dst--20_04_03--cutflow_data--data--2012--md.root"}), trig_run1));
  procs.push_back(Process::MakeShared<Baby_run2>("Data~2016", Process::Type::background, colors("run2"),
                                                set<string>({repofolder+"pre-0.9.0/Dst-cutflow_data/Dst--20_04_03--cutflow_data--data--2016--md.root"}), trig_run2));


  // Stripping cuts for Run 1 as a reference
  string strip1_k   = "(k_PIDK > 4) && (k_IPCHI2_OWNPV > 45) && (k_P > 2*1000) && (k_PT > 300) && (k_TRACK_GhostProb < 0.5)";
  string strip1_pi  = "(pi_P > 2*1000) && (pi_PT > 300) && (pi_IPCHI2_OWNPV > 45) && (pi_PIDK < 2) && (pi_TRACK_GhostProb < 0.5)";
  string strip1_d0  = "(k_PT+pi_PT > 1400) && (abs(d0_MM-1864.83) < 80) && (d0_ENDVERTEX_CHI2/d0_ENDVERTEX_NDOF < 4) && (d0_FDCHI2_OWNPV > 250) && (d0_DIRA_OWNPV > 0.9998)";
  string strip1_spi = "(spi_IPCHI2_OWNPV > 0) && (spi_TRACK_CHI2NDOF < 3) && (spi_TRACK_GhostProb < 0.25)";
  string strip1_dst = "(abs(dst_MM - 2010.26) < 125) && (dst_M - d0_M < 160) && (dst_ENDVERTEX_CHI2 / dst_ENDVERTEX_NDOF < 100)";
  string strip1_mu  = "(mu_IPCHI2_OWNPV > 45) && (mu_TRACK_GhostProb < 0.5) && (mu_PIDmu > 2) && (mu_P > 3*1000) && (mu_TRACK_CHI2NDOF < 3)";
  string strip1_b0  = "(0*1000 < b0_MM < 10*1000) && (b0_ENDVERTEX_CHI2 / b0_ENDVERTEX_NDOF < 6) && (b0_DIRA_OWNPV > 0.9995)";
  string strip_run1 = strip1_k+"&&"+strip1_pi+"&&"+strip1_spi+"&&"+strip1_mu+"&&"+strip1_d0+"&&"+strip1_dst+"&&"+strip1_b0;
  
  // Stripping cuts for Run 2 as a reference
  string strip2_k   = "(k_PIDK > 4) && (k_IPCHI2_OWNPV > 9) && (k_P > 2*1000) && (k_PT > 300) && (k_TRACK_GhostProb < 0.5)";
  string strip2_pi  = "(pi_P > 2*1000) && (pi_PT > 300) && (pi_IPCHI2_OWNPV > 9) && (pi_PIDK < 2) && (pi_TRACK_GhostProb < 0.5)";
  string strip2_d0  = "(k_PT+pi_PT > 2500) && (abs(d0_MM - 1864.83) < 80) && (d0_ENDVERTEX_CHI2 / d0_ENDVERTEX_NDOF < 4) && (d0_FDCHI2_OWNPV > 25) && (d0_DIRA_OWNPV > 0.999)";
  string strip2_spi = "(spi_IPCHI2_OWNPV > 0) && (spi_TRACK_CHI2NDOF < 3) && (spi_TRACK_GhostProb < 0.25)";
  string strip2_dst = "(abs(dst_MM - 2010.26) < 125) && (dst_M - d0_M < 160) && (dst_ENDVERTEX_CHI2 / dst_ENDVERTEX_NDOF < 100)";
  string strip2_mu  = "(mu_IPCHI2_OWNPV > 16) && (mu_TRACK_GhostProb < 0.5) && (mu_PIDmu > -200) && (mu_P > 3*1000) && (mu_TRACK_CHI2NDOF < 3)";
  string strip2_b0  = "(0*1000 < b0_MM < 10*1000) && (b0_ENDVERTEX_CHI2 / b0_ENDVERTEX_NDOF < 6) && (b0_DIRA_OWNPV > 0.999)";
  string strip_run2 = strip2_k+"&&"+strip2_pi+"&&"+strip2_spi+"&&"+strip2_mu+"&&"+strip2_d0+"&&"+strip2_dst+"&&"+strip2_b0;

  // Step 2 cuts

  // Missing HLT1 and HLT1_pT
  string step2_k   = "(Kplus_PT > 800) && (!Kplus_isMuon) && Kplus_IPCHI2_OWNPV > 45"; 
  // Missing HLT1 and HLT1_pT
  string step2_pi   = "(piminus0_PT > 800) && (!piminus0_isMuon) && piminus0_IPCHI2_OWNPV > 45";
  // Missing HLT2
  string step2_d0 = "D0_P>2 && D0_DIRA_OWNPV > 0.9998 && D0_FDCHI2_OWNPV > 250 && (D0_MM-1864.83) < 23.4 && (D0_MM-1864.83) > -23.4 && (Kplus_PT>1700 || piminus0_PT>1700)";

  // Missing eta, BDTmu, log10
  string step2_mu = "muplus_isMuon && muplus_PIDmu > 2 && muplus_PIDe < 1 && muplus_P < 100000 ";
  string step2_pis = "piminus_TRACK_GhostProb < 0.5";
  string step2_dsp = "(Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF) < 10 && (Dst_2010_minus_MM - D0_MM-145.43) < 2 &&  (Dst_2010_minus_MM - D0_MM-145.43) > -2";

  //Missing d_XY < 7 mm (FDCHI2_ORIVX/Y)
  string step2_b0 = "Y_ISOLATION_BDT < 0.15 && (Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF) < 10 && Y_MM<5280 && Y_DIRA_OWNPV>0.9995";
  


  PlotMaker pm;
  pm.Push<Table>("cutflow", vector<TableRow>{
  TableRow("Stripping","1", 0,0, "1"),
    TableRow("Kaon", step2_k,0,0, "1"),
    TableRow("Pion", step2_k+"&&"+step2_pi, 0,0, "1"),
    TableRow("$D^0 \\rightarrow K \\pi$", step2_k+"&&"+step2_pi+"&&"+step2_d0,0,0, "1"),
    TableRow("$\\mu$", step2_k+"&&"+step2_pi+"&&"+step2_d0+"&&"+step2_mu,0,0, "1"),
    TableRow("$\\pi_\\text{soft}$", step2_k+"&&"+step2_pi+"&&"+step2_d0+"&&"+step2_mu+"&&"+step2_pis,0,0, "1"),
    TableRow("$D^{*}$", step2_k+"&&"+step2_pi+"&&"+step2_d0+"&&"+step2_mu+"&&"+step2_pis+"&&"+step2_dsp,0,0, "1"),
    TableRow("$B^{0}$", step2_k+"&&"+step2_pi+"&&"+step2_d0+"&&"+step2_mu+"&&"+step2_pis+"&&"+step2_dsp+"&&"+step2_b0,0,0, "1"),

//    TableRow("$\\mu$ PID", muid,0,0, "1"),
//    TableRow("$\\text{IsoBDT}_{B} < 0.15$", iso,0,0, "1"),
//    TableRow("$m_{B} < 5280$", bmass,0,0, "1"),
//    TableRow("$p_{T}(K)>0.8$ GeV", kpt,0,0, "1"),
//    TableRow("$p_{T}(\\pi)>0.8$ GeV", pipt,0,0, "1"),
//    TableRow("$p_{T}(K)>1.7$ GeV or $p_{T}(\\pi)>1.7$ GeV", maxpt,0,0, "1"),

    },procs,0); // Pushing table
  
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
