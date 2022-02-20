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
    .Overflow(OverflowType::both);
  PlotOpt lin_lumi = log_lumi().YAxis(YAxisType::linear);
  PlotOpt log_shapes = log_lumi().Stack(StackType::shapes)
    .ShowBackgroundError(false).Bottom(BottomType::ratio);
  PlotOpt lin_shapes = log_shapes().YAxis(YAxisType::linear);
  PlotOpt log_lumi_info = log_lumi().Title(TitleType::info);
  PlotOpt lin_lumi_info = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt log_shapes_info = log_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_shapes_info = lin_shapes().Title(TitleType::info).Bottom(BottomType::ratio);
  PlotOpt lin_lumi_info_print = lin_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  PlotOpt log_lumi_info_print = log_lumi().Title(TitleType::info).Bottom(BottomType::ratio).PrintVals(true);
  
  vector<PlotOpt> linplot = {lin_shapes_info, lin_lumi_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string L0_run1 = "mu_L0Global_TIS && (b0_L0Global_TIS || dst_L0HadronDecision_TOS)";
  string HLT1_run1 = "k_Hlt1TrackAllL0Decision_TOS || pi_Hlt1TrackAllL0Decision_TOS";
  string HLT2_run1 = "d0_Hlt2CharmHadD02HH_D02KPiDecision_TOS";
  
  string L0_run2 = "b0_L0Global_TIS || dst_L0HadronDecision_TOS";
  string HLT1_run2 = "k_Hlt1TrackMVADecision_TOS || pi_Hlt1TrackMVADecision_TOS  || d0_Hlt1TwoTrackMVADecision_TOS";
  string HLT2_run2 = "d0_Hlt2XcMuXForTauB2XcMuDecision_TOS";
  
  string trig_run1 = L0_run1 + "&&(" + HLT1_run1 + ")&&"+HLT2_run1;
  string trig_run2 = L0_run2 + "&&(" + HLT1_run2 + ")&&"+HLT2_run2;

  // Custom NamedFunc
  NamedFunc mu_eta("mu_eta", [&](const Baby &b){
       return log((b.mu_P()+b.mu_PZ())/(b.mu_P()-b.mu_PZ()))/2.;
  });
  NamedFunc muk_log("muk_log", [&](const Baby &b){
       return log10(1-(b.mu_PX()*b.k_PX() + b.mu_PY()*b.k_PY() + b.mu_PZ()*b.k_PZ())/(b.mu_P()*b.k_P()));
  });
  NamedFunc mupi_log("mupi_log", [&](const Baby &b){
       return log10(1-(b.mu_PX()*b.pi_PX() + b.mu_PY()*b.pi_PY() + b.mu_PZ()*b.pi_PZ())/(b.mu_P()*b.pi_P()));
  });
  NamedFunc muspi_log("muspi_log", [&](const Baby &b){
       return log10(1-(b.mu_PX()*b.spi_PX() + b.mu_PY()*b.spi_PY() + b.mu_PZ()*b.spi_PZ())/(b.mu_P()*b.spi_P()));
  });
  NamedFunc log_ip("log_ip", [&](const Baby &b){
                               return log(b.d0_IP_OWNPV());
  });
  NamedFunc b0_dxy("b0_dxy", [&](const Baby &b){
                               return b.b0_FD_OWNPV()*sin(b.b0_FlightDir_Zangle());
  });

  //////////////////////// Step 2 cuts ////////////////////////
  NamedFunc step2_k   = "(k_PT > 800) && (!k_isMuon) && k_IPCHI2_OWNPV > 45"; 
  NamedFunc step2_pi   = "(pi_PT > 800) && (!pi_isMuon) && pi_IPCHI2_OWNPV > 45";
  NamedFunc step2_d0 = "d0_P>2000 && d0_FDCHI2_OWNPV > 250 && (d0_MM-1864.83) < 23.4 && (d0_MM-1864.83) > -23.4 && (k_PT>1700 || pi_PT>1700) && d0_IPCHI2_OWNPV > 9" && log_ip > -3.5;
  // Missing BDTmu
  NamedFunc step2_mu = "mu_isMuon && mu_PIDmu > 2 && mu_PIDe < 1 && mu_P < 100000 " && mu_eta > 1.7 && mu_eta < 5 && muk_log>-6.5 && mupi_log>-6.5 && muspi_log>-6.5;
  // Added the spi cut here since it does nothing
  NamedFunc step2_dsp = "spi_TRACK_GhostProb < 0.5 && (dst_ENDVERTEX_CHI2/dst_ENDVERTEX_NDOF) < 10 && (dst_MM - d0_MM-145.43) < 2 &&  (dst_MM - d0_MM-145.43) > -2";
  NamedFunc step2_b0 = "b0_ISOLATION_BDT < 0.15 && (b0_ENDVERTEX_CHI2/b0_ENDVERTEX_NDOF) < 6 && b0_MM<5280 && b0_DIRA_OWNPV>0.9995" && b0_dxy < 7;
  

  string repofolder = "ntuples/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run1_mc_b0>("D^{*+}#mu#nu 2012", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"0.9.3-production_for_validation/Dst_D0-mc/Dst_D0--21_01_30--mc--MC_2012_Beam4000GeV-2012-MagDown-Nu2.5-Pythia8_Sim08e_Digi13_Trig0x409f0045_Reco14a_Stripping20Filtered_11574020_DSTTAUNU.SAFESTRIPTRIG.DST.root"}), trig_run1 && step2_k && step2_pi && step2_d0 && step2_mu && step2_dsp && step2_b0));
  procs.push_back(Process::MakeShared<Baby_run2_mc_b0>("D^{*+}#mu#nu 2016", Process::Type::data, colors("data"),
                                                set<string>({repofolder+"0.9.4-trigger_emulation/Dst_D0-mc/Dst_D0--21_04_21--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11574021_D0TAUNU.SAFESTRIPTRIG.DST.root"}), trig_run2 && step2_k && step2_pi && step2_d0 && step2_mu && step2_dsp && step2_b0));



  PlotMaker pm;
  pm.Push<Hist1D>(Axis(40, 1, 6,mu_eta, "#eta(#mu)", {2.4, 4}), "1", procs, linplot).RatioTitle("2016", "2012");
  pm.Push<Hist1D>(Axis(80, -8, 12, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "1", procs, linplot).RatioTitle("2016", "2012");
  pm.Push<Hist1D>(Axis(80, -8, 12, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000<6", procs, linplot).RatioTitle("2016", "2012");
  pm.Push<Hist1D>(Axis(80, -8, 12, "FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>6", procs, linplot).RatioTitle("2016", "2012");
  
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
