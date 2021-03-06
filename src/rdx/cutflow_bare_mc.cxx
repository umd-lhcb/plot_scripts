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
  
  vector<PlotOpt> linplot = {lin_lumi_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string L0_run1 = "mu_L0Global_TIS && (b0_L0Global_TIS || d0_L0HadronDecision_TOS)";
  string HLT1_run1 = "(k_Hlt1TrackAllL0Decision_TOS || pi_Hlt1TrackAllL0Decision_TOS) && ((k_Hlt1TrackAllL0Decision_TOS && k_PT > 1700) || (pi_Hlt1TrackAllL0Decision_TOS && pi_PT > 1700))";
  string HLT2_run1 = "d0_Hlt2CharmHadD02HH_D02KPiDecision_TOS";
  
  string L0_run2 = "(b0_L0Global_TIS || d0_L0HadronDecision_TOS)";
  string HLT1_run2 = "(d0_Hlt1TrackMVALooseDecision_TOS  || d0_Hlt1TwoTrackMVADecision_TOS) && ((d0_Hlt1TwoTrackMVADecision_TOS || k_Hlt1TrackMVALooseDecision_TOS) && k_PT > 1700 || (d0_Hlt1TwoTrackMVADecision_TOS || pi_Hlt1TrackMVALooseDecision_TOS) && pi_PT > 1700)";
  string HLT2_run2 = "b0_Hlt2XcMuXForTauB2XcMuDecision_TOS";
  
  string trig_run1 = L0_run1 + "&& (" + HLT1_run1 + ") &&"+HLT2_run1;
  string trig_run2 = L0_run2 + "&& (" + HLT1_run2 + ") &&"+HLT2_run2;

  NamedFunc isNorm("isNorm", [&](const Baby &b){
       return abs(b.mu_MC_MOTHER_ID())==511 & abs(b.d0_MC_MOTHER_ID())==413 & abs(b.d0_MC_GD_MOTHER_ID())==511;
  });
  NamedFunc isSignal("isSignal", [&](const Baby &b){
       return abs(b.mu_MC_MOTHER_ID())==15 & abs(b.d0_MC_MOTHER_ID())==413 & abs(b.d0_MC_GD_MOTHER_ID())==511;
  });

  string repofolder = "ntuples/";
  vector<shared_ptr<Process> > procs_norm;
   procs_norm.push_back(Process::MakeShared<Baby_run1_bare>("D^{*+}#mu#nu~2011", Process::Type::background, colors("green"),
                                                 set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2011_Beam3500GeV-2011-MagDown-Nu2-Pythia8_Sim08h_Digi13_Trig0x40760037_Reco14c_Stripping20r1NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), isNorm));
   procs_norm.push_back(Process::MakeShared<Baby_run2_bare>("D^{*+}#mu#nu~2016", Process::Type::data, colors("darkblue"),
                                                 set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), isNorm));
  
  vector<shared_ptr<Process> > procs_sig;
  procs_sig.push_back(Process::MakeShared<Baby_run1_bare>("D^{*+}#tau#nu~2011", Process::Type::background, colors("green"),
                                               set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2011_Beam3500GeV-2011-MagDown-Nu2-Pythia8_Sim08h_Digi13_Trig0x40760037_Reco14c_Stripping20r1NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), isSignal));
  procs_sig.push_back(Process::MakeShared<Baby_run2_bare>("D^{*+}#tau#nu~2016", Process::Type::data, colors("darkblue"),
                                               set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), isSignal));
  
  vector<shared_ptr<Process> > procs_dss;
  procs_dss.push_back(Process::MakeShared<Baby_run1_bare>("D^{**}#mu#nu~2011", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2011_Beam3500GeV-2011-MagDown-Nu2-Pythia8_Sim08h_Digi13_Trig0x40760037_Reco14c_Stripping20r1NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), !isNorm && !isSignal));
  procs_dss.push_back(Process::MakeShared<Baby_run2_bare>("D^{**}#mu#nu~2016", Process::Type::data, colors("darkblue"),
                                                set<string>({repofolder+"0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), !isNorm && !isSignal));



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


  // Stripping cuts for Run 1
  string strip1_k   = "(k_PIDK > 4) && (k_IPCHI2_OWNPV > 45) && (k_P > 2*1000) && (k_PT > 300) && (k_TRACK_GhostProb < 0.5)";
  string strip1_pi  = "(pi_P > 2*1000) && (pi_PT > 300) && (pi_IPCHI2_OWNPV > 45) && (pi_PIDK < 2) && (pi_TRACK_GhostProb < 0.5)";
  string strip1_d0  = "(k_PT+pi_PT > 1400) && (d0_MM-1864.83) < 80 && (d0_MM-1864.83) > -80 && (d0_ENDVERTEX_CHI2/d0_ENDVERTEX_NDOF < 4) && (d0_FDCHI2_OWNPV > 250) && (d0_DIRA_OWNPV > 0.9998)";
  string strip1_spi = "(spi_IPCHI2_OWNPV > 0) && (spi_TRACK_CHI2NDOF < 3) && (spi_TRACK_GhostProb < 0.25)";
  string strip1_dst = "(dst_MM - 2010.26) < 125 && (dst_MM - 2010.26) > -125 && (dst_M - d0_M < 160) && (dst_ENDVERTEX_CHI2 / dst_ENDVERTEX_NDOF < 100)";
  string strip1_mu  = "(mu_IPCHI2_OWNPV > 45) && (mu_TRACK_GhostProb < 0.5) && (mu_PIDmu > 2) && (mu_P > 3*1000) && (mu_TRACK_CHI2NDOF < 3)";
  string strip1_b0  = "(b0_MM < 10*1000) && (b0_ENDVERTEX_CHI2 / b0_ENDVERTEX_NDOF < 6) && (b0_DIRA_OWNPV > 0.9995)";
  string strip_run1 = strip1_k+"&&"+strip1_pi+"&&"+strip1_spi+"&&"+strip1_mu+"&&"+strip1_d0+"&&"+strip1_dst+"&&"+strip1_b0;
  
  // Stripping cuts for Run 2
  string strip2_k   = "(k_PIDK > 4) && (k_IPCHI2_OWNPV > 9) && (k_P > 2*1000) && (k_PT > 300) && (k_TRACK_GhostProb < 0.5)";
  string strip2_pi  = "(pi_P > 2*1000) && (pi_PT > 300) && (pi_IPCHI2_OWNPV > 9) && (pi_PIDK < 2) && (pi_TRACK_GhostProb < 0.5)";
  string strip2_d0  = "(k_PT+pi_PT > 2500) && (d0_MM - 1864.83) < 80 && (d0_MM - 1864.83) > -80 && (d0_ENDVERTEX_CHI2 / d0_ENDVERTEX_NDOF < 4) && (d0_FDCHI2_OWNPV > 25) && (d0_DIRA_OWNPV > 0.999)";
  string strip2_spi = "(spi_IPCHI2_OWNPV > 0) && (spi_TRACK_CHI2NDOF < 3) && (spi_TRACK_GhostProb < 0.25)";
  string strip2_dst = "(dst_MM - 2010.26) < 125 && (dst_MM - 2010.26) > -125 && (dst_M - d0_M < 160) && (dst_ENDVERTEX_CHI2 / dst_ENDVERTEX_NDOF < 100)";
  string strip2_mu  = "(mu_IPCHI2_OWNPV > 16) && (mu_TRACK_GhostProb < 0.5) && (mu_PIDmu > -200) && (mu_P > 3*1000) && (mu_TRACK_CHI2NDOF < 3)";
  string strip2_b0  = "(b0_MM < 10*1000) && (b0_ENDVERTEX_CHI2 / b0_ENDVERTEX_NDOF < 6) && (b0_DIRA_OWNPV > 0.999)";
  string strip_run2 = strip2_k+"&&"+strip2_pi+"&&"+strip2_spi+"&&"+strip2_mu+"&&"+strip2_d0+"&&"+strip2_dst+"&&"+strip2_b0;

  ///////// Automatically appending cutflow cuts
  vector<NamedFunc> cuts1 = {L0_run1, HLT1_run1, HLT2_run1, strip_run1, "(k_PT > 800)  && k_P > 2000 && k_IPCHI2_OWNPV > 45 && k_TRACK_GhostProb < 0.5",
                            "(pi_PT > 800) && pi_P > 2000 && pi_IPCHI2_OWNPV > 45 && pi_TRACK_GhostProb < 0.5",
                            "d0_PT>2000 && (d0_ENDVERTEX_CHI2 / d0_ENDVERTEX_NDOF < 4) && d0_IPCHI2_OWNPV > 9  && (d0_DIRA_OWNPV > 0.9998) && d0_FDCHI2_OWNPV > 250 && (d0_M-1865.49) < 23.4 && (d0_M-1865.49) > -23.4"
                            && log_ip > -3.5,
                            "mu_P > 3000 && mu_P < 100000 && mu_IPCHI2_OWNPV > 45 && mu_TRACK_GhostProb < 0.5"
                            && mu_eta > 1.7 && mu_eta < 5 && muk_log>-6.5 && mupi_log>-6.5 && muspi_log>-6.5,
                            "spi_TRACK_GhostProb < 0.25 && (dst_ENDVERTEX_CHI2/dst_ENDVERTEX_NDOF) < 10 && (dst_M - d0_M-145.454) < 2 &&  (dst_M - d0_M-145.454) > -2",
                            "b0_ENDVERTEX_CHI2 < 24 && (b0_ENDVERTEX_CHI2/b0_ENDVERTEX_NDOF) < 6 && b0_M<5280 && b0_DIRA_OWNPV>0.9995 && b0_DISCARDMu_CHI2 <= 6" && b0_dxy < 7,
                            "b0_ISOLATION_BDT < 0.15",
                            "mu_isMuon && mu_PIDmu > 2 && mu_PIDe < 1 && (!k_isMuon) && (!pi_isMuon) "};

  bool doRun2Cuts = false;
  string trk_pT = "800", d0_dira = "0.9998", b0_dira = "0.9995", d0_fd = "250";
  if(doRun2Cuts) {
    trk_pT = "200"; d0_dira = "0.999"; b0_dira = "0.999"; d0_fd = "25";
  }
  vector<NamedFunc> cuts2 = {L0_run2, HLT1_run2, HLT2_run2, strip_run2, "(k_PT > "+trk_pT+")  && k_P > 2000 && k_IPCHI2_OWNPV > 45 && k_TRACK_GhostProb < 0.5",
                            "(pi_PT > "+trk_pT+") && pi_P > 2000 && pi_IPCHI2_OWNPV > 45 && pi_TRACK_GhostProb < 0.5",
                            "d0_PT>2000 && (d0_ENDVERTEX_CHI2 / d0_ENDVERTEX_NDOF < 4) && d0_IPCHI2_OWNPV > 9  && (d0_DIRA_OWNPV > "+d0_dira+") && d0_FDCHI2_OWNPV > "+d0_fd+" && (d0_M-1865.49) < 23.4 && (d0_M-1865.49) > -23.4"
                            && log_ip > -3.5,
                            "mu_P > 3000 && mu_P < 100000 && mu_IPCHI2_OWNPV > 45 && mu_TRACK_GhostProb < 0.5"
                            && mu_eta > 1.7 && mu_eta < 5 && muk_log>-6.5 && mupi_log>-6.5 && muspi_log>-6.5,
                            "spi_TRACK_GhostProb < 0.25 && (dst_ENDVERTEX_CHI2/dst_ENDVERTEX_NDOF) < 10 && (dst_M - d0_M-145.454) < 2 &&  (dst_M - d0_M-145.454) > -2",
                            "b0_ENDVERTEX_CHI2 < 24 && (b0_ENDVERTEX_CHI2/b0_ENDVERTEX_NDOF) < 6 && b0_M<5200 && b0_DIRA_OWNPV>"+b0_dira+" && b0_DISCARDMu_CHI2 <= 6" && b0_dxy < 7,
                            "b0_ISOLATION_BDT < 0.15",
                            "mu_isMuon && mu_PIDmu > 2 && mu_PIDe < 1 && (!k_isMuon) && (!pi_isMuon) "};
  
                            
  vector<string> rownames = {"L0", "HLT1", "HLT2", "Full strip.", "Kaon", "Pion", "$D^0 \\rightarrow K \\pi$","$\\mu$",
    "$D^{*+} \\rightarrow D^0 \\pi$", "$B^{0} \\rightarrow D^{*+} \\mu$", "ISO", "PID"};
  vector<TableRow> table_rows;
  NamedFunc fullcut1 = "1", fullcut2 = "1";
  for(size_t ind = 0; ind < cuts1.size(); ind++) {
    string title = (ind==0 ? rownames[ind] : "+ " + rownames[ind]);
    int lines = (ind==3 ? 1 : 0);
    fullcut1 = fullcut1 && cuts1[ind];
    fullcut2 = fullcut2 && cuts2[ind];
    table_rows.push_back(TableRow(title,{fullcut1, fullcut2}, 0,lines, "1"));
  }
  table_rows.push_back(TableRow("$q^2 < 6\\text{ GeV}^2$",{fullcut1 && "FitVar_q2/1000000 < 6", fullcut2 && "FitVar_q2/1000000 < 6"}, 1,0, "1"));
  table_rows.push_back(TableRow("$q^2 > 6\\text{ GeV}^2$",{fullcut1 && "FitVar_q2/1000000 > 6",fullcut2 && "FitVar_q2/1000000 > 6"}, 0,0, "1"));

  PlotMaker pm;
  pm.Push<Table>("cutflow", table_rows,procs_norm,0).TotColumn("Ratio").Tag("norm_");
  pm.Push<Table>("cutflow", table_rows,procs_sig,0).TotColumn("Ratio").Tag("sig_");
  pm.Push<Table>("cutflow", table_rows,procs_dss,0).TotColumn("Ratio").Tag("dss_");
  
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
