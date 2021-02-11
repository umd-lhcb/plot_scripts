// Program to compare run 1 MC between us/Phoebe and compare us run1/run2, to validate our reco script, using plot_scripts
// Created: Jan 17, 2021
// Last edited: Feb 9, 2021
// Note: relevant ntuples (see file names below) need to be downloaded, and the ntuple file paths in the processes
// should be edited to reflect their location. The ntuples are generated from dst files found (in DIRAC) at
// (for run 1) MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim08e/Digi13/Trig0x409f0045/Stripping20Filtered/11574020/DSTTAUNU.SAFESTRIPTRIG.DST
// (for run 2) MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09j/Trig0x6139160F/Reco16/Turbo03a/Filtered/11574021/D0TAUNU.SAFESTRIPTRIG.DST

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

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
#include "core/event_scan.hpp"
#include "core/hist1d.hpp"
#include "core/hist2d.hpp"
#include "core/plot_opt.hpp"

using namespace std;
using namespace PlotOptTypes;


int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  // User defined colors
  Palette colors("txt/colors.txt", "default");

  //////////////////////////// Define Plot Types

  PlotOpt lin_data_norm("txt/plot_styles.txt", "LHCbPaper");
  lin_data_norm.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none);//both);
  PlotOpt lin_sig_on_top = lin_data_norm().Stack(StackType::signal_on_top);

  vector<PlotOpt> plotnorm = {lin_data_norm};
  vector<PlotOpt> plotboth = {lin_data_norm, lin_sig_on_top};

  /////////////////////////// Define Cuts, Processes

  double dst_2010_m = 2010.26;
  //double d0_m = 1864.83;

  // Cuts are the same between the three ntuples, but individual NamedFuncs are created because of different variable names
  // Variables for Phoebe's ntuple should be defined in txt/variables/phoebe_step1_mc_Bd2DstMuNu
  NamedFunc cuts_phoebe("cuts_phoebe", [&](const Baby &b){
    // L0
    if ((b.Dst_2010_minus_L0HadronDecision_TOS() && !b.Y_L0Global_TIS()) && !(b.Dst_2010_minus_PT()>3700 && b.nSPDhits()<450)) return false;
    if (!(b.muplus_L0Global_TIS())) return false;
    // HLT2
    else if (!(b.Kplus_TRACK_CHI2NDOF()<3 && b.piminus0_TRACK_CHI2NDOF()<3 && b.Kplus_PT()>800 && b.piminus0_PT()>800 && b.Kplus_P()>5000 && b.piminus0_P()>5000 && (b.Kplus_PT()>1500 || b.piminus0_PT()>1500) && (b.Dst_2010_minus_ENDVERTEX_CHI2()/b.Dst_2010_minus_ENDVERTEX_NDOF())<10 && b.D0_PT()>2000
                && (b.Kplus_PT()+b.piminus0_PT())>2500 && 1830<b.D0_MM() && b.D0_MM()<1910 && b.D0_DIRA_OWNPV()>0.999)) return false;
    // Stripping
    else if (!(b.muplus_IPCHI2_OWNPV()>45 && b.muplus_TRACK_GhostProb()<0.5 && b.muplus_PIDmu()>2 && b.muplus_P()>3000 && b.piminus0_IPCHI2_OWNPV()>45 && b.Kplus_IPCHI2_OWNPV()>45 && b.Kplus_PIDK()>4 && b.piminus0_PIDK()<2 && b.piminus0_TRACK_GhostProb()<0.5 && b.Kplus_TRACK_GhostProb()<0.5
                && (b.D0_ENDVERTEX_CHI2()/b.D0_ENDVERTEX_NDOF())<4 && b.D0_FDCHI2_OWNPV()>250 && (b.Y_ENDVERTEX_CHI2()/b.Y_ENDVERTEX_NDOF())<6 && b.Y_DIRA_OWNPV()>0.9995)) return false;
    // Reco Script
    else if (!(b.piminus_TRACK_Type()==3 && b.Y_MM()<5400 && abs(b.Dst_2010_minus_MM()-dst_2010_m)<125 && abs(b.Dst_2010_minus_MM()-b.D0_MM())<160 && b.muplus_TRACK_CHI2NDOF()<3 && b.piminus_IPCHI2_OWNPV()>0 && b.piminus_TRACK_CHI2NDOF()<3 && b.piminus_TRACK_GhostProb()<0.25)) return false;
    else return true;
  });
  // Variables for our run 1 ntuple should be defined in txt/variables/run1_step1_mc_Bd2DstMuNu
  NamedFunc cuts_us_run1("cuts_us_run1", [&](const Baby &b){
    // L0
    if ((b.dst_L0HadronDecision_TOS() && !b.b0_L0Global_TIS()) && !(b.dst_PT()>3700 && b.nSPDhits()<450)) return false;
    if (!(b.mu_L0Global_TIS())) return false;
    // HLT2
    else if(!(b.k_TRACK_CHI2NDOF()<3 && b.pi_TRACK_CHI2NDOF()<3 && b.k_PT()>800 && b.pi_PT()>800 && b.k_P()>5000 && b.pi_P()>5000 && (b.k_PT()>1500 || b.pi_PT()>1500) && (b.dst_ENDVERTEX_CHI2()/b.dst_ENDVERTEX_NDOF())<10 && b.d0_PT()>2000 && (b.k_PT()+b.pi_PT())>2500 && 1830<b.d0_MM()
               && b.d0_MM()<1910 && b.d0_DIRA_OWNPV()>0.999)) return false;
    // Stripping
    else if (!(b.mu_IPCHI2_OWNPV()>45 && b.mu_TRACK_GhostProb()<0.5 && b.mu_PIDmu()>2 && b.mu_P()>3000 && b.pi_IPCHI2_OWNPV()>45 && b.k_IPCHI2_OWNPV()>45 && b.k_PIDK()>4 && b.pi_PIDK()<2 && b.pi_TRACK_GhostProb()<0.5 && b.k_TRACK_GhostProb()<0.5 && (b.d0_ENDVERTEX_CHI2()/b.d0_ENDVERTEX_NDOF())<4
               && b.d0_FDCHI2_OWNPV()>250 && (b.b0_ENDVERTEX_CHI2()/b.b0_ENDVERTEX_NDOF())<6 && b.b0_DIRA_OWNPV()>0.9995)) return false;
    // Reco Script
    else if (!(b.spi_TRACK_Type()==3 && b.b0_MM()<5400 && abs(b.dst_MM()-dst_2010_m)<125 && abs(b.dst_MM()-b.d0_MM())<160 && b.mu_TRACK_CHI2NDOF()<3 && b.spi_IPCHI2_OWNPV()>0 && b.spi_TRACK_CHI2NDOF()<3 && b.spi_TRACK_GhostProb()<0.25)) return false;
    else return true;
  });
  // Variables for our run 2 ntuple should be defined in txt/variables/run2_step1_mc_Bd2DstMuNu
  // In fact, the variable names here are identical to run 1; but there are still some different variables in the tree, so I'll keep this cut function separate...
  NamedFunc cuts_us_run2("cuts_us_run2", [&](const Baby &b){
    // L0
    if ((b.dst_L0HadronDecision_TOS() && !b.b0_L0Global_TIS()) && !(b.dst_PT()>3700 && b.nSPDhits()<450)) return false;
    if (!(b.mu_L0Global_TIS())) return false;
    // HLT2
    else if (!(b.k_TRACK_CHI2NDOF()<3 && b.pi_TRACK_CHI2NDOF()<3 && b.k_PT()>800 && b.pi_PT()>800 && b.k_P()>5000 && b.pi_P()>5000 && (b.k_PT()>1500 || b.pi_PT()>1500) && (b.dst_ENDVERTEX_CHI2()/b.dst_ENDVERTEX_NDOF())<10 && b.d0_PT()>2000 && (b.k_PT()+b.pi_PT())>2500 && 1830<b.d0_MM() && b.d0_MM()<1910 && b.d0_DIRA_OWNPV()>0.999)) return false;
    // Stripping
    else if (!(b.mu_IPCHI2_OWNPV()>45 && b.mu_TRACK_GhostProb()<0.5 && b.mu_PIDmu()>2 && b.mu_P()>3000 && b.pi_IPCHI2_OWNPV()>45 && b.k_IPCHI2_OWNPV()>45 && b.k_PIDK()>4 && b.pi_PIDK()<2 && b.pi_TRACK_GhostProb()<0.5 && b.k_TRACK_GhostProb()<0.5 && (b.d0_ENDVERTEX_CHI2()/b.d0_ENDVERTEX_NDOF())<4 && b.d0_FDCHI2_OWNPV()>250 && (b.b0_ENDVERTEX_CHI2()/b.b0_ENDVERTEX_NDOF())<6 && b.b0_DIRA_OWNPV()>0.9995)) return false;
    // Reco Script
    else if (!(b.spi_TRACK_Type()==3 && b.b0_MM()<5400 && abs(b.dst_MM()-dst_2010_m)<125 && abs(b.dst_MM()-b.d0_MM())<160 && b.mu_TRACK_CHI2NDOF()<3 && b.spi_IPCHI2_OWNPV()>0 && b.spi_TRACK_CHI2NDOF()<3 && b.spi_TRACK_GhostProb()<0.25)) return false;
    else return true;
  });


  string repofolder = "./";

  vector<shared_ptr<Process>> procs_us_phoebe;
  vector<shared_ptr<Process>> procs_us_run12;
  procs_us_phoebe.push_back(Process::MakeShared<Baby_phoebe_step1_mc_Bd2DstMuNu>("Phoebe", Process::Type::signal, colors("blue"),
                                                set<string>({repofolder+"../ref_ntuples/TFYCands_Bd2Dstmunu_e_Py8_MD.root"}), cuts_phoebe));
  procs_us_phoebe.push_back(Process::MakeShared<Baby_run1_step1_mc_Bd2DstMuNu>("Us (run1)", Process::Type::data, colors("data"),
                                                set<string>({repofolder+"ntuples/0.9.3-production_for_validation/Dst_D0-mc/Dst_D0--21_01_30--mc--MC_2012_Beam4000GeV-2012-MagDown-Nu2.5-Pythia8_Sim08e_Digi13_Trig0x409f0045_Reco14a_Stripping20Filtered_11574020_DSTTAUNU.SAFESTRIPTRIG.DST.root"}), cuts_us_run1));
  procs_us_run12.push_back(Process::MakeShared<Baby_run1_step1_mc_Bd2DstMuNu>("Run 1", Process::Type::data, colors("data"),
                                                set<string>({repofolder+"ntuples/0.9.3-production_for_validation/Dst_D0-mc/Dst_D0--21_01_30--mc--MC_2012_Beam4000GeV-2012-MagDown-Nu2.5-Pythia8_Sim08e_Digi13_Trig0x409f0045_Reco14a_Stripping20Filtered_11574020_DSTTAUNU.SAFESTRIPTRIG.DST.root"}), cuts_us_run1));
  procs_us_run12.push_back(Process::MakeShared<Baby_run2_step1_mc_Bd2DstMuNu>("Run 2", Process::Type::signal, colors("blue"),
                                                set<string>({repofolder+"../lhcb-ntuples-gen/run2-rdx/samples/Dst_D0--21_01_27--mc--Bd2DstMuNu--2016--md--py8-sim09j-dv45-subset.root"}), cuts_us_run2));


  ///////////////////////// Define Custom NamedFuncs for complicated variables

  NamedFunc dst_minus_d0_m_us("dst_minus_d0_m_us", [&](const Baby &b) {
    return (b.dst_MM()-b.d0_MM());
  });

  NamedFunc dst_minus_d0_m_phoebe("dst_minus_d0_m_phoebe", [&](const Baby &b) {
    return (b.Dst_2010_minus_MM()-b.D0_MM());
  });

  NamedFunc dst_endvxchi2_endvxndof_us("dst_endvxchi2_endvxndof_us", [&](const Baby &b) {
    return b.dst_ENDVERTEX_CHI2()/b.dst_ENDVERTEX_NDOF();
  });

  NamedFunc dst_endvxchi2_endvxndof_phoebe("dst_endvxchi2_endvxndof_phoebe", [&](const Baby &b) {
    return b.Dst_2010_minus_ENDVERTEX_CHI2()/b.Dst_2010_minus_ENDVERTEX_NDOF();
  });

  NamedFunc dst_l0hadtos_or_b0_l0tis_us("dst_l0hadtos_or_b0_l0tis_us", [&](const Baby &b) {
    return (b.dst_L0HadronDecision_TOS() || b.b0_L0Global_TIS());
  });

  NamedFunc dst_l0hadtos_or_b0_l0tis_phoebe("dst_l0hadtos_or_b0_l0tis_phoebe", [&](const Baby &b) {
    return (b.Dst_2010_minus_L0HadronDecision_TOS() || b.Y_L0Global_TIS());
  });

  NamedFunc nSPDhits_phoebe("nSPDhits_phoebe", [&](const Baby &b) {
    return b.nSPDhits();
  });


  ///////////////////// Make Plots

  PlotMaker pm;

  // Fit Variables (only for our scripts)
  pm.Push<Hist1D>(Axis(52,-3,10, "FitVar_Mmiss2/1e6","m_{miss}^{2} [GeV^{2}]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,2.5, "FitVar_El/1e3","E^{*}_{#mu} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(45,-3,12, "FitVar_q2/1e6","q^{2} [GeV^{2}]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");

  // Cut Variables (us)
  pm.Push<Hist1D>(Axis(60,-10,170, "k_PIDK","K PID K"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,300, "k_P/1e3","p_{K} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(40,0,0.08, "k_TRACK_GhostProb","K Track Ghost Prob"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,-170,10, "pi_PIDK","#pi PID K"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(40,0,0.08, "pi_TRACK_GhostProb","#pi Track Ghost Prob"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,300, "pi_P/1e3","p_{#pi} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(40,0,20, "mu_PIDmu","#mu PID #mu"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(2,0,2, "mu_L0Global_TIS","#mu L0 Global TIS"), "1", procs_us_run12, plotboth).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,300, "mu_P/1e3","p_{#mu} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,600, "d0_P/1e3","p_{D^{0}} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,600, "dst_P/1e3","p_{D^{*}} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(60,0,600, "b0_P/1e3","p_{B^{0}} [GeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(40,140,160, dst_minus_d0_m_us,"m_{D^{*}}-m_{D^{0}} (measured) [MeV]"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(40,0,4, "mu_TRACK_CHI2NDOF","#mu Track #chi^{2} ndof"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "dst_PT/1e3","D* p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "d0_PT/1e3","D p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "b0_PT/1e3","B p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "k_PT/1e3","K p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "pi_PT/1e3","#pi p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "mu_PT/1e3","#mu p_{T} (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,200, "b0_FDCHI2_OWNPV","B FD #chi^{2} (OWNPV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,200,600, "d0_FDCHI2_OWNPV","D FD #chi^{2} (OWNPV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,200, "dst_FDCHI2_OWNPV","D* FD #chi^{2} (OWNPV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(65,0,650, "nSPDhits","nSPD Hits"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,1.9,2.1, "dst_MM/1e3","m_{D*} (measured) (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,1.8,1.95, "d0_MM/1e3","m_{D} (measured) (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,5.2,5.4, "b0_MM/1e3","m_{B} (measured) (GeV)"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "b0_IPCHI2_OWNPV","B IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "d0_IPCHI2_OWNPV","D IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "dst_IPCHI2_OWNPV","D* IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "k_IPCHI2_OWNPV","K IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "mu_IPCHI2_OWNPV","#mu IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "pi_IPCHI2_OWNPV","#pi IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,100, "spi_IPCHI2_OWNPV","#pi_{s} IP #chi^{2}"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(36,0,12, dst_endvxchi2_endvxndof_us,"D^{*} ENDVERTEX_CHI2/ENDVERTEX_NDOF"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  // Also plot: endvertex_chi2/endvertex_ndof for B,D; maybe DIRA for B,D*,D ?

  // Others (us)
  pm.Push<Hist1D>(Axis(50,-3,2, "b0_ISOLATION_BDT","B^{0} Iso BDT"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,-3,2, "b0_ISOLATION_BDT2","B^{0} Iso BDT 2"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,-3,2, "b0_ISOLATION_BDT3","B^{0} Iso BDT 3"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,-3,2, "b0_ISOLATION_BDT4","B^{0} Iso BDT 4"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(50,0,5, "b0_ISOLATION_Type","B^{0} Iso Type"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(2,0,2, "b0_L0Global_TIS","B L0 Global TIS"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(2,0,2, "dst_L0HadronDecision_TOS","D^{*} L0 Had TOS"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(2,0,2, dst_l0hadtos_or_b0_l0tis_us, "B L0Global TIS or D* LOHadron TOS"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  pm.Push<Hist1D>(Axis(2,0,2, "mu_isMuon","#mu isMuon"), "1", procs_us_run12, plotnorm).RatioTitle("Run1", "Run2");
  // Also plot: k_isMuon, pi_isMuon ?

  // Note: Purposefully put Phoebe's process first, so filename will take on her (different) variable name-- if variable name the same, need custom NamedFunc...
  // Cut Variables (phoebe)
  pm.Push<Hist1D>(Axis(60,-10,170, {"Kplus_PIDK", "k_PIDK"},"K PID K"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,300, {"Kplus_P/1e3", "k_P/1e3"},"p_{K} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(40,0,0.08, {"Kplus_TRACK_GhostProb", "k_TRACK_GhostProb"},"K Track Ghost Prob"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,-170,10, {"piminus0_PIDK", "pi_PIDK"},"#pi PID K"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(40,0,0.08, {"piminus0_TRACK_GhostProb", "pi_TRACK_GhostProb"},"#pi Track Ghost Prob"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,300, {"piminus0_P/1e3", "pi_P/1e3"},"p_{#pi} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(40,0,20, {"muplus_PIDmu", "mu_PIDmu"},"#mu PID #mu"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(2,0,2, {"muplus_L0Global_TIS", "mu_L0Global_TIS"},"#mu L0 Global TIS"), "1", procs_us_phoebe, plotboth).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,300, {"muplus_P/1e3", "mu_P/1e3"},"p_{#mu} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,600, {"D0_P/1e3", "d0_P/1e3"},"p_{D^{0}} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,600, {"Dst_2010_minus_P/1e3", "dst_P/1e3"},"p_{D^{*}} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(60,0,600, {"Y_P/1e3", "b0_P/1e3"},"p_{B^{0}} [GeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(40,140,160, {dst_minus_d0_m_phoebe, dst_minus_d0_m_us},"m_{D^{*}}-m_{D^{0}} (measured) [MeV]"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(40,0,4, {"muplus_TRACK_CHI2NDOF", "mu_TRACK_CHI2NDOF"},"#mu Track #chi^{2} ndof"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Dst_2010_minus_PT/1e3", "dst_PT/1e3"},"D* p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"D0_PT/1e3", "d0_PT/1e3"},"D p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Y_PT/1e3", "b0_PT/1e3"},"B p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Kplus_PT/1e3", "k_PT/1e3"},"K p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"piminus0_PT/1e3", "pi_PT/1e3"},"#pi p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"muplus_PT/1e3", "mu_PT/1e3"},"#mu p_{T} (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,200, {"Y_FDCHI2_OWNPV", "b0_FDCHI2_OWNPV"},"B FD #chi^{2} (OWNPV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,200,600, {"D0_FDCHI2_OWNPV", "d0_FDCHI2_OWNPV"},"D FD #chi^{2} (OWNPV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,200, {"Dst_2010_minus_FDCHI2_OWNPV", "dst_FDCHI2_OWNPV"},"D* FD #chi^{2} (OWNPV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(65,0,650, {nSPDhits_phoebe, "nSPDhits"},"nSPD Hits"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,1.9,2.1, {"Dst_2010_minus_MM/1e3", "dst_MM/1e3"},"m_{D*} (measured) (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,1.8,1.95, {"D0_MM/1e3", "d0_MM/1e3"},"m_{D} (measured) (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,5.2,5.4, {"Y_MM/1e3", "b0_MM/1e3"},"m_{B} (measured) (GeV)"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Y_IPCHI2_OWNPV", "b0_IPCHI2_OWNPV"},"B IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"D0_IPCHI2_OWNPV", "d0_IPCHI2_OWNPV"},"D IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Dst_2010_minus_IPCHI2_OWNPV", "dst_IPCHI2_OWNPV"},"D* IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"Kplus_IPCHI2_OWNPV", "k_IPCHI2_OWNPV"},"K IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"muplus_IPCHI2_OWNPV", "mu_IPCHI2_OWNPV"},"#mu IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"piminus0_IPCHI2_OWNPV", "pi_IPCHI2_OWNPV"},"#pi IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,100, {"piminus_IPCHI2_OWNPV", "spi_IPCHI2_OWNPV"},"#pi_{s} IP #chi^{2}"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(36,0,12, {dst_endvxchi2_endvxndof_phoebe, dst_endvxchi2_endvxndof_us},"D^{*} ENDVERTEX_CHI2/ENDVERTEX_NDOF"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  // Also plot: endvertex_chi2/endvertex_ndof for B,D; maybe DIRA for B,D*,D ?

  // Others (phoebe)
  pm.Push<Hist1D>(Axis(50,-3,2, {"Y_ISOLATION_BDT","b0_ISOLATION_BDT"},"B^{0} Iso BDT"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,-3,2, {"Y_ISOLATION_BDT2","b0_ISOLATION_BDT2"},"B^{0} Iso BDT 2"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,-3,2, {"Y_ISOLATION_BDT3","b0_ISOLATION_BDT3"},"B^{0} Iso BDT 3"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,-3,2, {"Y_ISOLATION_BDT4","b0_ISOLATION_BDT4"},"B^{0} Iso BDT 4"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(50,0,5, {"Y_ISOLATION_Type","b0_ISOLATION_Type"},"B^{0} Iso Type"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(2,0,2, {"Y_L0Global_TIS", "b0_L0Global_TIS"},"B L0 Global TIS"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(2,0,2, {"Dst_2010_minus_L0HadronDecision_TOS", "dst_L0HadronDecision_TOS"},"D^{*} L0 Had TOS"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(2,0,2, {dst_l0hadtos_or_b0_l0tis_phoebe, dst_l0hadtos_or_b0_l0tis_us}, "B L0Global TIS or D* LOHadron TOS"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  pm.Push<Hist1D>(Axis(2,0,2, {"muplus_isMuon", "mu_isMuon"},"#mu isMuon"), "1", procs_us_phoebe, plotnorm).RatioTitle("Phoebe", "Us Run1");
  // Also plot: k_isMuon, pi_isMuon ?


  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
