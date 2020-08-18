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

  // Comparison between refit D* only and full refit
  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::both)
    .FileExtensions({"png"});

  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log);
  vector<PlotOpt> plottypes = {lin_lumi};
  vector<PlotOpt> plotboth  = {log_lumi,lin_lumi};

  string repofolder = "ntuples/0.9.1-dst_partial_refit/Dst_D0-cutflow_mc/";
  vector<shared_ptr<Process> > procs;

  procs.push_back(Process::MakeShared<Baby_rdx_std_step1>(
        "Run 2 cocktail: full refit", Process::Type::data, colors("data"),
        set<string>({repofolder+"Dst_D0--20_08_18--cutflow_mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), "1"));
  procs.push_back(Process::MakeShared<Baby_rdx_std_step1>(
        "Run 2 cocktail: refit D* only", Process::Type::signal, colors("blue"),
        set<string>({repofolder+"Dst_D0--20_08_18--cutflow_mc--refit_dst_only--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root"}), "1"));


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
       return log10(1-(b.mu_PX()*b.spi_PX() + b.mu_PY()*b.pi_PY() + b.mu_PZ()*b.pi_PZ())/(b.mu_P()*b.pi_P()));
  });
  NamedFunc log_ip("log_ip", [&](const Baby &b){
                               return log(b.d0_IP_OWNPV());
  });
  NamedFunc b0_dxy("b0_dxy", [&](const Baby &b){
                               return b.b0_FD_OWNPV()*sin(b.b0_FlightDir_Zangle());
  });

  // Automatically appending cutflow cuts
  vector<NamedFunc> cuts = {"1", "(k_PT > 800) && (!k_isMuon) && k_IPCHI2_OWNPV > 45",
                            "(pi_PT > 800) && (!pi_isMuon) && pi_IPCHI2_OWNPV > 45",
                            "d0_P>2000 && d0_FDCHI2_OWNPV > 250 && (d0_MM-1864.83) < 23.4 && (d0_MM-1864.83) > -23.4 && (k_PT>1700 || pi_PT>1700) && d0_IPCHI2_OWNPV > 9"
                            && log_ip > -3.5,
                            "mu_isMuon && mu_PIDmu > 2 && mu_PIDe < 1 && mu_P < 100000"
                            && mu_eta > 1.7 && mu_eta < 5 && muk_log>-6.5 && mupi_log>-6.5 && muspi_log>-6.5,
                            "spi_TRACK_GhostProb < 0.5 && (dst_ENDVERTEX_CHI2/dst_ENDVERTEX_NDOF) < 10 && (dst_MM - d0_MM-145.43) < 2 &&  (dst_MM - d0_MM-145.43) > -2",
                            "b0_ISOLATION_BDT < 0.15 && (b0_ENDVERTEX_CHI2/b0_ENDVERTEX_NDOF) < 6 && b0_MM<5280 && b0_DIRA_OWNPV>0.9995"};


  // MAKING PLOTS
  PlotMaker pm;

  // Comparison with DaVinci 36 and Yipeng's script
  pm.Push<Hist1D>(Axis(48,140,152, "dst_M - d0_M","unset",{160}), "1", procs, plottypes).RatioTitle("Full refit", "D* only");
  pm.Push<Hist1D>(Axis(48,140,152, "dst_MM - d0_MM","unset",{160}), "1", procs, plottypes).RatioTitle("Full refit","D* only");
  pm.Push<Hist1D>(Axis(35,-25,10, "dst_ENDVERTEX_X"), "1", procs, plotboth).RatioTitle("Full refit", "D* only");
  pm.Push<Hist1D>(Axis(25,0,6,"dst_ENDVERTEX_CHI2/dst_ENDVERTEX_NDOF","D*^{+} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{100}), "1", procs, plottypes).RatioTitle("Full refit", "D* only");
  pm.Push<Hist1D>(Axis(100,0,1500, "dst_FD_ORIVX"), "1", procs, plotboth).RatioTitle("Full refit","D* only");

  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
