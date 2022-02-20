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

  //// User defined colors
  Palette colors("txt/colors.txt", "default");

  ////////////////////////////////////// Comparisons between Yipeng and Phoebe's ntuples ///////////////////////////////////////

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::both);

  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log);
  vector<PlotOpt> plottypes = {lin_lumi};
  vector<PlotOpt> plotboth = {log_lumi,lin_lumi};

  NamedFunc totcut("(b_l0_tis || d0_l0had_tos) && mu_l0_tis && mu_pidb > 0 && dst_id_prod > 0 && id_prod > 0 && flag_sel_d0 && spi_gh_prob < 0.25 && dst_endvtx_chi2/dst_endvtx_ndof < 10 && (dst_m - d0_m - 145.454) < 2 && (dst_m - d0_m - 145.454) > -2 && mu_is_mu && (mu_p > 3 && mu_p < 100) && (mu_eta > 1.7 && mu_eta < 5)&& mu_pid_mu > 2 && mu_pid_e < 1 && mu_ip_chi2 > 45 && mu_gh_prob < 0.5 && b_discard_mu_chi2 <= 6 && b_endvtx_chi2 < 24 && b_endvtx_chi2/b_endvtx_ndof < 6 && b_fd_trans < 7 && b_dira > 0.9995 && b_m < 5280");
  
  string repofolder = "../lhcb-ntuples-gen/gen/";
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_step2>("Yipeng", Process::Type::data, colors("data"),
                                                set<string>({repofolder+"run1-Dst_D0-step2/Dst_D0--20_10_12--std--data--2011--md--step2.root"}),totcut && "iso_bdt < 0.2166 && q2>9000000"));
  procs.push_back(Process::MakeShared<Baby_phoebe_step1>("Phoebe", Process::Type::signal, colors("blue"),
                                                set<string>({repofolder+"run1-Dst-step2/Dst--20_09_16--std--data--2011--md--phoebe-step2.root"}), totcut && "iso_bdt < 0.15 && q2>9"));


  ///////////////////// MAKING PLOTS
  PlotMaker pm;
  
  // Fit variables
  pm.Push<Hist1D>(Axis(52,-3,10, "mm2","m_{miss}^{2} [GeV^{2}]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(50,0,2.5, "el","E^{*}_{#mu} [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(45,-3,12, {"q2/1000000","q2"},"q^{2} [GeV^{2}]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  // Stripping variables
  pm.Push<Hist1D>(Axis(100,0,100, "k_pid_k","unset",{4}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6000, "k_ip_chi2","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "k_gh_prob","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(105,-100,5, "pi_pid_k","unset",{2}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6000, "pi_ip_chi2","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "pi_gh_prob","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(160,1780,1940, "d0_m","unset",{1784.83, 1944.83}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(160,1780,1940, "d0_m","unset",{1784.83, 1944.83}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(84,0,4.2, "d0_endvtx_chi2/d0_endvtx_ndof","D^{0} endvtx_chi2/endvtx_ndof",{4}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,30000, "d0_fd_chi2","unset",{250}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(80,0.99978,1.00001, "d0_dira","unset",{0.9998}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,0.25, "spi_gh_prob","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(90,1920,2100, "dst_m","unset",{1885.26,2135.26}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(120,140,152, "dst_m - d0_m","unset",{160}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(120,0,6,"dst_endvtx_chi2/dst_endvtx_ndof","D*^{+} endvtx_chi2/endvtx_ndof",{100}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");


  pm.Push<Hist1D>(Axis(100,0,5000, "mu_ip_chi2","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "mu_gh_prob","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(140,0,14, "mu_pid_mu","unset",{2}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,2000,8000, "b_m","unset",{10000}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6, "b_endvtx_chi2/b_endvtx_ndof","unset",{6}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0.9995,1, "b_dira","unset",{0.9995}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");


  // Momentum
  pm.Push<Hist1D>(Axis(150,0,150, "k_p/1000", "p(K^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(150,0,150, "pi_p/1000", "p(#pi^{-}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,100, "mu_p/1000", "p(#mu^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");


  // Transverse momentum
  pm.Push<Hist1D>(Axis(48,0,12, "k_pt/1000", "p_{T}(K^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(48,0,12, "pi_pt/1000", "p_{T}(#pi^{-}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,10, "mu_pz/1000", "p_{z}(#mu^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(250,0,25, "d0_pt/1000", "p_{T}(D^{0}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(350,0,35, "b_pt/1000", "p_{T}(D*^{+}#mu) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
