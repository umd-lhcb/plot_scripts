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

  string repofolder = "ntuples/";
  string yipeng = "pre-0.9.0/Dst-std/Dst--19_09_05--std--data--2012--md.root";
  string phoebe = "ref-rdx-run1/Dst-std/Dst--19_09_05--std--data--2012--md--phoebe.root";

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_run1_old>("Yipeng 2012 MD",Process::Type::data, colors("data"),
                                                              set<string>({repofolder+yipeng}), "1"));
  procs.push_back(Process::MakeShared<Baby_phoebe_step1_dst>("Phoebe 2012 MD",Process::Type::signal, colors("blue"),
                                                             set<string>({repofolder+phoebe}), "piminus_TRACK_Type==3"));

  string repofolderDV = "ntuples/../run1-b2D0MuXB2DMuNuForTauMuLine/samples/";
  vector<shared_ptr<Process> > procsDV;
  procsDV.push_back(Process::MakeShared<Baby_run1>("DV45 2012 MD Yipeng",Process::Type::data, colors("data"),
                                                       set<string>({repofolderDV+"Dst--20_08_10--std--data--2012--md--dv45-subset-yipeng.root"}), "1"));
  procsDV.push_back(Process::MakeShared<Baby_run1_old_phoebe>("DV36 2012 MD Phoebe",Process::Type::signal, colors("blue"),
                                                       set<string>({repofolderDV+"Dst--20_08_07--std--data--2012--md--dv36-subset-phoebe.root"}), "1"));


  // Custom NamedFunc
  NamedFunc mu_eta("mu_eta", [&](const Baby &b){
       return log((b.muplus_P()+b.muplus_PZ())/(b.muplus_P()-b.muplus_PZ()))/2.;
  });
  NamedFunc muk_log("muk_log", [&](const Baby &b){
       return log10(1-(b.muplus_PX()*b.Kplus_PX() + b.muplus_PY()*b.Kplus_PY() + b.muplus_PZ()*b.Kplus_PZ())/(b.muplus_P()*b.Kplus_P()));
  });
  NamedFunc mupi_log("mupi_log", [&](const Baby &b){
       return log10(1-(b.muplus_PX()*b.piminus0_PX() + b.muplus_PY()*b.piminus0_PY() + b.muplus_PZ()*b.piminus0_PZ())/(b.muplus_P()*b.piminus0_P()));
  });
  NamedFunc muspi_log("muspi_log", [&](const Baby &b){
       return log10(1-(b.muplus_PX()*b.piminus_PX() + b.muplus_PY()*b.piminus_PY() + b.muplus_PZ()*b.piminus_PZ())/(b.muplus_P()*b.piminus_P()));
  });
  NamedFunc log_ip("log_ip", [&](const Baby &b){
                               return log(b.D0_IP_OWNPV());
  });
  NamedFunc b0_dxy("b0_dxy", [&](const Baby &b){
                               return b.Y_FD_OWNPV()*sin(b.Y_FlightDir_Zangle());
  });

  ///////// Automatically appending cutflow cuts
  vector<NamedFunc> cuts = {"1", "(Kplus_PT > 800) && (!Kplus_isMuon) && Kplus_IPCHI2_OWNPV > 45",
                            "(piminus0_PT > 800) && (!piminus0_isMuon) && piminus0_IPCHI2_OWNPV > 45",
                            "D0_P>2000 && D0_FDCHI2_OWNPV > 250 && (D0_MM-1864.83) < 23.4 && (D0_MM-1864.83) > -23.4 && (Kplus_PT>1700 || piminus0_PT>1700) && D0_IPCHI2_OWNPV > 9"
                            && log_ip > -3.5,
                            "muplus_isMuon && muplus_PIDmu > 2 && muplus_PIDe < 1 && muplus_P < 100000"
                            && mu_eta > 1.7 && mu_eta < 5 && muk_log>-6.5 && mupi_log>-6.5 && muspi_log>-6.5,
                            "piminus_TRACK_GhostProb < 0.5 && (Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF) < 10 && (Dst_2010_minus_MM - D0_MM-145.43) < 2 &&  (Dst_2010_minus_MM - D0_MM-145.43) > -2",
                            "Y_ISOLATION_BDT < 0.15 && (Y_ENDVERTEX_CHI2/Y_ENDVERTEX_NDOF) < 6 && Y_MM<5280 && Y_DIRA_OWNPV>0.9995"};


  ///////////////////// MAKING PLOTS
  PlotMaker pm;

  // Comparison with DaVinci 36 and Yipeng's script
  pm.Push<Hist1D>(Axis(48,140,152, "Dst_2010_minus_M - D0_M","unset",{160}), "1", procsDV, plottypes).RatioTitle("DV42","DV36").Tag("dv36");
  pm.Push<Hist1D>(Axis(48,140,152, "Dst_2010_minus_MM - D0_MM","unset",{160}), "1", procsDV, plottypes).RatioTitle("DV42","DV36").Tag("dv36");
  pm.Push<Hist1D>(Axis(35,-25,10, "Dst_2010_minus_ENDVERTEX_X"), "1", procsDV, plotboth).RatioTitle("Yipeng","Phoebe").Tag("dv36");
  pm.Push<Hist1D>(Axis(25,0,6,"Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF","D*^{+} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{100}), "1", procsDV, plottypes).RatioTitle("Yipeng","Phoebe").Tag("dv36");

  pm.Push<Hist1D>(Axis(48,140,152, "Dst_2010_minus_M - D0_M","unset",{160}), "1", procs, plottypes).RatioTitle("DV42","DV36").Tag("coarse");
  pm.Push<Hist1D>(Axis(48,140,152, "Dst_2010_minus_MM - D0_MM","unset",{160}), "1", procs, plottypes).RatioTitle("DV42","DV36").Tag("coarse");
  pm.Push<Hist1D>(Axis(35,-25,10, "Dst_2010_minus_ENDVERTEX_X"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe").Tag("coarse");
  pm.Push<Hist1D>(Axis(25,0,6,"Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF","D*^{+} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{100}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe").Tag("coarse");

  // Stripping variables
  pm.Push<Hist1D>(Axis(100,0,100, "Kplus_PIDK","unset",{4}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6000, "Kplus_IPCHI2_OWNPV","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "Kplus_TRACK_GhostProb","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(105,-100,5, "piminus0_PIDK","unset",{2}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6000, "piminus0_IPCHI2_OWNPV","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "piminus0_TRACK_GhostProb","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(160,1780,1940, "D0_M","unset",{1784.83, 1944.83}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(160,1780,1940, "D0_MM","unset",{1784.83, 1944.83}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(84,0,4.2, "D0_ENDVERTEX_CHI2/D0_ENDVERTEX_NDOF","D^{0} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{4}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,30000, "D0_FDCHI2_OWNPV","unset",{250}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(80,0.99978,1.00001, "D0_DIRA_OWNPV","unset",{0.9998}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,3, "piminus_TRACK_CHI2NDOF","unset",{3}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,200, "piminus_IPCHI2_OWNPV","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.25, "piminus_TRACK_GhostProb","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(90,1920,2100, "Dst_2010_minus_M","unset",{1885.26,2135.26}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(120,140,152, "Dst_2010_minus_M - D0_M","unset",{160}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(120,140,152, "Dst_2010_minus_MM - D0_MM","unset",{160}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(120,0,6,"Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF","D*^{+} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{100}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(200,-100,100, "Dst_2010_minus_ENDVERTEX_X"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(200,-100,100, "Dst_2010_minus_ENDVERTEX_Y"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(200,-1500,1000, "Dst_2010_minus_ENDVERTEX_Z"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(25,-1,1, "Dst_2010_minus_DIRA_ORIVX"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,25, "Dst_2010_minus_MMERR"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,1000, "Dst_2010_minus_FD_OWNPV"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(25,-1,1, "Dst_2010_minus_DIRA_OWNPV"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,1500, "Dst_2010_minus_FD_ORIVX"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,-15,15, "Dst_2010_minus_ORIVX_X"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,-15,15, "Dst_2010_minus_ORIVX_Y"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,-200,450, "Dst_2010_minus_ORIVX_Z"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0.35,0.85, "Dst_2010_minus_OWNPV_X"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,-200,200, "Dst_2010_minus_OWNPV_Z"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6, "Dst_2010_minus_IP_OWNPV"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,-50,50, "Dst_2010_minus_OWNPV_Z - D0_OWNPV_Z"), "1", procs, plotboth).RatioTitle("Yipeng","Phoebe");


  pm.Push<Hist1D>(Axis(100,0,5000, "muplus_IPCHI2_OWNPV","unset",{45}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,0.2, "muplus_TRACK_GhostProb","unset",{0.5}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(140,0,14, "muplus_PIDmu","unset",{2}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,3, "muplus_TRACK_CHI2NDOF","unset",{3}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,2000,8000, "Y_M","unset",{10000}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,2000,8000, "Y_MM","unset",{10000}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0,6, "Y_ENDVERTEX_CHI2/Y_ENDVERTEX_NDOF","unset",{6}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(100,0.9995,1, "Y_DIRA_OWNPV","unset",{0.9995}), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");


  // Momentum
  pm.Push<Hist1D>(Axis(150,0,150, "Kplus_P/1000", "p(K^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(150,0,150, "piminus0_P/1000", "p(#pi^{-}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(60,0,30, "piminus_P/1000", "p(#pi^{-}_{slow}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,100, "muplus_P/1000", "p(#mu^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(250,0,250, "D0_P/1000", "p(D^{0}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(250,0,250, "Dst_2010_minus_P/1000", "p(D*^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(350,0,350, "Y_P/1000", "p(D*^{+}#mu) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");


  // Transverse momentum
  pm.Push<Hist1D>(Axis(48,0,12, "Kplus_PT/1000", "p_{T}(K^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(48,0,12, "piminus0_PT/1000", "p_{T}(#pi^{-}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(50,0,2, "piminus_PT/1000", "p_{T}(#pi^{-}_{slow}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");

  pm.Push<Hist1D>(Axis(100,0,10, "muplus_PT/1000", "p_{T}(#mu^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(250,0,25, "D0_PT/1000", "p_{T}(D^{0}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(250,0,25, "Dst_2010_minus_PT/1000", "p_{T}(D*^{+}) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");
  pm.Push<Hist1D>(Axis(350,0,35, "Y_PT/1000", "p_{T}(D*^{+}#mu) [GeV]"), "1", procs, plottypes).RatioTitle("Yipeng","Phoebe");



  // Impact of cuts
  NamedFunc fullcut = "1";
  for(size_t ind = 0; ind < cuts.size(); ind++) {
    fullcut = fullcut && cuts[ind];
    pm.Push<Hist1D>(Axis(50,0,250, "D0_P/1000", "p(D^{0}) [GeV]"), fullcut, procs, plottypes).RatioTitle("Yipeng","Phoebe").Tag("cuts");
    pm.Push<Hist1D>(Axis(120,140,152, "Dst_2010_minus_M - D0_M","unset",{160}), fullcut, procs, plottypes).RatioTitle("Yipeng","Phoebe").Tag("cuts");
    pm.Push<Hist1D>(Axis(120,0,6,"Dst_2010_minus_ENDVERTEX_CHI2/Dst_2010_minus_ENDVERTEX_NDOF","D*^{+} ENDVERTEX_CHI2/ENDVERTEX_NDOF",{100}),
                    fullcut, procs, plottypes).RatioTitle("Yipeng","Phoebe").Tag("cuts");
  }

  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
