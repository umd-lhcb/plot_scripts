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

  //// Typically you would only have one PlotMaker so that you do not run on the same ntuples
  //// several times, but separated them here to make the examples modular and more clear
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// 1D plot mmiss sig/norm/D** ///////////////////////////////////////

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::data)
    .Bottom(BottomType::pull)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::none);
  PlotOpt log_lumi = lin_lumi().YAxis(YAxisType::log).Title(TitleType::preliminary).Bottom(BottomType::ratio).Overflow(OverflowType::both);
  PlotOpt lin_shapes = lin_lumi().Stack(StackType::shapes).Title(TitleType::info).Bottom(BottomType::off);
  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes);
  
  vector<PlotOpt> plottypes = {lin_lumi, log_lumi, lin_shapes, lin_lumi_shapes};

  //////// NamedFuncs to select signal, normalization, and D**
  NamedFunc is_dsptau("is_dsptau",[&](const Baby &b){
      return (abs(b.mu_MC_MOTHER_ID())==15 && abs(b.d0_MC_MOTHER_ID())==413 && abs(b.d0_MC_GD_MOTHER_ID())==511);
  });
  NamedFunc is_dspmu("is_dspmu",[&](const Baby &b){
      return (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==413 && abs(b.d0_MC_GD_MOTHER_ID())==511);
  });
  NamedFunc is_dss("is_dss",[&](const Baby &b){
    return ((abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==10411 && abs(b.d0_MC_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==10413 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==20413 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511) |
            (abs(b.mu_MC_MOTHER_ID())==511 && ((abs(b.d0_MC_MOTHER_ID())==415 && abs(b.d0_MC_GD_MOTHER_ID())==511) |
                                              (abs(b.d0_MC_GD_MOTHER_ID())==415 && abs(b.d0_MC_GD_GD_MOTHER_ID())==511))));
  });
  //////// Signal, normalization, D** processes for mmiss2 plot
  string repofolder = "ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  vector<shared_ptr<Process> > procs_mm;
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("Data (actually MC)",Process::Type::data, colors("data"),
                                                      set<string>({repofolder+run2bare}), "1"));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D*^{+} #tau #nu", Process::Type::signal, colors("green"),
                                                      set<string>({repofolder+run2bare}), is_dsptau));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D*^{+} #mu #nu", Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), is_dspmu));
  procs_mm.push_back(Process::MakeShared<Baby_run2_bare>("B #rightarrow D** #mu #nu", Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), is_dss));

  PlotMaker pm_mm;
  // Missing mass (in GeV^2, so need to divide by 1e6)
  pm_mm.Push<Hist1D>(Axis(75, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>8",
                  procs_mm, plottypes);
  pm_mm.Push<Hist1D>(Axis(75, -5, 10,"FitVar_Mmiss2/1000000", "m_{miss}^{2} [GeV^{2}]"), "FitVar_q2/1000000>8",
                  procs_mm, plottypes).Tag("example").TopRight("#font[82]{TopRight} label").RatioTitle("Fake data","MC");
  pm_mm.Push<Table>("pie", vector<TableRow>{
      TableRow("All events","1", 0,1, "1"),
        TableRow("$m_\\text{miss}^2 > 3\\text{ GeV}^2$",  "FitVar_Mmiss2/1000000 > 3",0,0, "1"),
        TableRow("BDT$_\\text{iso}<0.15$",  "FitVar_Mmiss2/1000000 > 3 && b0_ISOLATION_BDT < 0.15",0,0, "1"),
        },procs_mm, true, true, true, true, true, true).Precision(1); // Pushing table and pie charts
  pm_mm.Push<EventScan>("eventscan", !is_dsptau && !is_dspmu && !is_dss && "FitVar_Mmiss2/1000000 > 8", 
                        vector<NamedFunc>{"runNumber", "eventNumber", "mu_MC_MOTHER_ID", "d0_MC_MOTHER_ID",
                                            "d0_MC_GD_MOTHER_ID", "d0_MC_GD_GD_MOTHER_ID"}, procs_mm).Precision(10);


  pm_mm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to  

  ////////////////////////////////////// 1D plot mmiss sig/norm/D** ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// Slow pion pT and resolutions /////////////////////////////////////

  PlotOpt lin_lumi_spi = lin_lumi().Bottom(BottomType::off).Title(TitleType::simulation);
  PlotOpt log_lumi_spi = lin_lumi_spi().YAxis(YAxisType::log);
  vector<PlotOpt> plottypes_spi = {lin_lumi_spi, log_lumi_spi};

  //////// Processes for slow pion plots
  NamedFunc spi_mom_is_dsp("spi_mom_is_dsp",
                          [&](const Baby &b){
                            return ((b.spi_MC_MOTHER_ID()==413) || (b.spi_MC_MOTHER_ID()==-413));
  });
  vector<shared_ptr<Process> > procs_low_spi;
  procs_low_spi.push_back(Process::MakeShared<Baby_run2_bare>("250 < p_{T}^{reco}(#pi_{slow}) < 500 MeV",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), "spi_PT>250 && spi_PT<500 && spi_TRUEPT>0"));

  vector<shared_ptr<Process> > procs_comp_spi;
  procs_comp_spi.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(#pi_{slow}) > 300 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), "spi_TRUEPT>300"));
  procs_comp_spi.push_back(Process::MakeShared<Baby_run2_bare>("275 < p_{T}^{true}(#pi_{slow}) < 300 MeV",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), "spi_TRUEPT>275&&spi_TRUEPT<300"));
  procs_comp_spi.push_back(Process::MakeShared<Baby_run2_bare>("250 < p_{T}^{true}(#pi_{slow}) < 275 MeV",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), "spi_TRUEPT>250&&spi_TRUEPT<275"));
  procs_comp_spi.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(#pi_{slow}) < 250 MeV",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "spi_TRUEPT<250"));

 PlotMaker pm_spi;
  // Slow pion momentum and resolution
  pm_spi.Push<Hist1D>(Axis(100, 0, 500,"spi_PT", "p_{T}^{reco}(#pi_{slow}) [MeV]",{300}), "spi_TRUEPT>0",
                  procs_comp_spi,plottypes_spi).TopRight("13 TeV").RightLabel({"Slow pion p_{T}"});
  pm_spi.Push<Hist1D>(Axis(50, -0.5, 0.5,"(spi_TRUEPT-spi_PT)/spi_PT", "(p_{T}^{true}(#pi_{slow}) - p_{T}^{reco}(#pi_{slow}))/p_{T}^{reco}(#pi_{slow})",{-25/300., 50/300.}),
                  "1", procs_low_spi, plottypes_spi).TopRight("13 TeV").LeftLabel({"Slow pion", "p_{T} resolution"});
  pm_spi.MakePlots(1);  // The "1" is the luminosity to rescale the bkg to  

  ////////////////////////////////////// Slow pion pT and resolutions /////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// 2D scatter plot mu angles ////////////////////////////////////////

  PlotOpt style("txt/plot_styles.txt", "Scatter");
  vector<PlotOpt> scattertype = {style().Stack(StackType::signal_on_top).Title(TitleType::data)};
  //////// Processes for scatter plot
  auto mup_high = Process::MakeShared<Baby_run2_bare>("p^{reco}(#mu) > 100 GeV", Process::Type::data, kBlack,
                                                      set<string>({repofolder+run2bare}), "mu_P>100000");
  auto mupt_high = Process::MakeShared<Baby_run2_bare>("p^{reco}_{T}(#mu) > 8 GeV", Process::Type::data, kBlue,
                                                       set<string>({repofolder+run2bare}), "mu_PT>8000");
  auto all_mu = Process::MakeShared<Baby_run2_bare>("MC", Process::Type::background, kBlack,
                                                    set<string>({repofolder+run2bare}), "1");
  mup_high->SetMarkerStyle(20); mup_high->SetMarkerSize(0.4);
  mupt_high->SetMarkerStyle(21);mupt_high->SetMarkerSize(0.4);
  
  vector<shared_ptr<Process> > procs_mu = {mup_high, mupt_high, all_mu};

  PlotMaker pm_mu;
  pm_mu.Push<Hist2D>(Axis(55, -0.6, 0.5, "mu_TRUEP_X/mu_TRUEP_Z", "p_{x}^{true}/p_{z}^{true}(#mu)", {-0.38, 0.38}),
                  Axis(38, -0.38, 0.38, "mu_TRUEP_Y/mu_TRUEP_Z", "p_{y}^{true}/p_{z}^{true}(#mu)", {-0.28, 0.28}),
                  "1", procs_mu, scattertype).TopRight("");
  pm_mu.MakePlots(1);  // The "1" is the luminosity to rescale the bkg to  

  ////////////////////////////////////// 2D scatter plot mu angles ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

