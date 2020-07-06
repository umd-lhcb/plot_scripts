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
  log_lumi.Title(TitleType::simulation)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm)
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
  
  vector<PlotOpt> linplot = {lin_lumi, log_lumi};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string repofolder = "/Users/manuelf/code/lhcb-ntuples-gen/ntuples/";
  string run2bare = "0.9.0-cutflow/Dst-cutflow_mc/Dst--20_06_05--cutflow_mc--bare--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09b_Trig0x6138160F_Reco16_Turbo03_Stripping26NoPrescalingFlagged_11874091_ALLSTREAMS.DST.root";
  
  PlotMaker pm;

  /////////////////////////// Cut on pi, K pT ////////////////////////////////

  NamedFunc spi_mom_is_dsp("spi_mom_is_dsp",
                          [&](const Baby &b){
                            return ((b.spi_MC_MOTHER_ID()==413) || (b.spi_MC_MOTHER_ID()==-413));
  });
  
  vector<shared_ptr<Process> > procs_all_spi;
  procs_all_spi.push_back(Process::MakeShared<Baby_run2_bare>("Non-truthmatched #pi_{slow}",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "spi_PT>250 && spi_PT<500"));

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

  pm.Push<Hist1D>(Axis(100, 0, 500,"spi_PT", "p_{T}^{reco}(#pi_{slow}) [MeV]",{300}), "spi_TRUEPT>0", procs_comp_spi, linplot).TopRight("13 TeV");
  //pm.Push<Hist1D>(Axis(100, 0, 500,"spi_PT", "p_{T}^{reco}(#pi_{slow}) [MeV]",{300}), "spi_TRUEPT==0", procs_all_spi, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(60, -0.5, 0.5,"(spi_TRUEPT-spi_PT)/spi_PT", "(p_{T}^{true}(#pi_{slow}) - p_{T}^{reco}(#pi_{slow}))/p_{T}^{reco}(#pi_{slow})",{-25/300.}),
                  "1", procs_low_spi, linplot).TopRight("13 TeV");

  pm.Push<Table>("spi", vector<TableRow>{
      TableRow("$250 < p_{T}^\\text{reco}(\\pi_\\text{slow}) < $500 MeV","1", 0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 25 MeV/300 MeV",  "(spi_TRUEPT-spi_PT)/spi_PT < -25/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 50 MeV/300 MeV",  "(spi_TRUEPT-spi_PT)/spi_PT < -50/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 75 MeV/300 MeV",  "(spi_TRUEPT-spi_PT)/spi_PT < -75/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 100 MeV/300 MeV", "(spi_TRUEPT-spi_PT)/spi_PT < -100/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 150 MeV/300 MeV", "(spi_TRUEPT-spi_PT)/spi_PT < -150/300.",0,0, "1"),
    },procs_low_spi,0); // Pushing table

  ///////////// Kaon
  vector<shared_ptr<Process> > procs_low_k;
  procs_low_k.push_back(Process::MakeShared<Baby_run2_bare>("250 < p_{T}^{reco}(K) < 500 MeV",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), "k_PT>250 && k_PT<500"));

  vector<shared_ptr<Process> > procs_comp_k;
  procs_comp_k.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(K) > 300 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT>300"));
  procs_comp_k.push_back(Process::MakeShared<Baby_run2_bare>("275 < p_{T}^{true}(K) < 300 MeV",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT>275&&k_TRUEPT<300"));
  procs_comp_k.push_back(Process::MakeShared<Baby_run2_bare>("250 < p_{T}^{true}(K) < 275 MeV",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT>250&&k_TRUEPT<275"));
  procs_comp_k.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(K) < 250 MeV",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT<250"));

  pm.Push<Hist1D>(Axis(100, 200, 700,"k_PT", "p_{T}^{reco}(K) [MeV]",{300}), "1", procs_comp_k, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(60, -0.5, 0.5,"(k_TRUEPT-k_PT)/k_PT", "(p_{T}^{true}(K) - p_{T}^{reco}(K))/p_{T}^{reco}(K)",{-25/300.}),
                  "1", procs_low_k, linplot).TopRight("13 TeV");

  pm.Push<Table>("k", vector<TableRow>{
      TableRow("$250 < p_{T}^\\text{reco}(K) < $500 MeV","1", 0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 25 MeV/300 MeV",  "(k_TRUEPT-k_PT)/k_PT < -25/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 50 MeV/300 MeV",  "(k_TRUEPT-k_PT)/k_PT < -50/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 75 MeV/300 MeV",  "(k_TRUEPT-k_PT)/k_PT < -75/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 100 MeV/300 MeV", "(k_TRUEPT-k_PT)/k_PT < -100/300.",0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 150 MeV/300 MeV", "(k_TRUEPT-k_PT)/k_PT < -150/300.",0,0, "1"),
    },procs_low_k,0); // Pushing table


  ////////////////////////// Kaon + Pion pT//////////////////////////
  NamedFunc kpi_mom_is_d0("kpi_mom_is_d0",
                          [&](const Baby &b){
                            return (abs(b.pi_MC_MOTHER_ID())==421 && abs(b.k_MC_MOTHER_ID())==421);
  });
  NamedFunc kpi_pt_resol("kpi_pt_resol",
                          [&](const Baby &b){
                            double kpi_pt = b.k_PT()+b.pi_PT(), kpi_pt_tru = b.k_TRUEPT()+b.pi_TRUEPT();
                            return (kpi_pt_tru - kpi_pt)/kpi_pt;
  });
  vector<shared_ptr<Process> > procs_low_kpi;
  procs_low_kpi.push_back(Process::MakeShared<Baby_run2_bare>("2000 < p_{T}^{reco}(K)+p_{T}^{reco}(#pi) < 3000 MeV",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), "k_PT+pi_PT>2000 && k_PT+pi_PT<3000" && kpi_mom_is_d0));

  vector<shared_ptr<Process> > procs_comp_kpi;
  procs_comp_kpi.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(K)+p_{T}^{true}(#pi) > 2500 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT+pi_TRUEPT>2500"));
  procs_comp_kpi.push_back(Process::MakeShared<Baby_run2_bare>("2450 < p_{T}^{true}(K)+p_{T}^{true}(#pi) < 2500 MeV",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT+pi_TRUEPT>2450&&k_TRUEPT+pi_TRUEPT<2500"));
  procs_comp_kpi.push_back(Process::MakeShared<Baby_run2_bare>("2400 < p_{T}^{true}(K)+p_{T}^{true}(#pi) < 2450 MeV",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT+pi_TRUEPT>2400&&k_TRUEPT+pi_TRUEPT<2450"));
  procs_comp_kpi.push_back(Process::MakeShared<Baby_run2_bare>("p_{T}^{true}(K)+p_{T}^{true}(#pi) < 2400 MeV",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "k_TRUEPT+pi_TRUEPT<2400"));

  pm.Push<Hist1D>(Axis(40, 2000, 3000,"k_PT+pi_PT", "p_{T}^{reco}(K)+p_{T}^{reco}(#pi) [MeV]",{2500}), kpi_mom_is_d0, procs_comp_kpi, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(160, 0, 4000,"k_PT+pi_PT", "p_{T}^{reco}(K)+p_{T}^{reco}(#pi) [MeV]",{2500}), kpi_mom_is_d0, procs_comp_kpi, linplot).TopRight("13 TeV").Tag("full");
  pm.Push<Hist1D>(Axis(50, -0.25, 0.25,"(k_TRUEPT+pi_TRUEPT-k_PT-pi_PT)/(k_PT+pi_PT)", "p_{T}(K)+p_{T}(#pi) resolution (true-reco)/reco",{-50/2500.}),
                  kpi_mom_is_d0, procs_low_kpi, linplot).TopRight("13 TeV");

  pm.Push<Table>("kpi", vector<TableRow>{
      TableRow("$2000 < p_{T}^\\text{reco}(K)+p_{T}^\\text{reco}(\\pi) < $3000 MeV","1", 0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 25 MeV/2500 MeV", kpi_pt_resol < -25/2500.,0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 50 MeV/2500 MeV", kpi_pt_resol < -50/2500.,0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 75 MeV/2500 MeV", kpi_pt_resol < -75/2500.,0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 100 MeV/2500 MeV", kpi_pt_resol < -100/2500.,0,0, "1"),
        TableRow("$p_{T}^{res}$ worse than 150 MeV/2500 MeV", kpi_pt_resol < -150/2500.,0,0, "1"),
    },procs_low_kpi,0); // Pushing table
  

  /////////////////////////// Cut on muon P ////////////////////////////////
  vector<shared_ptr<Process> > procs_low_mu;
  procs_low_mu.push_back(Process::MakeShared<Baby_run2_bare>("3000 < p^{reco}(#mu) < 4000 MeV",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), "mu_P<4000"));

  vector<shared_ptr<Process> > procs_comp_mu;
  procs_comp_mu.push_back(Process::MakeShared<Baby_run2_bare>("p^{true}(#mu) > 3000 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), "mu_TRUEP_E>3000"));
  procs_comp_mu.push_back(Process::MakeShared<Baby_run2_bare>("2950 < p^{true}(#mu) < 3000 MeV",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), "mu_TRUEP_E>2950&&mu_TRUEP_E<3000"));
  procs_comp_mu.push_back(Process::MakeShared<Baby_run2_bare>("2925 < p^{true}(#mu) < 2950 MeV",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), "mu_TRUEP_E>2925&&mu_TRUEP_E<2950"));
  procs_comp_mu.push_back(Process::MakeShared<Baby_run2_bare>("p^{true}(#mu) < 2925 MeV",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), "mu_TRUEP_E<2925"));


  pm.Push<Hist1D>(Axis(50, 2500, 5000,"mu_P", "p^{reco}(#mu) [MeV]",{3000.}), "1", procs_comp_mu, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(250, 0, 25000,"mu_P", "p^{reco}(#mu) [MeV]",{3000.}), "1", procs_comp_mu, linplot).TopRight("13 TeV").Tag("mufull");
  pm.Push<Hist1D>(Axis(60, -0.25, 0.25,"(mu_TRUEP_E-mu_P)/mu_P", "(p^{true}(#mu) - p^{reco}(#mu))/p^{reco}(#mu)",{-50/3000.}),
                  "1", procs_low_mu, linplot).TopRight("13 TeV");

  pm.Push<Table>("mu", vector<TableRow>{
      TableRow("3000 < p^\\text{reco}(\\mu) < $4000 MeV","1", 0,0, "1"),
        TableRow("$p^{res}$ worse than 25 MeV/3000 MeV",  "(mu_TRUEP_E-mu_P)/mu_P < -25/3000.",0,0, "1"),
        TableRow("$p^{res}$ worse than 50 MeV/3000 MeV",  "(mu_TRUEP_E-mu_P)/mu_P < -50/3000.",0,0, "1"),
        TableRow("$p^{res}$ worse than 75 MeV/3000 MeV",  "(mu_TRUEP_E-mu_P)/mu_P < -75/3000.",0,0, "1"),
        TableRow("$p^{res}$ worse than 100 MeV/3000 MeV", "(mu_TRUEP_E-mu_P)/mu_P < -100/3000.",0,0, "1"),
    },procs_low_mu,0); // Pushing table


  /////////////////////////// Cut on muon angle ////////////////////////////////
  NamedFunc mu_thetax("mu_thetax",
                          [&](const Baby &b){
                            return fabs(b.mu_PX()/b.mu_PZ());
  });
  NamedFunc mu_thetax_tru("mu_thetax_tru",
                          [&](const Baby &b){
                            return fabs(b.mu_TRUEP_X()/b.mu_TRUEP_Z());
  });
  NamedFunc mu_thetay("mu_thetay",
                          [&](const Baby &b){
                            return fabs(b.mu_TRUEP_Y()/b.mu_TRUEP_Z());
  });
  NamedFunc mu_thetay_tru("mu_thetay_tru",
                          [&](const Baby &b){
                            return fabs(b.mu_TRUEP_Y()/b.mu_TRUEP_Z());
  });

  vector<shared_ptr<Process> > procs_edgex_mu;
  procs_edgex_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{x}^{reco}(#mu) > 0.25",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), mu_thetax>0.25));

  vector<shared_ptr<Process> > procs_comp_thetax_mu;
  procs_comp_thetax_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{x}^{true}(#mu) < 0.36 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), mu_thetax_tru<0.36));
  procs_comp_thetax_mu.push_back(Process::MakeShared<Baby_run2_bare>("0.36 < #theta_{x}^{true}(#mu) < 0.37",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), mu_thetax_tru>0.36 && mu_thetax_tru<0.37));
  procs_comp_thetax_mu.push_back(Process::MakeShared<Baby_run2_bare>("0.37 < #theta_{x}^{true}(#mu) < 0.38",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), mu_thetax_tru>0.37 && mu_thetax_tru<0.38));
  procs_comp_thetax_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{x}^{true}(#mu) > 0.38",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), mu_thetax_tru>0.38));

  vector<shared_ptr<Process> > procs_edgey_mu;
  procs_edgey_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{y}^{reco}(#mu) > 0.25",
                                                      Process::Type::background, colors("purple"),
                                                      set<string>({repofolder+run2bare}), mu_thetay>0.25));

  vector<shared_ptr<Process> > procs_comp_thetay_mu;
  procs_comp_thetay_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{y}^{true}(#mu) < 0.26 MeV",
                                                      Process::Type::background, colors("green"),
                                                      set<string>({repofolder+run2bare}), mu_thetay_tru<0.26));
  procs_comp_thetay_mu.push_back(Process::MakeShared<Baby_run2_bare>("0.26 < #theta_{y}^{true}(#mu) < 0.27",
                                                      Process::Type::background, colors("blue"),
                                                      set<string>({repofolder+run2bare}), mu_thetay_tru>0.26 && mu_thetay_tru<0.27));
  procs_comp_thetay_mu.push_back(Process::MakeShared<Baby_run2_bare>("0.27 < #theta_{y}^{true}(#mu) < 0.28",
                                                      Process::Type::background, colors("yellow"),
                                                      set<string>({repofolder+run2bare}), mu_thetay_tru>0.27 && mu_thetay_tru<0.28));
  procs_comp_thetay_mu.push_back(Process::MakeShared<Baby_run2_bare>("#theta_{y}^{true}(#mu) > 0.28",
                                                      Process::Type::background, colors("red"),
                                                      set<string>({repofolder+run2bare}), mu_thetay_tru>0.28));


  NamedFunc mu_eta("mu_eta", [&](const Baby &b){
       return log((b.muplus_P()+b.muplus_PZ())/(b.muplus_P()-b.muplus_PZ()))/2.;
  });

  
  pm.Push<Hist1D>(Axis(80, 0, 0.4, mu_thetax, "|p_{x}^{reco}(#mu)/p_{z}^{reco}(#mu)|", {0.361}), "1", procs_comp_thetax_mu, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(62, 0, 0.31, mu_thetay, "|p_{y}^{reco}(#mu)/p_{z}^{reco}(#mu)|"), "1", procs_comp_thetay_mu, linplot).TopRight("13 TeV");
  pm.Push<Hist1D>(Axis(100, -0.05, 0.05,(mu_thetax_tru - mu_thetax)/mu_thetax, "|p_{x}^{reco}(#mu)/p_{z}^{reco}(#mu)| resolution (true-reco)/reco",{0.01/0.37}),
                  "1", procs_edgex_mu, linplot).TopRight("13 TeV");


  
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
