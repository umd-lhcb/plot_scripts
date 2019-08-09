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

  /*PlotOpt log_lumi("txt/plot_styles.txt", "CMSPaper");
  log_lumi.Title(TitleType::data)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::log)
    .Stack(StackType::data_norm).LegendColumns(3);
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
  *///vector<PlotOpt> linplot = {lin_shapes_info, lin_lumi_info};
  
  PlotOpt q2_plots("txt/plot_styles.txt", "PRLPaper");
  q2_plots.Title(TitleType::info)
    //.Bottom(BottomType::diff)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none)
    .LegendColumns(4)
    .ShowBackgroundError(false)
    .LeftMargin(0.25)
    .YTitleOffset(2.0);
    //.LegendDensity(0.01);

  vector<PlotOpt> LHCb_plots = {q2_plots};


  //vector<PlotOpt> linplot = {lin_lumi_info};
  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  /*vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_std>("test1", Process::Type::background, colors("tt_1l"),
                                                set<string>({"ntuples/tree.root"}), "1"));
  procs.push_back(Process::MakeShared<Baby_std>("test2", Process::Type::background, colors("qcd"),
                                                set<string>({"ntuples/tree.root"}), "1"));
  */

  //vector<shared_ptr<Process> > procs_phoebe;
  //procs_phoebe.push_back(Process::MakeShared<Baby_phoebe>("test3", Process::Type::background, colors("tt_1l"),
  //                                              set<string>({"ntuples/phoebe.root"}), "1"));
  //procs_phoebe.push_back(Process::MakeShared<Baby_phoebe>("test4", Process::Type::background, colors("qcd"),
  //                                              set<string>({"ntuples/phoebe.root"}), "1"));

//COMMENT THIS OUT BECAUSE I WOULD HAVE TO CHANGE PRETTY MUCH ENTIRE INFRASTRUCTURE OF THE CODE FOR IT TO NOT ERROR
  //procs_LHC.push_back(Process::MakeShared<Baby_MC>("comb", Process::Type::background, colors("ttv"),
  //                                                 set<string>({"ntuples/data.root"}),
  //                                                 "((DstIDprod > 0 && IDprod < 0 &&  ((muPID < 1.) || (muPID > 0. && totWeight == 1.)) && isData > 0. && TMath::Abs(totWeight) <= 1.) && flagDstSB==0.)*totWeight"));


//PLOTTING PHOEBE'S NTUPLE DATA AND MC
  /*auto taunuProc = Process::MakeShared<Baby_MC>("B2Dtaunu", Process::Type::background, colors("t1tttt"),
                                                set<string>({"ntuples/MC.root"}),
                                                "(isData == 0. &&  flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511)*FFweight*mcWeight");
  
  //add processes for B->D*H_c(->lvX')X decay
  auto HXProc = Process::MakeShared<Baby_MC>("h_dDDmu", Process::Type::background, colors("wjets"),
                                             set<string>({"ntuples/MC.root"}),
                                             "(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && (Btype==511))*mcWeight");
  auto HXProc2 = Process::MakeShared<Baby_MC>("h_uDDmu", Process::Type::background, colors("wjets"),
                                             set<string>({"ntuples/MC.root"}),
                                             "(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && Btype==521 && flagTauonicD < 1.)*mcWeight");
  auto HXProc3 = Process::MakeShared<Baby_MC>("h_sDDmu", Process::Type::background, colors("wjets"),
                                             set<string>({"ntuples/MC.root"}),
                                             "(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && Btype==531 && flagTauonicD < 1.)*mcWeight");

  //add processes for B->D**lv decay
  auto lnuProc = Process::MakeShared<Baby_MC>("h_D1", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 10413)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc2 = Process::MakeShared<Baby_MC>("h_D2", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 415)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc3 = Process::MakeShared<Baby_MC>("h_D1p", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 20413)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc4 = Process::MakeShared<Baby_MC>("h_uD1", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==521 && Dststtype == 10423)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc5 = Process::MakeShared<Baby_MC>("h_uD2", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==521 && Dststtype == 425)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc6 = Process::MakeShared<Baby_MC>("h_uD1p", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "((((mm_mom <= 250. || Dststtype == Dst_2010_minus_MC_MOTHER_ID || Dststtype == -1*Dst_2010_minus_MC_MOTHER_ID))) && isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==521 && Dststtype == 20423)) && muPID == 1. && !ishigher)*FFweight*mcWeight");
  auto lnuProc7 = Process::MakeShared<Baby_MC>("h_Ds2", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "(isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==531 && Dststtype == 435)) && muPID == 1.)*FFweight*mcWeight");
  auto lnuProc8 = Process::MakeShared<Baby_MC>("h_Ds1p", Process::Type::background, colors("purple"),
                                              set<string>({"ntuples/MC.root"}),
                                              "(isData == 0. &&  (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==531 && Dststtype == 10433)) && muPID == 1.)*FFweight*mcWeight");

  auto munuProc = Process::MakeShared<Baby_MC>("B2Dmunu", Process::Type::background, colors("blue"),
                                               set<string>({"ntuples/MC.root"}),
                                               "(isData == 0. &&  flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0)*FFweight*mcWeight");
  auto combProc = Process::MakeShared<Baby_data>("comb", Process::Type::background, colors("single_t"),
                                               set<string>({"ntuples/data.root"}),
                                               "((DstIDprod > 0 && IDprod < 0 &&  ((muPID < 1.) || (muPID > 0. && totWeight == 1.)) && isData > 0. && ((totWeight <= 1.) && (totWeight >= -1.))) && flagDstSB==0.)*totWeight");
  auto misProc = Process::MakeShared<Baby_data>("misIDmu", Process::Type::background, colors("qcd"),
                                                set<string>({"ntuples/data.root"}),
                                                "((DstIDprod > 0 && IDprod > 0 && muPID < 1.) &&  isData > 0. && flagDstSB==0.)*totWeight");

//plot using alternate missing ID, h_misIDalt
//  auto misProc = Process::MakeShared<Baby_data>("misIDmu_alt", Process::Type::background, colors("qcd"),
//                                                set<string>({"ntuples/data.root"}),
//                                                "((DstIDprod > 0 && IDprod > 0 && muPID < 1.) &&  isData > 0. && flagDstSB==0.)*(((f_pi > 0 && f_pi < 1)*f_pi*pi2mu)+((f_k > 0 && f_k < 1)*f_k*k2mu)+((f_p > 0 && f_p < 1)*f_p*p2mu))");


  auto dataProc = Process::MakeShared<Baby_data>("raw_data", Process::Type::background, kBlack,
                                                 set<string>({"ntuples/data.root"}),
                                                 "(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0. && Y_M < 5280 && ((iso_BDT < 0.15) || (ishigher && keepme && isData==0.))");

  dataProc->SetMarkerStyle(20);
  dataProc->SetMarkerSize(1);

  //vector<shared_ptr<Process> > procs_LHC = {dataProc, taunuProc, HXProc, HXProc2, HXProc3, lnuProc, lnuProc2, lnuProc3, lnuProc4, lnuProc5, lnuProc6,
  //                                          lnuProc7, lnuProc8, munuProc, combProc, misProc};
*/

//PLOTTING OF PHOEBE'S DATA **ONLY** (excludes one of the four data files as this is included above as dataProc)
/*
  auto data2Proc = Process::MakeShared<Baby_data>("raw_data2", Process::Type::background, kBlack,
                                                  set<string>({"ntuples/data_MU_2011.root"}),
                                                  "(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0. && Y_M < 5280 && ((iso_BDT < 0.15) || (ishigher && keepme && isData==0.))");
  data2Proc->SetMarkerStyle(20);
  data2Proc->SetMarkerSize(1);

  auto data3Proc = Process::MakeShared<Baby_data>("raw_data3", Process::Type::background, kBlack,
                                                  set<string>({"ntuples/data_MD_2012.root"}),
                                                  "(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0. && Y_M < 5280 && ((iso_BDT < 0.15) || (ishigher && keepme && isData==0.))");
  data3Proc->SetMarkerStyle(20);
  data3Proc->SetMarkerSize(1);

  auto data4Proc = Process::MakeShared<Baby_data>("raw_data4", Process::Type::background, kBlack,
                                                  set<string>({"ntuples/data_MU_2012.root"}),
                                                  "(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0. && Y_M < 5280 && ((iso_BDT < 0.15) || (ishigher && keepme && isData==0.))");
  data4Proc->SetMarkerStyle(20);
  data4Proc->SetMarkerSize(1);

  //vector<shared_ptr<Process> > procs_LHC = {dataProc, data2Proc, data3Proc, data4Proc};
*/

  //PLOTTING OF OUR NTUPLES (generated by Yipeng)

  auto sypdataProc = Process::MakeShared<Baby_yipeng_data_MD>("syp_data_MD", Process::Type::background, kBlack,
                                                              set<string>({"ntuples/yipeng_data_MD.root"}),
                                                              "(Y_ISOLATION_BDT < 0.15)");
  sypdataProc->SetMarkerStyle(20);
  sypdataProc->SetMarkerSize(1);

  vector<shared_ptr<Process> > procs_LHC = {sypdataProc};



  PlotMaker pm;

  /*  
  pm.Push<Hist1D>(Axis(20, 0, 30,"y_pt/1000", "#Upsilon p_{T} (GeV)"), "1", procs, linplot);
  pm.Push<Hist1D>(Axis(20, 0, 300,"y_pz/1000", "#Upsilon p_{Z} (GeV)"), "1", procs, linplot);
  */

  //pm.Push<Hist1D>(Axis(20, 0, 60,"D0_PT/1000", "D_{0} p_{T} (GeV)"), "1", procs_phoebe, linplot);
  //pm.Push<Hist1D>(Axis(20, 0, 50,"K_PT/1000", "K p_{T} (GeV)"), "1", procs_phoebe, linplot);

//PLOTTING FOR PHOEBE'S FILES (data and MC)
/*
  pm.Push<Hist1D>(Axis(40, -2, 10, "m_nu1", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "q2 >= -400000 && q2 < 2850000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2, 10, "m_nu1", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "q2 >= 2850000 && q2 < 6100000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2, 10, "m_nu1", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "q2 >= 6100000 && q2 < 9350000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2, 10, "m_nu1", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "q2 >= 9350000 && q2 <= 12600000", procs_LHC, LHCb_plots);

  pm.Push<Hist1D>(Axis(32, 100, 2500, "El", "E_{u}^{*} (MeV)"),
                  "q2 >= -400000 && q2 < 2850000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "El", "E_{u}^{*} (MeV)"),
                  "q2 >= 2850000 && q2 < 6100000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "El", "E_{u}^{*} (MeV)"),
                  "q2 >= 6100000 && q2 < 9350000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "El", "E_{u}^{*} (MeV)"),
                  "q2 >= 9350000 && q2 <= 12600000", procs_LHC, LHCb_plots);
*/

//PLOTTING FOR OUR FILES
  pm.Push<Hist1D>(Axis(40, -2000000, 10000000, "FitVar_Mmiss2", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "FitVar_q2 >= -400000 && FitVar_q2 < 2850000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2000000, 10000000, "FitVar_Mmiss2", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "FitVar_q2 >= 2850000 && FitVar_q2 < 6100000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2000000, 10000000, "FitVar_Mmiss2", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "FitVar_q2 >= 6100000 && FitVar_q2 < 9350000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(40, -2000000, 10000000, "FitVar_Mmiss2", "m_{miss}^{2} (GeV^{2}/c^{4})"),
                  "FitVar_q2 >= 9350000 && FitVar_q2 <= 12600000", procs_LHC, LHCb_plots);

  pm.Push<Hist1D>(Axis(32, 100, 2500, "FitVar_El", "E_{u}^{*} (MeV)"),
                  "FitVar_q2 >= -400000 && FitVar_q2 < 2850000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "FitVar_El", "E_{u}^{*} (MeV)"),
                  "FitVar_q2 >= 2850000 && FitVar_q2 < 6100000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "FitVar_El", "E_{u}^{*} (MeV)"),
                  "FitVar_q2 >= 6100000 && FitVar_q2 < 9350000", procs_LHC, LHCb_plots);
  pm.Push<Hist1D>(Axis(32, 100, 2500, "FitVar_El", "E_{u}^{*} (MeV)"),
                  "FitVar_q2 >= 9350000 && FitVar_q2 <= 12600000", procs_LHC, LHCb_plots);


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
