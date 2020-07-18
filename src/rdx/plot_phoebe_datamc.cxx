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

enum class eventType {data, dsptau, dspmu, dss, dd, unknown};

eventType getType(double isData, double DstIDprod, double IDprod, float muPID, double flagDstSB,
                  double flagtaumu, double JustDst, double DstOk, int Btype, double flagBmu,
                  int Y_BKGCAT){
  if((DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0.)
    return eventType::data;
  else if(isData == 0. &&  flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511)
    return eventType::dsptau;
  else if(isData == 0. &&  flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0)
    return eventType::dspmu;
  
  return eventType::unknown;
}


int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  //// User defined colors
  Palette colors("txt/colors.txt", "default");

  ////////////////////////////////////// Phoebe's plots ///////////////////////////////////////

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none);
  
  vector<PlotOpt> plottypes = {lin_lumi};

  //////// Signal, normalization, D** processes for mmiss2 plot
  string repofolder = "ntuples/";
  repofolder = "/Volumes/manuelf_t5/lhcb-ntuples-gen/ntuples/";
  string ntuplefile = "ref-rdx-run1/Dst-mix/Dst--20_07_02--mix--all--2011-2012--md-mu--phoebe.root";
  
  NamedFunc event_type("event_type",[&](const Baby &b){
                                      return static_cast<double>(getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                                         b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT()));
                                    });

  NamedFunc weight("weight",[&](const Baby &b){
                              eventType type = getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                       b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT());
                              if(type == eventType::data) return 1.;
                              else return b.FFweight()*b.mcWeight();
                            });

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Data",Process::Type::data, colors("data"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::data)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D*^{+} #tau #nu",Process::Type::signal, colors("dsptau"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::dsptau)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D*^{+} #mu #nu",Process::Type::background, colors("dspmu"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::dspmu)));
  
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(260, -3, 10,"m_nu1", "m_{miss}^{2} [GeV^{2}]"), "1", procs, plottypes).Weight(weight);
  pm.MakePlots(0.5); // The "1" is the luminosity to rescale the bkg to  


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

