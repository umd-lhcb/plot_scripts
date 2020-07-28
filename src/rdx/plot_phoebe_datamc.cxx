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

// ** I could add ddtau_u, ddtau_d, ddmu_u, ddmu_d, ddmu_s types: should these all be differentiated? No ddtau_s??
// Maybe, at least, you should differentiate between ddtau and ddmu (using flagTauonicD)... I will do this
// and I'll call the types with mu dxx (this is based on table 1 pg 5 of ANA note) and those with tau dd... I
// may change this in the future. I've DELETED THE DSS TYPE MANUEL CREATED, replacing it with dss_mu and dss_tau.
// This may change, too...
enum class eventType {data, dsptau, dspmu, dss_mu, dss_tau, dd, dxx, unknown}; // ADD ANY OTHER EVENT TYPES HERE

eventType getType(double isData, double DstIDprod, double IDprod, float muPID, double flagDstSB,
                  double flagtaumu, double JustDst, double DstOk, int Btype, double flagBmu,
                  int Y_BKGCAT, Double_t flagTauonicD, Double_t flagDoubleD){
  if((DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0.)
    return eventType::data;
  else if(isData == 0. &&  flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511)
    return eventType::dsptau;
  else if(isData == 0. &&  flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0)
    return eventType::dspmu;
  // ADD CONDITIONS FOR DSS, DD, AND ANY OTHER EVENT TYPES YOU NEED TO LOOK AT
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0.
          && (Btype==511 || Btype==521))
    return eventType::dd;
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1.
          && (Btype==511 || Btype==521 || Btype==531))
    return eventType::dxx;
  /* Examples of type criteria you could set if being more specific about DD/DXX' events. Note that Phoebe
  seems to not check flagDoubleD for the ddtau types, and she seemingly doesn't consider a ddtau_s...
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && Btype==531 && flagTauonicD < 1.)
    return eventType::ddmu_s;
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0. && Btype==511) // no flagDoubleD??
    return eventType::ddtau_d;
  */

  // add conditionals for D**...

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
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none);

  vector<PlotOpt> plottypes = {lin_lumi}; // If needed, create more plot types and include in this vector; then when
                                          // provided as argument for histograms, all plot types should be created
                                          // as a separate pdf

  //////// Signal, normalization, D** processes for mmiss2 plot
  string repofolder = "ntuples/";
  string ntuplefile = "ref-rdx-run1/Dst-mix/Dst--20_07_02--mix--all--2011-2012--md-mu--phoebe.root";

  NamedFunc event_type("event_type",[&](const Baby &b){
                                      return static_cast<double>(getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                                         b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
                                                                         b.flagTauonicD(), b.flagDoubleD()));
                                    });

  NamedFunc weight("weight",[&](const Baby &b){
                              eventType type = getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                       b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
                                                       b.flagTauonicD(), b.flagDoubleD());
                              if(type == eventType::data) return 1.;
                              else return b.FFweight()*b.mcWeight(); // does this have to be modified?
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
  // ADD ANY OTHER PROCESSES HERE, eg. MC B->D** tau nu, ASSUMEDLY AS BACKGROUND PROCESS TYPES
  // this first process might have to be modified? also the ones following it too...
  /*procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Unknown process",Process::Type::background, colors("purple"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::unknown)));*/
  // Are reaction descriptions for DD and DXX' modes correct? I'm basing on table 1 pg 5 of ANA note
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D*^{+}(D_{s}^{-} #rightarrow #tau #nu) X",
                                                       Process::Type::background, colors("dd"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::dd)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D*(X_{c} #rightarrow #mu #nu X') X",
                                                       Process::Type::background, colors("orange"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<double>(eventType::dxx)));
  // include process (actually, 2 I think) for D** feed down



  NamedFunc isocuts = "1"; // ADD HERE CUTS FOR ISO SAMPLE, FOR EXAMPLE. (does this have to be modified?)
  // CREATE NAMEDFUNC FOR OTHER TYPES OF CUTS
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(40, -2, 10, "m_nu1", "m_{miss}^{2} [GeV^{2}]"), isocuts, procs, plottypes,
                  vector<NamedFunc>({weight}));
  pm.Push<Hist1D>(Axis(40, -3, 13, "q2/1e6", "q^{2} [GeV^{2}]"), isocuts, procs, plottypes,
                  vector<NamedFunc>({weight})); // might have to be modified...
  pm.Push<Hist1D>(Axis(40, -0.2, 3.8, "El/1e3", "E_{l} [GeV]"), isocuts, procs, plottypes,
                  vector<NamedFunc>({weight})); // might have to be modified...
  // PUSH OTHER HISTOGRAMS TO PLOTMAKER HERE, WITH APPROPRIATE SELECTION CUTS
  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
