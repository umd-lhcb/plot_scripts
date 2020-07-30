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


// Somewhat arbitrary groupings (may want to separate strange decays from other dss...)
enum class eventType {data, dsptau, dspmu, dss, dd, misID_plus_comb, unknown}; // add any other event types here


// Obviously, this getType func could be more efficient/concise, but I think it's useful to explicitly single
// out every process for now
eventType getType(Double_t isData, Double_t DstIDprod, Double_t IDprod, Float_t muPID, Double_t flagDstSB,
                  Double_t flagtaumu, Double_t JustDst, Double_t DstOk, Int_t Btype, Double_t flagBmu,
                  Int_t Y_BKGCAT, Double_t flagTauonicD, Double_t flagDoubleD, Bool_t ishigher,
                  Int_t Dststtype, Int_t Dst_2010_minus_MC_MOTHER_ID, Double_t mm_mom){
  // data
  if(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0. && isData > 0.)
    return eventType::data;
  // signal B0 -> D* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511)
    return eventType::dsptau;
  // normalization B0 -> D* mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0)
    return eventType::dspmu;
  // B0 -> D1 mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D1 -> D* pi pi... ***NOTE*** if you include, exclude, or negate this condition, you don't get exactly
  // analagous conditions for selecting these reactions but with tau instead of mu... is this a mistake? For now,
  // I'm going to negate this condition and put it in the cuts for the muonic D** processes, but this may change.)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 10413
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D2* mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D2* -> D* pi pi)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 415
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D1' mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D1' -> D* pi pi)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 20413
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1 mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D1 -> D* pi pi)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 10423
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D2* mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D2* -> D* pi pi)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 425
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1' mu nu (Phoebe also specifies that, if mm_mom>250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID),
  // then D1' -> D* pi pi)
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 20423
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D1 tau nu ***NOTE*** no more !ishigher? No mm_mom cut? Also, more critically, if these tau D** processes
  // are included, the histograms really get thrown off (weights off other processes dramatically effected, and
  // many of these tauonic D** processes get counted, and only in a single bin... I'm going to comment them out
  // for now, but you should do some experimenting/investigating
  /*else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 10413
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D2* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 415
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D1' tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 20413
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1 tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 10423
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D2* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 425
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1' tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 20423
     && muPID == 1. && !ishigher && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;*/
  // B0s -> Ds1' mu nu ***NOTE*** Manuel said these aren't dss... they seem to be in Phoebe's code, and from
  // what I can tell this is how they're categorized in the ANA note, too. Are Ds1' and Ds2* not D**?
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==531 && Dststtype == 435
     && muPID == 1.)
    return eventType::dss;
  // B0s -> Ds2* mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==531 && Dststtype == 10433
     && muPID == 1.)
    return eventType::dss;
  // B0 -> D* [Dq -> mu nu X'] X
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==511)
    return eventType::dd;
  // B- -> D* [Dq -> mu nu X'] X
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==521)
    return eventType::dd;
  // B0s -> D* [Dq -> mu nu X'] X ***NOTE*** this process is not mentioned in table 18; there is also no
  // corresponding process for tau
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==531)
    return eventType::dd;
  // B0 -> D* [Ds -> tau nu] X ***NOTE*** why no more flagDoubleD condition?
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0. && Btype==511)
    return eventType::dd;
  // B- -> D* [Ds -> tau nu] X
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0. && Btype==521)
    return eventType::dd;

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

  //////// Signal, normalization, D** and DD processes for plots
  string repofolder = "ntuples/";
  string ntuplefile = "ref-rdx-run1/Dst-mix/Dst--20_07_02--mix--all--2011-2012--md-mu--phoebe.root";


  NamedFunc event_type("event_type",[&](const Baby &b){
                                      return static_cast<Double_t>(getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                                         b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
                                                                         b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(),
                                                                         b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()));
                                    });

  NamedFunc weight("weight",[&](const Baby &b){
                              eventType type = getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                                       b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
                                                       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(),
                                                       b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom());
                              if(type == eventType::data) return 1.;
                              else return b.FFweight()*b.mcWeight(); // this should be modified to account for weights
                                                                     // changing in Phoebe's code... probably better to
                                                                     // define individual weight functions, actually
                            });


  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Data",Process::Type::data, colors("data"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::data)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* #tau #nu",Process::Type::signal, colors("dsptau"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dsptau)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* #mu #nu",Process::Type::background, colors("dspmu"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dspmu)));
  // Don't include unknown processes? They are concentrated heavily at large values (of mmiss2, El, q2), and they affect
  // the weights of the other (Monte Carlo) processes
  /*procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Unknown process",Process::Type::background, colors("purple"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::unknown)));*/
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D** (#mu/#tau) #nu",Process::Type::background, colors("dss"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dss)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* D X",
                                                       Process::Type::background, colors("dd"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dd)));



  NamedFunc isocuts = "1"; // ADD HERE CUTS FOR ISO SAMPLE, FOR EXAMPLE. (does this have to be modified?)
  // CREATE NAMEDFUNCS FOR OTHER TYPES OF CUTS


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
