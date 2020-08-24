/* Overall, the goal of this script is to translate Phoebe's code redoHistos_Dst in
RDRDsrRun1AnalaysisPreservation to work using the plotscripts packages. The plan (for now) is to
essentially copy her code and run over a run 1 ntuple (with data and MC). I'll try to mark the places
where it's unclear to me if Phoebe's code actually matches up with the ANA note, or anywhere where I
think work needs to be done for correctness' sake. I generally make comments like this with an "edit"
or, if I know what I want to do/investigate, a "TODO". */

/* Lines of Phoebe's code (at least partially) not implemented from Process function (including commented out code):
1100-1136, 1145-1171, 1172-1173, 1180-1192, 1197-1203, 1205-1254, 1260-1261,1263-1268, 1270, 1274, 1276-1279, 1281,
1282, 1287-1292, 1366-1379, 1380-1381
Still TODO lines 1382-1937, 2900/2989-2990 (have to check which of misID, comb, comb_D* event is, then add totweight
condition accordingly)... skip all other lines (mark them above... including lines with histogram filling you
didn't use/altered)
Also, I assume use_uBDT and calibTrig are true, and all the other variables in lines 35-49 are false
*/


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


// Groupings chosen to match Phoebe's plots in chapter 8
// Note dss is really dss_light, decaying to only 1 pion, while dss_2s decays to two pions
// Except for signal and normalization, I don't separate semimuonic and semitauonic modes
enum class eventType {data, dsptau, dspmu, dss, dss_str, dd, misID_plus_comb, dss_2s, unknown};


// Obviously, this getType func could be more efficient/concise, but I think it's useful to explicitly single
// out every process for now
// Note there is no D** (-> D* pi pi) tau nu component. Assumedly, it would be small if it were included, I think?
// Most channels get contributions from both B0 and B-, but not the B0s->Ds**mu channels (assumedly the B-
// would have the same shape, so it wouldn't matter for the fit, I think?). These strange channels also don't
// include Ds**tau, assumedly for the same reason that the contribution would be small.
eventType getType(Double_t isData, Double_t DstIDprod, Double_t IDprod, Float_t muPID, Double_t flagDstSB,
                  Double_t flagtaumu, Double_t JustDst, Double_t DstOk, Int_t Btype, Double_t flagBmu,
                  Int_t Y_BKGCAT, Double_t flagTauonicD, Double_t flagDoubleD, Bool_t ishigher,
                  Int_t Dststtype, Int_t Dst_2010_minus_MC_MOTHER_ID, Double_t mm_mom){
  // data
  if(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0. && isData > 0.)
    return eventType::data;
  // combinatorial background (part of data, taken as background). Phoebe's h_comb
  /* TODO take totweight conditionals into consideration when making corresponding weight NamedFunc
    To do this, insert at end of totweight NamedFunc,
    if ((muPID > 1. && (muPID < 0. || totWeight != 1)) || TMath::Abs(totWeight) > 1)
      totWeight = 0; // this says to cut the event
    Comment this process out until appropriate weight func is created... uncommented for testing
  */
  else if(DstIDprod > 0 && IDprod < 0 && isData > 0. && flagDstSB==0.)
    return eventType::misID_plus_comb;
  // Combinatorial D* (part of data, taken as background). Phoebe's h_SB
  // Does not include h_comb_usb (combinatorial B)
  /* TODO in totweight NamedFunc, to mimic line 2990, insert
    if (totWeight != 1) totWeight *= -1;
    Comment this process out until appropriate weight func is created... uncommented for testing
  */
  else if(DstIDprod > 0 && IDprod > 0 && isData > 0. && flagDstSB==1.)
    return eventType::misID_plus_comb;
  // Misidentified PID (part of data, but taken as background MC). Phoebe's h_misID
  else if(DstIDprod > 0 && IDprod > 0 && muPID < 1. && isData > 0. && flagDstSB==0.)
    return eventType::misID_plus_comb;

  // signal B0 -> D* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511)
    return eventType::dsptau;
  // normalization B0 -> D* mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0)
    return eventType::dspmu;

  // B -> D** (-> D* pi pi) mu nu, either B0 or B-, D** is one of heavy (2S, not 1P) states, Phoebe's h_D2Smu
  else if(isData == 0. && flagBmu > 0. && JustDst < 1.  && DstOk > 0. && ishigher && muPID == 1.)
    return eventType::dss_2s;

  // B0 -> D1 (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 10413
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D2* (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 415
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D1' (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 20413
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1 (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 10423
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D2* (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 425
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1' (-> D* pi) mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 20423
     && muPID == 1. && !ishigher && mm_mom < 250. && Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;

  // B0 -> D1 tau nu   edit? no more !ishigher?
  // Notably, there aren't many of these D**tau events--about a few thousand
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 10413
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D2* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 415
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B0 -> D1' tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==511 && Dststtype == 20413
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1 tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 10423
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D2* tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 425
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;
  // B- -> D1' tau nu
  else if(isData == 0. &&  flagtaumu > 0. && JustDst < 1. && DstOk > 0. && Btype==521 && Dststtype == 20423
     && muPID == 1. && Dststtype == TMath::Abs(Dst_2010_minus_MC_MOTHER_ID))
    return eventType::dss;

  // B0s -> Ds1' mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==531 && Dststtype == 435
     && muPID == 1.)
    return eventType::dss_str;
  // B0s -> Ds2* mu nu
  else if(isData == 0. &&  flagBmu > 0. && JustDst < 1. && DstOk > 0. && Btype==531 && Dststtype == 10433
     && muPID == 1.)
    return eventType::dss_str;

  // B0 -> D* [Dq -> mu nu X'] X
  // TODO note the weights in Phoebe's code for these components!! Line 2827, 2847, 2865
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==511)
    return eventType::dd;
  // B- -> D* [Dq -> mu nu X'] X
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==521)
    return eventType::dd;
  // B0s -> D* [Dq -> mu nu X'] X ***NOTE*** this process is not mentioned in table 18; there is also no
  // corresponding process for tau... I'm going to comment this out, for now. edit?
  // else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && Btype==531)
  //   return eventType::dd;
  // B0 -> D* [Ds -> tau nu] X ***NOTE*** why no more flagDoubleD condition? Look at def of this
  // variable and flagTauonicD in AddB.C... edit?
  // TODO note the weights in Phoebe's code for these components!! Line 2877, 2889
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0. && Btype==511)
    return eventType::dd;
  // B- -> D* [Ds -> tau nu] X
  else if(isData == 0. &&  DstOk > 0. && muPID == 1. && flagTauonicD > 0. && Btype==521)
    return eventType::dd;

  return eventType::unknown;
}

// Function is potentially also obsolete, like many weight calculations below...
Double_t GetFinalMCWeight(Double_t w_mc) { // run after all other MC weight calculations are done
  if (w_mc > 50) w_mc = 0;
  if (w_mc > 5) w_mc = 5;
  if (w_mc < 0) w_mc = 0;
  if (isnan(w_mc)) w_mc = 0; // not sure why cmath isnan works but TMath::IsNaN doesn't? Not worth looking into.
  return w_mc;
}

Bool_t PassesBasicCuts(const Baby &b) { // event must pass all these cuts and also other, more specific, cuts
  if ((b.isData() > 0 && b.muPID() > 0 && b.Y_M() < 5280 &&
    (//DOCAVAR > DOCAmax
    //|| DOCAVAR < DOCAmin
    b.m_nu1() < -2.0 // mm_low
    || b.m_nu1() > 10. // mm_high
    || b.El() > 2500 // El_high
    || b.El() < 100 // El_low
    || b.q2()*1e-6 > 12.6 // q2_high
    || b.q2()*1e-6 < -0.4)) // q2_low
    || b.Y_M() > 5280 || (b.isData() > 0 && b.muVeto()) || b.dxy() > 7 || !(b.Hlt1() && b.Hlt2()) || (!b.YTIS() && !b.YTOS())
    || !(((b.Hlt1TAL0K() && b.K_PT() > 1700) || (b.Hlt1TAL0pi() && b.pi_PT() > 1700)))
    || (b.isData() > 0 && b.muPID() > 0 && b.BDTmu() < 0.25)) {
      return false;
  }
  else if (b.Y_M() > 5280) return false;
  else if (b.isData() > 0 && b.muVeto()) return false;
  else if (b.dxy() > 7) return false;
  else if (!(b.Hlt1() && b.Hlt2())) return false;
  else if (!b.YTOS() && !b.YTIS()) return false;
  else if (!((b.Hlt1TAL0K() && b.K_PT() > 1700) || (b.Hlt1TAL0pi() && b.pi_PT() > 1700))) return false;
  else if (b.isData() > 0 && b.muPID() > 0 && b.BDTmu() < 0.25) return false;   // ASSUMES use_uBDT IS TRUE
  //else if (b.isData() > 0 && b.muPID() > 0 && b.BDTmu() > 0.25) return false    // ASSUMES use_notuBDT IS TRUE
  else if (b.mu_P() > 100e3) return false;
  // edit? This next line slightly different from Phoebe's code; I think it's ok, because I'm not plotting the
  // data Dst sideband like Phoebe is... when this next line is included, data yields look ok (at least for ISO)
  else if (TMath::Abs(b.Dst_M()-b.D0_M()-145.454) > 2.) return false;
  // **NOTE** Be careful of line 2305 for flagDstSB when filling components! (I think?) TODO

  return true;
}




int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  ///////////////////////// BEGIN PHOEBE PLOTS //////////////////////////////

  //// User defined colors
  Palette colors("txt/colors.txt", "default");

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::info)
    .Bottom(BottomType::ratio)
    .YAxis(YAxisType::linear)
    .Stack(StackType::data_norm)
    .Overflow(OverflowType::none);

  vector<PlotOpt> plottypes = {lin_lumi}; // If needed, create more plot types and include in this vector; then when
                                          // provided as argument for histograms, all plot types should be created
                                          // as a separate pdf

  // ntuple contains run 1 data, MC signal, MC normalization, MC D** and MC DD processes for plots
  string repofolder = "ntuples/";
  string ntuplefile = "ref-rdx-run1/Dst-mix/Dst--20_07_02--mix--all--2011-2012--md-mu--phoebe.root";


  /////////////////// Event type and weights for processes/plots /////////////////////////

  NamedFunc event_type("event_type",[&](const Baby &b){
    return static_cast<Double_t>(getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(), b.flagDstSB(),
                                       b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
                                       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(),
                                       b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()));
  });



/* All of these weight funcs are obsolete (they aren't fully implemented from Phoebe's code, either,
as I stopped working on them). I'll leave them here, commented out, but doubt I'll look at them again.

  // Weights used depends on the sample we're looking at, so make a weight func for each sample

  NamedFunc w_iso("w_iso",[&](const Baby &b){
    Double_t w_mc=b.mcWeight();
    Double_t w_ff=b.FFweight();

    w_mc = GetFinalMCWeight(w_mc);
    // for some reason, for D** tau nu processes, FFweight isn't used anymore... Note that flagtaumu and flagBmu
    // should be mutually exclusive, so for D** events, test if there's a tau, and set w_ff=1 if there is
    if (b.flagtaumu() > 0. && eventType::dss == getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(),
       b.flagDstSB(), b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(), b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()))
      w_ff=1;

    return w_mc*w_ff;
  });

  NamedFunc w_dd("w_dd",[&](const Baby &b){ // lines 1294-1325
    Double_t w_mc=b.mcWeight();
    Double_t w_ff=b.FFweight();

    Double_t iso_BDT = b.iso_BDT();
    Double_t iso_NNkw = b.iso_NNkw();
    if (iso_NNkw > -1 && iso_BDT > 0.15) {
      iso_NNkw = iso_NNkw*(b.iso_Type()==3);
      Double_t iso_NNkw2 = b.iso_NNkw2()*(b.iso_Type2()==3);
      Double_t kidweight = iso_NNkw;
      Double_t iso_BDT2 = b.iso_BDT2();
      if(iso_NNkw2 > -1 && iso_BDT2 > -2) {
        Double_t iso_NNkw3 = b.iso_NNkw3()*(b.iso_Type3()==3);
        Double_t iso_BDT3 = b.iso_BDT3();
        if(iso_NNkw3 > -1 && iso_BDT3 > -2) {
          kidweight=iso_NNkw+iso_NNkw2+iso_NNkw3-iso_NNkw*iso_NNkw2-iso_NNkw*iso_NNkw3-iso_NNkw2*iso_NNkw3+iso_NNkw*iso_NNkw2*iso_NNkw3;
        }
        else {
          kidweight=iso_NNkw+iso_NNkw2-iso_NNkw*iso_NNkw2;
        }
      }
      w_mc*=kidweight;
    }

    w_mc = GetFinalMCWeight(w_mc);
    // for some reason, for D** tau nu processes, FFweight isn't used anymore... Note that flagtaumu and flagBmu
    // should be mutually exclusive, so for D** events, test if there's a tau, and set w_ff=1 if there is
    if (b.flagtaumu() > 0. && eventType::dss == getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(),
       b.flagDstSB(), b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(), b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()))
      w_ff=1;

    return w_mc*w_ff;
  });

  // D** sample not explicitly mentioned in Phoebe's code: copy weights for 1OS or 2OS? Maybe see ANA pg 47?
  NamedFunc w_dss("w_dss",[&](const Baby &b){
    Double_t w_mc=b.mcWeight();
    Double_t w_ff=b.FFweight();

    Double_t iso_NNkw = b.iso_NNkw();
    if(iso_NNkw > -1) { // I won't worry about "thecut==true", since if it isn't true, this event will be cut anyway
      //iso_NNkw=iso_NNkw*(b.iso_Type()==3);
      double kidweight=(1-iso_NNkw);
      w_mc*=kidweight;
    }

    w_mc = GetFinalMCWeight(w_mc);
    // for some reason, for D** tau nu processes, FFweight isn't used anymore... Note that flagtaumu and flagBmu
    // should be mutually exclusive, so for D** events, test if there's a tau, and set w_ff=1 if there is
    if (b.flagtaumu() > 0. && eventType::dss == getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(),
       b.flagDstSB(), b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(), b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()))
      w_ff=1;

    return w_mc*w_ff;
  });

  NamedFunc w_2os("w_2os",[&](const Baby &b){
    Double_t w_mc=b.mcWeight();
    Double_t w_ff=b.FFweight();

    Double_t iso_NNkw = b.iso_NNkw();
    if(iso_NNkw > -1) { // I won't worry about "thecut==true", since if it isn't true, this event will be cut anyway
      iso_NNkw=iso_NNkw*(b.iso_Type()==3);
      Double_t iso_NNkw2 = b.iso_NNkw2();
      iso_NNkw2=iso_NNkw2*(b.iso_Type2()==3);
      double kidweight=(1-iso_NNkw);
      kidweight=1-(iso_NNkw+iso_NNkw2-iso_NNkw*iso_NNkw2);
      w_mc*=kidweight;
    }

    w_mc = GetFinalMCWeight(w_mc);
    // for some reason, for D** tau nu processes, FFweight isn't used anymore... Note that flagtaumu and flagBmu
    // should be mutually exclusive, so for D** events, test if there's a tau, and set w_ff=1 if there is
    if (b.flagtaumu() > 0. && eventType::dss == getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(),
       b.flagDstSB(), b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(), b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()))
      w_ff=1;

    return w_mc*w_ff;
  });


  NamedFunc w_1os("w_1os",[&](const Baby &b){
    Double_t w_mc=b.mcWeight();
    Double_t w_ff=b.FFweight();

    Double_t iso_NNkw = b.iso_NNkw();
    if(iso_NNkw > -1) { // I won't worry about "thecut==true", since if it isn't true, this event will be cut anyway
      //iso_NNkw=iso_NNkw*(b.iso_Type()==3);
      double kidweight=(1-iso_NNkw);
      w_mc*=kidweight;
    }

    w_mc = GetFinalMCWeight(w_mc);
    // for some reason, for D** tau nu processes, FFweight isn't used anymore... Note that flagtaumu and flagBmu
    // should be mutually exclusive, so for D** events, test if there's a tau, and set w_ff=1 if there is
    if (b.flagtaumu() > 0. && eventType::dss == getType(b.isData(), b.DstIDprod(), b.IDprod(), b.muPID(),
       b.flagDstSB(), b.flagtaumu(), b.JustDst(), b.DstOk(), b.Btype(), b.flagBmu(), b.Y_BKGCAT(),
       b.flagTauonicD(), b.flagDoubleD(), b.ishigher(), b.Dststtype(), b.Dst_2010_minus_MC_MOTHER_ID(), b.mm_mom()))
      w_ff=1;

    return w_mc*w_ff;
  });

  // NamedFunc totweight() TODO -- for data components used in background, will have to be sure of process order and
                                // include this weight func in the right position in the vector given as argument
                                // to the plot pushed to plotmaker

*/


  ////////////////////////// Processes //////////////////////////////////

  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Data",Process::Type::data, colors("data"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::data)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Comb+MisID",Process::Type::background, colors("purple"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::misID_plus_comb)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* #tau #nu",Process::Type::signal, colors("dsptau"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dsptau)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* #mu #nu",Process::Type::background, colors("dspmu"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dspmu)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D** l #nu",Process::Type::background, colors("dss"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dss)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D_{H}** #mu #nu",Process::Type::background, kBlue,
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dss_2s)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D_{s}** #mu #nu",Process::Type::background, colors("orange"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dss_str)));
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("B #rightarrow D* D_{q,s} X", Process::Type::background, colors("dd"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::dd)));

  // Don't include unknown processes? They are concentrated heavily at large values (of mmiss2, El, q2), and there
  // are a lot of them...
  /*procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Unknown process",Process::Type::background, colors("purple"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       event_type == static_cast<Double_t>(eventType::unknown)));*/





  ///////////////////////// Selection cuts ////////////////////////
  // ISO cuts not explicitly mentioned in Phoebe's code
  NamedFunc isocuts("ISO", [&](const Baby &b){ // edit?
    if (!(b.ISOnum()==0)) return false; // singleCand, line 1282
    return PassesBasicCuts(b) && (b.iso_BDT() < 0.15 || (b.ishigher() && b.keepme() && b.isData()==0.));
  });

  // DD cuts starts at line 1294
  NamedFunc ddcuts("DD", [&](const Baby &b){ // edit? this differs (more complicated) from table 15 of ANA note (July 9)
    if (!(b.AntiISOnum()==0)) return false; // singleCand, line 1300
    if (b.iso_NNkw() > -1) return (b.iso_BDT() > 0.15);

    Double_t iso_NNk = b.iso_NNk();
    Double_t iso_NNk2 = b.iso_NNk2();
    Double_t iso_NNk3 = b.iso_NNk3();
    Float_t iso_P = b.iso_P();
    Float_t iso_PT = b.iso_PT();
    Float_t iso_P2 = b.iso_P2();
    Float_t iso_PT2 = b.iso_PT2();
    Float_t iso_P3 = b.iso_P3();
    Float_t iso_PT3 = b.iso_PT3();
    Double_t iso_BDT = b.iso_BDT();
    Double_t iso_BDT2 = b.iso_BDT2();
    Double_t iso_BDT3 = b.iso_BDT3();
    if (iso_BDT <= -1.1) iso_NNk = 0.;
    if (iso_BDT2 <= -1.1) iso_NNk2 = 0.;
    if (iso_BDT3 <= -1.1) iso_NNk3 = 0.;


    return PassesBasicCuts(b) && (b.iso_BDT() > 0.15 && (iso_NNk > 0.2 || iso_NNk2 > 0.2 || iso_NNk3 > 0.2)
            && TMath::Max(iso_P*(iso_PT > 150),TMath::Max(iso_P2*(iso_PT2 > 150)*(iso_BDT2 > -1.1),
            iso_P3*(iso_PT3 > 150)*(iso_BDT3 > -1.1))) > 5e3);
  });

  // 2OS cuts starts at line 1326
  NamedFunc twooscuts("2OS", [&](const Baby &b){ // edit? this differs (more complicated) from table 15 of ANA note (July 9)
                                                 // Actually, I think Phoebe's code is logically incorrect, too. Lines 1342
                                                 // and 1343 shouldn't be inside the "else" (as then they'll only be applied
                                                 // if "thecut" is already false)
     if (!(b.AntiISOnum()==0)) return false; // singleCand, line 1329

     Double_t iso_NNk = b.iso_NNk();
     Double_t iso_NNk2 = b.iso_NNk2();
     Float_t iso_P = b.iso_P();
     Float_t iso_PT = b.iso_PT();
     Float_t iso_P2 = b.iso_P2();
     Float_t iso_PT2 = b.iso_PT2();
     Double_t iso_BDT = b.iso_BDT();
     Double_t iso_BDT2 = b.iso_BDT2();
     Double_t iso_BDT3 = b.iso_BDT3();
     Float_t iso_CHARGE = b.iso_CHARGE();
     Float_t iso_CHARGE2 = b.iso_CHARGE2();


     return PassesBasicCuts(b) && (iso_BDT > 0.15 && iso_BDT2 > 0.15 && iso_BDT3 < 0.15 && iso_CHARGE != iso_CHARGE2
       && iso_CHARGE != 0 && iso_CHARGE2 != 0 && iso_CHARGE < 100
       && TMath::Max(iso_P*(iso_PT > 150),iso_P2*(iso_PT2 > 150)) > 5e3
       && iso_NNk < 0.2 && iso_NNk2 < 0.2);
  });

  // 1OS cuts starts at line 1346
  NamedFunc oneoscuts("1OS", [&](const Baby &b){ // edit? this differs (more complicated) from table 15 of ANA note (July 9)
                                                 // Actually, I think Phoebe's code is logically incorrect, too. Line 1363
                                                 // shouldn't be inside the "else" (as then it'll only be applied
                                                 // if "thecut" is already false)
     if (!(b.AntiISOnum()==0)) return false; // singleCand, line 1350

     Double_t iso_NNk = b.iso_NNk();
     Float_t iso_P = b.iso_P();
     Float_t iso_PT = b.iso_PT();
     Double_t iso_BDT = b.iso_BDT();
     Double_t iso_BDT2 = b.iso_BDT2();
     Float_t iso_CHARGE = b.iso_CHARGE(); // charge not mentioned in table 15?
     Int_t Dst_ID = b.Dst_ID();

     return PassesBasicCuts(b) && (iso_BDT > 0.15 && iso_BDT2 < 0.15 && iso_CHARGE*Dst_ID < 0 && iso_P > 5e3
       && iso_PT > 150 && iso_NNk < 0.2);
  });

  // D** cuts just the same as 1OS but with mass cut on D*+pi1... I actually think this is what is implemented
  // in Phoebe's code currently, so remove the mass cut above, and just copy the 1OS cut w the mass cut here. Also,
  // use the mass cut given in table 12, not in the code.
  NamedFunc dsscuts("DSS", [&](const Baby &b) { // edit?
    if (!(b.AntiISOnum()==0)) return false; // singleCand, line 1350

    Double_t iso_NNk = b.iso_NNk();
    Float_t iso_P = b.iso_P();
    Float_t iso_PT = b.iso_PT();
    Double_t iso_BDT = b.iso_BDT();
    Double_t iso_BDT2 = b.iso_BDT2();
    Float_t iso_CHARGE = b.iso_CHARGE(); // charge not mentioned in table 15?
    Int_t Dst_ID = b.Dst_ID();
    Double_t iso_DeltaM = b.iso_DeltaM();

    return PassesBasicCuts(b) && (iso_BDT > 0.15 && iso_BDT2 < 0.15 && iso_CHARGE*Dst_ID < 0 && (iso_DeltaM > 390
      && iso_DeltaM < 510) && iso_P > 5e3 && iso_PT > 150 && iso_NNk < 0.2);
  });


  /////////////////////// Plotting variables ////////////////////////

  NamedFunc mm2("mm2", [&](const Baby &b) {
    if (b.isData() > 0.) return b.m_nu1(); // data
    else return b.m_nu1smear()*1e-6; // MC
  });

  NamedFunc el("el", [&](const Baby &b) {
    if (b.isData() > 0.) return b.El()*1e-3; // data
    else return b.Elsmear()*1e-3; // MC
  });

  NamedFunc q2("q2", [&](const Baby &b) {
    if (b.isData() > 0.) return b.q2()*1e-6; // data
    else return b.q2smear()*1e-6; // MC
  });


  /////////////////////////// Plots /////////////////////////////

  PlotMaker pm;
  vector<NamedFunc> cuts{isocuts};//, ddcuts, twooscuts, oneoscuts, dsscuts};
  //vector<NamedFunc> MC_w{w_iso, w_dd, w_2os, w_1os, w_dss}; // should correspond to cut in same order as vector above
  for (Int_t i=0; i<static_cast<Int_t>(cuts.size()); i++) {
    pm.Push<Hist1D>(Axis(40, -2.0, 10.0, mm2, "m_{miss}^{2} [GeV^{2}]"), cuts[i], procs, plottypes,
                    vector<NamedFunc>({"1"}));
                    //vector<NamedFunc>({"1", "1", MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i]}));
    pm.Push<Hist1D>(Axis(4, -0.4, 12.6, q2, "q^{2} [GeV^{2}]"), cuts[i], procs, plottypes,
                    vector<NamedFunc>({"1"}));
                    //vector<NamedFunc>({"1", "1", MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i]}));
    pm.Push<Hist1D>(Axis(32, 0.1, 2.5, el, "E_{l} [GeV]"), cuts[i], procs, plottypes,
                    vector<NamedFunc>({"1"}));
                    //vector<NamedFunc>({"1", "1", MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i], MC_w[i]}));
  }
  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to



  //////////////////////// END PHOEBE PLOTS ///////////////////////////

  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
