// Code for making plots for run2 fullsim MC validation (using plot_scripts)

// This is going to be slightly ugly and repetitive code... but I'm not sure there's any way around that,
// and being ugly may make it the easiest to edit in this case, honestly.

// Created: Feb 22, 2021
// Last edited: Mar 10, 2021
// Note: relevant ntuples (see file names below) need to be downloaded, and the repofolder variable should
// be edited to reflect their location.

// Naming convention: if I don't specify Bd or Bu, it's Bd

#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <set>
#include <string>
#include <assert.h>

#include <getopt.h>

#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TLorentzVector.h"

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

// just to associate processes with tags used for filenames when looping
class proc_tag {
private:
  vector<shared_ptr<Process>> procs;
  string tag;
public:
  proc_tag(vector<shared_ptr<Process>> _procs, string _tag) { // no default constructor
    procs.insert(procs.end(),_procs.begin(),_procs.end());
    tag=_tag;
  }
  vector<shared_ptr<Process>> Processes() {return procs;}
  string FileTag() {return tag;}
  ~proc_tag() {}; // no internal memory leaks
};

/*
// types of decays
enum class eventType {test_type,unknown};

// returns decay type for MC reconstructed as B0 -> D*+ [-> D0 [-> K- pi+] spi+] mu-
eventType getEventTypeBd(const Baby &b) {
  ULong64_t totcands = b.totCandidates(); // just to check I can indeed access baby variables this way...
  totcands+=0;
  if (totcands>0) return eventType::test_type;
  else return eventType::unknown;
}

// returns decay type for MC reconstructed as B- -> D0 [-> K- pi+] mu-
eventType getEventTypeBu(const Baby &b) {
  ULong64_t totcands = b.totCandidates(); // just to check I can indeed access baby variables this way...
  totcands+=0;
  if (totcands>0) return eventType::test_type;
  else return eventType::unknown;
}
*/



int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  // User defined colors
  Palette colors("txt/colors.txt", "default");

  //////////////////////////// Define Plot Types

  PlotOpt lin_shapes("txt/plot_styles.txt", "LHCbPaper");
  lin_shapes.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes)
    .Overflow(OverflowType::none);//both);

  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes);
  PlotOpt log_lumi_shapes = lin_lumi_shapes.YAxis(YAxisType::log);

  vector<PlotOpt> plotshapes = {lin_shapes};
  vector<PlotOpt> plotlumi = {lin_lumi_shapes};
  vector<PlotOpt> plotloglumi = {log_lumi_shapes};
  vector<PlotOpt> plotboth = {lin_shapes, lin_lumi_shapes};

  /////////////////////////// Define Cuts

  // Truly, this is very ugly code, not at all an improvement over Phoebe's implementation. I'm not sure why I started this
  // way, but if I'm being a good programmer, I should really implement the cuts using the commented out functions
  // above and below. This may be a good idea when sclaing up to all samples, because it should make defining the cuts
  // faster.
  // Note: I'm worried about some of my selections not being exclusive from one another because I'm not implementing all
  // of Phoebe's selections. In this current implementation, the non-exclusivity won't be an issue, but to ensure it isn't
  // an issue if pursuing the (better) implementation, you should have the "getEventType" funcs return a vector of eventType
  // (I think NamedFunc can handle returning vectors, naively looking at the named_func.cpp source code; I don't think it
  // can handle sets, though).
  // Actually, too, though this current implementation is uglier/less portable, I actually think it should be faster.
  // If time starts to become an issue, it may be worth considering using this current implementation instead of the
  // "better" one.

/*
  NamedFunc get_event_type_bd("get_event_type_bd",[&](const Baby &b){
    return static_cast<Double_t>(getEventTypeBd(b));
    return s;
  });

  NamedFunc get_event_type_bu("get_event_type_bu",[&](const Baby &b){
    return static_cast<Double_t>(getEventTypeBu(b));
  });
*/

  // General Cuts function, to be applied to EVERY decay that is reconstructed using TupleB0
  NamedFunc cuts_bd("cuts_bd", [&](const Baby &b){
    // DstOk from Phoebe AddB.C line 2547
    if (!(b.dst_BKGCAT()==0 || (b.dst_BKGCAT()==50 && b.d0_BKGCAT()==50))) return false;
    // muPID from Phoebe AddB.C line 2540
    if (!(abs(b.mu_TRUEID())==13)) return false;
    return true;
  });

  // General Cuts function, to be applied to EVERY decay that is reconstructed using TupleBminus
  NamedFunc cuts_bu("cuts_bu", [&](const Baby &b){
    // mcut from Phoebe redoHistos_D0.C line 2101
    if (!(abs(b.d0_M()-1864.6)<23.4)) return false;
    // muPID from Phoebe AddB.C line 2540
    if (!(abs(b.mu_TRUEID())==13)) return false;
    return true;
  });

  NamedFunc selec_Bd2DstMuNu("selec_Bd2DstMuNu", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: flagBmu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype == 511 && Y_BKGCAT==0
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())==511)) return false;
    // Y_BKGCAT just my b0_BKGCAT
    if (!(b.b0_BKGCAT()==0)) return false;
    return true;
  });

  NamedFunc selec_Bd2DstTauNu("selec_Bd2DstTauNu", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: flagtaumu > 0. && JustDst > 0. && DstOk > 0. && muPID == 1. && Btype==511
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagtaumu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==15 && abs(b.mu_MC_GD_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())==511)) return false;
    return true;
  });

  // note D0* decay "missing" because D0* can't decay to D*+ with only one pion (two pions, it can, but those two pion decays should be separate plots)
  NamedFunc selec_Bd2DststMuNu_D1("selec_Bd2DststMuNu_D1", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 10413)) && muPID == 1. && !ishigher (I ignore ishigher)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==10413)) return false;
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D1p("selec_Bd2DststMuNu_D1p", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 20413)) && muPID == 1. && !ishigher (I ignore ishigher)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==20413)) return false;
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D2st("selec_Bd2DststMuNu_D2st", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511)) && muPID == 1. && !ishigher && Dststtype==415 (I ignore ishigher)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==415)) return false;
    return true;
  });

/* these are not properly implemented... these NamedFuncs may actually be fine (?), but you shouldn't be using the higher D** state decays for your input ntuple!
 NamedFunc selec_Bd2DststMuNu_D1_pipi("selec_Bd2DststMuNu_D1_pipi", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 10413)) && muPID == 1. && !ishigher && (mm_mom > 250. || Dststtype!=TMath::Abs(Dst_2010_minus_MC_MOTHER_ID)) (I ignore ishigher, also I think last condition is always false)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==10413)) return false; // my implementation (and I think Phoebe's too) ensures Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID)
    // mm_mom from Phoebe AddB.C line 3136-3141
    TLorentzVector dst_mom_p(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector dst_p(b.dst_TRUEP_X(),b.dst_TRUEP_Y(),b.dst_TRUEP_Z(),b.dst_TRUEP_E());
    double mm2_mom=(dst_mom_p-dst_p).M2();
    if (!(mm2_mom>0 && sqrt(mm2_mom)>250)) return false;
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D1p_pipi("selec_Bd2DststMuNu_D1p_pipi", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511 && Dststtype == 20413)) && muPID == 1. && !ishigher (I ignore ishigher)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==20413)) return false; // my implementation (and I think Phoebe's too) ensures Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID)
    // mm_mom from Phoebe AddB.C line 3136-3141
    TLorentzVector dst_mom_p(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector dst_p(b.dst_TRUEP_X(),b.dst_TRUEP_Y(),b.dst_TRUEP_Z(),b.dst_TRUEP_E());
    double mm2_mom=(dst_mom_p-dst_p).M2();
    if (!(mm2_mom>0 && sqrt(mm2_mom)>250)) return false;
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D2st_pipi("selec_Bd2DststMuNu_D2st_pipi", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==511)) && muPID == 1. && !ishigher && Dststtype==415 (I ignore ishigher)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID())==415)) return false; // my implementation (and I think Phoebe's too) ensures Dststtype==TMath::Abs(Dst_2010_minus_MC_MOTHER_ID)
    // mm_mom from Phoebe AddB.C line 3136-3141
    TLorentzVector dst_mom_p(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector dst_p(b.dst_TRUEP_X(),b.dst_TRUEP_Y(),b.dst_TRUEP_Z(),b.dst_TRUEP_E());
    double mm2_mom=(dst_mom_p-dst_p).M2();
    if (!(mm2_mom>0 && sqrt(mm2_mom)>250)) return false;
    return true;
  });
*/
  NamedFunc selec_Bd2DststMuNu_higher_pipi_Dst2S("selec_Bd2DststMuNu_higher_pipi_Dst2S", [&](const Baby &b){
     // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && (ishigher) && muPID == 1 (I ignore ishigher)
     // I also add a requirement that D* mom be D*(2S) and gdmom be B0
     // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst from Phoebe AddB.C line 2814 (simplified for me)
     if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
     // my own requirements here...
     if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID()==100413))) return false;
     return true;
  });

  NamedFunc selec_Bd2DststMuNu_higher_pipi_D2S("selec_Bd2DststMuNu_higher_pipi_D2S", [&](const Baby &b){
     // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && (ishigher) && muPID == 1 (I ignore ishigher)
     // I also add a requirement that D* mom be D(2S) and gdmom be B0
     // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst from Phoebe AddB.C line 2814 (simplified for me)
     if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
     // my own requirements here...
     if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID()==100411))) return false;
     return true;
   });

   NamedFunc selec_Bd2DststMuNu_higher_pipi_D2750("selec_Bd2DststMuNu_higher_pipi_D2750", [&](const Baby &b){
     // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && (ishigher) && muPID == 1 (I ignore ishigher)
     // I also add a requirement that D* mom be D(2750) and gdmom be B0
     // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst from Phoebe AddB.C line 2814 (simplified for me)
     if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
     // my own requirements here...
     if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID()==415))) return false; // I'd really like to have "ishigher"/"CocktailHigher" here to differentiate from D2*...
     return true;
   });

   NamedFunc selec_Bd2DststMuNu_higher_pipi_D3000("selec_Bd2DststMuNu_higher_pipi_D3000", [&](const Baby &b){
     // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && (ishigher) && muPID == 1 (I ignore ishigher)
     // I also add a requirement that D* mom be D(3000) and gdmom be B0
     // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst from Phoebe AddB.C line 2814 (simplified for me)
     if (!(abs(b.dst_MC_MOTHER_ID())!=511)) return false;
     // my own requirements here...
     if (!(abs(b.dst_MC_GD_MOTHER_ID())==511 && abs(b.dst_MC_MOTHER_ID()==10413))) return false; // I'd really like to have "ishigher"/"CocktailHigher" here to differentiate from D1...
     return true;
  });

  NamedFunc selec_Bs2DststMuNu_Ds1p("selec_Bs2DststMuNu_Ds1p", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==531 && Dststtype == 10433)) && muPID == 1
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==531 || abs(b.dst_MC_GD_MOTHER_ID())==531 || abs(b.dst_MC_GD_GD_MOTHER_ID())==531)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==531)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=531)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==531 && abs(b.dst_MC_MOTHER_ID())==10433)) return false;
    return true;
  });

  NamedFunc selec_Bs2DststMuNu_Ds2st("selec_Bs2DststMuNu_Ds2st", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: (flagBmu > 0.) && JustDst < 1.  && DstOk > 0. && ((Btype==531 && Dststtype == 435)) && muPID == 1
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==531 || abs(b.dst_MC_GD_MOTHER_ID())==531 || abs(b.dst_MC_GD_GD_MOTHER_ID())==531)) return false;
    // flagBmu from Phoebe's AddB.C lines 2726-2754... I simplified a little, I don't think is totally equivalent to Phoebe (mine should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==531)) return false;
    // justDst from Phoebe AddB.C line 2814 (simplified for me)
    if (!(abs(b.dst_MC_MOTHER_ID())!=531)) return false;
    // Dststtype from Phoebe AddB.C 2814-2864... I simplified a bit, I don't think is totally equivalent to Phoebe (mine should be a looser selection)
    if (!(abs(b.dst_MC_GD_MOTHER_ID())==531 && abs(b.dst_MC_MOTHER_ID())==435)) return false;
    return true;
  });

  NamedFunc selec_Bd2DstDXMuNu("selec_Bd2DstDXMuNu", [&](const Baby &b){
    // Phoebe redoHistos_Dst.C: DstOk > 0. && muPID == 1. && flagDoubleD > 0. && flagTauonicD < 1. && (Btype==511)
    // Btype defined in Phoebe's AddB.C lines 2726-2754, then redefined in 2814-2864: quite confusing, and I've simplified here (should be equivalent for my purposes)
    if (!(abs(b.dst_MC_MOTHER_ID())==511 || abs(b.dst_MC_GD_MOTHER_ID())==511 || abs(b.dst_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagDoubleD is defined in Phoebe's AddB.C on lines 2726-2759: it's confusing, and I'll only partially implement it (my selection should be looser)
    if (!((abs(b.dst_MC_MOTHER_ID())==511 && (abs(b.mu_MC_GD_GD_MOTHER_ID())==20433 || abs(b.mu_MC_GD_GD_MOTHER_ID())==10433)) || (abs(b.mu_MC_MOTHER_ID())==411 || abs(b.mu_MC_MOTHER_ID())==421 || abs(b.mu_MC_MOTHER_ID())==431))) return false;
    // flagTauonicD is defined in Phoebe's AddB.C line 2791
    if (abs(b.mu_MC_MOTHER_ID())==15 && abs(b.mu_MC_GD_MOTHER_ID())==431) return false;
    return true;
  });



  NamedFunc selec_Bu2D0MuNu("selec_Bu2D0MuNu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: flagBmu > 0. && JustDst > 0. && muPID==1. && Y_BKGCAT==0 && Btype==521 && mcut (I think "justdst" should really be called "justd0" here, but whatever)
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if (!(abs(b.d0_MC_MOTHER_ID())==521)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==521)) return false;
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==521 || abs(b.d0_MC_GD_MOTHER_ID())==521 || abs(b.d0_MC_GD_GD_MOTHER_ID())==521)) return false;
    // Y_BKGCAT just my b_BKGCAT
    if (!(b.b_BKGCAT()==0)) return false;
    return true;
  });

  NamedFunc selec_Bu2D0TauNu("selec_Bu2D0TauNu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: flagtaumu > 0. && JustDst > 0. && muPID == 1. && Y_BKGCAT==10 && mcut (I think "justdst" should really be called "justd0" here, but whatever. Also, I expect this line is supposed to have Btype==521 too,so I'll implement this)
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if (!(abs(b.d0_MC_MOTHER_ID())==521)) return false;
    // flagtamu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==15 && abs(b.mu_MC_GD_MOTHER_ID())==521)) return false;
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==521 || abs(b.d0_MC_GD_MOTHER_ID())==521 || abs(b.d0_MC_GD_GD_MOTHER_ID())==521)) return false;
    // Y_BKGCAT just my b_BKGCAT
    if (!(b.b_BKGCAT()==10)) return false;  // remove?
    return true;
  });

  NamedFunc selec_Bu2DstMuNu("selec_Bu2DstMuNu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: flagBmu > 0. && muPID == 1. && Y_BKGCAT==5 && Btype==521 && mcut
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==521)) return false;
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==521 || abs(b.d0_MC_GD_MOTHER_ID())==521 || abs(b.d0_MC_GD_GD_MOTHER_ID())==521)) return false;
    // Y_BKGCAT redefined in Phoebe AddD0B.C line 2579- I'll implement something equivalent for my purposes
    if (!(abs(b.d0_MC_MOTHER_ID())==423)) return false;
    return true;
  });

  NamedFunc selec_Bu2DstTauNu("selec_Bu2DstTauNu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: flagtaumu > 0. && muPID == 1. && Y_BKGCAT==15 && mcut && Btype==521
    // flagtamu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==15 && abs(b.mu_MC_GD_MOTHER_ID())==521)) return false;
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==521 || abs(b.d0_MC_GD_MOTHER_ID())==521 || abs(b.d0_MC_GD_GD_MOTHER_ID())==521)) return false;
    // Y_BKGCAT redefined in Phoebe AddD0B.C line 2580- I'll implement something equivalent for my purposes
    if (!(abs(b.d0_MC_MOTHER_ID())==423)) return false;
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D0st_bu("selec_Bd2DststMuNu_D0st_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==511 && Dststtype == 10411)) && muPID == 1. && mcut && !ishigher (I ignore ishigher)
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==511 || abs(b.d0_MC_GD_MOTHER_ID())==511 || abs(b.d0_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!(abs(b.d0_MC_GD_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==10411)) return false; // D0* only decays to D0 (with only one pion)
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D1_bu("selec_Bd2DststMuNu_D1_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==511 && Dststtype == 10413)) && muPID == 1. && mcut && !ishigher) (I ignore ishigher)
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==511 || abs(b.d0_MC_GD_MOTHER_ID())==511 || abs(b.d0_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==10413)) return false; // D1 only decays to D* (with only one pion)
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D1p_bu("selec_Bd2DststMuNu_D1p_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==511 && Dststtype == 20413)) && muPID == 1. && mcut && !ishigher) (I ignore ishigher)
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==511 || abs(b.d0_MC_GD_MOTHER_ID())==511 || abs(b.d0_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==20413)) return false; // D1' only decays to D* (with only one pion)
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_D2st_bu("selec_Bd2DststMuNu_D2st_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==511 && Dststtype == 415)) && muPID == 1. && mcut &&  !ishigher (I ignore ishigher)
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==511 || abs(b.d0_MC_GD_MOTHER_ID())==511 || abs(b.d0_MC_GD_GD_MOTHER_ID())==511)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!((abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID())==415) || (abs(b.d0_MC_GD_MOTHER_ID())==511 && abs(b.d0_MC_MOTHER_ID())==415))) return false; // D2* can decay to D* or D0 (even with only one pion)
    return true;
  });

  NamedFunc selec_Bd2DststMuNu_higher_pipi_Dst2S_bu("selec_Bd2DststMuNu_higher_pipi_Dst2S_bu", [&](const Baby &b){
     // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  (ishigher) && muPID == 1 && mcut && (TMath::Abs(D0_MC_MOTHER_ID)==413 || TMath::Abs(D0_MC_MOTHER_ID)==423) (I ignore ishigher)
     // I also add a requirement that D gdmom be D*(2S) and gdgdmom be B0
     // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
     if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
     // implement requirement from Phoebe that D0 mom be D* (adjusted for my purposes)
     if (!(abs(b.d0_MC_MOTHER_ID())==413)) return false;
     // my own requirements here...
     if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID()==100413))) return false;
     return true;
  });

  NamedFunc selec_Bd2DststMuNu_higher_pipi_D2S_bu("selec_Bd2DststMuNu_higher_pipi_D2S_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  (ishigher) && muPID == 1 && mcut && (TMath::Abs(D0_MC_MOTHER_ID)==413 || TMath::Abs(D0_MC_MOTHER_ID)==423) (I ignore ishigher)
     // I also add a requirement that D gdmom be D(2S) and gdgdmom be B0
     // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
     if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
     // implement requirement from Phoebe that D0 mom be D* (adjusted for my purposes)
     if (!(abs(b.d0_MC_MOTHER_ID())==413)) return false;
     // my own requirements here...
     if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID()==100411))) return false;
     return true;
   });

   NamedFunc selec_Bd2DststMuNu_higher_pipi_D2750_bu("selec_Bd2DststMuNu_higher_pipi_D2750_bu", [&](const Baby &b){
     // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  (ishigher) && muPID == 1 && mcut && (TMath::Abs(D0_MC_MOTHER_ID)==413 || TMath::Abs(D0_MC_MOTHER_ID)==423) (I ignore ishigher)
     // I also add a requirement that D gdmom be D(2750) and gdgdmom be B0
     // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
     if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
     // implement requirement from Phoebe that D0 mom be D* (adjusted for my purposes)
     if (!(abs(b.d0_MC_MOTHER_ID())==413)) return false;
     // my own requirements here...
     if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID()==415))) return false; // I'd really like to have "ishigher"/"CocktailHigher" here to differentiate from D2*...
     return true;
   });

   NamedFunc selec_Bd2DststMuNu_higher_pipi_D3000_bu("selec_Bd2DststMuNu_higher_pipi_D3000_bu", [&](const Baby &b){
     // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  (ishigher) && muPID == 1 && mcut && (TMath::Abs(D0_MC_MOTHER_ID)==413 || TMath::Abs(D0_MC_MOTHER_ID)==423) (I ignore ishigher)
     // I also add a requirement that D gdmom be D(3000) and gdgdmom be B0
     // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
     if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==511)) return false;
     // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
     if ((abs(b.d0_MC_MOTHER_ID())==511)) return false;
     // implement requirement from Phoebe that D0 mom be D* (adjusted for my purposes)
     if (!(abs(b.d0_MC_MOTHER_ID())==413)) return false;
     // my own requirements here...
     if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==511 && abs(b.d0_MC_GD_MOTHER_ID()==10413))) return false; // I'd really like to have "ishigher"/"CocktailHigher" here to differentiate from D1...
     return true;
  });

  NamedFunc selec_Bs2DststMuNu_Ds1p_bu("selec_Bs2DststMuNu_Ds1p_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==531 && Dststtype == 10433)) && muPID == 1. && mcut
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==531 || abs(b.d0_MC_GD_MOTHER_ID())==531 || abs(b.d0_MC_GD_GD_MOTHER_ID())==531)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==531)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==531)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==531 && abs(b.d0_MC_GD_MOTHER_ID())==10433)) return false;
    return true;
  });

  NamedFunc selec_Bs2DststMuNu_Ds2st_bu("selec_Bs2DststMuNu_Ds2st_bu", [&](const Baby &b){
    // Phoebe redoHistos_D0.C: (flagBmu > 0.) && JustDst < 1.  &&  ((Btype==531 && Dststtype == 435)) && muPID == 1. && mcut
    // Btype is defined in Phoebe AddD0B_temp.C lines 2120-2137, then re-defined in lines 2148-2168: this is too complicated, so I've simplified (not equivalent to Phoebe and I might actually be TIGHTER, but I don't expect it to be significant)
    if (!(abs(b.d0_MC_MOTHER_ID())==531 || abs(b.d0_MC_GD_MOTHER_ID())==531 || abs(b.d0_MC_GD_GD_MOTHER_ID())==531)) return false;
    // flagBmu is from Phoebe AddD0B_temp.C lines 2120-2137: it's a bit complicated, so I've simplified (my selection should be looser)
    if (!(abs(b.mu_TRUEID())==13 && abs(b.mu_MC_MOTHER_ID())==531)) return false;
    // justDst defined in Phoebe AddD0B_temp.C line 2148 (I've simplified it for my purposes)
    if ((abs(b.d0_MC_MOTHER_ID())==531)) return false;
    // Dststtype from Phoebe AddD0B_temp.C 2154-2198... I simplified a bit, I don't think is totally equivalent to Phoebe... not sure if I'm looser or not
    if (!(abs(b.d0_MC_GD_GD_MOTHER_ID())==531 && abs(b.d0_MC_GD_MOTHER_ID())==435)) return false;
    return true;
  });



  /////////////////////////// Define Processes

  string repofolder = "ntuples/0.9.3-production_for_validation/Dst_D0-mc/";

  // first define one process per sample, then combine processes into new process as desired (I want to keep the vector structure for indiv processes in case I decide to plot indiv processes- PlotMaker requires a vector as input, I think)
  vector<shared_ptr<Process>> proc_Bd2DstMuNu;
  vector<shared_ptr<Process>> proc_Bd2DstTauNu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1p;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D2st;
  /*vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1_pipi;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1p_pipi;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D2st_pipi;*/
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_Dst2S;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D2S;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D2750;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D3000;
  vector<shared_ptr<Process>> proc_Bs2DststMuNu_Ds1p;
  vector<shared_ptr<Process>> proc_Bs2DststMuNu_Ds2st;
  vector<shared_ptr<Process>> proc_Bd2DstDXMuNu;

  vector<shared_ptr<Process>> proc_Bu2D0MuNu;
  vector<shared_ptr<Process>> proc_Bu2D0TauNu;
  vector<shared_ptr<Process>> proc_Bu2DstMuNu;
  vector<shared_ptr<Process>> proc_Bu2DstTauNu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D0st_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D1p_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_D2st_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_Dst2S_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D2S_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D2750_bu;
  vector<shared_ptr<Process>> proc_Bd2DststMuNu_higher_pipi_D3000_bu;
  vector<shared_ptr<Process>> proc_Bs2DststMuNu_Ds1p_bu;
  vector<shared_ptr<Process>> proc_Bs2DststMuNu_Ds2st_bu;


  proc_Bd2DstMuNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*} #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11574021_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DstMuNu));
  proc_Bd2DstTauNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*} #tau #nu", Process::Type::signal, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11574011_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DstTauNu));
  proc_Bd2DststMuNu_D1.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D_{1}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D1));
  proc_Bd2DststMuNu_D1p.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D'_{1}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D1p));
  proc_Bd2DststMuNu_D2st.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D_{2}^{*}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D2st));
/*
  proc_Bd2DststMuNu_D1_pipi.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D_{1}[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D1_pipi));
  proc_Bd2DststMuNu_D1p_pipi.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D'_{1}[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D1p_pipi));
  proc_Bd2DststMuNu_D2st_pipi.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D_{2}^{*}[#rightarrow D^{*}#pi#pi||D^{*0}] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_D2st_pipi));
*/
  proc_Bd2DststMuNu_higher_pipi_Dst2S.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*}(2S)[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_higher_pipi_Dst2S));
  proc_Bd2DststMuNu_higher_pipi_D2S.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D(2S)[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_higher_pipi_D2S));
  proc_Bd2DststMuNu_higher_pipi_D2750.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D(2750)[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_higher_pipi_D2750));
  proc_Bd2DststMuNu_higher_pipi_D3000.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D(3000)[#rightarrow D^{*}#pi#pi] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DststMuNu_higher_pipi_D3000));
  proc_Bs2DststMuNu_Ds1p.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D'_{s1}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_13674000_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bs2DststMuNu_Ds1p));
  proc_Bs2DststMuNu_Ds2st.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D_{s2}^{*}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_13674000_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bs2DststMuNu_Ds2st));
  proc_Bd2DstDXMuNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*} X_{c}[#rightarrow #mu #nu X'] X", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11894610_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DstDXMuNu));

  proc_Bu2D0MuNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{-} #rightarrow D^{0} #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_12573012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bu2D0MuNu));
  proc_Bu2D0TauNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{-} #rightarrow D^{0} #tau #nu", Process::Type::signal, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_12573001_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bu2D0TauNu));
  proc_Bu2DstMuNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{-} #rightarrow D^{*0} #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_12773410_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bu2DstMuNu));
  proc_Bu2DstTauNu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{-} #rightarrow D^{*0} #tau #nu", Process::Type::signal, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_12773400_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bu2DstTauNu));
  proc_Bd2DststMuNu_D0st_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D_{0}^{*}[#rightarrow D^{0}] #mu #nu", Process::Type::background, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_D0st_bu));
  proc_Bd2DststMuNu_D1_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D_{1}[#rightarrow D^{*}||D_{0}^{*}] #mu #nu", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_D1_bu));
  proc_Bd2DststMuNu_D1p_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D'_{1}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_D1p_bu));
  proc_Bd2DststMuNu_D2st_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D_{2}^{*}[#rightarrow D^{*}||D_{0}^{*}||D^{0}] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11874430_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_D2st_bu));
  proc_Bd2DststMuNu_higher_pipi_Dst2S_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D^{*}(2S)[#rightarrow D^{*}[#rightarrow D^{0}#pi]#pi#pi] #mu #nu", Process::Type::background, colors("red"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_higher_pipi_Dst2S_bu));
  proc_Bd2DststMuNu_higher_pipi_D2S_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D(2S)[#rightarrow D^{*}[#rightarrow D^{0}#pi]#pi#pi] #mu #nu", Process::Type::background, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_higher_pipi_D2S_bu));
  proc_Bd2DststMuNu_higher_pipi_D2750_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D(2750)[#rightarrow D^{*}[#rightarrow D^{0}#pi]#pi#pi] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_higher_pipi_D2750_bu));
  proc_Bd2DststMuNu_higher_pipi_D3000_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D(3000)[#rightarrow D^{*}[#rightarrow D^{0}#pi]#pi#pi] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11676012_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bd2DststMuNu_higher_pipi_D3000_bu));
  proc_Bs2DststMuNu_Ds1p_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D'_{s1}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("orange"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_13674000_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bs2DststMuNu_Ds1p_bu));
  proc_Bs2DststMuNu_Ds2st_bu.push_back(Process::MakeShared<Baby_run2_fullsim_Bu>("B^{0} #rightarrow D_{s2}^{*}[#rightarrow D^{*}] #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_15--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_13674000_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bu&&selec_Bs2DststMuNu_Ds2st_bu));


  // define multiple-process processes
  vector<shared_ptr<Process>> procs_Bd2DstMuNu_Bd2DstTauNu(proc_Bd2DstMuNu.begin(),proc_Bd2DstMuNu.end());
  procs_Bd2DstMuNu_Bd2DstTauNu.insert(procs_Bd2DstMuNu_Bd2DstTauNu.end(),proc_Bd2DstTauNu.begin(),proc_Bd2DstTauNu.end());
  vector<shared_ptr<Process>> procs_Bu2D0MuNu_Bu2D0TauNu(proc_Bu2D0MuNu.begin(),proc_Bu2D0MuNu.end());
  procs_Bu2D0MuNu_Bu2D0TauNu.insert(procs_Bu2D0MuNu_Bu2D0TauNu.end(),proc_Bu2D0TauNu.begin(),proc_Bu2D0TauNu.end());
  vector<shared_ptr<Process>> procs_Bu2DstMuNu_Bu2DstTauNu(proc_Bu2DstMuNu.begin(),proc_Bu2DstMuNu.end());
  procs_Bu2DstMuNu_Bu2DstTauNu.insert(procs_Bu2DstMuNu_Bu2DstTauNu.end(),proc_Bu2DstTauNu.begin(),proc_Bu2DstTauNu.end());
  vector<shared_ptr<Process>> procs_Bd2DststMuNu(proc_Bd2DststMuNu_D1.begin(),proc_Bd2DststMuNu_D1.end());
  procs_Bd2DststMuNu.insert(procs_Bd2DststMuNu.end(),proc_Bd2DststMuNu_D1p.begin(),proc_Bd2DststMuNu_D1p.end());
  procs_Bd2DststMuNu.insert(procs_Bd2DststMuNu.end(),proc_Bd2DststMuNu_D2st.begin(),proc_Bd2DststMuNu_D2st.end());
  /*vector<shared_ptr<Process>> procs_Bd2DststMuNu_pipi(proc_Bd2DststMuNu_D1_pipi.begin(),proc_Bd2DststMuNu_D1_pipi.end());
  procs_Bd2DststMuNu_pipi.insert(procs_Bd2DststMuNu_pipi.end(),proc_Bd2DststMuNu_D1p_pipi.begin(),proc_Bd2DststMuNu_D1p_pipi.end());
  procs_Bd2DststMuNu_pipi.insert(procs_Bd2DststMuNu_pipi.end(),proc_Bd2DststMuNu_D2st_pipi.begin(),proc_Bd2DststMuNu_D2st_pipi.end());*/
  vector<shared_ptr<Process>> procs_Bd2DststMuNu_higher_pipi(proc_Bd2DststMuNu_higher_pipi_Dst2S.begin(),proc_Bd2DststMuNu_higher_pipi_Dst2S.end());
  procs_Bd2DststMuNu_higher_pipi.insert(procs_Bd2DststMuNu_higher_pipi.end(),proc_Bd2DststMuNu_higher_pipi_D2S.begin(),proc_Bd2DststMuNu_higher_pipi_D2S.end());
  procs_Bd2DststMuNu_higher_pipi.insert(procs_Bd2DststMuNu_higher_pipi.end(),proc_Bd2DststMuNu_higher_pipi_D2750.begin(),proc_Bd2DststMuNu_higher_pipi_D2750.end());
  procs_Bd2DststMuNu_higher_pipi.insert(procs_Bd2DststMuNu_higher_pipi.end(),proc_Bd2DststMuNu_higher_pipi_D3000.begin(),proc_Bd2DststMuNu_higher_pipi_D3000.end());
  vector<shared_ptr<Process>> procs_Bs2DststMuNu(proc_Bs2DststMuNu_Ds1p.begin(),proc_Bs2DststMuNu_Ds1p.end());
  procs_Bs2DststMuNu.insert(procs_Bs2DststMuNu.end(),proc_Bs2DststMuNu_Ds2st.begin(),proc_Bs2DststMuNu_Ds2st.end());
  vector<shared_ptr<Process>> procs_Bd2DststMuNu_bu(proc_Bd2DststMuNu_D0st_bu.begin(),proc_Bd2DststMuNu_D0st_bu.end());
  procs_Bd2DststMuNu_bu.insert(procs_Bd2DststMuNu_bu.end(),proc_Bd2DststMuNu_D1_bu.begin(),proc_Bd2DststMuNu_D1_bu.end());
  procs_Bd2DststMuNu_bu.insert(procs_Bd2DststMuNu_bu.end(),proc_Bd2DststMuNu_D1p_bu.begin(),proc_Bd2DststMuNu_D1p_bu.end());
  procs_Bd2DststMuNu_bu.insert(procs_Bd2DststMuNu_bu.end(),proc_Bd2DststMuNu_D2st_bu.begin(),proc_Bd2DststMuNu_D2st_bu.end());
  vector<shared_ptr<Process>> procs_Bd2DststMuNu_higher_pipi_bu(proc_Bd2DststMuNu_higher_pipi_Dst2S_bu.begin(),proc_Bd2DststMuNu_higher_pipi_Dst2S_bu.end());
  procs_Bd2DststMuNu_higher_pipi_bu.insert(procs_Bd2DststMuNu_higher_pipi_bu.end(),proc_Bd2DststMuNu_higher_pipi_D2S_bu.begin(),proc_Bd2DststMuNu_higher_pipi_D2S_bu.end());
  procs_Bd2DststMuNu_higher_pipi_bu.insert(procs_Bd2DststMuNu_higher_pipi_bu.end(),proc_Bd2DststMuNu_higher_pipi_D2750_bu.begin(),proc_Bd2DststMuNu_higher_pipi_D2750_bu.end());
  procs_Bd2DststMuNu_higher_pipi_bu.insert(procs_Bd2DststMuNu_higher_pipi_bu.end(),proc_Bd2DststMuNu_higher_pipi_D3000_bu.begin(),proc_Bd2DststMuNu_higher_pipi_D3000_bu.end());
  vector<shared_ptr<Process>> procs_Bs2DststMuNu_bu(proc_Bs2DststMuNu_Ds1p_bu.begin(),proc_Bs2DststMuNu_Ds1p_bu.end());
  procs_Bs2DststMuNu_bu.insert(procs_Bs2DststMuNu_bu.end(),proc_Bs2DststMuNu_Ds2st_bu.begin(),proc_Bs2DststMuNu_Ds2st_bu.end());


  ///////////////////////// Define Custom NamedFuncs

  NamedFunc true_q2_Bd2Dst("true_q2_Bd2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Dst(b.dst_TRUEP_X(),b.dst_TRUEP_Y(),b.dst_TRUEP_Z(),b.dst_TRUEP_E());
    return (B-Dst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bd2Dst("true_mm2_Bd2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Dst(b.dst_TRUEP_X(),b.dst_TRUEP_Y(),b.dst_TRUEP_Z(),b.dst_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bd2Dst("true_el_Bd2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bd2Dstst("true_q2_Bd2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==10413||abs(b.dst_MC_MOTHER_ID())==20413||abs(b.dst_MC_MOTHER_ID())==415);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bd2Dstst("true_mm2_Bd2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==10413||abs(b.dst_MC_MOTHER_ID())==20413||abs(b.dst_MC_MOTHER_ID())==415);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bd2Dstst("true_el_Bd2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bd2Dstst_higher_pipi("true_q2_Bd2Dstst_higher_pipi", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==100413||abs(b.dst_MC_MOTHER_ID())==100411||abs(b.dst_MC_MOTHER_ID())==415||abs(b.dst_MC_MOTHER_ID())==10413);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bd2Dstst_higher_pipi("true_mm2_Bd2Dstst_higher_pipi", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==100413||abs(b.dst_MC_MOTHER_ID())==100411||abs(b.dst_MC_MOTHER_ID())==415||abs(b.dst_MC_MOTHER_ID())==10413);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bd2Dstst_higher_pipi("true_el_Bd2Dstst_higher_pipi", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bs2Dstst("true_q2_Bs2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==10433||abs(b.dst_MC_MOTHER_ID())==435);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bs2Dstst("true_mm2_Bs2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    assert(abs(b.dst_MC_MOTHER_ID())==10433||abs(b.dst_MC_MOTHER_ID())==435);
    TLorentzVector Dstst(b.dst_MC_MOTHER_TRUEPX(),b.dst_MC_MOTHER_TRUEPY(),b.dst_MC_MOTHER_TRUEPZ(),b.dst_MC_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bs2Dstst("true_el_Bs2Dstst", [&](const Baby &b) {
    TLorentzVector B(b.b0_TRUEP_X(),b.b0_TRUEP_Y(),b.b0_TRUEP_Z(),b.b0_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });



  NamedFunc true_q2_Bu2D0("true_q2_Bu2D0", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector D0(b.d0_TRUEP_X(),b.d0_TRUEP_Y(),b.d0_TRUEP_Z(),b.d0_TRUEP_E());
    return (B-D0).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bu2D0("true_mm2_Bu2D0", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector D0(b.d0_TRUEP_X(),b.d0_TRUEP_Y(),b.d0_TRUEP_Z(),b.d0_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-D0-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bu2D0("true_el_Bu2D0", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bu2Dst("true_q2_Bu2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_MOTHER_ID())==423); // require d0 mom to be D*0
    TLorentzVector Dst(b.d0_MC_MOTHER_TRUEPX(),b.d0_MC_MOTHER_TRUEPY(),b.d0_MC_MOTHER_TRUEPZ(),b.d0_MC_MOTHER_TRUEPE());
    return (B-Dst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bu2Dst("true_mm2_Bu2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_MOTHER_ID())==423); // require d0 mom to be D*0
    TLorentzVector Dst(b.d0_MC_MOTHER_TRUEPX(),b.d0_MC_MOTHER_TRUEPY(),b.d0_MC_MOTHER_TRUEPZ(),b.d0_MC_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bu2Dst("true_el_Bu2Dst", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bd2Dstst_bu("true_q2_Bd2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_MOTHER_ID())==10411||abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==20413||abs(b.d0_MC_MOTHER_ID())==415||abs(b.d0_MC_GD_MOTHER_ID())==415);
    TLorentzVector Dstst;
    bool mom_cond=false;
    if (abs(b.d0_MC_MOTHER_ID())==10411 || abs(b.d0_MC_MOTHER_ID())==415) {
      Dstst.SetPxPyPzE(b.d0_MC_MOTHER_TRUEPX(),b.d0_MC_MOTHER_TRUEPY(),b.d0_MC_MOTHER_TRUEPZ(),b.d0_MC_MOTHER_TRUEPE());
      mom_cond=true;
      if (abs(b.d0_MC_MOTHER_ID())==10411 && (abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==415)) mom_cond=false; // this is ok, and the Dstst will get reset correctly in next if statement
    }
    if (abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==20413||abs(b.d0_MC_GD_MOTHER_ID())==415) {
      assert(!mom_cond);
      Dstst.SetPxPyPzE(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    }
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bd2Dstst_bu("true_mm2_Bd2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_MOTHER_ID())==10411||abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==20413||abs(b.d0_MC_MOTHER_ID())==415||abs(b.d0_MC_GD_MOTHER_ID())==415);
    TLorentzVector Dstst;
    bool mom_cond=false;
    if (abs(b.d0_MC_MOTHER_ID())==10411 || abs(b.d0_MC_MOTHER_ID())==415) {
      Dstst.SetPxPyPzE(b.d0_MC_MOTHER_TRUEPX(),b.d0_MC_MOTHER_TRUEPY(),b.d0_MC_MOTHER_TRUEPZ(),b.d0_MC_MOTHER_TRUEPE());
      mom_cond=true;
      if (abs(b.d0_MC_MOTHER_ID())==10411 && (abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==415)) mom_cond=false; // this is ok, and the Dstst will get reset correctly in next if statement
    }
    if (abs(b.d0_MC_GD_MOTHER_ID())==10413||abs(b.d0_MC_GD_MOTHER_ID())==20413||abs(b.d0_MC_GD_MOTHER_ID())==415) {
      assert(!mom_cond);
      Dstst.SetPxPyPzE(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    }
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bd2Dstst_bu("true_el_Bd2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bd2Dstst_higher_pipi_bu("true_q2_Bd2Dstst_higher_pipi_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_GD_MOTHER_ID())==100413||abs(b.d0_MC_GD_MOTHER_ID())==100411||abs(b.d0_MC_GD_MOTHER_ID())==415||abs(b.d0_MC_GD_MOTHER_ID())==10413);
    TLorentzVector Dstst(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bd2Dstst_higher_pipi_bu("true_mm2_Bd2Dstst_higher_pipi_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_GD_MOTHER_ID())==100413||abs(b.d0_MC_GD_MOTHER_ID())==100411||abs(b.d0_MC_GD_MOTHER_ID())==415||abs(b.d0_MC_GD_MOTHER_ID())==10413);
    TLorentzVector Dstst(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bd2Dstst_higher_pipi_bu("true_el_Bd2Dstst_higher_pipi_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });

  NamedFunc true_q2_Bs2Dstst_bu("true_q2_Bs2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_GD_MOTHER_ID())==10433||abs(b.d0_MC_GD_MOTHER_ID())==435);
    TLorentzVector Dstst(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    return (B-Dstst).M2()/1e6; // in GeV^2
  });

  NamedFunc true_mm2_Bs2Dstst_bu("true_mm2_Bs2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    assert(abs(b.d0_MC_GD_MOTHER_ID())==10433||abs(b.d0_MC_GD_MOTHER_ID())==435);
    TLorentzVector Dstst(b.d0_MC_GD_MOTHER_TRUEPX(),b.d0_MC_GD_MOTHER_TRUEPY(),b.d0_MC_GD_MOTHER_TRUEPZ(),b.d0_MC_GD_MOTHER_TRUEPE());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    return (B-Dstst-Mu).M2()/1e6; // in GeV^2
  });

  NamedFunc true_el_Bs2Dstst_bu("true_el_Bs2Dstst_bu", [&](const Baby &b) {
    TLorentzVector B(b.b_TRUEP_X(),b.b_TRUEP_Y(),b.b_TRUEP_Z(),b.b_TRUEP_E());
    TLorentzVector Mu(b.mu_TRUEP_X(),b.mu_TRUEP_Y(),b.mu_TRUEP_Z(),b.mu_TRUEP_E());
    Mu.Boost(-B.BoostVector());
    return Mu.E()/1e3; // in GeV
  });


  NamedFunc fit_q2("fit_q2", [&](const Baby &b) {
    return b.FitVar_q2()/1e6; // in GeV^2
  });

  NamedFunc fit_mm2("fit_mm2", [&](const Baby &b) {
    return b.FitVar_Mmiss2()/1e6; // in GeV^2
  });

  NamedFunc fit_el("fit_el", [&](const Baby &b) {
    return b.FitVar_El()/1e3; // in GeV
  });


  NamedFunc b0_trueid("b0_trueid", [&](const Baby &b){
    return abs(b.b0_TRUEID());
  });

  NamedFunc dst_trueid("dst_trueid", [&](const Baby &b){
    return abs(b.dst_TRUEID());
  });

  NamedFunc mu_trueid("mu_trueid", [&](const Baby &b){
    return abs(b.mu_TRUEID());
  });

  NamedFunc d0_trueid("d0_trueid", [&](const Baby &b){
    return abs(b.d0_TRUEID());
  });

  NamedFunc spi_trueid("spi_trueid", [&](const Baby &b){
    return abs(b.spi_TRUEID());
  });

  NamedFunc k_trueid("k_trueid", [&](const Baby &b){
    return abs(b.k_TRUEID());
  });

  NamedFunc pi_trueid("pi_trueid", [&](const Baby &b){
    return abs(b.pi_TRUEID());
  });

  NamedFunc b0_mom_id("b0_mom_id", [&](const Baby &b){
    return abs(b.b0_MC_MOTHER_ID());
  });

  NamedFunc dst_mom_id("dst_mom_id", [&](const Baby &b){
    return abs(b.dst_MC_MOTHER_ID());
  });

  NamedFunc mu_mom_id("mu_mom_id", [&](const Baby &b){
    return abs(b.mu_MC_MOTHER_ID());
  });

  NamedFunc d0_mom_id("d0_mom_id", [&](const Baby &b){
    return abs(b.d0_MC_MOTHER_ID());
  });

  NamedFunc k_mom_id("k_mom_id", [&](const Baby &b){
    return abs(b.k_MC_MOTHER_ID());
  });

  NamedFunc pi_mom_id("pi_mom_id", [&](const Baby &b){
    return abs(b.pi_MC_MOTHER_ID());
  });

  NamedFunc spi_mom_id("spi_mom_id", [&](const Baby &b){
    return abs(b.spi_MC_MOTHER_ID());
  });

  NamedFunc b_trueid("b_trueid", [&](const Baby &b){
    return abs(b.b_TRUEID());
  });

  NamedFunc b_mom_id("b_mom_id", [&](const Baby &b){
    return abs(b.b_MC_MOTHER_ID());
  });



  ///////////////////// Make Plots

  PlotMaker pm;

  // Convention: end the tag with _bd if reconstructed as B0 -> D*+ [-> D0 [-> K- pi+] spi+] mu-, _bu if reconstructed as B- -> D0 [-> K- pi+] mu-
  vector<proc_tag> vec_proc_tag;
  proc_tag Bd2DstMuNu_Bd2DstTauNu(procs_Bd2DstMuNu_Bd2DstTauNu, "Bd2DstMuTauNu_bd"); vec_proc_tag.push_back(Bd2DstMuNu_Bd2DstTauNu);
  proc_tag Bu2D0MuNu_Bu2D0TauNu(procs_Bu2D0MuNu_Bu2D0TauNu, "Bu2D0MuTauNu_bu"); vec_proc_tag.push_back(Bu2D0MuNu_Bu2D0TauNu);
  proc_tag Bu2DstMuNu_Bu2DstTauNu(procs_Bu2DstMuNu_Bu2DstTauNu, "Bu2DstMuTauNu_bu"); vec_proc_tag.push_back(Bu2DstMuNu_Bu2DstTauNu);
  proc_tag Bd2DststMuNu(procs_Bd2DststMuNu, "Bd2DststMuNu_bd"); vec_proc_tag.push_back(Bd2DststMuNu);
  //proc_tag Bd2DststMuNu_pipi(procs_Bd2DststMuNu_pipi, "Bd2DststMuNu_pipi_bd"); vec_proc_tag.push_back(Bd2DststMuNu_pipi);
  proc_tag Bd2DststMuNu_higher_pipi(procs_Bd2DststMuNu_higher_pipi, "Bd2DststMuNu_higher_pipi_bd"); vec_proc_tag.push_back(Bd2DststMuNu_higher_pipi);
  proc_tag Bs2DststMuNu(procs_Bs2DststMuNu, "Bs2DststMuNu_bd"); vec_proc_tag.push_back(Bs2DststMuNu);
  proc_tag Bd2DstDXMuNu(proc_Bd2DstDXMuNu, "Bd2DstDXMuNu_bd"); vec_proc_tag.push_back(Bd2DstDXMuNu);
  proc_tag Bd2DststMuNu_bu(procs_Bd2DststMuNu_bu, "Bd2DststMuNu_bu"); vec_proc_tag.push_back(Bd2DststMuNu_bu);
  proc_tag Bd2DststMuNu_higher_pipi_bu(procs_Bd2DststMuNu_higher_pipi_bu, "Bd2DststMuNu_higher_pipi_bu"); vec_proc_tag.push_back(Bd2DststMuNu_higher_pipi_bu);
  proc_tag Bs2DststMuNu_bu(procs_Bs2DststMuNu_bu, "Bs2DststMuNu_bu"); vec_proc_tag.push_back(Bs2DststMuNu_bu);

  for (vector<proc_tag>::iterator pt = vec_proc_tag.begin(); pt != vec_proc_tag.end(); pt++) {
    vector<shared_ptr<Process>> processes = pt->Processes();
    string filetag = pt->FileTag();
    Double_t d_dst_mom_upper;
    if (filetag.find("higher") != string::npos) d_dst_mom_upper=118000;
    else d_dst_mom_upper=22000;
    if (filetag.find("_bd") != string::npos) { // reconstructed as B0 -> D*+ [-> D0 [-> K- pi+] spi+] mu-
      // ID plots
      pm.Push<Hist1D>(Axis(240,-2000,22000, b0_trueid,"|B TrueID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(240,-2000,22000, dst_trueid,"|D^{*} TrueID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(210,-200,4000, spi_trueid,"|#pi_{s} TrueID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(240,-2000,22000, b0_mom_id,"|B Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(240,-2000,d_dst_mom_upper, dst_mom_id,"|D^{*} Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(240,-2000,22000, spi_mom_id,"|#pi_{s} Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
      // q2, mm2, El
      if (filetag=="Bd2DstMuTauNu_bd") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dst, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dst, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dst, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bd2DststMuNu_bd") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dstst, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dstst, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dstst, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bd2DststMuNu_higher_pipi_bd") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dstst_higher_pipi, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dstst_higher_pipi, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dstst_higher_pipi, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bs2DststMuNu_bd") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bs2Dstst, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bs2Dstst, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bs2Dstst, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
    }
    else if (filetag.find("_bu") != string::npos) { // reconstructed as B- -> D0 [-> K- pi+] mu-
      // ID plots
      pm.Push<Hist1D>(Axis(240,-2000,22000, b_trueid,"|B^{-} TrueID|"), "1", processes, plotloglumi).Tag(filetag);
      pm.Push<Hist1D>(Axis(240,-2000,22000, b_mom_id,"|B^{-} Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
      // q2, mm2, El
      if (filetag=="Bu2D0MuTauNu_bu") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bu2D0, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bu2D0, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bu2D0, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bu2DstMuTauNu_bu") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bu2Dst, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bu2Dst, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bu2Dst, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bd2DststMuNu_bu") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dstst_bu, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dstst_bu, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dstst_bu, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bd2DststMuNu_higher_pipi_bu") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dstst_higher_pipi_bu, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dstst_higher_pipi_bu, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dstst_higher_pipi_bu, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
      if (filetag=="Bs2DststMuNu_bu") {
        pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bs2Dstst_bu, "True m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bs2Dstst_bu, "True q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
        pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bs2Dstst_bu, "True E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
      }
    }
    // Everyone should have these ID plots
    pm.Push<Hist1D>(Axis(300,-20,580, mu_trueid,"|#mu TrueID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(250,-20,480, d0_trueid,"|D TrueID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(210,-200,4000, pi_trueid,"|#pi TrueID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(185,-20,350, k_trueid,"|K TrueID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(240,-2000,22000, mu_mom_id,"|#mu Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(240,-2000,d_dst_mom_upper, d0_mom_id,"|D Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(250,-20,480, k_mom_id,"|K Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
    pm.Push<Hist1D>(Axis(250,-20,480, pi_mom_id,"|#pi Mom ID|"), "1", processes, plotloglumi).Tag(filetag);
    // Everyone should have these fitvars
    pm.Push<Hist1D>(Axis(40,-2,10, fit_mm2, "m_{miss}^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
    pm.Push<Hist1D>(Axis(45,-3,12, fit_q2, "q^{2} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
    pm.Push<Hist1D>(Axis(32,0.1,2.5, fit_el, "E_{#mu}^{*} [GeV^{2}]"), "1", processes, plotshapes).Tag(filetag);
  }


  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
