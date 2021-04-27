// Code for making plots for run2 trackeronly MC validation (using plot_scripts) using D* normalization mode

// Created: Apr 26, 2021
// Last edited: Apr 26, 2021
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


int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  // Start measuring time
  time_t begtime, endtime;
  time(&begtime);

  // User defined colors
  Palette colors("txt/colors.txt", "default");

  //////////////////////////////////////////////////// Define Plot Types ///////////////////////////////////////////

  PlotOpt lin_shapes("txt/plot_styles.txt", "LHCbPaper");
  lin_shapes.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::shapes)
    .Overflow(OverflowType::none);//both);

  PlotOpt lin_lumi_shapes = lin_shapes().Stack(StackType::lumi_shapes).Overflow(OverflowType::both);
  PlotOpt log_lumi_shapes = lin_lumi_shapes().YAxis(YAxisType::log);
  PlotOpt lin_shapes_ratio = lin_shapes().Bottom(BottomType::ratio);

  vector<PlotOpt> plotshapes = {lin_shapes};
  vector<PlotOpt> plotshapesratio = {lin_shapes_ratio};
  vector<PlotOpt> plotlumi = {lin_lumi_shapes};
  vector<PlotOpt> plotloglumi = {log_lumi_shapes};
  vector<PlotOpt> plotboth = {lin_shapes, lin_lumi_shapes};

  ///////////////////////////////////////////////// Define Cuts //////////////////////////////////////////////////

  // Note: Phoebe has some basic cuts that she applies to all templates; see, for example, my plot_phoebe_datamc.cxx script. I am not going to include these, instead
  // opting to implement a subset of her cuts. In particular, I will focus on cuts that help define the specific decay for the templates, and for the fit samples, I
  // will focus on cuts that make the samples distinct from one another.

  // General Cuts function, to be applied to EVERY decay that is reconstructed using TupleB0
  NamedFunc cuts_bd("cuts_bd", [&](const Baby &b){
    // DstOk from Phoebe AddB.C line 2547
    if (!(b.dst_BKGCAT()==0 || (b.dst_BKGCAT()==50 && b.d0_BKGCAT()==50))) return false;
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


  // ISO cuts mentioned at line 1284
  NamedFunc isocuts("ISO", [&](const Baby &b){
    return (b.b0_ISOLATION_BDT() < 0.15); // simplified slightly from Phoebe's line 1284 (she's doing something with the higher D** states here... not sure why)
  });

  // DD cuts starts at line 1294
  // Note: our ntuples seem to be missing variables "iso_NNkw" that Phoebe makes use of... I'm not sure what these are (I think the
  // iso_NNk in general are an output from a trained neural network ["NN"] that give a probability that a track is a kaon ["k"];
  // in this case, the numbers [iso_NNl1,2,3] seem to refer to different tracks selected by the isolation ["iso"] of the B), but
  // I have no choice but to ignore these variables in my selections.
  NamedFunc ddcuts("DD", [&](const Baby &b){ // this differs (more complicated) from table 15 of ANA note (July 9)
    // if (b.iso_NNkw() > -1) return (b.b0_ISOLATION_BDT() > 0.15); // based on lines 1303-1320

    Double_t iso_NNk = b.b0_ISOLATION_NNk();
    Double_t iso_NNk2 = b.b0_ISOLATION_NNk2();
    Double_t iso_NNk3 = b.b0_ISOLATION_NNk3();
    Float_t iso_P = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY()+b.b0_ISOLATION_PZ()*b.b0_ISOLATION_PZ());
    Float_t iso_PT = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY());
    Float_t iso_P2 = sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2()+b.b0_ISOLATION_PZ2()*b.b0_ISOLATION_PZ2());
    Float_t iso_PT2 = sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2());
    Float_t iso_P3 = sqrt(b.b0_ISOLATION_PX3()*b.b0_ISOLATION_PX3()+b.b0_ISOLATION_PY3()*b.b0_ISOLATION_PY3()+b.b0_ISOLATION_PZ3()*b.b0_ISOLATION_PZ3());
    Float_t iso_PT3 = sqrt(b.b0_ISOLATION_PX3()*b.b0_ISOLATION_PX3()+b.b0_ISOLATION_PY3()*b.b0_ISOLATION_PY3());
    Double_t iso_BDT = b.b0_ISOLATION_BDT();
    Double_t iso_BDT2 = b.b0_ISOLATION_BDT2();
    Double_t iso_BDT3 = b.b0_ISOLATION_BDT3();
    if (iso_BDT <= -1.1) iso_NNk = 0.;
    if (iso_BDT2 <= -1.1) iso_NNk2 = 0.;
    if (iso_BDT3 <= -1.1) iso_NNk3 = 0.;

    return (b.b0_ISOLATION_BDT() > 0.15 && (iso_NNk > 0.2 || iso_NNk2 > 0.2 || iso_NNk3 > 0.2)
            && TMath::Max(iso_P*(iso_PT > 150),TMath::Max(iso_P2*(iso_PT2 > 150)*(iso_BDT2 > -1.1),
                          iso_P3*(iso_PT3 > 150)*(iso_BDT3 > -1.1))) > 5e3);
  });

  // 2OS cuts starts at line 1326 (again, I'll ignore anything with iso_NNkw variables)
  NamedFunc twooscuts("2OS", [&](const Baby &b){ // this differs (more complicated) from table 15 of ANA note (July 9)
                                                 // Actually, I think Phoebe's code is logically incorrect, too. Lines 1342
                                                 // and 1343 shouldn't be inside the "else" (as then they'll only be applied
                                                 // if "thecut" is already false)
     Double_t iso_NNk = b.b0_ISOLATION_NNk();
     Double_t iso_NNk2 = b.b0_ISOLATION_NNk2();
     Float_t iso_P = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY()+b.b0_ISOLATION_PZ()*b.b0_ISOLATION_PZ());
     Float_t iso_PT = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY());
     Float_t iso_P2 = sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2()+b.b0_ISOLATION_PZ2()*b.b0_ISOLATION_PZ2());
     Float_t iso_PT2 = sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2());
     Double_t iso_BDT = b.b0_ISOLATION_BDT();
     Double_t iso_BDT2 = b.b0_ISOLATION_BDT2();
     Double_t iso_BDT3 = b.b0_ISOLATION_BDT3();
     Float_t iso_CHARGE = b.b0_ISOLATION_CHARGE();
     Float_t iso_CHARGE2 = b.b0_ISOLATION_CHARGE2();

     return (iso_BDT > 0.15 && iso_BDT2 > 0.15 && iso_BDT3 < 0.15 && iso_CHARGE != iso_CHARGE2
             && iso_CHARGE != 0 && iso_CHARGE2 != 0 && iso_CHARGE < 100
             && TMath::Max(iso_P*(iso_PT > 150),iso_P2*(iso_PT2 > 150)) > 5e3
             && iso_NNk < 0.2 && iso_NNk2 < 0.2);
  });

  // 1OS cuts starts at line 1346 (again, I'll ignore anything with iso_NNkw variables)
  NamedFunc oneoscuts("1OS", [&](const Baby &b){ // this differs (more complicated) from table 15 of ANA note (July 9)
                                                 // Actually, I think Phoebe's code is logically incorrect, too. Line 1363
                                                 // shouldn't be inside the "else" (as then it'll only be applied
                                                 // if "thecut" is already false)
     Double_t iso_NNk = b.b0_ISOLATION_NNk();
     Float_t iso_P = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY()+b.b0_ISOLATION_PZ()*b.b0_ISOLATION_PZ());
     Float_t iso_PT = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY());
     Double_t iso_BDT = b.b0_ISOLATION_BDT();
     Double_t iso_BDT2 = b.b0_ISOLATION_BDT2();
     Float_t iso_CHARGE = b.b0_ISOLATION_CHARGE();
     Int_t Dst_ID = b.dst_ID();

     return (iso_BDT > 0.15 && iso_BDT2 < 0.15 && iso_CHARGE*Dst_ID < 0 && iso_P > 5e3
             && iso_PT > 150 && iso_NNk < 0.2);
  });

  // D** cuts just the same as 1OS but with mass cut on D*+pi1... I actually think this is what is implemented
  // in Phoebe's code currently, so remove the mass cut above for 1OS, and just copy the 1OS cut w the mass cut here. Also,
  // use the mass cut given in table 12, not in the code.
  NamedFunc dsscuts("DSS", [&](const Baby &b) {
    Double_t iso_NNk = b.b0_ISOLATION_NNk();
    Float_t iso_P = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY()+b.b0_ISOLATION_PZ()*b.b0_ISOLATION_PZ());
    Float_t iso_PT = sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY());
    Double_t iso_BDT = b.b0_ISOLATION_BDT();
    Double_t iso_BDT2 = b.b0_ISOLATION_BDT2();
    Float_t iso_CHARGE = b.b0_ISOLATION_CHARGE();
    Int_t Dst_ID = b.dst_ID();
    TLorentzVector B_iso(b.b0_ISOLATION_PX(),b.b0_ISOLATION_PY(),b.b0_ISOLATION_PZ(),b.b0_ISOLATION_PE());
    TLorentzVector Dst(b.dst_PX(),b.dst_PY(),b.dst_PZ(),b.dst_PE());
    Double_t iso_DeltaM = (B_iso+Dst).M()-b.dst_M();

    return (iso_BDT > 0.15 && iso_BDT2 < 0.15 && iso_CHARGE*Dst_ID < 0 && (iso_DeltaM > 390
            && iso_DeltaM < 510) && iso_P > 5e3 && iso_PT > 150 && iso_NNk < 0.2);
  });


  ////////////////////////////////////////////////// Define Processes ////////////////////////////////////////////////

  string repofolder = "ntuples/0.9.3-production_for_validation/Dst_D0-mc/";

  // first define one process per sample, then combine processes into new process as desired (I want to keep the vector structure for indiv processes in case I decide to plot indiv processes- PlotMaker requires a vector as input, I think)
  vector<shared_ptr<Process>> proc_Bd2DstMuNu_fullsim;
  vector<shared_ptr<Process>> proc_Bd2DstMuNu_trackeronly;


  proc_Bd2DstMuNu_fullsim.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*} #mu #nu", Process::Type::background, colors("blue"),
                                                set<string>({repofolder+"Dst_D0--21_02_14--mc--MC_2016_Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8_Sim09j_Trig0x6139160F_Reco16_Turbo03a_Filtered_11574021_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DstMuNu));
  proc_Bd2DstMuNu_trackeronly.push_back(Process::MakeShared<Baby_run2_fullsim_Bd>("B^{0} #rightarrow D^{*} #mu #nu (tra. only)", Process::Type::data, colors("green"),
                                                set<string>({repofolder+"Dst_D0--21_03_10--mc--tracker_only--MC_2016_Beam6500GeV-2016-MagDown-TrackerOnly-Nu1.6-25ns-Pythia8_Sim09j_Reco16_Filtered_11574021_D0TAUNU.SAFESTRIPTRIG.DST.root"}), cuts_bd&&selec_Bd2DstMuNu));

  // // define multiple-process processes
  vector<shared_ptr<Process>> procs_Bd2DstMuNu_trackeronly_fullsim(proc_Bd2DstMuNu_fullsim.begin(),proc_Bd2DstMuNu_fullsim.end());
  procs_Bd2DstMuNu_trackeronly_fullsim.insert(procs_Bd2DstMuNu_trackeronly_fullsim.end(),proc_Bd2DstMuNu_trackeronly.begin(),proc_Bd2DstMuNu_trackeronly.end());

  ///////////////////////////////////////////// Define Custom NamedFuncs ////////////////////////////////////////////

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

  NamedFunc fit_q2("fit_q2", [&](const Baby &b) {
    return b.FitVar_q2()/1e6; // in GeV^2
  });

  NamedFunc fit_mm2("fit_mm2", [&](const Baby &b) {
    return b.FitVar_Mmiss2()/1e6; // in GeV^2
  });

  NamedFunc fit_el("fit_el", [&](const Baby &b) {
    return b.FitVar_El()/1e3; // in GeV
  });

  NamedFunc iso_P("iso_P", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY()+b.b0_ISOLATION_PZ()*b.b0_ISOLATION_PZ());
  });

  NamedFunc iso_PT("iso_PT", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX()*b.b0_ISOLATION_PX()+b.b0_ISOLATION_PY()*b.b0_ISOLATION_PY());
  });

  NamedFunc iso_P2("iso_P2", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2()+b.b0_ISOLATION_PZ2()*b.b0_ISOLATION_PZ2());
  });

  NamedFunc iso_PT2("iso_PT2", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX2()*b.b0_ISOLATION_PX2()+b.b0_ISOLATION_PY2()*b.b0_ISOLATION_PY2());
  });

  NamedFunc iso_P3("iso_P3", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX3()*b.b0_ISOLATION_PX3()+b.b0_ISOLATION_PY3()*b.b0_ISOLATION_PY3()+b.b0_ISOLATION_PZ3()*b.b0_ISOLATION_PZ3());
  });

  NamedFunc iso_PT3("iso_PT3", [&](const Baby &b) {
    return sqrt(b.b0_ISOLATION_PX3()*b.b0_ISOLATION_PX3()+b.b0_ISOLATION_PY3()*b.b0_ISOLATION_PY3());
  });

  NamedFunc iso_CHARGE_times_dst_ID("iso_CHARGE_times_dst_ID", [&](const Baby &b) {
    return b.b0_ISOLATION_CHARGE()*b.dst_ID()<0;
  });

  NamedFunc iso_DeltaM("iso_DeltaM", [&](const Baby &b) {
    TLorentzVector B_iso(b.b0_ISOLATION_PX(),b.b0_ISOLATION_PY(),b.b0_ISOLATION_PZ(),b.b0_ISOLATION_PE());
    TLorentzVector Dst(b.dst_PX(),b.dst_PY(),b.dst_PZ(),b.dst_PE());
    return (B_iso+Dst).M()-b.dst_M();
  });


  ///////////////////////////////////////////// Make Plots /////////////////////////////////////////////////

  PlotMaker pm;

  // Convention: end the tag with _bd if reconstructed as B0 -> D*+ [-> D0 [-> K- pi+] spi+] mu-, _bu if reconstructed as B- -> D0 [-> K- pi+] mu-
  vector<proc_tag> vec_proc_tag;
  proc_tag Bd2DstMuNu_trackeronly_fullsim(procs_Bd2DstMuNu_trackeronly_fullsim, "Bd2DstMuNu_trackeronly_fullsim_bd"); vec_proc_tag.push_back(Bd2DstMuNu_trackeronly_fullsim);

  // list the fit samples you want to have fit variables plotted for
  vector<NamedFunc> fitsample_selections{"1", isocuts, ddcuts, twooscuts, oneoscuts, dsscuts};

  for (vector<proc_tag>::iterator pt = vec_proc_tag.begin(); pt != vec_proc_tag.end(); pt++) {
    vector<shared_ptr<Process>> processes = pt->Processes();
    string filetag = pt->FileTag();
    if (filetag.find("_bd") != string::npos) { // reconstructed as B0 -> D*+ [-> D0 [-> K- pi+] spi+] mu-
      // q2, mm2, El
      if (filetag=="Bd2DstMuNu_trackeronly_fullsim_bd") {
        for (Int_t i=0; i<static_cast<Int_t>(fitsample_selections.size()); i++) {
          pm.Push<Hist1D>(Axis(40,-2,10, true_mm2_Bd2Dst, "True m_{miss}^{2} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
          pm.Push<Hist1D>(Axis(45,-3,12, true_q2_Bd2Dst, "True q^{2} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
          pm.Push<Hist1D>(Axis(32,0.1,2.5, true_el_Bd2Dst, "True E_{#mu}^{*} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
        }
      }
      // variables used for fit samples
      pm.Push<Hist1D>(Axis(30,-2,1, "b0_ISOLATION_BDT", "iso_BDT"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(30,-2,1, "b0_ISOLATION_BDT2", "iso_BDT2"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(30,-2,1, "b0_ISOLATION_BDT3", "iso_BDT3"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,-1,1, "b0_ISOLATION_NNk", "iso_NNk"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,-1,1, "b0_ISOLATION_NNk2", "iso_NNk2"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,-1,1, "b0_ISOLATION_NNk3", "iso_NNk3"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,-5,5, "b0_ISOLATION_CHARGE", "iso_CHARGE"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,-5,5, "b0_ISOLATION_CHARGE2", "iso_CHARGE2"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,0,40000, iso_P, "iso_P (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(50,0,1000, iso_PT, "iso_PT (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,0,40000, iso_P2, "iso_P2 (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(50,0,1000, iso_PT2, "iso_PT2 (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(40,0,40000, iso_P3, "iso_P3 (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(50,0,1000, iso_PT3, "iso_PT3 (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(2,0,2, iso_CHARGE_times_dst_ID, "Iso track and D^{*} have opp. charge"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      pm.Push<Hist1D>(Axis(60,-100,2000, iso_DeltaM, "iso_DeltaM (MeV)"), "1", processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      // Everyone should have these fitvars
      for (Int_t i=0; i<static_cast<Int_t>(fitsample_selections.size()); i++) {
        pm.Push<Hist1D>(Axis(40,-2,10, fit_mm2, "m_{miss}^{2} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
        pm.Push<Hist1D>(Axis(45,-3,12, fit_q2, "q^{2} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
        pm.Push<Hist1D>(Axis(32,0.1,2.5, fit_el, "E_{#mu}^{*} [GeV^{2}]"), fitsample_selections[i], processes, plotshapesratio).Tag(filetag).RatioTitle("TrackerOnly","FullSim");
      }
    }
    else if (filetag.find("_bu") != string::npos) { // reconstructed as B- -> D0 [-> K- pi+] mu-
      continue; // edit this if using TupleBminus at all
    }
  }


  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
