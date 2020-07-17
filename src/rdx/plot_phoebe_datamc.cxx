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

  ////////////////////////////////////// Phoebe's plots ///////////////////////////////////////

  PlotOpt lin_lumi("txt/plot_styles.txt", "LHCbPaper");
  lin_lumi.Title(TitleType::info)
    .Bottom(BottomType::off)
    .YAxis(YAxisType::linear)
    .Stack(StackType::signal_on_top)
    .Overflow(OverflowType::none);
  
  vector<PlotOpt> plottypes = {lin_lumi};

  //////// Signal, normalization, D** processes for mmiss2 plot
  string repofolder = "ntuples/";
  repofolder = "/Volumes/manuelf_t5/lhcb-ntuples-gen/ntuples/";
  string ntuplefile = "ref-rdx-run1/Dst-mix/Dst--20_07_02--mix--all--2011-2012--md-mu--phoebe.root";
  
  vector<shared_ptr<Process> > procs;
  procs.push_back(Process::MakeShared<Baby_phoebe_dsp>("Data",Process::Type::data, colors("data"),
                                                       set<string>({repofolder+ntuplefile}),
                                                       "(DstIDprod > 0 && IDprod > 0 && muPID > 0. && flagDstSB==0.) && isData > 0."));
  
  PlotMaker pm;
  pm.Push<Hist1D>(Axis(1000, -2, 10,"m_nu1", "m_{miss}^{2} [GeV^{2}]"), "1",
                  procs, plottypes);
  pm.MakePlots(1); // The "1" is the luminosity to rescale the bkg to  


  time(&endtime);
  cout<<endl<<"Making plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

