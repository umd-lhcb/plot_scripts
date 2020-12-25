#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <getopt.h>

#include "TStyle.h"
#include "TError.h" // Controls error level reporting
#include "TColor.h" // Controls error level reporting
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "core/baby.hpp"
#include "core/named_func.hpp"

using namespace std;

void GetOptions(int argc, char *argv[]);


namespace{
  float lumi = 4.3;
  string example = "search";
}

struct sevent{
  int run;
  long event, entry;
};

template<class BabyType>
vector<sevent> getUniqueEvents(BabyType &baby, NamedFunc &cut, vector<sevent> &repeated);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  string repofolder = "../lhcb-ntuples-gen/gen/";
  Baby_step2 ybaby(std::set<std::string>{repofolder+"run1-Dst_D0-step2/Dst_D0--20_10_12--std--data--2011--md--step2.root"});
  auto activator = ybaby.Activate();
  //Baby_phoebe_step2 pbaby(std::set<std::string>{repofolder+"run1-Dst-step2/Dst--20_07_02--mix--data--2011--md--phoebe-step2.root"});
  Baby_phoebe_step1 pbaby(std::set<std::string>{repofolder+"run1-Dst-step2/Dst--20_09_16--std--data--2011--md--phoebe-step2.root"});
  auto pactivator = pbaby.Activate();

  NamedFunc basecuts("(b_l0_tis || d0_l0had_tos) && mu_l0_tis && mu_pidb > 0 && dst_id_prod > 0 && id_prod > 0");
  NamedFunc pheebcuts(basecuts);

  vector<sevent> prep, yrep;
  vector<sevent> pevents = getUniqueEvents(pbaby, pheebcuts, prep);
  vector<sevent> yevents = getUniqueEvents(ybaby, basecuts, yrep);

  cout<<"Found "<<pevents.size()<<" unique and "<<prep.size()<<" repeated events in Phoebe"<<endl;
  cout<<"Found "<<yevents.size()<<" unique and "<<yrep.size()<<" repeated events in Yipeng"<<endl;
  time(&endtime);
  cout<<endl<<"Getting events took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;

  vector<TString> types({"pcommon", "ycommon", "pdiff", "ydiff"});
  vector<vector<double> > binnings({});

  vector<NamedFunc> variables({"mm2", "el", "q2", "dst_m", "iso_bdt", "d0_m", "mu_p", "b_endvtx_chi2",
      "d0_ip_chi2", "d0_pt"});
  binnings.push_back({104, -3, 10});
  binnings.push_back({100, 0, 2.5});
  binnings.push_back({90, -3, 12});
  binnings.push_back({170, 1980, 2050});
  binnings.push_back({40, -2, 1});
  binnings.push_back({140, 1830, 1910});
  binnings.push_back({100, 0, 120});
  binnings.push_back({100, 0, 30});
  binnings.push_back({100, 0, 1000});
  binnings.push_back({100, 0, 20});
 
  TH1D *histos[4][20];
  for(unsigned ivar=0; ivar<variables.size(); ivar++){
    for(unsigned type=0; type<4; type++){
      TString var(variables[ivar].Name());
      histos[type][ivar] = new TH1D(var+"_"+types[type], var, static_cast<int>(binnings[ivar][0]),  binnings[ivar][1], binnings[ivar][2]);
      histos[type][ivar]->SetMarkerStyle(8);
      histos[type][ivar]->SetFillColor(0);
      histos[type][ivar]->SetLineWidth(2);
     if(type%2 == 0) {
        histos[type][ivar]->SetLineColor(kRed+1);
        histos[type][ivar]->SetMarkerColor(1);
      } else {
        histos[type][ivar]->SetLineColor(4);
        histos[type][ivar]->SetMarkerColor(4);
      }
      if(type>=2) histos[type][ivar]->SetLineStyle(2);
    }
  } // For over histos

  bool common;
  int ncommon = 0;
  for(unsigned ind=0; ind<pevents.size(); ind++){
    common = false;
    for(unsigned yind=0; yind<yevents.size(); yind++){
       if(yevents[yind].run == pevents[ind].run){
        if(yevents[yind].event == pevents[ind].event){
          ncommon++;
          pbaby.GetEntry(pevents[ind].entry);
          ybaby.GetEntry(yevents[yind].entry);
          for(unsigned ivar=0; ivar<variables.size(); ivar++){
            histos[0][ivar]->Fill(variables[ivar].GetScalar(pbaby)); // Filling Phoebe's common events
            // Filling Yipeng's common events
            if(variables[ivar].Name() == "q2") histos[1][ivar]->Fill(variables[ivar].GetScalar(ybaby)/1000000); 
            else histos[1][ivar]->Fill(variables[ivar].GetScalar(ybaby)); 
          }
          common = true;
          break;
        } // if is common
      }
    }//For yevents

    if(!common){
      pbaby.GetEntry(pevents[ind].entry);
      for(unsigned ivar=0; ivar<variables.size(); ivar++)
        histos[2][ivar]->Fill(variables[ivar].GetScalar(pbaby)); // Filling Phoebe's non-common events
    }
  }//For pevents
  // Filling Yipeng's non-common events
  for(unsigned ind=0; ind<yevents.size(); ind++){
    common = false;
    for(unsigned pind=0; pind<pevents.size(); pind++){
       if(pevents[pind].run == yevents[ind].run){
        if(pevents[pind].event == yevents[ind].event){
          common = true;
          break;
        } // if is common
      }
    }//For pevents

    if(!common){
      ybaby.GetEntry(yevents[ind].entry);
      for(unsigned ivar=0; ivar<variables.size(); ivar++)
        // Filling Yipeng's non-common events
        if(variables[ivar].Name() == "q2") histos[3][ivar]->Fill(variables[ivar].GetScalar(ybaby)/1000000); 
        else histos[3][ivar]->Fill(variables[ivar].GetScalar(ybaby)); 
    }
  }//For yevents
     
  gStyle->SetOptStat(0);              // No Stats box
  gStyle->SetPadTickX(1);             // No ticks at the right
  gStyle->SetPadTickY(1);             // No ticks at the top
  TCanvas can("can","histos",1200,800);

  vector<TString> hnames({"Common Phoebe", "Common Yipeng", "Diff. Phoebe", "Diff. Yipeng"});
  vector<TString> loptions({"f", "lep", "l", "l"});
  for(unsigned ivar=0; ivar<variables.size(); ivar++){
    can.cd();
    histos[0][ivar]->SetFillColor(kRed-10);
    histos[0][ivar]->SetLineColor(kRed-10);
    TLegend leg(0.55, 0.6, 0.65, 0.87);
    leg.SetTextSize(0.042); leg.SetTextFont(132); leg.SetFillColor(0); 
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    TLegend leg2(0.12, 0.65, 0.17, 0.87);
    leg2.SetTextSize(0.042); leg2.SetTextFont(132); leg2.SetFillColor(0); 
    leg2.SetFillStyle(0); leg2.SetBorderSize(0);
    float ymax=0;
    for(unsigned type=0; type<4; type++){
      histos[type][ivar]->SetBinContent(1, histos[type][ivar]->GetBinContent(0)+histos[type][ivar]->GetBinContent(1));
      histos[type][ivar]->SetBinContent(histos[type][ivar]->GetNbinsX(),
                                        histos[type][ivar]->GetBinContent(histos[type][ivar]->GetNbinsX())+histos[type][ivar]->GetBinContent(histos[type][ivar]->GetNbinsX()+1));
      histos[type][ivar]->Scale(1/histos[type][ivar]->Integral());
      if(histos[type][ivar]->GetMaximum() > ymax) ymax = histos[type][ivar]->GetMaximum();
      TString options="hist";
      if(type==1) options = "e same";
      if(type>1) options = "hist same";
      histos[type][ivar]->Draw(options);
      TString lname = hnames[type]; lname += " ("; lname += histos[type][ivar]->GetEntries()-2; lname += ")";
      leg.AddEntry(histos[type][ivar], lname, loptions[type]);
      leg2.AddEntry(histos[type][ivar], lname, loptions[type]);
    }
    if(variables[ivar].Name() == "q2" || variables[ivar].Name() == "el" || variables[ivar].Name() == "iso_bdt") leg2.Draw();
    else leg.Draw();
    histos[0][ivar]->SetMaximum(ymax*1.1);
    TString pname("plots/common_"); pname += variables[ivar].Name(); pname += ".pdf";
    can.SaveAs(pname);
    cout<<" open "<<pname<<endl;
  } // variables
  
  cout<<endl<<"Phoebe has "<<pevents.size()<<" unique events, Yipeng "<<yevents.size()
      <<". In common "<<ncommon<<endl;
  time(&endtime);
  cout<<endl<<"Finding unique events took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

template<class BabyType>
vector<sevent> getUniqueEvents(BabyType &baby, NamedFunc &cut, vector<sevent> &repeated){
  vector<sevent> events;

  int run=0;
  long nentries = baby.GetEntries(), event=0;
  //nentries = 25000;
  bool unique;
  for(long entry = 0; entry < nentries; ++entry){
    baby.GetEntry(entry);
    if(!cut.GetScalar(baby)) continue;
    
    run = baby.runNumber();
    event = baby.eventNumber();
    unique = true;
    for(unsigned ind=0; ind<events.size(); ind++){
      if(run == events[ind].run){
        if(event == events[ind].event){
          unique = false;
          break;
        }
      }
    } // for loop
    sevent newevent;
    newevent.run = run; newevent.event = event; newevent.entry = entry;
    if(unique) events.push_back(newevent);
    else repeated.push_back(newevent);
  } // Loop over entries

  return events;
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
