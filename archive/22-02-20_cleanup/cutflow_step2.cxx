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

  Palette colors("txt/colors.txt", "default");

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// Defining processes //////////////////////////////////////////

  string basecuts = "(b_l0_tis || d0_l0had_tos) && mu_l0_tis && mu_pidb > 0 && dst_id_prod > 0 && id_prod > 0";
  string repofolder = "../lhcb-ntuples-gen/gen/";
  vector<shared_ptr<Process> > procs;
  // procs.push_back(Process::MakeShared<Baby_phoebe_step2>("Phoebe", Process::Type::background, colors("run2"),
  //                                               set<string>({repofolder+"run1-Dst-step2/Dst--20_07_02--mix--data--2011--md--phoebe-step2.root"}),
  //                                                 "flag_2011 && polarity < 0 && spi_trk_type == 3"));
  procs.push_back(Process::MakeShared<Baby_phoebe_step1>("Phoebe", Process::Type::background, colors("run2"),
                                                set<string>({repofolder+"run1-Dst-step2/Dst--20_09_16--std--data--2011--md--phoebe-step2.root"}), "1"));
  procs.push_back(Process::MakeShared<Baby_step2>("Yipeng", Process::Type::background, colors("run1"),
                                                set<string>({repofolder+"run1-Dst_D0-step2/Dst_D0--20_10_12--std--data--2011--md--step2.root"}),basecuts));

  ///////// Automatically appending cutflow cuts
  vector<TString> scuts = {"1",
      "spi_gh_prob < 0.25", 
      "dst_endvtx_chi2/dst_endvtx_ndof < 10", 
      "(dst_m - d0_m - 145.454) < 2 && (dst_m - d0_m - 145.454) > -2", 
      "mu_is_mu", 
      "(mu_p > 3 && mu_p < 100)", 
      "(mu_eta > 1.7 && mu_eta < 5)", 
      "mu_pid_mu > 2", 
      "mu_pid_e < 1", 
      "mu_ip_chi2 > 45",   
      "mu_gh_prob < 0.5", 
      "b_discard_mu_chi2 <= 6",  
      "b_endvtx_chi2 < 24",  
      "b_endvtx_chi2/b_endvtx_ndof < 6", 
      "b_fd_trans < 7", 
      "b_dira > 0.9995", 
      "b_m < 5280 ", 
      "iso_bdt < 0.15"
  };
  vector<NamedFunc> cuts = {"1",
      "spi_gh_prob < 0.25", 
      "dst_endvtx_chi2/dst_endvtx_ndof < 10", 
      "(dst_m - d0_m - 145.454) < 2 && (dst_m - d0_m - 145.454) > -2", 
      "mu_is_mu", 
      "(mu_p > 3 && mu_p < 100)", 
      "(mu_eta > 1.7 && mu_eta < 5)", 
      "mu_pid_mu > 2", 
      "mu_pid_e < 1", 
      "mu_ip_chi2 > 45",   
      "mu_gh_prob < 0.5", 
      "b_discard_mu_chi2 <= 6",  
      "b_endvtx_chi2 < 24",  
      "b_endvtx_chi2/b_endvtx_ndof < 6", 
      "b_fd_trans < 7", 
      "b_dira > 0.9995", 
      "b_m < 5280 ", 
      "iso_bdt < 0.15"
  };
  vector<string> rownames = {"2011 data MD",
      "$spi\\_gh\\_prob < 0.25$", 
      "$dst\\_endvtx\\_chi2/dst\\_endvtx\\_ndof < 10$", 
      "$(dst\\_m - d0\\_m - 145.454) < 2 \\&\\& (dst\\_m - d0\\_m - 145.454) > -2$", 
      "$mu\\_is\\_mu$", 
      "$(mu\\_p > 3 \\&\\& mu\\_p < 100)$", 
      "$(mu\\_eta > 1.7 \\&\\& mu\\_eta < 5)$", 
      "$mu\\_pid\\_mu > 2$", 
      "$mu\\_pid\\_e < 1$", 
      "$mu\\_ip\\_chi2 > 45$",   
      "$mu\\_gh\\_prob < 0.5$", 
      "$b0\\_discard\\_mu\\_chi2 <= 6$",  
      "$b0\\_endvtx\\_chi2 < 24$",  
      "$b0\\_endvtx\\_chi2/b0\\_endvtx\\_ndof < 6$", 
      "$b0\\_fd\\_trans < 7$", 
      "$b0\\_dira > 0.9995$", 
      "$b0\\_m < 5280 $", 
      "$iso\\_bdt < 0.15$"
  };
  
  // vector<NamedFunc> cuts = {"1",
  //                           "k_pt > 0.5", "pi_pt > 0.5", "k_pt+pi_pt > 1.4",
  //                           "(k_hlt1ta_tos || pi_hlt1ta_tos)",
  //                           "k_p > 2", "pi_p > 2", "k_ip_chi2 > 45", "pi_ip_chi2 > 45",
  //                           "k_pid_k > 4", "pi_pid_k < 2", "!mu_veto",
  //                           "k_gh_prob < 0.5", "pi_gh_prob < 0.5",
  //                           "d0_pt > 2", 
  //                           "d0_hlt2charmhad_tos", 
  //                           "d0_endvtx_chi2/d0_endvtx_ndof < 4", 
  //                           "d0_ip_chi2 > 9", 
  //                           "d0_dira > 0.9998", 
  //                           "d0_fd_chi2 > 250", 
  //                           "(d0_m - 1865.49) < 23.4 && (-d0_m + 1865.49) < 23.4"
  // };
  
  // vector<string> rownames = {"2011 data MD",
  //                            "$k-pt > 0.5$", "$pi-pt > 0.5$", "$k-pt+pi-pt > 1.4$",
  //                            "$(k-hlt1ta-tos || pi-hlt1ta-tos)$",
  //                            "$k-p > 2$", "$pi-p > 2$", "$k-ip-chi2 > 45$", "$pi-ip-chi2 > 45$",
  //                            "$k-pid-k > 4$", "$pi-pid-k < 2$", "$!mu-veto$",
  //                            "$k-gh-prob < 0.5$", "$pi-gh-prob < 0.5$",
  //                            "$d0-pt > 2$", 
  //                            "$d0-hlt2charmhad-tos$", 
  //                            "$d0-endvtx-chi2/d0-endvtx-ndof < 4$", 
  //                            "$d0-ip-chi2 > 9$", 
  //                            "$d0-dira > 0.9998$", 
  //                            "$d0-fd-chi2 > 250$", 
  //                            "$(d0-m - 1865.49) < 23.4\\&\\& (-d0-m + 1865.49) < 23.4$"
  // };

  vector<TableRow> table_rows;
  NamedFunc fullcut = "1";
  TString scut ="";
  for(size_t ind = 0; ind < cuts.size(); ind++) {
    string title = (ind==0 ? rownames[ind] : "+ " + rownames[ind]);
    int lines = (ind==0 ? 1 : 0);
    fullcut = fullcut && cuts[ind];
    
    title = rownames[ind];
    fullcut = cuts[ind];
    scut += " && "; scut += scuts[ind];
    
    table_rows.push_back(TableRow(title,fullcut, 0,lines, "1"));
  }
  cout<<scut<<endl;
  PlotMaker pm;
  pm.Push<Table>("cutflow", table_rows,procs,0).TotColumn("Ratio"); // Pushing table
  
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
