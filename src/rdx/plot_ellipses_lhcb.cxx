#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TCanvas.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TEllipse.h"

#include "rdx/results.hpp"
#include "core/utilities.hpp"

using namespace std;

void getEllipse(float s1, float s2, float corr, float &r1, float &r2, float &angle){
  float sum = pow(s1,2)+pow(s2,2);
  float diff = pow(s1,2)-pow(s2,2);
  float root = pow(s1,4) - 2*pow(s1,2)*pow(s2,2) + 4*pow(s1,2)*pow(corr,2)*pow(s2,2) + pow(s2,4);
  if(root<0) cout<<"Bad root "<<root<<endl;
  root = sqrt(root);

  if(corr==0) angle = 90;
  else angle = -180/3.141593*atan((diff-root)/(2*corr*s1*s2)); // Eigenvector for r1
  r1 = sqrt((sum+root)/2.); // Largest eigenvalue of the cov. matrix
  r2 = sqrt((sum-root)/2.); // Smallest eigenvalue of the cov. matrix
}

int main(){

  bool isWide = true;

  //////// Results
  // int cBabar = kGreen, cLHCb = kCyan+1, cBelleHT = kBlue, cBelleST = kOrange;
  // int cAverage = kRed+1, cSM = kGreen+2;

  // Numbers from HFLAV Spring 2019
  // https://hflav-eos.web.cern.ch/hflav-eos/semi/spring19/html/RDsDsstar/RDRDs.html
  Results RD_SM ("SM pred.", 0.299, {0.003}, {0}, kBlue-4);
  Results RDs_SM("SM pred.", 0.258, {0.005}, {0}, kViolet+2);

  Results RD_HFLAV ("HFLAV aver.", 0.340, {0.027}, {0.013}, kRed, -0.38);
  Results RDs_HFLAV("HFLAV aver.", 0.295, {0.011}, {0.008}, kRed, -0.38);

  Results RD_BABAR ("BABAR (HT)", 0.440, {0.058}, {0.042}, kGreen+2, -0.31);
  Results RDs_BABAR("BABAR (HT)", 0.332, {0.024}, {0.018}, kGreen+2, -0.31);

  Results RD_BelleHT ("Belle (HT)", 0.375, {0.064}, {0.026}, kOrange-7, -0.50);
  Results RDs_BelleHT("Belle (HT)", 0.293, {0.038}, {0.015}, kBlue, -0.50);

  Results RD_BelleST("Belle (ST)", 0.307, {0.037}, {0.016}, kOrange, -0.51);
  Results RDs_BelleST("Belle (ST)", 0.283, {0.018}, {0.014}, kOrange, -0.51);

  Results RDs_Bellepi("Belle (#pi/#rho)", 0.270, {0.035}, {0.027}, kOrange+1);

  float syst_rd  = sqrt(pow(4.45,2) + pow(4.47,2) + pow(2.5,2) + pow(2,2))*0.01, stat_rd = 5.82*0.01; //Add. syst, MC stats, DD, Mult. syst
  float syst_rds = sqrt(pow(1.86,2) + pow(1.69,2) + pow(0.9,2) + pow(0.9,2))*0.01, stat_rds = 1.8*0.01; //Add. syst, MC stats, DD, Mult. syst
  Results RD_LHCb("LHCb (#mu)", 0.34, {stat_rd}, {syst_rd}, kRed-7, -0.38);
  Results RDs_LHCb("LHCb (#mu)", 0.336, {stat_rds}, {syst_rds}, kCyan-4, -0.38);
  
  Results RDs_LHCb1("LHCb (#mu)", 0.336, {0.027}, {0.030}, kCyan-2);
  
  Results RD_LHCb2("LHCb2 (#mu)", 0.34, {0.09*0.34}, {0.09*0.34}, kRed+2, -0.38);
  Results RDs_LHCb2("LHCb2 (#mu)", 0.336, {0.03*0.336}, {0.03*0.336}, kCyan-2, -0.38);

  Results RDs_LHCb3pi("LHCb (#pi#pi#pi)", 0.280, {0.018}, {0.029}, kMagenta+2);

  if(isWide){
    RD_BABAR.name = "BaBar, PRL #font[62]{109}, 101802 (2012)";
    RD_BelleHT.name = "Belle, PRD #font[62]{92}, 072014 (2015)";
    //RDs_BelleST.name = "Belle, PRD #font[62]{94}, 072007 (2016)";
    RD_BelleST.name = "Belle, PRL #font[62]{124}, 161803 (2020)";
    RDs_BelleST.name = "Belle, PRL #font[62]{124}, 161803 (2020)";
    RDs_Bellepi.name = "Belle, PRL #font[62]{118}, 211801 (2017)";
    RD_LHCb.name = "LHCb muonic Run 1 projection";
    RD_LHCb2.name = "LHCb muonic Run 2 projection";
    RDs_LHCb1.name = "LHCb, PRL #font[62]{115}, 111803 (2015)";
    RDs_LHCb3pi.name = "LHCb, PRL #font[62]{120}, 171802 (2018)";
    RD_HFLAV.name = "HFLAV average Spring 2019";
    RD_SM.name = "SM predictions";
  }

  vector<vector<Results> > results;
  results.push_back({RD_BABAR,  RDs_BABAR});
  results.push_back({RDs_LHCb1});
  results.push_back({RD_BelleHT,  RDs_BelleHT});
  results.push_back({RDs_Bellepi});
  results.push_back({RDs_LHCb3pi});
  results.push_back({RD_BelleST, RDs_BelleST});
  //results.push_back({RD_HFLAV,  RDs_HFLAV});
  results.push_back({RD_SM,    RDs_SM});
  results.push_back({RD_LHCb, RDs_LHCb});
  results.push_back({RD_LHCb2, RDs_LHCb2});


  //////// Setting style
  float axisTextSize = 0.065;
  float lMargin = 0.15, rMargin = 0.05;
  float bMargin = 0.17, tMargin = 0.05;
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(bMargin);
  gStyle->SetPadTopMargin(tMargin);
  gStyle->SetPadLeftMargin(lMargin);
  gStyle->SetPadRightMargin(rMargin);
  gStyle->SetNdivisions(907, "xy");   
  gStyle->SetTitleOffset(1.17,"x");
  gStyle->SetTitleOffset(1.17,"y");     
  gStyle->SetTitleFont(132,"xyz");          // Set the all 2 axes title font
  gStyle->SetLabelFont(132,"xyz");          // Set the all 2 axes label font
  gStyle->SetTextFont(132);                // Set global text font
  gStyle->SetPadTickX(1);             // Ticks at the top
  gStyle->SetPadTickY(1);             // Ticks at the right

  if(isWide) {
    gStyle->SetTitleOffset(0.8,"y");
    axisTextSize *= 1.1;
  }
  gStyle->SetTitleSize(axisTextSize,"xy");     // Set the 2 axes title size
  gStyle->SetLabelSize(axisTextSize/1.1,"xy");     // Set the 2 axes label size


  //////// Creating canvas and base histogram
  int cW = 500, cH = 480;
  if(isWide) cW = 800;
  TCanvas can("can","", cW, cH);
  //TCanvas can("can","", 500, 410);
  float minX=0.20, maxX=0.54, minY=0.21, maxY=0.49;
  if(isWide){
    minX = 0.19; maxX = 0.54; minY = 0.207; maxY = 0.52;
  }
  TH1D histo("histo", "", 10, minX, maxX);
  histo.SetMinimum(minY);
  histo.SetMaximum(maxY);
  histo.GetXaxis()->CenterTitle(true);
  histo.GetYaxis()->CenterTitle(true);
  histo.SetXTitle("R(D)");
  histo.SetYTitle("R(D#lower[-.1]{*})");
  histo.Draw("axis");
  

  double Nsingle = 0; // Number of 1D measurements

  //////// Drawing results
  float alpha = 0.04;
  int lWidth = 2;
  TEllipse ellipse; 
  TBox box; box.SetLineWidth(0);
  TLine line;  line.SetLineWidth(lWidth+1);
  float maxR, minR, angle;
  for(auto &result : results){
    if(result.size()>1){
      float rd = result[0].value, rds = result[1].value;
      float erd = result[0].errUp(), erds = result[1].errUp();
      getEllipse(erd, erds, result[0].correl, maxR, minR, angle);
      if(result[0].name.Contains("aver") || result[0].name.Contains("SM")) 
	ellipse.SetFillColorAlpha(result[0].color, alpha+0.32);
      else if(result[0].name.Contains("muonic")) 
	ellipse.SetFillColorAlpha(result[0].color, .8);
      else ellipse.SetFillColorAlpha(result[0].color, alpha);
      ellipse.SetLineColor(result[0].color); ellipse.SetLineWidth(lWidth);
      ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360, angle, "c");
      if(result[0].name.Contains("SM")) ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360*10, angle, "c");
      //cout<<"RDs = "<<rds<<" +- "<<erds/rds*100<<endl;
      //cout<<"RD  = "<<rd<<" +- "<<erd/rd*100<<endl<<endl;

    } else {
      float dx = 0.003;
      float rd = 0.23+0.01*Nsingle;
      float rds = result[0].value;
      float rdsUp = result[0].value + result[0].errUp();
      float rdsDown = result[0].value - result[0].errDown();
      //cout<<"RDs = "<<rds<<" +- "<<(rdsUp-rds)/rds*100<<endl<<endl;

      ellipse.SetLineColor(result[0].color); ellipse.SetLineWidth(lWidth+2);
      ellipse.DrawEllipse(rd, rds, 0.0023, 0.0023, 0, 360, 0, "c");
      line.SetLineWidth(lWidth);
      line.SetLineColor(result[0].color);
      line.DrawLine(rd, rdsUp, rd, rdsDown);
      line.DrawLine(rd-dx, rdsUp, rd+dx, rdsUp);
      line.DrawLine(rd-dx, rdsDown, rd+dx, rdsDown);
      // // Drawing bands for RD* results only
      // box.SetFillColorAlpha(result[0].color, alpha*2);
      // box.DrawBox(minX, rdsDown, maxX, rdsUp);
      // line.SetLineWidth(1);
      // line.DrawLine(minX, rdsUp, maxX, rdsUp);
      // line.DrawLine(minX, rdsDown, maxX, rdsDown);
      Nsingle++;
    }
  } // Loop over results

  // if(isWide){
  //   Results RD_HFLAV4 ("noplot", RD_HFLAV.value, {static_cast<float>(3.56 * RD_HFLAV.errUp())}, {0}, kRed, RD_HFLAV.correl);
  //   Results RDs_HFLAV4("noplot", RDs_HFLAV.value, {static_cast<float>(3.56* RDs_HFLAV.errUp())}, {0}, kRed, RD_HFLAV.correl);
  //   vector<Results>  result({RD_HFLAV4,  RDs_HFLAV4});

  //   float rd = result[0].value, rds = result[1].value;
  //   float erd = result[0].errUp(), erds = result[1].errUp();
  //   getEllipse(erd, erds, result[0].correl, maxR, minR, angle);
  //   ellipse.SetFillStyle(0);
  //   ellipse.SetLineWidth(1);
  //   ellipse.SetLineStyle(2);
  //   ellipse.SetLineColor(result[0].color); 
  //   //ellipse.DrawEllipse(rd, rds, maxR, minR, 0, 360, angle, "c");
  //   TLatex label;  //label.SetNDC(kTRUE);
  //   label.SetTextAlign(23); label.SetTextSize(0.06);
  //   label.SetTextColor(kRed+1);
  //   //label.DrawLatex(0.3, 0.365, "#font[62]{3.1#sigma}");

  // }

  //////// Legend
  int Ncols = 3; if(isWide) Ncols = 2;
  double legX(lMargin+0.05), legY = 1-tMargin-0.04, legSingle = 0.06, leg2X = legX;
  double legW = 0.25*Ncols, legH = legSingle*(8+Ncols-1)/Ncols, legFontSize = 0.045;
  if(results.size()<=2) legW = 0.217*Ncols;
  if(results.size()<=3) {
    legH = legSingle*(8+Ncols-1)/Ncols/3;   
  } else if(results.size()<=6) legH = legSingle*(8+Ncols-1)/Ncols/3*2;
  if(isWide) {
    Ncols = 1;
    legX -= 0.01;
    leg2X = legX+0.37;
    legW = 0.1;
    legH = legSingle * results.size()/1.8;
    legFontSize = 0.045;
  } 
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(legFontSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetNColumns(Ncols);
  TLegend leg2(leg);
  leg2.SetX1(leg2X); leg2.SetX2(leg2X+legW); 
  vector<TH1D*> hLeg;
  if(!isWide){
    for(size_t ind=0; ind<results.size(); ind++){
      TString hname = "hLeg"+to_string(ind);
      hLeg.push_back(new TH1D(hname,"",10,0,10));
      hLeg.back()->SetLineWidth(lWidth);
      hLeg.back()->SetLineColor(results[ind][0].color);
      hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
      leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
    }
    leg.Draw();
  }else {
    for(size_t ind=0; ind<results.size(); ind++){
      TString hname = "hLeg"+to_string(ind);
      hLeg.push_back(new TH1D(hname,"",10,0,10));
      hLeg.back()->SetLineWidth(lWidth);
      hLeg.back()->SetLineColor(results[ind][0].color);
      hLeg.back()->SetFillColorAlpha(results[ind][0].color, alpha);
      if(ind%2==0) leg.AddEntry((hLeg.back()), results[ind][0].name, "lf");
      else leg2.AddEntry((hLeg.back()), results[ind][0].name, "lf");
    }
    leg.Draw(); leg2.Draw();
  }

  histo.Draw("same axis");


  TString plotName = "plots/rdx_ellipses.pdf";
  can.SaveAs(plotName);

  return 0;
}

