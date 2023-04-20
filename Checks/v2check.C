/*#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"
*/
using namespace std;


  //initialize some common variables
  double pi = 3.14159265359;
  double lbound = 10.0;
  double rbound = 190.0;
  double eventsperjet = 49.5;
  bool cent = false;
  bool semi = false;
  bool fine = false;
  bool coar = false;
  bool r02 = false;
  bool r04 = false;
  int centcolor = kViolet+5;
  int semicolor = kCyan-4;

  //declare functions
  void setcolor(TH1F* h, int kcolor);
  void setcoloralpha(TH1F* h, int kcolor, double trans);


int ci = 1756; // color index
TColor *bluecolor = new TColor(ci, 0, 125, 197);
vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};
TH1D *spacer = new TH1D("spacer", "spacer", 10, 20, 120);
TH1D *spacer2 = (TH1D*)spacer->Clone("spacer2");

//=======================================================================
//=======================================================================
//========================= main function ===============================
//=======================================================================
//=======================================================================

void v2check (int iteration, string filename, int planestatus = 0) {

  const char* r02file = filename.c_str();
  const char* r04file = filename.c_str();
  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

      //set custom colors
      int kCaitieBlue, kCaitieDarkBlue;
        kCaitieBlue = TColor::GetColor("#6f73c3");
        kCaitieDarkBlue = TColor::GetColor("#1c2280");
      int kCaitieGreen, kCaitieDarkGreen;
        kCaitieGreen = TColor::GetColor("#CCFFE4");
        kCaitieDarkGreen = TColor::GetColor("#066C67");
      int kCaitiePink, kCaitieLightPink;
        kCaitiePink = TColor::GetColor("#CA267A");
        kCaitieLightPink = TColor::GetColor("#F42ECF");
      int kSunBlue, kSunPurple, kSunPink, kSunOrange, kSunGreen, kDarkSunGreen;
        kSunBlue = TColor::GetColor("#006a80"); //blue
        kSunPurple = TColor::GetColor("#2100a3");  //purple
        kSunPink = TColor::GetColor("#b8239c");//pink
        kSunOrange = TColor::GetColor("#f8842d"); //e86400
        kSunGreen = TColor::GetColor("#9FD064"); //green
        kDarkSunGreen = TColor::GetColor("#67AD2B");

   //==================================================
   //======== 1. get file tag for legend 
   //==================================================
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = "_";
   for (int i = 0; i < filename.length(); i++)    {
       const char *a = &r02file[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
       else ss << r02file[i];
       }

  
   //========================================================
   //====== 2. establish centrality from unfolding filename
   //========================================================  
    string tags;
       vector<string> words(0);
       while (ss >> tags)   words.push_back(tags);
       stringstream s2;
    const char* ccent = "010";
    const char* scent3050 = "3050";
    const char* scent3035 = "3035";
    const char* scent3540 = "3540";
    const char* scent4045 = "4045";
    const char* scent4550 = "4550";
    const char* fin = "fine";
    const char* coa = "coarse";
    const char* r020 = "R02";
    const char* r040 = "R04";
      for (int i = 2; i < words.size() - 1; i++) {
        s2 << words[i] << " ";
        const char *newword = words[i].c_str();
        if (strncmp(newword, ccent, 3) == 0)   cent = true;
        if (strncmp(newword, scent3050, 3) == 0)   semi = true;
        if (strncmp(newword, scent3035, 3) == 0 || strncmp(newword, scent3540, 3) == 0)   semi = true;
        if (strncmp(newword, scent4045, 3) == 0 || strncmp(newword, scent4550, 3) == 0)   semi = true;
        if (strncmp(newword, fin, 3) == 0)     fine = true;
        if (strncmp(newword, coa, 3) == 0)     coar = true;
        if (strncmp(newword, r020, 3) == 0)     r02 = true;
        if (strncmp(newword, r040, 3) == 0)     r04 = true;
      } 
    double centl, centr;
     if (cent == true)  { 
       centl = 0.0;
       centr = 10.0;
      }
     if (semi == true)  {
       centl = 30.0;
       centr = 50.0;
      } 
    if (fine == false && coar == false)   fine = true;
    int R;
    if (r02 == true)  R = 2;
    if (r04 == true)  R = 4;





  //======================================================
  //========= 4. retrieve unfolded spectra from file
  //=======================================================
    TFile *lefilel = new TFile(r02file);
    TFile *lefileh = new TFile(r04file);

    //Get 2D unfolded hists for low and high q2
    vector<TH2F*> unfoldhistsl(0); 
    vector<TH2F*> unfoldhistsh(0);
    for (int i = 1; i <= 15; i++)  {
        stringstream unfoldName;
        unfoldName << "Bayesian_Unfolded_" << i;
        unfoldhistsl.push_back((TH2F*)lefilel->Get(unfoldName.str().c_str()));
        unfoldhistsh.push_back((TH2F*)lefileh->Get(unfoldName.str().c_str()));
    }

    //Convert 2D hists into 1D projections
    vector<TH1F*> unfoldHistsLIP(0); 
    vector<TH1F*> unfoldHistsLOP(0); 
    vector<TH1F*> unfoldHistsLow(0); 
    vector<TH1F*> unfoldHistsHIP(0);
    vector<TH1F*> unfoldHistsHOP(0); 
    vector<TH1F*> unfoldHistsHigh(0); 
    for (int i = 0; i < 15; i++) {
       stringstream sslip, sslop, sslow, sship, sshop, sshigh;
       sslip << "tempHistLIP" << i;
       sslop << "tempHistLOP" << i;
       sslow << "tempHistLow" << i;
       sship << "tempHistHIP" << i;
       sshop << "tempHistHOP" << i;
       sshigh << "tempHistHigh" << i;
       TH2F* tempHistLIP = (TH2F*)unfoldhistsl[i]->Clone(sslip.str().c_str());
       TH2F* tempHistLOP = (TH2F*)unfoldhistsl[i]->Clone(sslop.str().c_str());
       TH2F* tempHistLow = (TH2F*)unfoldhistsl[i]->Clone(sslow.str().c_str());
       TH2F* tempHistHIP = (TH2F*)unfoldhistsh[i]->Clone(sship.str().c_str());
       TH2F* tempHistHOP = (TH2F*)unfoldhistsh[i]->Clone(sshop.str().c_str());
       TH2F* tempHistHigh = (TH2F*)unfoldhistsh[i]->Clone(sshigh.str().c_str());
       tempHistLOP->GetYaxis()->SetRangeUser(0.0, sqrt(2)/2.0);
       tempHistLIP->GetYaxis()->SetRangeUser(sqrt(2)/2.0, 1.0);
       tempHistLow->GetYaxis()->SetRangeUser(0.0, 1.0);
       tempHistHOP->GetYaxis()->SetRangeUser(0.0, sqrt(2)/2.0);
       tempHistHIP->GetYaxis()->SetRangeUser(sqrt(2)/2.0, 1.0);
       tempHistHigh->GetYaxis()->SetRangeUser(0.0, 1.0);
       unfoldHistsLIP.push_back((TH1F*)tempHistLIP->ProjectionX());
       unfoldHistsLOP.push_back((TH1F*)tempHistLOP->ProjectionX());
       unfoldHistsLow.push_back((TH1F*)tempHistLow->ProjectionX());
       unfoldHistsHIP.push_back((TH1F*)tempHistHIP->ProjectionX());
       unfoldHistsHOP.push_back((TH1F*)tempHistHOP->ProjectionX());
       unfoldHistsHigh.push_back((TH1F*)tempHistHigh->ProjectionX());
    }


  //=======================================================
  //========= 5. retrieve efficiencies from file ==========
  //=======================================================
    
        //Get Kinematic Efficiency
        //File Low
        TH1F *kinpre1 = (TH1F*)lefilel->Get("KinPre");
        TH1F *kinpos1 = (TH1F*)lefilel->Get("KinPos");
             TH1F *effic1 = (TH1F*)kinpos1->Clone("effic1");
             effic1->Divide(kinpre1);
        //File High
        TH1F *kinpre2 = (TH1F*)lefileh->Get("KinPre");
        TH1F *kinpos2 = (TH1F*)lefileh->Get("KinPos");
             TH1F *effic2 = (TH1F*)kinpos2->Clone("effic2");
             effic2->Divide(kinpre2);

       //Get Reconstruction Efficiency
       TFile *recfile = new TFile("ChargedJetRecEfficiencies.root");
       stringstream recEfficName;
       recEfficName << "RecEff_R0" << R << "0_5GeV";
       TH1D  *receffic = (TH1D*)recfile->Get(recEfficName.str().c_str());
       

  //===================================================
  //============= Rebin PbPb for propper edges ========
  //===================================================
      //vector<double> pTedgerebin = {20, 30, 40, 50, 60, 70, 85, 100};
      vector<double> pTedgerebin;
          //if (fine == true)     pTedgerebin = {30, 35, 40, 50, 60, 70, 85, 100, 120, 140};
          //if (coar == true)     pTedgerebin  = {30, 35, 40, 50, 65, 80, 100, 120, 140};
            if (R == 2) pTedgerebin = {35, 40, 50, 60, 80, 100, 120};
            if (R == 4) pTedgerebin = {40, 50, 60, 80, 100, 120};
      TH1F* efficiency1 = (TH1F*)effic1->Rebin(pTedgerebin.size()-1, "efficiency1", pTedgerebin.data());
      TH1F* efficiency2 = (TH1F*)effic2->Rebin(pTedgerebin.size()-1, "efficiency2", pTedgerebin.data());
      TH1F* recefficiency = (TH1F*)receffic->Rebin(pTedgerebin.size()-1, "recefficiency", pTedgerebin.data());
      TH1F* yieldLow = (TH1F*)unfoldHistsLow[iteration]->Rebin(pTedgerebin.size()-1, "yieldLow", pTedgerebin.data());
            yieldLow->Sumw2();
      TH1F* yieldLowInPlane = (TH1F*)unfoldHistsLIP[iteration]->Rebin(pTedgerebin.size()-1, "yieldLowInPlane", pTedgerebin.data());
            yieldLowInPlane->Sumw2();
      TH1F* yieldLowOutPlane = (TH1F*)unfoldHistsLOP[iteration]->Rebin(pTedgerebin.size()-1, "yieldLowOutPlane", pTedgerebin.data());
            yieldLowOutPlane->Sumw2();
      TH1F* yieldHigh = (TH1F*)unfoldHistsHigh[iteration]->Rebin(pTedgerebin.size()-1, "yieldHigh", pTedgerebin.data());
            yieldHigh->Sumw2();
      TH1F* yieldHighInPlane = (TH1F*)unfoldHistsHIP[iteration]->Rebin(pTedgerebin.size()-1, "yieldHighInPlane", pTedgerebin.data());
            yieldHighInPlane->Sumw2();
      TH1F* yieldHighOutPlane = (TH1F*)unfoldHistsHOP[iteration]->Rebin(pTedgerebin.size()-1, "yieldHighOutPlane", pTedgerebin.data());
            yieldHighOutPlane->Sumw2();
  


//==============================================================
  //========= 6. apply scaling for bin width + N_events ========
  //============================================================

        double TPCfid = 2*(0.9-(0.1*R));

        //Unfolded Spectra 
        for (int i = 0; i < unfoldHistsLIP.size(); i++)    {
        unfoldHistsLIP[i]->Scale(1.0, "width");
        //unfoldHistsLIP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsLOP[i]->Scale(1.0, "width");
        //unfoldHistsLOP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsLow[i]->Scale(1.0, "width");
        //unfoldHistsLow[i]->Scale(1.0/((double)N_event*TPCfid));
        }


        for (int i = 0; i < unfoldHistsHIP.size(); i++)    {
        unfoldHistsHIP[i]->Scale(1.0, "width");
        //unfoldHistsHIP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsHOP[i]->Scale(1.0, "width");
        //unfoldHistsHOP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsHigh[i]->Scale(1.0, "width");
        //unfoldHistsHigh[i]->Scale(1.0/((double)N_event*TPCfid));
        }

      //grab fractional stat errors before applying efficiencies
      vector<double> statYL (0);
      vector<double> statYLIP (0);
      vector<double> statYLOP (0);
      vector<double> statYH (0);
      vector<double> statYHIP (0);
      vector<double> statYHOP (0);
      for (int i = 1; i <= yieldLow->GetNbinsX(); i++) {
          statYL.push_back(yieldLow->GetBinError(i)/yieldLow->GetBinContent(i));
          statYLIP.push_back(yieldLowInPlane->GetBinError(i)/yieldLowInPlane->GetBinContent(i));
          statYLOP.push_back(yieldLowOutPlane->GetBinError(i)/yieldLowOutPlane->GetBinContent(i));
          statYH.push_back(yieldHigh->GetBinError(i)/yieldHigh->GetBinContent(i));
          statYHIP.push_back(yieldHighInPlane->GetBinError(i)/yieldHighInPlane->GetBinContent(i));
          statYHOP.push_back(yieldHighOutPlane->GetBinError(i)/yieldHighOutPlane->GetBinContent(i));
     }

       
       //kinematic efficiency + reconstruction efficiency
       yieldLow->Divide(efficiency1);
       yieldLow->Divide(recefficiency);
       yieldLowInPlane->Divide(efficiency1);
       yieldLowInPlane->Divide(recefficiency);
       yieldLowOutPlane->Divide(efficiency1);
       yieldLowOutPlane->Divide(recefficiency);
       yieldHigh->Divide(efficiency2);
       yieldHigh->Divide(recefficiency);
       yieldHighInPlane->Divide(efficiency2);
       yieldHighInPlane->Divide(recefficiency);
       yieldHighOutPlane->Divide(efficiency2);
       yieldHighOutPlane->Divide(recefficiency);

 
  //=====================================================
  //============= Calculate the Yield Ratio =============
  //=====================================================
 
      //ratio high q2/lowq2
      TH1F* yieldratio =(TH1F*)yieldHigh->Clone("yieldratio");
      yieldratio->Divide(yieldLow);
        //out/in (high)
        TH1F* ratioOIHigh = (TH1F*)yieldHighOutPlane->Clone("ratioOIHigh");
        ratioOIHigh->Divide(yieldHighInPlane);
        //out/in (low)
        TH1F* ratioOILow = (TH1F*)yieldLowOutPlane->Clone("ratioOILow");
        ratioOILow->Divide(yieldLowInPlane);
      double lowedge;
      if (R == 2)  lowedge = 35.0;
      if (R == 4)  lowedge = 40.0;
      int startbin = yieldratio->FindBin(lowedge+0.01);
      
      yieldratio->SetLineColor(kCaitieDarkGreen);
      yieldratio->SetMarkerColor(kCaitieDarkGreen);
      yieldratio->SetMarkerStyle(20);
      yieldratio->SetMarkerSize(0.8);
      yieldratio->SetLineWidth(2);

      setcolor(ratioOILow, kSunPink);
      setcolor(ratioOIHigh, kCaitieDarkBlue);

      setcolor(yieldHighInPlane, kSunBlue);
      setcolor(yieldHighOutPlane, kSunOrange);
        yieldHighInPlane->SetMarkerStyle(21);
        yieldHighOutPlane->SetMarkerStyle(21);
      setcolor(yieldLowInPlane, kSunPink);
      setcolor(yieldLowOutPlane, kSunPurple);

      yieldLowInPlane->SetLineColorAlpha(kWhite, 0.0);
      yieldLowInPlane->SetMarkerColorAlpha(kWhite, 0.0);

  //===============================================================
  //============= Calculate Statistical Uncertainties =============
  //===============================================================
    
      //reset to original values after applying efficiencies
      //note: indices do not match on vectors and hists because vectors 0 indexed, but hists start with bin 1
      for (int i = 1; i <= yieldLow->GetNbinsX(); i++)  {
          //spectra (fractional unc * efficiency corrected values)
          yieldLow->SetBinError(i, statYL[i-1]*yieldLow->GetBinContent(i));
          yieldLowInPlane->SetBinError(i, statYLIP[i-1]*yieldLowInPlane->GetBinContent(i));
          yieldLowOutPlane->SetBinError(i, statYLOP[i-1]*yieldLowOutPlane->GetBinContent(i));
          yieldHigh->SetBinError(i, statYH[i-1]*yieldHigh->GetBinContent(i));
          yieldHighInPlane->SetBinError(i, statYHIP[i-1]*yieldHighInPlane->GetBinContent(i));
          yieldHighOutPlane->SetBinError(i, statYHOP[i-1]*yieldHighOutPlane->GetBinContent(i));
          //ratios (add fractional unc in quadrature, multiply by ratio value)
          yieldratio->SetBinError(i, sqrt(pow(statYL[i-1],2)+pow(statYH[i-1],2))*yieldratio->GetBinContent(i));
          ratioOILow->SetBinError(i, sqrt(pow(statYLIP[i-1],2)+pow(statYLOP[i-1],2))*ratioOILow->GetBinContent(i));
          ratioOIHigh->SetBinError(i, sqrt(pow(statYHIP[i-1],2)+pow(statYHOP[i-1],2))*ratioOIHigh->GetBinContent(i));
      }



  //===============================================================
  //============= Calculate Systematic Uncertainties ==============
  //===============================================================


      //Penguin
      //coordinates for TGraphErrors
      double bincenters3050[8] = {31., 36., 42., 52., 62., 73., 88., 108.5};
      double bincentersR02[6] = {37.5, 45.0, 55.0, 70.0, 90.0, 110.0};
      double bincentersR04[5] = {45.0, 55.0, 70.0, 90.0, 110.0};
      double bincenters3050_hi[8] = {31., 36., 42., 52., 62., 73., 88., 108.5};
      double bincenters3050_li[8] = {32., 37., 44., 54., 64., 76., 91., 109.5};
      double bincenters3050_lo[8] = {33., 38., 46., 56., 66., 79., 94., 110.5};
      double bincenters3050_ho[8] = {34., 39., 48., 58., 68., 82., 97., 111.5};

      //spectral values
      //out/in ratio
      double lowvalues3050[8] = {ratioOILow->GetBinContent(startbin), ratioOILow->GetBinContent(startbin+1),
                                ratioOILow->GetBinContent(startbin+2), ratioOILow->GetBinContent(startbin+3),
                                ratioOILow->GetBinContent(startbin+4), ratioOILow->GetBinContent(startbin+5),
                                ratioOILow->GetBinContent(startbin+6), ratioOILow->GetBinContent(startbin+7)};
      double highvalues3050[8] = {ratioOIHigh->GetBinContent(startbin), ratioOIHigh->GetBinContent(startbin+1), 
                                 ratioOIHigh->GetBinContent(startbin+2), ratioOIHigh->GetBinContent(startbin+3),
                                 ratioOIHigh->GetBinContent(startbin+4), ratioOIHigh->GetBinContent(startbin+5),
                                 ratioOIHigh->GetBinContent(startbin+6), ratioOIHigh->GetBinContent(startbin+7)};
      //high/low ratio
      double nomlo[8] = {yieldLow->GetBinContent(startbin), yieldLow->GetBinContent(startbin+1), yieldLow->GetBinContent(startbin+2), 
                        yieldLow->GetBinContent(startbin+3), yieldLow->GetBinContent(startbin+4), yieldLow->GetBinContent(startbin+5),
                        yieldLow->GetBinContent(startbin+6), yieldLow->GetBinContent(startbin+7)};
      double nomhi[8] = {yieldHigh->GetBinContent(startbin), yieldHigh->GetBinContent(startbin+1), yieldHigh->GetBinContent(startbin+2), 
                        yieldHigh->GetBinContent(startbin+3), yieldHigh->GetBinContent(startbin+4), yieldHigh->GetBinContent(startbin+5),
                        yieldHigh->GetBinContent(startbin+6), yieldHigh->GetBinContent(startbin+7)};
      double nomrat[8] = {yieldratio->GetBinContent(startbin), yieldratio->GetBinContent(startbin+1), yieldratio->GetBinContent(startbin+2),
                         yieldratio->GetBinContent(startbin+3), yieldratio->GetBinContent(startbin+4), yieldratio->GetBinContent(startbin+5),
                          yieldratio->GetBinContent(startbin+6), yieldratio->GetBinContent(startbin+8)};
      double nul3050[8] = {2.5, 2.5, 5.0, 5.0, 5.0, 7.5, 7.5, 10.0};
      double nul3050_crop[8] = {0.75, 0.75, 1.5, 1.5, 1.5, 2.25, 2.25, 3.};
      double nul3050_short[6] = {2.5, 5.0, 5.0, 5.0, 7.5, 7.5};
      double nulR02[6] = {2.5, 5.0, 5.0, 10.0, 10.0, 10.0};
      double nulR04[5] = {5.0, 5.0, 10.0, 10.0, 10.0};






  //======================================
  //========== Redmer's Data  ============
  //======================================
  //https://www.hepdata.net/record/ins1394678
  //Figure 4B
  double jetv2[7] = {0.0811, 0.0914, 0.0756, 0.0625, 0.0575, 0.0588, 0.0611};
  double shapesyst[7] = {0.0392, 0.0291, 0.0179, 0.0188, 0.0111, 0.0154, 0.0148};
  vector<double>  staterr = {0.00753, 0.0114, 0.0165, 0.0213, 0.0286, 0.0383, 0.0489};
  //vector<double>  jetv2 = {0.06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.04};
  double R2 = 0.75;
  //double R2 = 1.0;
  double k1 = (pi/4)*(1.0/R2);
  vector<double> k2(0);
  vector<double> k0(0);
  vector<double>  binning = {20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
  double redmerCent[7] = {25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0};
  double redmerNul[7] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  TH1D *hJetv2 = new TH1D("hJetv2", "hJetv2", binning.size()-1, binning.data());
  

  //reproduce paper figure
  double newR2 = 0.77;
  //double newR2 = 1.0;
  TH1D *paperFig = new TH1D("paperFig", "paperFig", binning.size()-1, binning.data());
  TH1D *myv2 = new TH1D("myv2", "myv2", pTedgerebin.size()-1, pTedgerebin.data());
  vector<double> vecMyv2(0);
  vector<double> myRat(0);
  for (int i = 0; i < binning.size()-1; i++)   {
      paperFig->Fill(binning[i]+0.01, jetv2[i]); 
      paperFig->SetBinError(i+1, staterr[i]);}
  for (int i = 0; i < pTedgerebin.size()-1; i++)  {
      double Nin = yieldLowInPlane->GetBinContent(yieldLowInPlane->FindBin(pTedgerebin[i]+0.01));
      double NinErr = yieldLowInPlane->GetBinError(yieldLowInPlane->FindBin(pTedgerebin[i]+0.01));
      double Nout = yieldLowOutPlane->GetBinContent(yieldLowOutPlane->FindBin(pTedgerebin[i]+0.01));
      double NoutErr = yieldLowOutPlane->GetBinError(yieldLowOutPlane->FindBin(pTedgerebin[i]+0.01));
             double numErr = sqrt((NinErr*NinErr)+(NoutErr*NoutErr))/(Nin-Nout);
             double denErr = sqrt((NinErr*NinErr)+(NoutErr*NoutErr))/(Nin+Nout);
      myRat.push_back(Nout/Nin);
      vecMyv2.push_back((pi/4.0)*(1.0/newR2)*(Nin-Nout)/(Nin+Nout));
      myv2->Fill(pTedgerebin[i]+0.01, vecMyv2[i]);
      myv2->SetBinError(i+1, myv2->GetBinContent(i+1)*sqrt((numErr*numErr)+(denErr*denErr)));}

      double v2caitie[6] = {vecMyv2[0], vecMyv2[1], vecMyv2[2], vecMyv2[3], vecMyv2[4], vecMyv2[5]};

  paperFig->SetLineColor(kBlack);
  paperFig->SetLineWidth(2);
  paperFig->SetMarkerColor(kBlack);
  paperFig->SetMarkerStyle(20);
  paperFig->SetMarkerSize(0.7);
 
  myv2->SetLineColor(kMagenta+2);
  myv2->SetLineWidth(2);
  myv2->SetMarkerColor(kMagenta+2);
  myv2->SetMarkerStyle(20);
  myv2->SetMarkerSize(0.7);

  //fill vector (then hist) with values for yield ratio
  for (int i = 0; i < binning.size()-1; i++) {
                    k2.push_back(k1/jetv2[i]);
                    k0.push_back((k2[i]-1.0)/(k2[i]+1.0));}
  for (int i = 0; i <= binning.size()-1; i++)   {
                    hJetv2->Fill(binning[i]+0.01, k0[i]);
                    hJetv2->SetBinError(i, 0.0001);}

  //fill vector with values for errors
  vector<double> l2(0), l0(0), h2(0), h0(0);
  vector<double> el0(0), el2(0), eh0(0), eh2(0), e0(0);
  for (int i = 0; i < binning.size()-1; i++)   { 
      l2.push_back(k1/(jetv2[i]-shapesyst[i]));
      h2.push_back(k1/(jetv2[i]+shapesyst[i]));
      el2.push_back(k1/(jetv2[i]-staterr[i]));
      eh2.push_back(k1/(jetv2[i]+staterr[i]));
      //if ((jetv2[i] + shapesyst[i]) > 1.047) h2.push_back(k1/1.04);
      l0.push_back((l2[i]-1.0)/(l2[i]+1.0));
      h0.push_back((h2[i]-1.0)/(h2[i]+1.0));
      el0.push_back((el2[i]-1.0)/(el2[i]+1.0));
      eh0.push_back((eh2[i]-1.0)/(eh2[i]+1.0));
      //cout << "el0: " << el0[i] << ", eh0: " << eh0[i] <<"\n";
        // cout << "k0: " << k0[i] <<"\n";
         double el = el0[i]-k0[i];
         double eh = k0[i]-eh0[i];
         //cout << "	el: " << el << ", eh: " << eh <<"\n";
      e0.push_back((el+eh)/2.0);
      hJetv2->SetBinError(i+1, e0[i]);
 }
      double v2bin[7] = {25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0};
      double v2nul[7] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
      double v2sysl[7] = {l0[0]-k0[0], l0[1]-k0[1], l0[2]-k0[2], l0[3]-k0[3], l0[4]-k0[4], l0[5]-k0[5], l0[6]-k0[6]};
      double v2sysh[7] = {k0[0]-h0[0], k0[1]-h0[1], k0[2]-h0[2], k0[3]-h0[3], abs(h0[4]-k0[4]), abs(h0[5]-k0[5]), abs(h0[6]-k0[6])};
      double v2nom[7] = {k0[0], k0[1], k0[2], k0[3], k0[4], k0[5], k0[6]};

      auto systsh = new TGraphAsymmErrors(8, v2bin, v2nom, v2nul, v2nul, v2sysh, v2sysl);
           systsh->SetMarkerSize(0.6);
           systsh->SetMarkerStyle(20);
           systsh->SetFillColorAlpha(kCaitieBlue, 0.15);
           systsh->SetMarkerColor(kCaitieDarkBlue);
           systsh->SetLineColor(kCaitieDarkBlue);



  hJetv2->SetLineColor(kCaitieDarkBlue);
  hJetv2->SetLineWidth(2);
  hJetv2->SetMarkerStyle(20);
  hJetv2->SetMarkerColor(kCaitieDarkBlue);
  hJetv2->SetMarkerSize(0.6);

  //=================================================================
  //============= Create TGraph Errors for Systematics ==============
  //=================================================================
      //Hedgehog
      int nbins = 0;
      const double *bincenters, *nul;
      if (R == 2) {
         bincenters = (const double*)bincentersR02;
         nul = (const double*)nulR02;
         nbins = 6;}
      if (R == 4) {
         bincenters = (const double*)bincentersR04;
         nul = (const double*) nulR04;
         nbins = 5;}

      auto sys_redmer = new TGraphErrors(7, redmerCent, jetv2, redmerNul, shapesyst);
      sys_redmer->SetMarkerSize(0.8);
      sys_redmer->SetMarkerStyle(20);
      sys_redmer->SetFillColorAlpha(kBlack, 0.0);
      sys_redmer->SetMarkerColor(kBlack);
      sys_redmer->SetLineColor(kBlack);

  //====================================
  //========== create lines ============
  //====================================
  TLine *line = new TLine(20.0, 1, 100.0, 1);
  line->SetLineStyle(2);


  //======================================
  //========== create legends ============
  //======================================
   auto datadetails = new TLegend(0.55, 0.75, 0.85, 0.85);
        datadetails->SetTextSize(0.035);
        datadetails->SetBorderSize(0);
        datadetails->AddEntry(paperFig, "2.76 TeV, PLB", "plf");
        datadetails->AddEntry(myv2, "5.02 TeV, Work in Progress", "plf");



  //====================================
  //========== draw histograms =========
  //====================================


 
  TCanvas *c2 = new TCanvas("Plot 2", "Plot 2", 800, 600);
  c2->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
          pad2->SetBottomMargin(0.1); // Upper and lower plot are joined
          pad2->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad2->SetRightMargin(0.05); // Upper and lower plot are joined
          pad2->SetTopMargin(0.05); // Upper and lower plot are joined
          pad2->Draw();
          pad2->cd();               // pad1 becomes the current pad   
         //pad1->SetLogx();
         //pad1->SetLogy();
             spacer2->SetTitle("");
             spacer2->SetMinimum(-0.1);
             spacer2->SetMaximum(0.2);
             spacer2->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             spacer2->GetXaxis()->SetTitleSize(0.035);
             spacer2->GetYaxis()->SetTitleSize(0.04);
             spacer2->GetXaxis()->SetRangeUser(20.0, 120.0);
             spacer2->SetYTitle("v2");
         spacer2->Draw("same");
         paperFig->Draw("same");
         sys_redmer->Draw("ezp 5 same");
         myv2->Draw("same");
         datadetails->Draw("same");
 }



void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}

void setcoloralpha(TH1F* h, int kcolor, double trans) {
  h->SetLineColorAlpha(kcolor, trans);
  h->SetMarkerColorAlpha(kcolor, trans);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}
 
