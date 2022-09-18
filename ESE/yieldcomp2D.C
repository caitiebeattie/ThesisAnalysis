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
TH1D *spacer = new TH1D("spacer", "spacer", 10, 30, 140);


//=======================================================================
//=======================================================================
//========================= main function ===============================
//=======================================================================
//=======================================================================

void yieldcomp2D (int iteration, string lofilename, string hifilename, int planestatus = 0) {

  const char* lofile = lofilename.c_str();
  const char* hifile = hifilename.c_str();
  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

      //set custom colors
      int kCaitieBlue, kCaitieDarkBlue;
        kCaitieBlue = TColor::GetColor("#00c6c9");
        kCaitieDarkBlue = TColor::GetColor("#00576A");
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
        kSunOrange = TColor::GetColor("#e86400"); //e86400
        kSunGreen = TColor::GetColor("#9FD064"); //green
        kDarkSunGreen = TColor::GetColor("#67AD2B");

   //==================================================
   //======== 1. get file tag for legend 
   //==================================================
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = "_";
   for (int i = 0; i < lofilename.length(); i++)    {
       const char *a = &lofile[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
       else ss << lofile[i];
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



  //==============================================================
  //========= 3. retrieve info from data files used to unfold
  //==============================================================

      //Get File
      TFile *specfile;
      if (r02 == true) specfile = new TFile("AnalysisResults8144.root");
      if (r04 == true) specfile = new TFile("AnalysisResults8119.root");
   
      //Get N_event for scaling
      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      stringstream tlistName;
      tlistName << "AliAnalysisTaskJetExtractor_Jet_AKTChargedR0" << R << "0_tracks_pT0150_pt_scheme_Rho_Jet_histos";
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get(tlistName.str().c_str()); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral();  //(8.16.22 - used to divide by a factor 5. Not sure why, but it's been removed)
      cout << "N_event: " << N_event <<"\n";




  //======================================================
  //========= 4. retrieve unfolded spectra from file
  //=======================================================
    TFile *lefilel = new TFile(lofile);
    TFile *lefileh = new TFile(hifile);

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
       tempHistLOP->GetYaxis()->SetRangeUser(0.0, sqrt(1)/2.0);
       tempHistLIP->GetYaxis()->SetRangeUser(sqrt(3)/2.0, 1.0);
       tempHistLow->GetYaxis()->SetRangeUser(0.0, 1.0);
       tempHistHOP->GetYaxis()->SetRangeUser(0.0, sqrt(1)/2.0);
       tempHistHIP->GetYaxis()->SetRangeUser(sqrt(3)/2.0, 1.0);
       tempHistHigh->GetYaxis()->SetRangeUser(0.0, 1.0);
       unfoldHistsLIP.push_back((TH1F*)tempHistLIP->ProjectionX());
       unfoldHistsLOP.push_back((TH1F*)tempHistLOP->ProjectionX());
       unfoldHistsLow.push_back((TH1F*)tempHistLow->ProjectionX());
       unfoldHistsHIP.push_back((TH1F*)tempHistHIP->ProjectionX());
       unfoldHistsHOP.push_back((TH1F*)tempHistHOP->ProjectionX());
       unfoldHistsHigh.push_back((TH1F*)tempHistHigh->ProjectionX());
             //unfoldHistsLIP[i]->Sumw2();
             //unfoldHistsLOP[i]->Sumw2();
             //unfoldHistsLow[i]->Sumw2();
             //unfoldHistsHIP[i]->Sumw2();
             //unfoldHistsHOP[i]->Sumw2();
             //unfoldHistsHigh[i]->Sumw2();
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
       TFile *recfile = new TFile("../Unfolding/ChargedJetRecEfficiencies.root");
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
          if (R == 2)   pTedgerebin  = {35, 40, 50, 65, 80, 100, 120};
          if (R == 4)   pTedgerebin  = {50, 65, 80, 100, 120};
      TH1F* efficiency1 = (TH1F*)effic1->Rebin(pTedgerebin.size()-1, "efficiency1", pTedgerebin.data());
      TH1F* efficiency2 = (TH1F*)effic2->Rebin(pTedgerebin.size()-1, "efficiency2", pTedgerebin.data());
      TH1F* recefficiency = (TH1F*)receffic->Rebin(pTedgerebin.size()-1, "recefficiency", pTedgerebin.data());
            recefficiency->Scale(1.0, "width");
           
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
 
      // 0.3 = fraction of events considered with ESE cuts
       yieldLow->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");
       yieldLowInPlane->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");
       yieldLowOutPlane->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");
       yieldHigh->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");
       yieldHighInPlane->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");
       yieldHighOutPlane->Scale(1.0/((double)N_event*TPCfid*(3./10.)), "width");


        //Unfolded Spectra 
        for (int i = 0; i < unfoldHistsLIP.size(); i++)    {
        unfoldHistsLIP[i]->Scale(1.0, "width");
        unfoldHistsLIP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsLOP[i]->Scale(1.0, "width");
        unfoldHistsLOP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsLow[i]->Scale(1.0, "width");
        unfoldHistsLow[i]->Scale(1.0/((double)N_event*TPCfid));
        }


        for (int i = 0; i < unfoldHistsHIP.size(); i++)    {
        unfoldHistsHIP[i]->Scale(1.0, "width");
        unfoldHistsHIP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsHOP[i]->Scale(1.0, "width");
        unfoldHistsHOP[i]->Scale(1.0/((double)N_event*TPCfid));
        unfoldHistsHigh[i]->Scale(1.0, "width");
        unfoldHistsHigh[i]->Scale(1.0/((double)N_event*TPCfid));
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
      if (R == 4)  lowedge = 50.0;
      int startbin = yieldratio->FindBin(lowedge+0.01);

      
      yieldratio->SetLineColor(kCaitieDarkGreen);
      yieldratio->SetMarkerColor(kCaitieDarkGreen);
      yieldratio->SetMarkerStyle(20);
      yieldratio->SetMarkerSize(0.8);
      yieldratio->SetLineWidth(2);

      setcolor(ratioOILow, kCaitiePink);
      setcolor(ratioOIHigh, kCaitieDarkBlue);

      setcolor(yieldHighInPlane, kSunBlue);
      setcolor(yieldHighOutPlane, kSunOrange);
        yieldHighInPlane->SetMarkerStyle(21);
        yieldHighOutPlane->SetMarkerStyle(21);

      setcolor(yieldLowInPlane, kSunPurple);
   
      setcolor(yieldLowOutPlane, kSunPink);

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
      double bincentersR02[5] = {37.5, 45.0, 57.5, 72.5, 90.0};
      double bincentersR04[3] = {57.5, 72.5, 90.0};
      double bincenters3050_hi[3] = {54.5, 69.5, 87.};//, 52., 62., 73., 88., 108.5};
        double exl_hi[3] = {4.5, 4.5, 7.};
        double exr_hi[3] = {10.5, 10.5, 13.};
      double bincenters3050_li[3] = {56.5, 71.5, 89.};//, 54., 64., 76., 91., 109.5};
        double exl_li[3] = {6.5, 6.5, 9.};
        double exr_li[3] = {8.5, 8.5, 11.};
      double bincenters3050_lo[3] = {58.5, 73.5, 91.};//, 56., 66., 79., 94., 110.5};
        double exl_lo[3] = {8.5, 8.5, 11.};
        double exr_lo[3] = {6.5, 6.5, 9.};
      double bincenters3050_ho[3] = {60.5, 75.5, 93.};//, 58., 68., 82., 97., 111.5};
        double exl_ho[3] = {10.5, 10.5, 13.};
        double exr_ho[3] = {4.5, 4.5, 7.};

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
      //4 sorted spectra
      double nomloin[8] = {yieldLowInPlane->GetBinContent(startbin), yieldLowInPlane->GetBinContent(startbin+1),
                          yieldLowInPlane->GetBinContent(startbin+2), yieldLowInPlane->GetBinContent(startbin+3),
                          yieldLowInPlane->GetBinContent(startbin+4), yieldLowInPlane->GetBinContent(startbin+5),
                          yieldLowInPlane->GetBinContent(startbin+6), yieldLowInPlane->GetBinContent(startbin+7)};
      double nomloout[8] = {yieldLowOutPlane->GetBinContent(startbin), yieldLowOutPlane->GetBinContent(startbin+1), 
                           yieldLowOutPlane->GetBinContent(startbin+2), yieldLowOutPlane->GetBinContent(startbin+3),
                           yieldLowOutPlane->GetBinContent(startbin+4), yieldLowOutPlane->GetBinContent(startbin+5),
                           yieldLowOutPlane->GetBinContent(startbin+6), yieldLowOutPlane->GetBinContent(startbin+7)};
      double nomhiout[8] = {yieldHighOutPlane->GetBinContent(startbin), yieldHighOutPlane->GetBinContent(startbin+1),
                            yieldHighOutPlane->GetBinContent(startbin+2), yieldHighOutPlane->GetBinContent(startbin+3), 
                            yieldHighOutPlane->GetBinContent(startbin+4), yieldHighOutPlane->GetBinContent(startbin+5),
                            yieldHighOutPlane->GetBinContent(startbin+6), yieldHighOutPlane->GetBinContent(startbin+7)};
      double nomhiin[8] = {yieldHighInPlane->GetBinContent(startbin), yieldHighInPlane->GetBinContent(startbin+1),
                           yieldHighInPlane->GetBinContent(startbin+2), yieldHighInPlane->GetBinContent(startbin+3),
                           yieldHighInPlane->GetBinContent(startbin+4), yieldHighInPlane->GetBinContent(startbin+5),
                           yieldHighInPlane->GetBinContent(startbin+6), yieldHighInPlane->GetBinContent(startbin+7)};
      double nul3050[8] = {2.5, 2.5, 5.0, 5.0, 5.0, 7.5, 7.5, 10.0};
      double nul3050_crop[8] = {0.75, 0.75, 1.5, 1.5, 1.5, 2.25, 2.25, 3.};
      double nul3050_short[6] = {2.5, 5.0, 5.0, 5.0, 7.5, 7.5};
      double nulR02[5] = {2.5, 5.0, 7.5, 7.5, 10.0};
      double nulR04[3] = {7.5, 7.5, 10.0};



      //------------ 
      //yield errors
      //-------------
     string date;
        if (R == 2) date = "Aug15";
        if (R == 4) date = "Sep9";
     TFile *eFileLo, *eFileHi, *eFileLoOut, *eFileLoIn, *eFileHiOut, *eFileHiIn;
      eFileLo = new TFile("../Systematics/Systematics3050_lo.root");
        eFileHi = new TFile("../Systematics/Systematics3050_hi.root");
        eFileLoOut = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_loout_%s.root", R, date.c_str()));
        eFileLoIn = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_loin_%s.root", R, date.c_str()));
        eFileHiOut = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_hiout_%s.root", R, date.c_str()));
        eFileHiIn = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_hiin_%s.root", R, date.c_str()));

      //get hists from file
      TH1F *errorsLo = (TH1F*)eFileLo->Get("h_err"); 
        TH1F *errorsHi = (TH1F*)eFileHi->Get("h_err"); 
        TH1F *errorsLoOut = (TH1F*)eFileLoOut->Get("h_err"); 
        TH1F *errorsLoIn = (TH1F*)eFileLoIn->Get("h_err"); 
        TH1F *errorsHighOut = (TH1F*)eFileHiOut->Get("h_err"); 
        TH1F *errorsHighIn = (TH1F*)eFileHiIn->Get("h_err"); 
      int startbinerror = errorsLo->FindBin(lowedge);  //was 30.01, changed 8.16.2022

      //get fractional errors
      double syst_FracLow[8] = {errorsLo->GetBinContent(startbinerror),  errorsLo->GetBinContent(startbinerror+1), 
                               errorsLo->GetBinContent(startbinerror+2), errorsLo->GetBinContent(startbinerror+3), 
                               errorsLo->GetBinContent(startbinerror+4), errorsLo->GetBinContent(startbinerror+5),
                               errorsLo->GetBinContent(startbinerror+6), errorsLo->GetBinContent(startbinerror+7)};
      double syst_FracHigh[8] = {errorsHi->GetBinContent(startbinerror),   errorsLo->GetBinContent(startbinerror+1),
                                errorsLo->GetBinContent(startbinerror+2),  errorsLo->GetBinContent(startbinerror+3),
                                errorsLo->GetBinContent(startbinerror+4),  errorsLo->GetBinContent(startbinerror+5),
                                errorsLo->GetBinContent(startbinerror+6), errorsLo->GetBinContent(startbinerror+7)};
      double syst_FracLowIn[8] = {errorsLoIn->GetBinContent(startbinerror),  errorsLoIn->GetBinContent(startbinerror+1),
                                 errorsLoIn->GetBinContent(startbinerror+2), errorsLoIn->GetBinContent(startbinerror+3),
                                 errorsLoIn->GetBinContent(startbinerror+4), errorsLoIn->GetBinContent(startbinerror+5),
                                 errorsLoIn->GetBinContent(startbinerror+6), errorsLoIn->GetBinContent(startbinerror+7)};
      double syst_FracLowOut[8] = {errorsLoOut->GetBinContent(startbinerror),   errorsLoOut->GetBinContent(startbinerror+1),
                                   errorsLoOut->GetBinContent(startbinerror+2), errorsLoOut->GetBinContent(startbinerror+3),
                                   errorsLoOut->GetBinContent(startbinerror+4), errorsLoOut->GetBinContent(startbinerror+5),
                                   errorsLoOut->GetBinContent(startbinerror+6), errorsLoOut->GetBinContent(startbinerror+7)};
      double syst_FracHighIn[8] = {errorsHighIn->GetBinContent(startbinerror),  errorsHighIn->GetBinContent(startbinerror+1),
                                   errorsHighIn->GetBinContent(startbinerror+2), errorsHighIn->GetBinContent(startbinerror+3),
                                   errorsHighIn->GetBinContent(startbinerror+4), errorsHighIn->GetBinContent(startbinerror+5),
                                   errorsHighIn->GetBinContent(startbinerror+6), errorsHighIn->GetBinContent(startbinerror+7)};
      double syst_FracHighOut[8] = {errorsHighOut->GetBinContent(startbinerror), errorsHighOut->GetBinContent(startbinerror+1),
                                    errorsHighOut->GetBinContent(startbinerror+2), errorsHighOut->GetBinContent(startbinerror+3),
                                    errorsHighOut->GetBinContent(startbinerror+4), errorsHighOut->GetBinContent(startbinerror+5),
                                    errorsHighOut->GetBinContent(startbinerror+6), errorsHighOut->GetBinContent(startbinerror+7)};
      //get absolute errors
      double syslo[8] = {yieldLow->GetBinContent(startbin)*syst_FracLow[0], yieldLow->GetBinContent(startbin+1)*syst_FracLow[1], 
                        yieldLow->GetBinContent(startbin+2)*syst_FracLow[2], yieldLow->GetBinContent(startbin+3)*syst_FracLow[3], 
                        yieldLow->GetBinContent(startbin+4)*syst_FracLow[4], yieldLow->GetBinContent(startbin+5)*syst_FracLow[5],
                        yieldLow->GetBinContent(startbin+6)*syst_FracLow[6], yieldLow->GetBinContent(startbin+7)*syst_FracLow[7]};
      double syshi[8] = {yieldHigh->GetBinContent(startbin)*syst_FracHigh[0], yieldHigh->GetBinContent(startbin+1)*syst_FracHigh[1],
                        yieldHigh->GetBinContent(startbin+2)*syst_FracHigh[2], yieldHigh->GetBinContent(startbin+3)*syst_FracHigh[3],
                        yieldHigh->GetBinContent(startbin+4)*syst_FracHigh[4], yieldHigh->GetBinContent(startbin+5)*syst_FracHigh[5],
                        yieldHigh->GetBinContent(startbin+6)*syst_FracHigh[6], yieldHigh->GetBinContent(startbin+7)*syst_FracHigh[7]};
      double sysloout[8] = {yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[0],
                            yieldLowOutPlane->GetBinContent(startbin+1)*syst_FracLowOut[1],
                            yieldLowOutPlane->GetBinContent(startbin+2)*syst_FracLowOut[2], 
                            yieldLowOutPlane->GetBinContent(startbin+3)*syst_FracLowOut[3], 
                            yieldLowOutPlane->GetBinContent(startbin+4)*syst_FracLowOut[4], 
                            yieldLowOutPlane->GetBinContent(startbin+5)*syst_FracLowOut[5],
                            yieldLowOutPlane->GetBinContent(startbin+6)*syst_FracLowOut[6],
                            yieldLowOutPlane->GetBinContent(startbin+7)*syst_FracLowOut[7]};
     double syshiout[8] = {yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[0], 
                           yieldHighOutPlane->GetBinContent(startbin+1)*syst_FracHighOut[1],
                           yieldHighOutPlane->GetBinContent(startbin+2)*syst_FracHighOut[2], 
                           yieldHighOutPlane->GetBinContent(startbin+3)*syst_FracHighOut[3],
                           yieldHighOutPlane->GetBinContent(startbin+4)*syst_FracHighOut[4],
                           yieldHighOutPlane->GetBinContent(startbin+5)*syst_FracHighOut[5], 
                           yieldHighOutPlane->GetBinContent(startbin+6)*syst_FracHighOut[6],
                           yieldHighOutPlane->GetBinContent(startbin+7)*syst_FracHighOut[7]};
     double sysloin[8] = {yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[0], yieldLowInPlane->GetBinContent(startbin+1)*syst_FracLowIn[1],
                          yieldLowInPlane->GetBinContent(startbin+2)*syst_FracLowIn[2], yieldLowInPlane->GetBinContent(startbin+3)*syst_FracLowIn[3],
                          yieldLowInPlane->GetBinContent(startbin+4)*syst_FracLowIn[4], yieldLowInPlane->GetBinContent(startbin+5)*syst_FracLowIn[5],
                          yieldLowInPlane->GetBinContent(startbin+6)*syst_FracLowIn[6], yieldLowInPlane->GetBinContent(startbin+7)*syst_FracLowIn[7]};
     double syshiin[8] = {yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[0], 
                          yieldHighInPlane->GetBinContent(startbin+1)*syst_FracHighIn[1],
                          yieldHighInPlane->GetBinContent(startbin+2)*syst_FracHighIn[2],
                          yieldHighInPlane->GetBinContent(startbin+3)*syst_FracHighIn[3],
                          yieldHighInPlane->GetBinContent(startbin+4)*syst_FracHighIn[4],
                          yieldHighInPlane->GetBinContent(startbin+5)*syst_FracHighIn[5],
                          yieldHighInPlane->GetBinContent(startbin+6)*syst_FracHighIn[6], 
                          yieldHighInPlane->GetBinContent(startbin+7)*syst_FracHighIn[7]};



      //------------ 
      //ratio errors
      //-------------
      TFile *eFileHiLo, *eFileHiInOut, *eFileLoInOut, *eFileAllInOut, *eFileHiInOutUp, *eFileLoInOutUp;
      eFileHiLo = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_hilo_%s.root", R, date.c_str()));
      eFileHiInOut = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_hiinout_%s.root", R, date.c_str()));
      eFileLoInOut = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_loinout_%s.root",R, date.c_str()));
      eFileAllInOut = new TFile("../Systematics/Systematics3050_allinout.root");
      eFileHiInOutUp = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_hiinoutup_%s.root", R, date.c_str()));
      eFileLoInOutUp = new TFile(Form("../Systematics/Systematics3050_R0%.0d0_loinoutup_%s.root", R, date.c_str()));
      TH1F *errorsHiLo = (TH1F*)eFileHiLo->Get("h_err"); 
      TH1F *errorsHiInOut = (TH1F*)eFileHiInOut->Get("h_err"); 
      TH1F *errorsLoInOut = (TH1F*)eFileLoInOut->Get("h_err");
      TH1F *errorsAllInOut = (TH1F*)eFileAllInOut->Get("h_err");
      TH1F *errorsHiInOutUp = (TH1F*)eFileHiInOutUp->Get("h_err"); 
      TH1F *errorsLoInOutUp = (TH1F*)eFileLoInOutUp->Get("h_err"); 
 
      
     //Low q2, out/in
     double syslooutin[5] = {lowvalues3050[0]*errorsLoInOut->GetBinContent(startbinerror), 
                             lowvalues3050[1]*errorsLoInOut->GetBinContent(startbinerror+1),
                             lowvalues3050[2]*errorsLoInOut->GetBinContent(startbinerror+2),
                             lowvalues3050[3]*errorsLoInOut->GetBinContent(startbinerror+3),
                             lowvalues3050[4]*errorsLoInOut->GetBinContent(startbinerror+4)};
     double syslooutinup[5] = {lowvalues3050[0]*errorsLoInOut->GetBinContent(startbinerror), 
                             lowvalues3050[1]*errorsLoInOutUp->GetBinContent(startbinerror+1),
                             lowvalues3050[2]*errorsLoInOutUp->GetBinContent(startbinerror+2),
                             lowvalues3050[3]*errorsLoInOutUp->GetBinContent(startbinerror+3),
                             lowvalues3050[4]*errorsLoInOutUp->GetBinContent(startbinerror+4)};
     //High q2, out/in 
     double syshioutin[5] = {highvalues3050[0]*errorsHiInOut->GetBinContent(startbinerror),  
                             highvalues3050[1]*errorsHiInOut->GetBinContent(startbinerror+1),
                             highvalues3050[2]*errorsHiInOut->GetBinContent(startbinerror+2),
                             highvalues3050[3]*errorsHiInOut->GetBinContent(startbinerror+3), 
                             highvalues3050[4]*errorsHiInOut->GetBinContent(startbinerror+4)};
     double syshioutinup[5] = {highvalues3050[0]*errorsHiInOut->GetBinContent(startbinerror),  
                             highvalues3050[1]*errorsHiInOutUp->GetBinContent(startbinerror+1),
                             highvalues3050[2]*errorsHiInOutUp->GetBinContent(startbinerror+2),
                             highvalues3050[3]*errorsHiInOutUp->GetBinContent(startbinerror+3), 
                             highvalues3050[4]*errorsHiInOutUp->GetBinContent(startbinerror+4)};
 
     //High/low q2
     double sysraerrors[7] = {1.23456, 0.183724, 0.0961714, 0.0676367, 0.0438335, 0.0414292, 0};
     double sysra[8] = {yieldratio->GetBinContent(startbin)*errorsHiLo->GetBinContent(startbinerror), 
                        yieldratio->GetBinContent(startbin+1)*errorsHiLo->GetBinContent(startbinerror+1), 
                        yieldratio->GetBinContent(startbin+2)*errorsHiLo->GetBinContent(startbinerror+2), 
                        yieldratio->GetBinContent(startbin+3)*errorsHiLo->GetBinContent(startbinerror+3),
                        yieldratio->GetBinContent(startbin+4)*errorsHiLo->GetBinContent(startbinerror+4), 
                        yieldratio->GetBinContent(startbin+5)*errorsHiLo->GetBinContent(startbinerror+5), 
                        yieldratio->GetBinContent(startbin+6)*errorsHiLo->GetBinContent(startbinerror+6),
                        yieldratio->GetBinContent(startbin+7)*errorsHiLo->GetBinContent(startbinerror+7)};
  

  //=================================================================
  //============= Create TGraph Errors for Systematics ==============
  //=================================================================
      //Hedgehog
      int nbins = 0;
      const double *bincenters, *nul;
      if (R == 2) {
         bincenters = (const double*)bincentersR02;
         nul = (const double*)nulR02;
         nbins = 5;}
      if (R == 4) {
         bincenters = (const double*)bincentersR04;
         nul = (const double*) nulR04;
         nbins = 3;}

      //ESE ratios
      auto systlo = new TGraphErrors(8, bincenters3050, nomlo, nul3050, syslo);
           systlo->SetMarkerSize(0.6);
           systlo->SetMarkerStyle(20);
           systlo->SetFillColorAlpha(kCaitieBlue, 0.15);
           systlo->SetMarkerColor(kCaitieBlue);
           systlo->SetLineColor(kCaitieBlue);
      auto systhi = new TGraphErrors(8, bincenters3050, nomhi, nul3050, syshi);
           systhi->SetMarkerSize(0.6);
           systhi->SetMarkerStyle(20);
           systhi->SetFillColorAlpha(kCaitiePink, 0.2);
           systhi->SetMarkerColor(kCaitiePink);
           systhi->SetLineColor(kCaitiePink);

      //ratio of high q2/low q2
      auto systra = new TGraphErrors(nbins, bincenters, nomrat, nul, sysra);
           systra->SetMarkerSize(0.8);
           systra->SetMarkerStyle(20);
           systra->SetFillColorAlpha(kCaitieGreen, 0.9);
           systra->SetMarkerColor(kCaitieDarkGreen);
           systra->SetLineColor(kCaitieDarkGreen);
           systra->SetLineWidth(2);

      //out/in
      /*auto systlooutin = new TGraphErrors(6, bincenters3050_short, lowvalues3050s, nul3050_short, syslooutin);
      systlooutin->SetMarkerSize(0.9);
      systlooutin->SetMarkerStyle(20);
      systlooutin->SetFillColorAlpha(kCaitieLightPink, 0.2);
      systlooutin->SetMarkerColor(kCaitiePink);
      systlooutin->SetLineColor(kCaitiePink);*/
      auto systlooutinAsym = new TGraphAsymmErrors(nbins, bincenters, lowvalues3050, nul, nul, syslooutin, syslooutinup); 
      systlooutinAsym->SetMarkerSize(0.9);
      systlooutinAsym->SetMarkerStyle(20);
      systlooutinAsym->SetFillColorAlpha(kCaitieLightPink, 0.2);
      systlooutinAsym->SetMarkerColor(kCaitiePink);
      systlooutinAsym->SetLineColor(kCaitiePink);

      /*auto systhioutin = new TGraphErrors(6, bincenters3050_short, highvalues3050s, nul3050_short, syshioutin);
      systhioutin->SetMarkerSize(0.9);
      systhioutin->SetMarkerStyle(21);
      systhioutin->SetFillColorAlpha(kCaitieBlue, 0.3);    
      systhioutin->SetMarkerColor(kCaitieDarkBlue);           
      systhioutin->SetLineColor(kCaitieDarkBlue);*/           
      auto systhioutinAsym = new TGraphAsymmErrors(nbins, bincenters, highvalues3050, nul, nul, syshioutin, syshioutinup);
      systhioutinAsym->SetMarkerSize(0.9);
      systhioutinAsym->SetMarkerStyle(21);
      systhioutinAsym->SetFillColorAlpha(kCaitieBlue, 0.3);   
      systhioutinAsym->SetMarkerColor(kCaitieDarkBlue);        
      systhioutinAsym->SetLineColor(kCaitieDarkBlue);       


      //4 yield spectra
      //auto sys_yieldLowIn = new TGraphAsymmErrors(nbins, bincenters3050_li, nomloin, exl_li, exr_li, sysloin, sysloin);
      auto sys_yieldLowIn = new TGraphErrors(nbins, bincenters, nomloin, nul, sysloin);
      sys_yieldLowIn->SetMarkerSize(1.0);
      sys_yieldLowIn->SetMarkerStyle(20);
      sys_yieldLowIn->SetFillColorAlpha(kSunPurple, 0.0);
      sys_yieldLowIn->SetLineWidth(2);
      sys_yieldLowIn->SetMarkerColor(kSunPurple);
      sys_yieldLowIn->SetLineColor(kSunPurple);

      //auto sys_yieldLowOut = new TGraphAsymmErrors(nbins, bincenters3050_lo, nomloout, exl_lo, exr_lo, sysloout, sysloout);
      auto sys_yieldLowOut = new TGraphErrors(nbins, bincenters, nomloout, nul, sysloout);
      sys_yieldLowOut->SetMarkerSize(1.0);
      sys_yieldLowOut->SetMarkerStyle(20);
      sys_yieldLowOut->SetFillColorAlpha(kSunPink, 0.0);
      sys_yieldLowOut->SetMarkerColor(kSunPink);
      sys_yieldLowOut->SetLineColor(kSunPink);
      sys_yieldLowOut->SetLineWidth(2);

      //auto sys_yieldHighIn = new TGraphAsymmErrors(nbins, bincenters3050_hi, nomhiin, exl_hi, exr_hi, syshiin, syshiin);
      auto sys_yieldHighIn = new TGraphErrors(nbins, bincenters, nomhiin, nul, syshiin);
      sys_yieldHighIn->SetMarkerSize(1.0);
      sys_yieldHighIn->SetMarkerStyle(21);
      sys_yieldHighIn->SetFillColorAlpha(kCaitieBlue, 0.4);
      sys_yieldHighIn->SetMarkerColor(kSunBlue);
      sys_yieldHighIn->SetLineColor(kSunBlue);
    
      //auto sys_yieldHighOut = new TGraphAsymmErrors(nbins, bincenters3050_ho, nomhiout, exl_ho, exr_ho, syshiout, syshiout);
      auto sys_yieldHighOut = new TGraphErrors(nbins, bincenters, nomhiout, nul, syshiout);
      sys_yieldHighOut->SetMarkerSize(1.0);
      sys_yieldHighOut->SetMarkerStyle(21);
      sys_yieldHighOut->SetFillColorAlpha(kSunOrange, 0.4);
      sys_yieldHighOut->SetLineWidth(2);
      sys_yieldHighOut->SetMarkerColor(kSunOrange);
      sys_yieldHighOut->SetLineColor(kSunOrange);


  //====================================
  //========== Polynomials =============
  //====================================
  //fit polynomial
  TF1* linFit1 = new TF1("linFit1", "pol0", 35,100);
  linFit1->SetRange(35, 100);
  ratioOILow->Fit(linFit1, "N R I", "");
  linFit1->SetLineColor(kCaitieDarkBlue); 
 
  TF1* linFit2 = new TF1("linFit2", "pol0", 35,100);
  linFit2->SetRange(35, 100);
  ratioOIHigh->Fit(linFit2, "N R I", "");
  linFit2->SetLineColor(kViolet);

  //====================================
  //========== create lines ============
  //====================================
  TLine *line = new TLine(lbound, 1, rbound, 1);
  line->SetLineStyle(2);
  TLine *sline = new TLine(30.0, 1, 100.0, 1);
  sline->SetLineStyle(2);
  TLine *mline = new TLine(35.0, 1, 100.0, 1);
  mline->SetLineStyle(2);

  TLine *dcut = new TLine(lbound, 0.9, rbound, 0.9);
  dcut->SetLineStyle(3);
  TLine *dcut2 = new TLine (lbound, 0.8, rbound, 0.8);
  dcut2->SetLineStyle(3);
  TLine *dcut3 = new TLine(lbound, 0.7, rbound, 0.7);
  dcut3->SetLineStyle(3);

  TLine *ucut = new TLine(lbound, 1.1, rbound, 1.1);
  ucut->SetLineStyle(3);
  TLine *ucut2 = new TLine(lbound, 1.2, rbound, 1.2);
  ucut2->SetLineStyle(3);
  TLine *ucut3 = new TLine(lbound, 1.3, rbound, 1.3);
  ucut3->SetLineStyle(3);


  //======================================
  //========== create legends ============
  //======================================
   auto fileleg = new TLegend(0.6, 0.8, 0.85, 0.85);
        fileleg->SetTextSize(0.055);
        fileleg->SetBorderSize(0);
        fileleg->AddEntry((TObject*)0, s2.str().c_str(), "");
   auto gen = new TLegend(0.12, 0.65, 0.3, 0.9);
        gen->SetTextSize(0.04);
        gen->SetBorderSize(0);
        gen->SetFillColorAlpha(kWhite, 0.0);
        gen->AddEntry((TObject*)0, "Work in Progress", "");
        if (cent == true)    gen->AddEntry((TObject*)0, Form("#sqrt{#it{s}_{NN}} = 5.02 TeV %.0f#font[122]{-}%.0f%% Pb#font[122]{-}Pb", centl, centr), "");
        if (semi == true)    gen->AddEntry((TObject*)0, Form("%.0f#font[122]{-}%.0f%% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", centl, centr), "");
        gen->AddEntry((TObject*)0, Form("Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.%.0d", R), "");
        gen->AddEntry((TObject*)0, Form("#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.%.0d", 9-R), "");
        //gen->AddEntry((TObject*)0, "Work in Progress", "");
   auto gen2 = new TLegend(0.25, 0.7, 0.7, 0.9);
        gen2->SetTextSize(0.035);
        gen2->SetBorderSize(0);
        gen2->AddEntry((TObject*)0, "Work in Progess", "");
        if (cent == true)    gen2->AddEntry((TObject*)0, Form("%.0f#font[122]{-}%.0f%% Pb#font[122]{-}Pb, 5.02 TeV", centl, centr), "");
        if (semi == true)    gen2->AddEntry((TObject*)0, Form("%.0f#font[122]{-}%.0f%% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV", centl, centr), "");
        gen2->AddEntry((TObject*)0, Form("Charged-particle jets, anti-#it{k}_{T}, #it{R} = 0.%.0d", R), "");
        gen2->AddEntry((TObject*)0, Form("#it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.%.0d", 9-R), "");
        //gen2->AddEntry((TObject*)0, "Work in Progress", "");
   auto ali = new TLegend(0.15, 0.6, 0.3, 0.85);
        ali->SetTextSize(0.05);
        ali->SetBorderSize(0);
        ali->AddEntry((TObject*)0, "Work in Progress", "");
        ali->AddEntry((TObject*)0, Form("#sqrt{s_{NN}} = 5.02 TeV, Pb#font[122]{-}Pb %.0f#font[122]{-}%.0f%%", centl, centr), "");
        ali->AddEntry((TObject*)0, Form("Charged Jets, anti-#it{k}_{T}, #it{R} = 0.%.0d, #it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.%.0d", R, 9-R), "");
        //ali->AddEntry((TObject*)0, "Work in Progress", "");
   auto datadetails = new TLegend(0.65, 0.75, 0.85, 0.9);
        datadetails->SetTextSize(0.035);
        datadetails->SetBorderSize(0);
        datadetails->AddEntry(systhioutinAsym, "30\% highest #it{q}_{2}^{V0M}", "plf");
        datadetails->AddEntry(systlooutinAsym, "10-40\% lowest #it{q}_{2}^{V0M}", "plf");

   auto esedetails = new TLegend(0.55, 0.15, 0.75, 0.3);
        esedetails->SetTextSize(0.04);
        esedetails->SetBorderSize(0);
        esedetails->AddEntry(systra, "large (70\%-100\%)/small (10\%-40\%) #it{q}_{2}^{V0M}", "plf");

   auto specleg = new TLegend(0.6, 0.4, 0.85, 0.525);
        specleg->SetTextSize(0.055);
        specleg->SetBorderSize(0);
        specleg->AddEntry(systhi, "30\% highest #it{q}_{2}", "plf");
        specleg->AddEntry(systlo, "30\% lowest #it{q}_{2}", "plf");  

   auto yieldLeg = new TLegend(0.2, 0.15, 0.4, 0.3);
        yieldLeg->SetTextSize(0.03);
        yieldLeg->SetBorderSize(0);
        yieldLeg->AddEntry(sys_yieldHighIn, "30\% large #it{q}_{2}, in-plane", "plf");
        yieldLeg->AddEntry(sys_yieldHighOut, "30\% large #it{q}_{2}, out-of-plane", "plf");
        yieldLeg->AddEntry(sys_yieldLowIn, "10-40\% small #it{q}_{2}, in-plane", "plf");
        yieldLeg->AddEntry(sys_yieldLowOut, "10-40\% small #it{q}_{2}, out-of-plane", "plf");


  //====================================
  //========== draw histograms =========
  //====================================
  //In-Plane/Out-of-Plane
  TCanvas *c1 = new TCanvas("QM Plot 1", lofile, 800, 600);
  c1->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad1->SetRightMargin(0.05); // Upper and lower plot are joined
          pad1->SetTopMargin(0.05); // Upper and lower plot are joined
          pad1->Draw();
          pad1->cd();               // pad1 becomes the current pad   
         //pad1->SetLogx();
         //pad1->SetLogy();
             spacer->SetTitle("");
             spacer->SetMinimum(0.01);
             spacer->SetMaximum(1.99);
             spacer->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             spacer->GetXaxis()->SetTitleSize(0.035);
             spacer->GetYaxis()->SetTitleSize(0.04);
             spacer->GetXaxis()->SetRangeUser(35.0, 140.0);
             spacer->SetYTitle("Ratio (out-of-plane/in-plane)");
         spacer->Draw("same");
         ratioOILow->GetXaxis()->SetRangeUser(ratioOILow->GetBinLowEdge(startbin), 120.0);
         ratioOIHigh->GetXaxis()->SetRangeUser(ratioOIHigh->GetBinLowEdge(startbin), 120.0);
         ratioOILow->Draw("same");
         ratioOIHigh->Draw("same");
         systlooutinAsym->Draw("ezp 5 same");
         systhioutinAsym->GetXaxis()->SetRangeUser(35.0, 120.0);
         systhioutinAsym->Draw("ezp 5 same");
         ratioOILow->Draw("same");
         ratioOIHigh->Draw("same");
     //caitie->Draw("ezp 2 same");
   gen->Draw("same");
   datadetails->Draw("same");
   sline->Draw("same");
   //linFit1->Draw("same");
   //linFit2->Draw("same");
       cout <<"\n";
       cout << "Fractional Statistical Error\n";
       cout << "pT Range      YieldLowInPlane      YieldHighInPlane     YieldLowOutPlane      YieldHighOutPlane\n"; 
       cout << "-------------------------------------------------------------------------------------------------\n";
   for (int i = 0; i < ratioOILow->GetNbinsX(); i++)  {
       cout << ratioOILow->GetBinLowEdge(i+1) << "-" << ratioOILow->GetBinLowEdge(i+2) <<" GeV	" 
            << yieldLowInPlane->GetBinError(i+1)/yieldLowInPlane->GetBinContent(i+1) << "		" 
            << yieldHighInPlane->GetBinError(i+1)/yieldHighInPlane->GetBinContent(i+1)<< "		" 
            << yieldLowOutPlane->GetBinError(i+1)/yieldLowOutPlane->GetBinContent(i+1)<< "		"
            << yieldHighOutPlane->GetBinError(i+1)/yieldHighOutPlane->GetBinContent(i+1)<<"\n"; 
       }
       cout << "-------------------------------------------------------------------------------------------------\n";
       cout << "pT Range	RatioLow	RatioHigh\n";
       cout << "-------------------------------------------------------------------------------------------------\n";
   for (int i = 0; i < ratioOILow->GetNbinsX(); i++)  {
       cout << ratioOILow->GetBinLowEdge(i+1) << "-" << ratioOILow->GetBinLowEdge(i+2) <<" GeV  "
            << ratioOILow->GetBinError(i+1)/ratioOILow->GetBinContent(i+1) << "	"
            << ratioOIHigh->GetBinError(i+1)/ratioOIHigh->GetBinContent(i+1)<<"\n"; 
       }

  TCanvas *c2 = new TCanvas("QM Plot 2", lofile, 800, 600);
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
             yieldratio->SetTitle("");
             yieldratio->SetMinimum(0.01);
             yieldratio->SetMaximum(2.39);
             yieldratio->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             yieldratio->GetXaxis()->SetTitleOffset(1.15);
             yieldratio->GetXaxis()->SetTitleSize(0.04);
             yieldratio->GetYaxis()->SetTitleSize(0.04);
             yieldratio->GetXaxis()->SetRangeUser(yieldratio->GetBinLowEdge(startbin), 120.0);
             spacer->SetYTitle("Ratio (large/small #it{q}_{2})");
         spacer->Draw("same");
         yieldratio->Draw("same");
         systra->Draw("ezp 5 same");
         yieldratio->Draw("same");
     //caitie->Draw("ezp 2 same");
   gen->Draw("same");
   esedetails->Draw("same");
   sline->Draw("same");



   stringstream foutname;
   foutname << "ESEplots_R0" << R << ".root";
   TFile *fout=new TFile (foutname.str().c_str(), "RECREATE");
   fout->cd();
   yieldratio->Write();
   systra->SetName("systematics");
   systra->Write();

/*

  //Overall ESE
  TCanvas *c2 = new TCanvas("QM Plot 2", "raw style", 600, 500);
  c2->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad2a = new TPad("pad2a", "", 0.0, 0.4, 1.0, 1.0);
          pad2a->SetBottomMargin(0.0); // Upper and lower plot are joined
          pad2a->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad2a->Draw();
       TPad *pad2b = new TPad("pad2b", "", 0.0, 0.05, 1.0, 0.4);
          pad2b->SetBottomMargin(0.2); // Upper and lower plot are joined
          pad2b->SetTopMargin(0.0); // Upper and lower plot are joined
          pad2b->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad2b->Draw();
          pad2a->cd();               // pad1 becomes the current pad   
         //pad1a->SetLogx();
          pad2a->SetLogy();
       //setcoloralpha(yieldLow, kBlack, 0.0);
       //setcoloralpha(yieldHigh, kBlack, 0.0);
       yieldLow->SetMarkerColor(kCaitieBlue);
       yieldHigh->SetMarkerColor(kViolet);
       yieldLow->SetMaximum(yieldLow->GetBinContent(2)*50);
       yieldLow->SetTitle("");
       yieldLow->GetXaxis()->SetRangeUser(35.0, 100.0);
       yieldLow->GetYaxis()->SetLabelSize(0.05);
       yieldLow->Draw("same");
       yieldHigh->Draw("same");
       systlo->Draw("ezp 2 same");
       specleg->Draw("same");
       ali->Draw("same");
     c2->cd();
     pad2b->cd();
       yieldratio->SetTitle("");
       yieldratio->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       yieldratio->SetYTitle("Ratio (high/low)");
       yieldratio->GetXaxis()->SetLabelSize(0.08);
       yieldratio->GetYaxis()->SetLabelSize(0.08);
       yieldratio->GetYaxis()->SetTitleSize(0.09);
       yieldratio->GetXaxis()->SetTitleSize(0.09);
       yieldratio->SetMinimum(0.01);
       yieldratio->SetMaximum(2.49); 
       yieldratio->GetXaxis()->SetRangeUser(35.0, 100.0);
       yieldratio->Draw("same");
       systra->Draw("ezp 3 same");
       yieldratio->Draw("same");
       sline->Draw("same"); 
*/
 

  //Yields
  TCanvas *c3 = new TCanvas("Spectra", "Spectra", 600, 700);
  c3->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad3 = new TPad("pad2a", "", 0.0, 0.05, 1.0, 1.0);
          pad3->SetBottomMargin(0.1); // Upper and lower plot are joined
          pad3->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad3->SetRightMargin(0.05); // Upper and lower plot are joined
          pad3->SetTopMargin(0.05); // Upper and lower plot are joined
          pad3->Draw();
          pad3->cd();
          pad3->SetLogy();
       yieldLowInPlane->GetXaxis()->SetRangeUser(35.0, 100.0);
       yieldLowInPlane->SetTitle("");
       yieldLowInPlane->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       yieldLowInPlane->GetXaxis()->SetTitleOffset(1.2);
       yieldLowInPlane->SetYTitle("#frac{1}{N_{event}}#frac{d^{2}N}{d#it{p}_{T}d#it{#eta}}");
       yieldLowInPlane->SetMaximum(2*pow(10,-4));
       yieldLowInPlane->SetMinimum(8*pow(10,-8));
       yieldLowInPlane->Draw("same");
       //yieldHighInPlane->Draw("same");
       //yieldHighOutPlane->Draw("same");
       //yieldLowOutPlane->Draw("same");
       yieldLeg->Draw("same");
       sys_yieldHighIn->Draw("ezp 2 same");
       sys_yieldHighOut->Draw("ezp 2 same");
       sys_yieldLowOut->Draw("ezp 5 same");
       sys_yieldLowIn->Draw("ezp 5 same");
       yieldHighInPlane->Draw("same ex0");
       yieldHighOutPlane->Draw("same ex0");
       yieldLowInPlane->Draw("same ex0");
       yieldLowOutPlane->Draw("same ex0");
       gen2->Draw("same");

   stringstream foutname2;
   foutname2 << "spectraplots_R0" << R << ".root";
   TFile *fout2=new TFile (foutname2.str().c_str(), "RECREATE");
   fout2->cd();
   yieldLow->SetName("yieldLow");
   yieldLow->Write();
   yieldHigh->SetName("yieldHigh");
   yieldHigh->Write();
   yieldLowInPlane->SetName("yieldLowInPlane");
   yieldLowInPlane->Write();
   yieldLowOutPlane->SetName("yieldLowOutPlane");
   yieldLowOutPlane->Write();
   yieldHighInPlane->SetName("yieldHighInPlane");
   yieldHighInPlane->Write();
   yieldHighOutPlane->SetName("yieldHighOutPlane");
   yieldHighOutPlane->Write();
  

  TCanvas *c4 = new TCanvas("rec", "rec", 500, 500);
  c4->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
       pad4->Draw();
       pad4->cd();
       recefficiency->Draw("same");

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
 
