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

double pi = 3.14159265359;
double lbound = 10.0;
double rbound = 190.0;
double eventsperjet = 49.5;

void setcolor(TH1F* h, int kcolor);
void setcoloralpha(TH1F* h, int kcolor, double trans);
void updateR (double &phiscale, double &etascale, double jetR) {
    phiscale = (2*pi)/(1.92 - 2*jetR);
    etascale = (1.4 - 2*jetR);
}

bool cent = false;
bool semi = false;


int centcolor = kViolet+5;
int semicolor = kCyan-4;

int ci = 1756; // color index
TColor *bluecolor = new TColor(ci, 0, 125, 197);

vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};



//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void yieldcomp2D (int iteration, string lofilename, string hifilename, int planestatus = 0) {

  const char* lofile = lofilename.c_str();
  const char* hifile = hifilename.c_str();
  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);
      int kCaitieBlue, kCaitieDarkBlue;
      kCaitieBlue = TColor::GetColor("#13C2E5");
      kCaitieDarkBlue = TColor::GetColor("#0496B4");

      int kCaitieGreen, kCaitieDarkGreen;
      kCaitieGreen = TColor::GetColor("#A5F6AE");
      kCaitieDarkGreen = TColor::GetColor("#199C28");

      int kCaitiePink;
      kCaitiePink = TColor::GetColor("#FF12D8");

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
      for (int i = 2; i < words.size() - 1; i++) {
        s2 << words[i] << " ";
        const char *newword = words[i].c_str();
        if (strncmp(newword, ccent, 3) == 0)   cent = true;
        if (strncmp(newword, scent3050, 3) == 0)   semi = true;
        if (strncmp(newword, scent3035, 3) == 0 || strncmp(newword, scent3540, 3) == 0)   semi = true;
        if (strncmp(newword, scent4045, 3) == 0 || strncmp(newword, scent4550, 3) == 0)   semi = true;
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
  //==============================================================
  //========= 3. retrieve info from data files used to unfold
  //==============================================================

      //Get File
      TFile *specfile = new TFile("AnalysisResults7812.root");
   
      //Get N_event for scaling
      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral()/5.0;
      cout << "N_event: " << N_event <<"\n";


      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)specfile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2V0C_EPV0C");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                             //number percentile we're on
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      int centdiff = (int)centr-(int)centl;
      int leftc, rightc;
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < centdiff; j++) {
        leftc = (int)centl + j; 
        rightc = (int)centl + j + 1;    
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);
        TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();
        q2hist->Rebin(100);
        long double percentilemarker = (long double)q2hist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        for (int i = 0; i <= 100; i++)   {
          percentileticker += q2hist->GetBinContent(i);         //increment the percentileticker
          if (percentileticker >= percentilen*percentilemarker)  {
              percentiles.push_back(q2hist->GetBinCenter(i));
              percentilen++;
          }
        if (i == 100)  {
            allPercent.push_back(percentiles);
            percentiles.resize(0);
            break;}             
        }}

  //======================================================
  //========= 4. retrieve unfolded spectra from file
  //=======================================================
    TFile *lefilel = new TFile(lofile);
    TFile *lefileh = new TFile(hifile);
      double phifrac, etafrac;
      updateR(phifrac, etafrac, 0.4);

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
    }


  //==============================================================
  //========= 3. Pull statistical error info
  //==============================================================
      //Get Statistical Uncertainty
      vector<double> stat_lo_OP(7);
      vector<double> stat_lo_IP(7);
      vector<double> stat_lo(7);
      vector<double> stat_hi_OP(7);
      vector<double> stat_hi_IP(7);
      vector<double> stat_hi(7);
      vector<vector<double>> stat_vector(0);
      
      for (int i = 0; i < 7; i++)   stat_lo_OP[i] = unfoldHistsLOP[iteration]->GetBinError(i+unfoldHistsLOP[iteration]->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_lo_IP[i] = unfoldHistsLIP[iteration]->GetBinError(i+unfoldHistsLIP[iteration]->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_lo[i] = unfoldHistsLow[iteration]->GetBinError(i+unfoldHistsLow[iteration]->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_hi_OP[i] = unfoldHistsHOP[iteration]->GetBinError(i+unfoldHistsHOP[iteration]->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_hi_IP[i] = unfoldHistsHIP[iteration]->GetBinError(i+unfoldHistsHIP[iteration]->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_hi[i] = unfoldHistsHigh[iteration]->GetBinError(i+unfoldHistsHigh[iteration]->FindBin(25.01));
      
      stat_vector.push_back(stat_lo_OP);
      stat_vector.push_back(stat_lo_IP);
      stat_vector.push_back(stat_lo);
      stat_vector.push_back(stat_hi_OP);
      stat_vector.push_back(stat_hi_IP);
      stat_vector.push_back(stat_hi);
      
      //convert to fractional error
      for (int z = 0; z < 6; z++) {
        for (int i = 0; i < 7; i++)        {
          stat_vector[z][i] = (1.0/(double)stat_vector[z][i])*pow(stat_vector[z][i], 0.5);
        }
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
             effic2->Divide(kinpos2);

       //Get Reconstruction Efficiency
       TFile *recfile = new TFile("../Unfolding/ChargedJetRecEfficiencies.root");
       TH1D  *receffic = (TH1D*)recfile->Get("RecEff_R020_5GeV");


  //=============================================
  //============= Rebin PbPb to match pp ========
  //=============================================
      vector<double> pTedgerebin = {25, 30, 40, 50, 60, 70, 85, 100};
      TH1F* efficiency1 = (TH1F*)effic1->Rebin(pTedgerebin.size()-1, "efficiency1", pTedgerebin.data());
          efficiency1->Scale(1.0, "width");
      TH1F* efficiency2 = (TH1F*)effic2->Rebin(pTedgerebin.size()-1, "efficiency2", pTedgerebin.data());
          efficiency2->Scale(1.0, "width");
      TH1F* recefficiency = (TH1F*)receffic->Rebin(pTedgerebin.size()-1, "recefficiency", pTedgerebin.data());
          recefficiency->Scale(1.0, "width");
      TH1F* yieldLow = (TH1F*)unfoldHistsLow[iteration]->Rebin(pTedgerebin.size()-1, "yieldLow", pTedgerebin.data());
      TH1F* yieldLowInPlane = (TH1F*)unfoldHistsLIP[iteration]->Rebin(pTedgerebin.size()-1, "yieldLowInPlane", pTedgerebin.data());
      TH1F* yieldLowOutPlane = (TH1F*)unfoldHistsLOP[iteration]->Rebin(pTedgerebin.size()-1, "yieldLowOutPlane", pTedgerebin.data());
      TH1F* yieldHigh = (TH1F*)unfoldHistsHigh[iteration]->Rebin(pTedgerebin.size()-1, "yieldHigh", pTedgerebin.data());
      TH1F* yieldHighInPlane = (TH1F*)unfoldHistsHIP[iteration]->Rebin(pTedgerebin.size()-1, "yieldHighInPlane", pTedgerebin.data());
      TH1F* yieldHighOutPlane = (TH1F*)unfoldHistsHOP[iteration]->Rebin(pTedgerebin.size()-1, "yieldHighOutPlane", pTedgerebin.data());
  


//============================================================
  //========= 6. apply scaling for bin width + N_events ========
  //============================================================
        //Unfolded Spectra 
        for (int i = 0; i < unfoldHistsLIP.size(); i++)    {
        unfoldHistsLIP[i]->Scale(1.0, "width");
        //unfoldhistsl[i]->Scale(phifrac/etafrac);
        unfoldHistsLIP[i]->Scale(1.0/((double)N_event*1.4));
        unfoldHistsLOP[i]->Scale(1.0, "width");
        unfoldHistsLOP[i]->Scale(1.0/((double)N_event*1.4));
        unfoldHistsLow[i]->Scale(1.0, "width");
        unfoldHistsLow[i]->Scale(1.0/((double)N_event*1.4));
        }


        for (int i = 0; i < unfoldHistsHIP.size(); i++)    {
        unfoldHistsHIP[i]->Scale(1.0, "width");
        //unfoldhistsl[i]->Scale(phifrac/etafrac);
        unfoldHistsHIP[i]->Scale(1.0/((double)N_event*1.4));
        unfoldHistsHOP[i]->Scale(1.0, "width");
        unfoldHistsHOP[i]->Scale(1.0/((double)N_event*1.4));
        unfoldHistsHigh[i]->Scale(1.0, "width");
        unfoldHistsHigh[i]->Scale(1.0/((double)N_event*1.4));
        }


       yieldLow->Scale(1.0/((double)N_event*1.4), "width");
       yieldLowInPlane->Scale(1.0/((double)N_event*1.4), "width");
       yieldLowOutPlane->Scale(1.0/((double)N_event*1.4), "width");
       yieldHigh->Scale(1.0/((double)N_event*1.4), "width");
       yieldHighInPlane->Scale(1.0/((double)N_event*1.4), "width");
       yieldHighOutPlane->Scale(1.0/((double)N_event*1.4), "width");

       yieldLow->Divide(efficiency1);
       yieldLow->Divide(recefficiency);
       yieldHigh->Divide(efficiency2);
       yieldHigh->Divide(recefficiency);

 
  //=====================================================
  //============= Calculate the Yield Ratio =============
  //=====================================================
 
      TH1F* yieldratio =(TH1F*)yieldHigh->Clone("yieldratio");
      yieldratio->Divide(yieldLow);
        TH1F* ratioOIHigh = (TH1F*)yieldHighOutPlane->Clone("ratioOIHigh");
        ratioOIHigh->Divide(yieldHighInPlane);
        TH1F* ratioOILow = (TH1F*)yieldLowOutPlane->Clone("ratioOILow");
        ratioOILow->Divide(yieldLowInPlane);

      double lowedge;
      if (cent == true)  lowedge = 60.0;
      if (semi == true)  lowedge = 30.0;
      int startbin = yieldratio->FindBin(lowedge+0.01);
      yieldratio->GetXaxis()->SetRangeUser(lowedge, 100.0);
      ratioOIHigh->GetXaxis()->SetRangeUser(lowedge, 100.0);
      ratioOILow->GetXaxis()->SetRangeUser(lowedge, 100.0);

      
      yieldratio->SetLineColorAlpha(kCaitieDarkGreen, 0.0);
      yieldratio->SetMarkerColor(kCaitieDarkGreen);
      yieldratio->SetMarkerStyle(20);
      yieldratio->SetMarkerSize(0.6);

      setcolor(ratioOILow, kCaitieDarkBlue);
      setcolor(ratioOIHigh, kViolet-4);

      setcolor(yieldHighInPlane, kRed);
      setcolor(yieldHighOutPlane, kOrange);
      setcolor(yieldLowInPlane, kGreen);
      setcolor(yieldLowOutPlane, kBlue);


  //===============================================================
  //============= Calculate Statistical Uncertainties =============
  //===============================================================
     
      //get spectra errors
      //fractional
      vector<double> rat_vector(7);
      for (int i = 0; i < 7; i++)   rat_vector[i] = pow(pow(stat_vector[2][i],2) + pow(stat_vector[5][i],2), 0.5);
      //absolute
      for (int i = 0; i < 7; i++)   rat_vector[i] = rat_vector[i]*yieldratio->GetBinContent(i+yieldratio->FindBin(25.01)); 
      //absolute
      for (int i = 0; i < 7; i++)   stat_vector[2][i] = stat_vector[2][i]*yieldLow->GetBinContent(i+yieldLow->FindBin(25.01));   //low q2
      for (int i = 0; i < 7; i++)   stat_vector[5][i] = stat_vector[5][i]*yieldHigh->GetBinContent(i+yieldHigh->FindBin(25.01));   //high q2
      //sum out- and in-plane statistical errors in quadrature
      vector<double> stat_low_OI(7); //stat_lo/hi[0] -> 25-30 GeV
      vector<double> stat_high_OI(7);
      for (int i = 0; i < 7; i++)   stat_low_OI[i] = pow(pow(stat_vector[0][i],2) + pow(stat_vector[1][i],2), 0.5);
      for (int i = 0; i < 7; i++)   stat_high_OI[i] = pow(pow(stat_vector[3][i],2) + pow(stat_vector[4][i],2), 0.5);
      //multiply fractional errors by yield ratio to get absolute errors
      vector<double> stat_low(7);
      vector<double> stat_high(7);
      for (int i = 0; i < 7; i++)   stat_low[i] = stat_low[i]*ratioOILow->GetBinContent(i+ratioOILow->FindBin(25.01));
      for (int i = 0; i < 7; i++)   stat_high[i] = stat_high[i]*ratioOIHigh->GetBinContent(i+ratioOIHigh->FindBin(25.01));


      //create TGraphErrors for statistical errors
      double bincenters010[3] = {65.0, 77.5, 92.5};
      double bincenters3050[6] = {35.0, 45.0, 55.0, 65.0, 77.5, 92.5};
      double lowvalues010[3] = {ratioOILow->GetBinContent(startbin), ratioOILow->GetBinContent(startbin+1), ratioOILow->GetBinContent(startbin+2)};
      double lowvalues3050[6] = {ratioOILow->GetBinContent(startbin), ratioOILow->GetBinContent(startbin+1), ratioOILow->GetBinContent(startbin+2), ratioOILow->GetBinContent(startbin+3), ratioOILow->GetBinContent(startbin+4), ratioOILow->GetBinContent(startbin+5)};
      double highvalues3050[6] = {ratioOIHigh->GetBinContent(startbin), ratioOIHigh->GetBinContent(startbin+1), ratioOIHigh->GetBinContent(startbin+2), ratioOIHigh->GetBinContent(startbin+3), ratioOIHigh->GetBinContent(startbin+4), ratioOIHigh->GetBinContent(startbin+5)};
      double nomlo[6] = {yieldLow->GetBinContent(startbin), yieldLow->GetBinContent(startbin+1), yieldLow->GetBinContent(startbin+2), yieldLow->GetBinContent(startbin+3), yieldLow->GetBinContent(startbin+4), yieldLow->GetBinContent(startbin+5)};
      double nomhi[6] = {yieldHigh->GetBinContent(startbin), yieldHigh->GetBinContent(startbin+1), yieldHigh->GetBinContent(startbin+2), yieldHigh->GetBinContent(startbin+3), yieldHigh->GetBinContent(startbin+4), yieldHigh->GetBinContent(startbin+5)};
      double nomrat[6] = {yieldratio->GetBinContent(startbin), yieldratio->GetBinContent(startbin+1), yieldratio->GetBinContent(startbin+2), yieldratio->GetBinContent(startbin+3), yieldratio->GetBinContent(startbin+4), yieldratio->GetBinContent(startbin+5)};
      double nul010[3] = {5.0, 7.5, 7.5};
      double nul3050[6] = {5.0, 5.0, 5.0, 5.0, 7.5, 7.5};
      //convert stat_comb to correct type
      double stat_low010[3] = {stat_low[2], stat_low[3], stat_low[4]};
      double stat_low3050[6] = {stat_low[1], stat_low[2], stat_low[3], stat_low[4], stat_low[5], stat_low[6]};
      double statlo3050[6] = {stat_vector[2][1], stat_vector[2][2], stat_vector[2][3], stat_vector[2][4], stat_vector[2][5], stat_vector[2][6]};
      double stathi3050[6] = {stat_vector[5][1], stat_vector[5][2], stat_vector[5][3], stat_vector[5][4], stat_vector[5][5], stat_vector[5][6]};
      double statra3050[6] = {rat_vector[1], rat_vector[2], rat_vector[3], rat_vector[4], rat_vector[5], rat_vector[6]};
      //input to TGraphErrors
      TGraphErrors *statlo, *stathi, *statlio, *stathio, *stat, *statra;
      if (cent == true) stat = new TGraphErrors(3, bincenters010, lowvalues010, nul010, stat_low010);
      if (semi == true) stat = new TGraphErrors(6, bincenters3050, lowvalues3050, nul3050, stat_low3050);
           stat->SetMarkerStyle(20);
           stat->SetMarkerSize(0.7);
           stat->SetMarkerColor(kBlack);
           stat->SetLineColor(kGreen);
      statlo = new TGraphErrors(6, bincenters3050, nomlo, nul3050, statlo3050);
           statlo->SetLineColor(kCaitieDarkBlue);
           statlo->SetMarkerColor(kCaitieDarkBlue);
           statlo->SetLineWidth(1);
           statlo->SetMarkerSize(0.6);
           statlo->SetMarkerStyle(20);
      stathi = new TGraphErrors(6, bincenters3050, nomhi, nul3050, stathi3050);
           stathi->SetLineColor(kViolet);
           stathi->SetMarkerColor(kViolet);
           stathi->SetLineWidth(1);
           stathi->SetMarkerSize(0.6);
           stathi->SetMarkerStyle(20);
      statra = new TGraphErrors(6, bincenters3050, nomrat, nul3050, statra3050);
           statra->SetLineColor(kCaitieDarkGreen);


  //===============================================================
  //============= Calculate Systematic Uncertainties ==============
  //===============================================================
      //------------ 
      //yield errors
      //-------------
      TFile *eFileLo, *eFileHi, *eFileLoOut, *eFileLoIn, *eFileHiOut, *eFileHiIn;
      eFileLo = new TFile("../Systematics/Systematics3050_lo.root");
        eFileHi = new TFile("../Systematics/Systematics3050_hi.root");
        eFileLoOut = new TFile("../Systematics/Systematics3050_loout.root");
        eFileLoIn = new TFile("../Systematics/Systematics3050_loin.root");
        eFileHiOut = new TFile("../Systematics/Systematics3050_hiout.root");
        eFileHiIn = new TFile("../Systematics/Systematics3050_hiin.root");
      //get hists from file
      TH1F *errorsLo = (TH1F*)eFileLo->Get("h_err"); 
        TH1F *errorsHi = (TH1F*)eFileHi->Get("h_err"); 
        TH1F *errorsLoOut = (TH1F*)eFileLoOut->Get("h_err"); 
        TH1F *errorsLoIn = (TH1F*)eFileHLoIn->Get("h_err"); 
        TH1F *errorsHiOut = (TH1F*)eFileHiOut->Get("h_err"); 
        TH1F *errorsHiIn = (TH1F*)eFileHiIn->Get("h_err"); 
      int startbinerror = errorsLo->FindBin(30.01);
      //get fractional errors
      double syst_FracLow[6] = {errorsLo->GetBinContent(startbinerror),  errorsLo->GetBinContent(startbinerror+1), 
                               errorsLo->GetBinContent(startbinerror+2), errorsLo->GetBinContent(startbinerror+3), 
                               errorsLo->GetBinContent(startbinerror+4), errorsLo->GetBinContent(startbinerror+5)};
      double syst_FracHigh[6] = {errorsHi->GetBinContent(startbinerror),   errorsLo->GetBinContent(startbinerror+1),
                                errorsLo->GetBinContent(startbinerror+2),  errorsLo->GetBinContent(startbinerror+3),
                                errorsLo->GetBinContent(startbinerror+4),  errorsLo->GetBinContent(startbinerror+5)};
      double syst_FracLowIn[6] = {errorsLoIn->GetBinContent(startbinerror),  errorsLoIn->GetBinContent(startbinerror+1),
                                 errorsLoIn->GetBinContent(startbinerror+2), errorsLoIn->GetBinContent(startbinerror+3),
                                 errorsLoIn->GetBinContent(startbinerror+4), errorsLoIn->GetBinContent(startbinerror+5)};
      double syst_FracLowOut[6] = {errorsLoOut->GetBinContent(startbinerror),   errorsLoOut->GetBinContent(startbinerror+1),
                                   errorsLoOut->GetBinContent(startbinerror+2), errorsLoOut->GetBinContent(startbinerror+3),
                                   errorsLoOut->GetBinContent(startbinerror+4), errorsLoOut->GetBinContent(startbinerror+5)};
      double syst_FracHighIn[6] = {errorsHighIn->GetBinContent(startbinerror),  errorsHighIn->GetBinContent(startbinerror+1),
                                   errorsHighIn->GetBinContent(startbinerror+2), errorsHighIn->GetBinContent(startbinerror+3),
                                   errorsHighIn->GetBinContent(startbinerror+4), errorsHighIn->GetBinContent(startbinerror+5)};
      double syst_FracHighOut[6] = {errorsHighOut->GetBinContent(startbinerror), errorsHighOut->GetBinContent(startbinerror+1),
                                    errorsHighOut->GetBinContent(startbinerror+2), errorsHighOut->GetBinContent(startbinerror+3),
                                    errorsHighOut->GetBinContent(startbinerror+4), errorsHighOut->GetBinContent(startbinerror+5)};
      //get absolute errors
      double syslo[6] = {yieldLow->GetBinContent(startbin)*syst_FracLow[0], yieldLow->GetBinContent(startbin+1)*syst_FracLow[1], 
                        yieldLow->GetBinContent(startbin+2)*syst_FracLow[2], yieldLow->GetBinContent(startbin+3)*syst_FracLow[3], 
                        yieldLow->GetBinContent(startbin+4)*syst_FracLow[4], yieldLow->GetBinContent(startbin+5)*syst_FracLow[5]};
      double syshi[6] = {yieldHigh->GetBinContent(startbin)*syst_FracHigh[0], yieldHigh->GetBinContent(startbin+1)*syst_FracHigh[1],
                        yieldHigh->GetBinContent(startbin+2)*syst_FracHigh[2], yieldHigh->GetBinContent(startbin+3)*syst_FracHigh[3],
                        yieldHigh->GetBinContent(startbin+4)*syst_FracHigh[4], yieldHigh->GetBinContent(startbin+5)*syst_FracHigh[5]};
     double sysloout[6] = {yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[0], yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[1],
                           yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[2], yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[3], 
                           yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[4], yieldLowOutPlane->GetBinContent(startbin)*syst_FracLowOut[5]};
     double syshiout[6] = {yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[0], yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[1],
                           yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[2], yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[3],
                           yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[4], yieldHighOutPlane->GetBinContent(startbin)*syst_FracHighOut[5]};
     double sysloin[6] = {yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[0], yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[1],
                          yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[2], yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[3],
                          yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[4], yieldLowInPlane->GetBinContent(startbin)*syst_FracLowIn[5]};
     double syshiin[6] = {yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[0], yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[1],
                          yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[2], yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[3],
                          yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[4], yieldHighInPlane->GetBinContent(startbin)*syst_FracHighIn[5]};



     //ratio errors
     double syslooutin[6] = {lowvalues3050[0]*1.159, lowvalues3050[1]*0.185797, lowvalues3050[2]*0.106635, lowvalues3050[3]*0.089712, lowvalues3050[4]*0.110993, lowvalues3050[5]*0.0845498};
     double syshioutin[6] = {highvalues3050[0]*0.722617,  highvalues3050[1]*0.164787, highvalues3050[2]*0.0986985, highvalues3050[3]*0.0689082, highvalues3050[4]*0.0588923, highvalues3050[5]*0.0695121};    
 

  double sysraerrors[6] = {1.23456, 0.183724, 0.0961714, 0.0676367, 0.0438335, 0.0414292};

     double sysra[6] = {yieldratio->GetBinContent(startbin)*sysraerrors[0], yieldratio->GetBinContent(startbin+1)*sysraerrors[1], yieldratio->GetBinContent(startbin+2)*sysraerrors[2], yieldratio->GetBinContent(startbin+3)*sysraerrors[3], yieldratio->GetBinContent(startbin+4)*sysraerrors[4], yieldratio->GetBinContent(startbin+5)*sysraerrors[5]};
  

      //ctraye TGraph Errors for systematic errors
      auto systlo = new TGraphErrors(6, bincenters3050, nomlo, nul3050, syslo);
      //auto systhi = new TGraphErrors(5, bincenters, binvalues, nul, sys3050);
      systlo->SetMarkerSize(0.6);
      systlo->SetMarkerStyle(20);
      systlo->SetFillColorAlpha(kCaitieBlue, 0.15);
      systlo->SetMarkerColor(kCaitieDarkBlue);
      systlo->SetLineColor(kCaitieDarkBlue);
      //syst3050->SetLineColor(kBlack);
      //syst3050->SetFillStyle(0);   

      auto systhi = new TGraphErrors(6, bincenters3050, nomhi, nul3050, syshi);
      systhi->SetMarkerSize(0.6);
      systhi->SetMarkerStyle(20);
      systhi->SetFillColorAlpha(kCaitiePink, 0.2);
      systhi->SetMarkerColor(kViolet);
      systhi->SetLineColor(kViolet);


      auto systra = new TGraphErrors(6, bincenters3050, nomrat, nul3050, sysra);
      systra->SetMarkerSize(0.6);
      systra->SetMarkerStyle(20);
      systra->SetFillColorAlpha(kCaitieGreen, 0.5);
      systra->SetMarkerColor(kCaitieDarkGreen);
      systra->SetLineColor(kCaitieDarkGreen);
 
      auto systlooutin = new TGraphErrors(6, bincenters3050, lowvalues3050, nul3050, syslooutin);
      systlooutin->SetMarkerSize(0.6);
      systlooutin->SetMarkerStyle(20);
      systlooutin->SetFillColorAlpha(kCaitieBlue, 0.15);
      systlooutin->SetMarkerColor(kCaitieDarkBlue);
      systlooutin->SetLineColor(kCaitieDarkBlue);

    
 
      auto systhioutin = new TGraphErrors(6, bincenters3050, highvalues3050, nul3050, syshioutin);
      systhioutin->SetMarkerSize(0.6);
      systhioutin->SetMarkerStyle(20);
      systhioutin->SetFillColorAlpha(kCaitiePink, 0.2);    //65
      systhioutin->SetMarkerColor(kViolet);           //65
      systhioutin->SetLineColor(kViolet);           //65


/*
      //auto caitie = new TGraphErrors(3, bincenters010, binvalues010, nul010, sys010);
      auto caitie = new TGraphErrors(5, bincenters, binvalues, nul, sys3050);
      caitie->SetMarkerSize(0.6);
      caitie->SetMarkerStyle(20);
      if (cent == true) {
      caitie->SetFillColorAlpha(centcolor, 0.3);
      caitie->SetMarkerColor(centcolor);
      caitie->SetLineColor(centcolor);}
      if (semi == true) {
      caitie->SetFillColorAlpha(kCyan-3, 0.3);
      caitie->SetMarkerColor(kCyan+3);
      caitie->SetLineColor(kCyan+3);}

*/



  //====================================
  //========== create lines ============
  //====================================
  TLine *line = new TLine(lbound, 1, rbound, 1);
  line->SetLineStyle(2);
  TLine *sline = new TLine(40.0, 1, 100.0, 1);
  sline->SetLineStyle(2);
  TLine *mline = new TLine(30.0, 1, 100.0, 1);
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
   auto gen = new TLegend(0.2, 0.7, 0.3, 0.85);
        gen->SetTextSize(0.035);
        gen->SetBorderSize(0);
        if (cent == true)    gen->AddEntry((TObject*)0, Form("ALICE Pb#font[122]{-}Pb, 5.02 TeV, %.0f#font[122]{-}%.0f%%, LHC18q_pass3", centl, centr), "");
        if (semi == true)    gen->AddEntry((TObject*)0, Form("ALICE Pb#font[122]{-}Pb, 5.02 TeV, %.0f#font[122]{-}%.0f%%, LHC18q_pass3", centl, centr), "");
        gen->AddEntry((TObject*)0, "Charged Jets, anti-#it{k}_{T}, R = 0.2, |#it{#eta}_{jet}| < 0.7", "");
        gen->AddEntry((TObject*)0, "Work in Progress", "");
   auto ali = new TLegend(0.15, 0.6, 0.3, 0.85);
        ali->SetTextSize(0.05);
        ali->SetBorderSize(0);
        ali->AddEntry((TObject*)0, Form("ALICE Pb#font[122]{-}Pb, 5.02 TeV, %.0f#font[122]{-}%.0f%%, LHC18q_pass3", centl, centr), "");
        ali->AddEntry((TObject*)0, "Charged Jets, anti-#it{k}_{T}, R = 0.2, #it{p}_{T}^{lead track} > 5 GeV/#it{c}, |#it{#eta}_{jet}| < 0.7", "");
        ali->AddEntry((TObject*)0, "Work in Progress", "");
   auto datadetails = new TLegend(0.55, 0.25, 0.75, 0.4);
        datadetails->SetTextSize(0.035);
        datadetails->SetBorderSize(0);
        datadetails->AddEntry(ratioOIHigh, "30\% Highest q_{2}", "pl");
        datadetails->AddEntry(ratioOILow, "30 \% Lowest q_{2}", "pl");

   auto specleg = new TLegend(0.6, 0.4, 0.85, 0.525);
        specleg->SetTextSize(0.055);
        specleg->SetBorderSize(0);
        specleg->AddEntry(stathi, "30\% Highest q_{2}", "plf");
        specleg->AddEntry(statlo, "30\% Lowest q_{2}", "plf");  

   auto yieldLeg = new TLegend(0.6, 0.6, 0.85, 0.85);
        yieldLeg->SetTextSize(0.04);
        yieldLeg->SetBorderSize(0);
        yieldLeg->AddEntry(yieldHighInPlane, "High q_{2}, in-plane", "pl");
        yieldLeg->AddEntry(yieldHighOutPlane, "High q_{2}, out-of-plane", "pl");
        yieldLeg->AddEntry(yieldLowInPlane, "Low q_{2}, in-plane", "pl");
        yieldLeg->AddEntry(yieldLowOutPlane, "Low q_{2}, out-of-plane", "pl");


  //====================================
  //========== draw histograms =========
  //====================================
  //In-Plane/Out-of-Plane
  TCanvas *c1 = new TCanvas("QM Plot 1", lofile, 600, 500);
  c1->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad1->Draw();
          pad1->cd();               // pad1 becomes the current pad   
         //pad1->SetLogx();
         //pad1->SetLogy();
         //unfoldhistsl[0]->SetMaximum(pow(10,-1));
         //unfoldhistsl[0]->SetTitle("");
         //for (int i = 0; i < 15; i++)         unfoldhistsl[i]->Draw("same");
             ratioOILow->SetTitle("");
             ratioOILow->SetMinimum(0.01);
             ratioOILow->SetMaximum(1.99);
             ratioOILow->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             ratioOILow->GetXaxis()->SetTitleSize(0.035);
             ratioOILow->GetYaxis()->SetTitleSize(0.04);
             ratioOILow->GetXaxis()->SetRangeUser(30.0, 100.0);
             ratioOILow->SetYTitle("Ratio (out-of-plane/in-plane)");
         //yieldratio->Draw("same");
         ratioOILow->Draw("same");
         ratioOIHigh->Draw("same");
         systlooutin->Draw("ezp 2 same");
         systhioutin->Draw("ezp 2 same");
         ratioOILow->Draw("same");
         ratioOIHigh->Draw("same");
    //stat->Draw("ez same");
     //caitie->Draw("ezp 2 same");
   gen->Draw("same");
   datadetails->Draw("same");
   sline->Draw("same");





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
       setcoloralpha(yieldLow, kBlack, 0.0);
       setcoloralpha(yieldHigh, kBlack, 0.0);
       yieldLow->SetMarkerColor(kCaitieBlue);
       yieldHigh->SetMarkerColor(kViolet);
       yieldLow->SetMaximum(yieldLow->GetBinContent(2)*50);
       yieldLow->SetTitle("");
       yieldLow->GetXaxis()->SetRangeUser(30.0, 100.0);
       yieldLow->GetYaxis()->SetLabelSize(0.05);
       yieldLow->Draw("same");
       yieldHigh->Draw("same");
       systlo->Draw("ezp 2 same");
       systhi->Draw("ezp 2 same");
       statlo->Draw("ez same");
       stathi->Draw("ez same");   
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
       yieldratio->GetXaxis()->SetRangeUser(30.0, 100.0);
       yieldratio->Draw("same");
       statra->Draw("ez same");
       systra->Draw("ezp 2 same");
       yieldratio->Draw("same");
       statra->Draw("ez same");
       sline->Draw("same"); 

 

  //Yields
  TCanvas *c3 = new TCanvas("Spectra", "Spectra", 600, 500);
  c3->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad3 = new TPad("pad2a", "", 0.0, 0.05, 1.0, 1.0);
          pad3->SetBottomMargin(0.1); // Upper and lower plot are joined
          pad3->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad3->Draw();
          pad3->cd();
          pad3->SetLogy();
       yieldHighInPlane->GetXaxis()->SetRangeUser(30.0, 100.0);
       yieldHighInPlane->SetTitle("");
       yieldHighInPlane->SetXTitle("#it{p}_{T} (GeV/#it{c})");
       yieldHighInPlane->SetMaximum(yieldHighInPlane->GetBinContent(1)*100);
       yieldHighInPlane->SetMinimum(yieldHighInPlane->GetBinContent(1)*pow(10,-3));
       yieldHighInPlane->Draw("same");
       yieldHighOutPlane->Draw("same");
       yieldLowInPlane->Draw("same");
       yieldLowOutPlane->Draw("same");
       yieldLeg->Draw("same");


   cout << "Draw Hist 2\n";
   cout << "	Low q2 Errors:\n";
   cout << "	30-40 GeV: " << statlo3050[0] << "\n";  
   cout << "	40-50 GeV: " << statlo3050[1] << "\n";  
   cout << "	50-60 GeV: " << statlo3050[2] << "\n";  
   cout << "	60-70 GeV: " << statlo3050[3] << "\n";  


 }


void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}

void setcoloralpha(TH1F* h, int kcolor, double trans) {
  h->SetLineColorAlpha(kcolor, trans);
  h->SetMarkerColorAlpha(kcolor, trans);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}
 
