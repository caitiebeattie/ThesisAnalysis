//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

/*
  RooSimplepTPbPb_split.cxx : Script to perform the split MC test on the area based method.
  Hannah Bossi <hannah.bossi@yale.edu>
  3/4/2020
 */


#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"

#include "TFile.h"
#include "TVectorD.h"

#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TPostScript.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
#include "TNtuple.h"
#include "TProfile.h"
#inlcude "TRandom3.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"
//#include "RooUnfoldTestHarness2D.h"
#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;
double pi = 3.1415926;
//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================



 long int preeta = 0;
 long int posteta = 0;
 long double totevents = 0.0;





TH2D* CorrelationHist (const TMatrixD& cov,const char* name, const char* title,
		       Double_t lo, Double_t hi,Double_t lon,Double_t hin)
{
  int nb = cov.GetNrows();
  int na = cov.GetNcols();
  cout<<nb<<" "<<na<<endl;
  TH2D* h= new TH2D (name, title, nb, 0, nb, na, 0, na);
  h->SetAxisRange (-1.0, 1.0, "Z");
  for(int i=0; i < na; i++)
  for(int j=0; j < nb; j++) {
  Double_t Viijj = cov(i,i)*cov(j,j);
      if (Viijj>0.0) h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
    }
  return h;
}

//==============================================================================
// Example Unfolding
//==============================================================================



void RooSimplepTPbPb_differential_Caitie(double centl, double centr, int ese = 2, int R = 2, TString cFiles2="files15.txt")
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  int difference=1;
  int Ppol=0;
  cout << "==================================== pick up the response matrix (ese) for background==========================" << endl;
  ///////////////////parameter setting
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  
  TRandom3* rand = new TRandom3();

   //***************************************************
 

  double truncation = 30.0;

  std::vector<double> kBinsUnfolded = {10.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0, 140.0, 190.0/*, 250.0*/};
  std::vector<double> kBinsMeasured = {truncation, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0};   //semi-central
  //std::vector<double> kBinsMeasured = {truncation, 50.0, 60, 70, 85, 100, 120};   //central
  

  vector<const char*> histnames = {"cent1", "cent2", "cent3", "cent4", "cent5", "cent6", "cent7", "cent8", "cent9", "cent10", "cent11", "cent12", "cent13", "cent14", "cent15", "cent16", "cent17", "cent18", "cent19", "cent20"};
  vector<const char*> histhi = {"hi1", "hi2", "hi3", "hi4", "hi5", "hi6", "hi7", "hi8", "hi9", "hi10", "hi11", "hi12", "hi13", "hi14", "hi15", "hi16", "hi17", "hi18", "hi19", "hi20"};
  vector<const char*> histlo = {"lo1", "lo2", "lo3", "lo4", "lo5", "lo6", "lo7", "lo8", "lo9", "lo10", "lo11", "lo12", "lo13", "lo14", "lo15", "lo16", "lo17", "lo18", "lo19", "lo20"};
  vector<TH1D*> hists(0);
  vector<TH1D*> histlos(0);
  

  for (int i = 0; i < 20; i++)  {
      hists.push_back(new TH1D(histnames[i], histnames[i], 4, 20.0, 100.0));
      histlos.push_back(new TH1D(histlo[i], histlo[i], 4, 20.0, 100.0));}
      

  //the raw correlation (data or psuedodata)
  TH1D *h1raw(0);
  h1raw = new TH1D("raw", "raw", kBinsMeasured.size()-1, kBinsMeasured.data());
  //h1raw = new TH1D("raw", "raw", 95.0, 25.0, 120.0);
  //TH1D *h1rawfull(0);
  //h1rawfull = new TH1D("rawfull", "rawfull", kBinsMeasured.size()-1, kBinsMeasured.data());

  //detector measure level (reco or raw MC)
  //TH1F *h1smeared(0);
  //h1smeared = new TH1F("smeared", "smeared",kBinsMeasured.size()-1, kBinsMeasured.data());
  //h1smeared = new TH1F("smeared", "smeared", 95.0, 25.0, 120.0);
  //true correlations with measured cuts
  TH1F *h1true(0);
  h1true = new TH1F("true", "true", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  //h1true = new TH1F("true", "true", 240.0, 10.0, 250.0);
  //TH1F *h1truefull(0);
  //h1truefull = new TH1F("truefull", "truefull", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  //full true correlation
  TH1F *h1fulleff(0);
  h1fulleff=new TH1F("truef", "truef", kBinsUnfolded.size()-1, kBinsUnfolded.data()); 
  //h1fulleff=new TH1F("truef", "truef", 240.0, 10.0, 250.0); 
  TH1F *mr = new TH1F("mr", "mr", 30, 0.0, 0.3);
  //Kinematic efficiency
  TH1F *hKinPre = new TH1F("KinPre", "KinPre", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  TH1F *hKinPos = new TH1F("KinPos", "KinPos", kBinsUnfolded.size()-1, kBinsUnfolded.data());
 

 
  //============================================
  //==========   Get Raw Spectra    ============
  //============================================
   //Raw Spectra from Trains
   //string rootfilename = "../Spectra/AnalysisResults7283.root";   // semi-central
    string rootfilename = "AnalysisResults7812.root";             // R=0.2, 30-50%, pass1
    const char *lerootfile = rootfilename.c_str(); 

    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);
    TTree *T = nullptr;
    if (R == 2) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
    if (R == 4) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_Jet", T);
    float jetpT, centra, qvec, ep;
    T->SetBranchAddress("Jet_Pt", &jetpT);
    T->SetBranchAddress("Event_Centrality", &centra);
    T->SetBranchAddress("Event_Q2Vector", &qvec);    
    T->SetBranchAddress("Jet_EPangle", &ep);    

      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2V0C_EPV0C");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                             //number percentile we're on
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      double centdiff = (centr-centl)/2.0;
      int leftc, rightc;
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < (double)centdiff; j++) {
        leftc = (int)centl + 2*j; 
        rightc = (int)centl + 2*j + 2;    
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);       //restrict 2D hist to current 1% centrality bin
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
        }
      }
  TH1D *perc20 = new TH1D("perc20", "perc20", centdiff, centl, centr);
  TH1D *perc50 = new TH1D("perc50", "perc50", centdiff, centl, centr);
  TH1D *perc80 = new TH1D("perc80", "perc80", centdiff, centl, centr);
  for (int i = 1; i <= centdiff; i++) {
      perc20->SetBinContent(i, allPercent[i-1][2]);
      perc50->SetBinContent(i, allPercent[i-1][4]);
      perc80->SetBinContent(i, allPercent[i-1][6]);
  }
  perc20->SetLineColor(kRed);
    perc20->SetLineWidth(2);
  perc80->SetLineColor(kRed);
    perc80->SetLineWidth(2);
  perc50->SetLineColor(kBlack);
    perc50->SetLineWidth(2);



    //Determine Extractor Bin Scaling
    const int extraBins = 9;
    double extraFrac[extraBins] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double extraEdge[extraBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0};
    
    //const int pbins = 18;
    //double pTedgepp[pbins + 1]  = {25, 30, 40, 50, 60, 70, 85, 100};
    TH1D* specraw = new TH1D("specraw", "specraw", kBinsMeasured.size()-1, kBinsMeasured.data());
    int entries = T->GetEntries();

    vector<long int>  centCheck(0);
    vector<long int>  centOG(0);
    for (int i = 0; i < 20; i++) {
            centCheck.push_back(0);
            centOG.push_back(0);}
    for (int i = 0; i < entries; i++)   {
       T->GetEntry(i);     
       if (centra < centl || centra > centr)    continue;  
       if (jetpT < truncation || jetpT > 120.0)       continue;
       centOG[(int)centra-30]++;     
       if (ese == 0)  {if (qvec > allPercent[(int)centra-30][1]) continue;} //20th percentile
       if (ese == 1)  {if (qvec < allPercent[(int)centra-30][7]) continue;} //80th percentile
       centCheck[(int)centra-30]++;
       //if (abs(cos(ep)) < sqrt(2)/2.0) continue;  //in-plane
       //if (abs(cos(ep)) > sqrt(2)/2.0) continue;  //out-plane
       int scaleindex = -1;
         if (jetpT >= extraEdge[0] && jetpT < extraEdge[1])  scaleindex = 0;
         if (jetpT >= extraEdge[1] && jetpT < extraEdge[2])  scaleindex = 1;
         if (jetpT >= extraEdge[2] && jetpT < extraEdge[3])  scaleindex = 2;
         if (jetpT >= extraEdge[3] && jetpT < extraEdge[4])  scaleindex = 3;
         if (jetpT >= extraEdge[4] && jetpT < extraEdge[5])  scaleindex = 4;
         if (jetpT >= extraEdge[5] && jetpT < extraEdge[6])  scaleindex = 5;
         if (jetpT >= extraEdge[6] && jetpT < extraEdge[7])  scaleindex = 6;
         if (jetpT >= extraEdge[7] && jetpT < extraEdge[8])  scaleindex = 7;
           if (scaleindex == 0)   continue;      
       specraw->Fill(jetpT, 1.0/extraFrac[scaleindex]);
       hists[(int)centra-30]->Fill(jetpT, 1.0/extraFrac[scaleindex]);

    }

    for (int i = 0; i < 20; i++)    cout <<"centPass/Total(" << i << "): " << centCheck[i] << "/" << centOG[i] << " = " << (double)centCheck[i]/(double)centOG[i] <<"\n";


     //for (int i = 0; i < 20; i++)  hists[i]->Divide(histlos[i]);

     //TH1D* specraw = (TH1D*)lespec->Rebin(kBinsMeasured.size()-1, "specraw", kBinsMeasured.data());

  int tot = 0;
  int small = 0;
  int large = 0;
      int largepT = 0;
  
  TH2F *hcovariance(0);
  hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);

 
   //h1smeared->Sumw2();
   h1true->Sumw2();
   h1raw->Sumw2();
   specraw->Sumw2();




  //============================================
  //=========   Fill Response Matrix     =======
  //============================================

   //branches in the tree that you need in this analysis
   float ptJet,ptPar, ptdet, centrality, matchradius, area, eta, q2, epangle;
  
   int nEv=0;; 
 
   ifstream infile2;
   infile2.open(cFiles2.Data());
   char filename2[300];
   //overall
   //RooUnfoldResponse responsefull;
   //responsefull.Setup(h1rawfull, h1truefull);
   //split test
   RooUnfoldResponse responsesplit;
   responsesplit.Setup(h1raw, h1true);
 
 
  

    while(infile2 >> filename2){
      int pthardbin=0;
      cout << "Filename: " << filename2 <<"\n";
      TFile *input=TFile::Open(filename2);
      TList *list=(TList*) input->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");
      TList *list2=(TList*) list->FindObject("EventCuts");
      //TH1D *hcent=(TH1D*)list2->FindObject("Centrality_selected");
      TProfile *hcross=(TProfile*)list->FindObject("fHistXsection");
      TH1D *htrials=(TH1D*)list->FindObject("fHistTrials");
      TH1D *hpthard=(TH1D*)list->FindObject("fHistPtHard");
      TH1D *hnevent=(TH1D*)list->FindObject("fHistEventCount");
      
    //determine pT hard scaling
    for(int i=1;i<=htrials->GetNbinsX();i++)   { if(htrials->GetBinContent(i)!=0) pthardbin = i;}
      double n_event = hnevent->GetBinContent(1);
        totevents += n_event;
      double xsect = hcross->GetBinContent(htrials->GetMean() + 1);
      double xscale = xsect*hcross->GetEntries();
      double trials = htrials->GetBinContent(htrials->GetMean() + 1);
      double pTHardscalefactor=xscale/(n_event*trials);  

      std::cout << "pT Hard Bin: " << pthardbin-1 << ", Scale Factor:  " << pTHardscalefactor << std::endl;
    TFile *input2=TFile::Open(filename2);
    // get the mc tree
    TTree *mc = nullptr;
    if (R == 2) mc=(TTree*)input2->Get("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_hybJet"); 
    if (R == 4) mc=(TTree*)input2->Get("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet"); 
    int nEv=mc->GetEntries(); 
    cout << "Tree Obtained\n";

    int numTracks;
    // set the branches of the mc tree
    mc->SetBranchAddress("Jet_Pt", &ptJet); 
    mc->SetBranchAddress("Jet_Eta", &eta);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &ptPar);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptdet);
    mc->SetBranchAddress("Jet_NumTracks", &numTracks);
    float jetTrackPt[400];
    //mc->SetBranchAddress("Jet_Track_Pt", &jetTrackPt);
    mc->SetBranchAddress("Event_Centrality", &centrality);
    mc->SetBranchAddress("Jet_Area", &area);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius); 
    mc->SetBranchAddress("Event_Q2Vector", &q2);
    mc->SetBranchAddress("Jet_EPangle", &epangle);
    
  
  int countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
      mc->GetEntry(iEntry);
      //****** Get the scaling factor
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      if(ptJet >= -20. && ptJet < 10.)        EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)    EBscale = 1.0/0.01;
      else if(ptJet >= 20. && ptJet < 30.)    EBscale = 1.0/0.03;
      else if(ptJet >= 30. && ptJet < 40.)    EBscale = 1.0/0.15;
      else if(ptJet >= 40. && ptJet < 60.)    EBscale = 1.0/0.20;
      else if(ptJet >= 60. && ptJet < 80.)    EBscale = 1.0/0.20;
      else if(ptJet >= 80. && ptJet < 100.)   EBscale = 1.0/0.20;
      else if(ptJet >= 100. && ptJet < 140.)  EBscale = 1.0/0.15;
      else if(ptJet >= 140. && ptJet < 200.)  EBscale = 1.0/0.5;
      scalefactor*=EBscale;
      // *********

      //if (abs(cos(epangle)) < sqrt(2)/2.0) continue;  //in-plane
      //if (abs(cos(epangle)) > sqrt(2)/2.0) continue;  //out-plane
      if (centrality < centl || centrality > centr)  continue;
      if (ese == 0)  {if (q2 > allPercent[(int)centrality-30][1]) continue;} //20th percentile
      if (ese == 1)  {if (q2 < allPercent[(int)centrality-30][7]) continue;} //80th percentile

      if (abs(eta) > (0.9-(double)R*0.1))  continue;  //TPC fiducial cut
      if(ptPar > 200.) continue;    
  
      
      hKinPre->Fill(ptPar, scalefactor);
      if (ptJet < truncation || ptJet > 120) continue;
      hKinPos->Fill(ptPar, scalefactor);
      if (ptPar < 10.0 || ptJet < 10.0 )   continue; 
      tot++;

	  h1true->Fill(ptPar, scalefactor);
	  h1raw->Fill(ptJet, scalefactor);
	  responsesplit.Fill(ptJet, ptPar, scalefactor);
    }}


    
    TH1F *htrueptd = (TH1F*) h1fulleff->Clone("trueptd");
    TH1F *htruept  = (TH1F*) h1fulleff->Clone("truept"); 


    long double avgevents = (long double)totevents/20.0;
    h1raw->Scale(avgevents);
    h1true->Scale(avgevents);
    //////////efficiencies done////////////////////////////////////


    stringstream ww;
    if (ese == 0)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_lo20_Dec6.root"; 
    if (ese == 1)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_hi20_Dec6.root"; 
    if (ese == 2)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_allq_Dec6.root"; 
    TFile *fout=new TFile (ww.str().c_str(),"RECREATE");
    fout->cd();
    h1raw->SetName("raw");
    h1raw->Write();
    //h1rawfull->Write();
    //h1smeared->SetName("smeared");
    //h1smeared->Write();
    htrueptd->Write();
    h1true->SetName("true");
    h1true->Write();
    specraw->SetName("specraw");
    specraw->Write();
    hKinPre->Write();
    hKinPos->Write();
    mr->Write(); 
   //h1truefull->Write();
    responsesplit.Write();
    //projection of the responsesplit onto the truth axis
    TH1D* htruth = (TH1D*)responsesplit.Htruth();
    htruth->SetName("htruth");
    htruth->Write();


   cout << "Preeta: " << preeta <<"\n";
   cout << "Posteta: " << posteta <<"\n";

    for(int jar = 1; jar < 16; jar++) {
      int iter = jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&responsesplit, specraw, iter, false);    // OR
      //RooUnfoldBayes unfold(&responsefull, h1rawfull, iter, false);    // OR
      TH1D* hunf = (TH1D*)unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH1* hfold = responsesplit.ApplyToTruth(hunf, "");

      TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
      htempUnf->SetName(Form("Bayesian_Unfolded_%d",iter));
      
      TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
      htempFold->SetName(Form("Bayesian_Folded_%d",iter));        

      htempUnf->Write();
      htempFold->Write();
	  
    }
	  
  TCanvas *c1 = new TCanvas("Q2 vs Cent", "Q2 vs Cent", 500, 500);
  c1->cd();
       gStyle->SetOptStat(0);
       //gStyle->SetPadTickX(1);
       //gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad1->SetRightMargin(0.15);
          pad1->Draw();
          pad1->cd();
          pad1->SetLogz();
       q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
       q2hist2D->GetXaxis()->SetTitle("Centrality (%)");
       q2hist2D->GetYaxis()->SetTitle("q_{2}^{V0C}");
       q2hist2D->Draw("same colz");
       perc20->Draw("same");
       perc80->Draw("same");
       //perc50->Draw("same");


  TCanvas *c2 = new TCanvas("centbinsp", "centbinsp", 1000, 500);
  c2->cd();
     gStyle->SetOptStat(0);
     TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
        pad2->Draw();
        pad2->cd();
     for (int i = 0; i < 20; i++)   {
         cout << "hist[i]->GetEntries(): " << hists[i]->GetEntries() <<"\n";
         hists[i]->SetLineColor(30+i); 
         hists[i]->Draw("same");
         }
}
#ifndef __CINT__
//int main () { RooSimplepTPbPb_data_Caitie(); return 0; }  // Main program when run stand-alone
#endif
