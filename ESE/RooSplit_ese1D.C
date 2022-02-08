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
  RooSimplepTPbPb_split.cxx : Script to perform 2D unfolding.
  Caitie Beattie <caitie.beattie@yale.edu>
  1/13/2022
 
  Class documentation can be found at https://hepunx.rl.ac.uk/~adye/software/unfold/htmldoc/RooUnfold.html
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

 long int preeta = 0;
 long int posteta = 0;
 long double totevents = 0.0;




//==============================================================================
// Unfolding
//==============================================================================

void RooSplit_ese1D(double centl, double centr, int ese, int R = 2, int epcut = 2, TString cFiles2="files14.txt")
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
 

  //========================================================================
  //=================   Determine Unfolding Parameters  ====================   
  //========================================================================
  
  //establish centrality trigger
  bool cent = false;
  bool semi = false;
  if (centl == 0.0  && centr == 10.0)   cent = true;
  if (centl == 30.0 && centr == 50.0)   semi = true;

  //establish unfolding range
  double truncation;
  if (cent == true)   truncation = 50.0;
  if (semi == true)   truncation = 30.0;




  //=============================================================================
  //=================   Initialize Histograms and Responses  ====================   
  //=============================================================================
  std::vector<double> kBinsUnfolded = {10.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0, 140.0, 190.0/*, 250.0*/};
  std::vector<double> kBinsMeasured = {truncation, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0};   //semi-central
  //std::vector<double> kBinsMeasured = {truncation, 50.0, 60, 70, 85, 100, 120};   //central
  std::vector<double> kBinsEPangle = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  //raw spectra from embedding (pseudo-data)
  TH1D *h1rawresponse = new TH1D("rawresponse", "rawresponse", kBinsMeasured.size()-1, kBinsMeasured.data());
  //true spectra from embedding
  TH1F *h1true = new TH1F("true", "true", kBinsUnfolded.size()-1, kBinsUnfolded.data());

  //raw spectra from data
  TH1D* h1rawdata = new TH1D("h1rawdata", "h1rawdata", kBinsMeasured.size()-1, kBinsMeasured.data());

  //Kinematic efficiency
  TH1F *hKinPre = new TH1F("KinPre", "KinPre", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  TH1F *hKinPos = new TH1F("KinPos", "KinPos", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  
  //Responses
  RooUnfoldResponse response1D;
  response1D.Setup(h1rawresponse, h1true);



  //=============================================================================
  //==========================   Get Raw Spectra  ===============================   
  //=============================================================================
   
    //Raw Spectra from Trains
    string rootfilename;
    if (semi == true)   rootfilename = "AnalysisResults7784.root";
    const char *lerootfile = rootfilename.c_str(); 

    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);      

      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2V0C_EPV0M");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrFullV0");
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
        }
      }



    //Determine Extractor Bin Scaling
    const int extraBins = 9;
    double extraFrac[extraBins] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double extraEdge[extraBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0};
    
  
  TH2F *hcovariance(0);
  hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);

 


  //==================================================================================
  //==========================   Get Embedded Spectra  ===============================   
  //==================================================================================

   //branches in the tree that you need in this analysis
   float ptJet,ptPar, ptdet, centrality, matchradius, area, eta, q2, epangle, epPar;
  
   int nEv=0;; 
 
   ifstream infile2;
   infile2.open(cFiles2.Data());
   char filename2[300];
   //RooUnfoldResponse responsefull;
   //responsefull.Setup(h1rawfull, h1truefull);
  

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

    // set the branches of the mc tree
    mc->SetBranchAddress("Jet_Pt", &ptJet); 
    mc->SetBranchAddress("Jet_Eta", &eta);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &ptPar);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptdet);
    mc->SetBranchAddress("Event_Centrality", &centrality);
    mc->SetBranchAddress("Jet_Area", &area);
    mc->SetBranchAddress("Event_Q2Vector", &q2);
    mc->SetBranchAddress("Jet_EPangle", &epangle);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_EPangle", &epPar);
    
  
  int countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
      mc->GetEntry(iEntry);
      //****** Get the scaling factor
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      if(ptJet >= -20. && ptJet < 10.)         EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)     EBscale = 1.0/0.01;
      else if(ptJet >= 20. && ptJet < 30.)     EBscale = 1.0/0.03;
      else if(ptJet >= 30. && ptJet < 40.)     EBscale = 1.0/0.15;
      else if(ptJet >= 40. && ptJet < 60.)     EBscale = 1.0/0.20;
      else if(ptJet >= 60. && ptJet < 80.)     EBscale = 1.0/0.20;
      else if(ptJet >= 80. && ptJet < 100.)    EBscale = 1.0/0.20;
      else if(ptJet >= 100. && ptJet < 140.)   EBscale = 1.0/0.12;
      else if(ptJet >= 140. && ptJet < 200.)   EBscale = 1.0/0.02;
      scalefactor*=EBscale;
      // *********

      if (centrality < centl || centrality > centr)  continue;
      if (epcut == 0)    {if (abs(cos(epangle)) > sqrt(2)/2.0) continue;}  //out-plane
      if (epcut == 1)    {if (abs(cos(epangle)) < sqrt(2)/2.0) continue;}  //in-plane
      if (ese == 0)           {if (q2 > allPercent[int(centrality-centl)][1])    continue;}  //20th percentile
      if (ese == 1)           {if (q2 < allPercent[int(centrality-centl)][7])    continue;}  //80th percentile

      preeta++;      
      if (abs(eta) > (0.9-(double)R*0.1))  continue;  //TPC fiducial cut
      posteta++;
      if(ptPar > 200.) continue;    
  
      
      hKinPre->Fill(ptPar, scalefactor);
      if (ptJet < truncation || ptJet > 120) continue;
      hKinPos->Fill(ptPar, scalefactor);
      if (ptPar < 10.0 || ptJet < 10.0 )   continue; 

      double split = rand->Rndm();
      if (split < 0.7)	{
	  h1rawresponse->Fill(ptJet, scalefactor);
	  response1D.Fill(ptJet, ptPar, scalefactor);
      }
      else  {
          h1rawdata->Fill(ptJet, scalefactor);
          h1true->Fill(ptPar, scalefactor);
      }
    }
    }

    long double avgevents = (long double)totevents/20.0;
    h1rawresponse->Scale(avgevents);
    h1rawdata->Scale(avgevents);
    h1true->Scale(avgevents);
    //////////efficiencies done////////////////////////////////////



  //==================================================================================
  //==========================   Setup the Outputs   =================================   
  //==================================================================================

    stringstream ww;
    if (ese == 0)  ww << "UnfoldingSplit_1D_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_lo20_Jan21.root"; 
    if (ese == 1)  ww << "UnfoldingSplit_1D_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_hi20_Jan21.root"; 
    if (ese == 2)  ww << "UnfoldingSplit_1D_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_allq_Jan21.root"; 
    TFile *fout=new TFile (ww.str().c_str(),"RECREATE");
    fout->cd();
    h1rawresponse->SetName("rawresponse1");
    h1rawresponse->Write();
    h1rawdata->SetName("rawdata1");
    h1rawdata->Write();
    h1true->SetName("true1");
    h1true->Write();
    hKinPre->Write();
    hKinPos->Write();
    response1D.Write();
    //projection of the response1D onto the truth axis
    TH1D* htruth = (TH1D*)response1D.Htruth();
    htruth->SetName("htruth");
    htruth->Write();



  //==================================================================================
  //==========================   Do the Unfolding ====================================   
  //==================================================================================

    for(int jar = 1; jar < 16; jar++) {
      int iter = jar;
      cout<< "Iteration: " << iter << endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&response1D, h1rawdata, iter, false);    // OR
      TH1D* hunf = (TH1D*)unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH1* hfold = response1D.ApplyToTruth(hunf, "");

      TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
      htempUnf->SetName(Form("Bayesian_Unfolded_%d",iter));
      
      TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
      htempFold->SetName(Form("Bayesian_Folded_%d",iter));        

      htempUnf->Write();
      htempFold->Write();
	  
    }
	  
}
#ifndef __CINT__
//int main () { RooSimplepTPbPb_data_Caitie(); return 0; }  // Main program when run stand-alone
#endif
