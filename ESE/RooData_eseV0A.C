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

//if R = 0.4, use files25.txt
//if R = 0.2, use files27.txt
void RooData_eseV0A(double centl, double centr, int ese, int R = 4, int epcut = 2, TString cFiles2="files29.txt")
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
  if (R == 2)   truncation = 20.0;
  if (R == 4)   truncation = 30.0;
  //if (R == 4)   truncation = 50.0;

  //establish detectors
  string epDet = "V0A";
  string q2Det = "V0C";


  //=============================================================================
  //=================   Initialize Histograms and Responses  ====================   
  //=============================================================================
  //std::vector<double> kBinsUnfolded = {10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 70.0, 90.0, 120.0, 140.0, 190.0, 250.0};
  std::vector<double> kBinsUnfolded = {10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 140.0, 190.0, 250.0};
  //std::vector<double> kBinsUnfolded = {10.0, 20.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0};
  //std::vector<double> kBinsUnfolded = {10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 85.0, 100.0, 120.0};
  std::vector<double> kBinsMeasured;
        if (R == 2)  kBinsMeasured = {truncation, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0};   //R=0.2
        //if (R == 2)  kBinsMeasured = {truncation, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 85.0, 100.0, 120.0};
        //if (R == 2)  kBinsMeasured = {truncation, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 85.0, 100.0};
        //if (R == 4)  kBinsMeasured = {truncation, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 85.0, 100.0};
        //if (R == 4)  kBinsMeasured = {truncation, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 85.0, 100.0, 120.0};
        if (R == 4)  kBinsMeasured = {truncation, 35.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0};   //R=0.2
  std::vector<double> kBinsEPangle = {0.0, sqrt(1)/2.0, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};

  //raw spectra from embedding (pseudo-data)
  TH1D *h1raw = new TH1D("raw1", "raw1", kBinsMeasured.size()-1, kBinsMeasured.data());
  TH2D *h2raw = new TH2D("raw2", "raw2", kBinsMeasured.size()-1, kBinsMeasured.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
  TH1D *h1rawep = new TH1D("h1rawep", "h1rawep", kBinsEPangle.size()-1, kBinsEPangle.data());
        h1raw->Sumw2();
        h2raw->Sumw2();
        h1rawep->Sumw2();
  //true spectra from embedding
  TH1F *h1true = new TH1F("true1", "true1", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  TH2F *h2true = new TH2F("true2", "true2", kBinsUnfolded.size()-1, kBinsUnfolded.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
        h1true->Sumw2();
        h2true->Sumw2();

  //raw spectra from data
  TH1D* specraw1 = new TH1D("specraw1", "specraw1", kBinsMeasured.size()-1, kBinsMeasured.data());
  TH2D* specraw2 = new TH2D("specraw2", "specraw2", kBinsMeasured.size()-1, kBinsMeasured.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
  TH1D* specrawep1 = new TH1D("specrawep1", "specrawep1", kBinsEPangle.size()-1, kBinsEPangle.data());
        specraw1->Sumw2();
        specraw2->Sumw2();
        specrawep1->Sumw2();

  //Kinematic efficiency
  TH1F *hKinPre = new TH1F("KinPre", "KinPre", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  TH1F *hKinPos = new TH1F("KinPos", "KinPos", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  TH2F *h2KinPre = new TH2F("KinPre2", "KinPre2", kBinsUnfolded.size()-1, kBinsUnfolded.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
  TH2F *h2KinPos = new TH2F("KinPos2", "KinPos2", kBinsUnfolded.size()-1, kBinsUnfolded.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
        hKinPre->Sumw2();
        hKinPos->Sumw2();  
        h2KinPre->Sumw2();
        h2KinPos->Sumw2();

  //Responses
  RooUnfoldResponse response1D;
  response1D.Setup(h1raw, h1true);
  RooUnfoldResponse response2D;
  response2D.Setup(h2raw, h2true);


  //=============================================================================
  //==========================   Get Raw Spectra  ===============================   
  //=============================================================================
   
    //Raw Spectra from Trains
    string rootfilename;
    if (R == 2)   rootfilename = "AnalysisResults8144.root";
    //if (R == 4)   rootfilename = "AnalysisResults8294.root";
    if (R == 4)   rootfilename = "AnalysisResults8119.root";
    //if (R == 4)   rootfilename = "AnalysisResults8316.root";
    const char *lerootfile = rootfilename.c_str(); 

    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);      

      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2_EP");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      TH2F *qAhist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      int centdiff = (int)centr-(int)centl;
      int leftc, rightc;
      int iRebin = q2hist2D->GetNbinsY()/q2hist2D->GetNbinsY();
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < centdiff; j++) {
        //reset percentile ticker
        long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
        int percentilen = 1;                                             //number percentile we're on
        //update centrality
        leftc = (int)centl + j; 
        rightc = (int)centl + j + 1;    
        cout << "Cent: " << leftc <<"\n";
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);
        //project slice onto q2 axis
        TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();
        q2hist->Rebin(iRebin);
        long double percentilemarker = (long double)q2hist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        //increment along q2 axis until percentile is surpassed
        for (int i = 0; i <= q2hist->GetNbinsX(); i++)   {
          percentileticker += q2hist->GetBinContent(i);         //increment the percentileticker
          if (percentileticker >= percentilen*percentilemarker)  {
              percentiles.push_back(q2hist->GetBinCenter(i));
              percentilen++;
          }
        if (i == q2hist->GetNbinsX())  {
            allPercent.push_back(percentiles);
            percentiles.resize(0);
            break;}             
        }
      }
      //calculate percentiles in 2nd detector
      vector<long double> percentilesA(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercentA(0);                      //vector that stores percentiles vectors for all centralities
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < centdiff; j++) {
        //reset percentile ticker
        long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
        int percentilen = 1;                                             //number percentile we're on
        leftc = (int)centl + j; 
        rightc = (int)centl + j + 1;    
        qAhist2D->GetXaxis()->SetRangeUser(leftc, rightc);
        TH1F *qAhist = (TH1F*)qAhist2D->ProjectionY();
        qAhist->Rebin(iRebin);
        long double percentilemarker = (long double)qAhist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        for (int i = 0; i <= qAhist->GetNbinsX(); i++)   {
          percentileticker += qAhist->GetBinContent(i);         //increment the percentileticker
          if (percentileticker >= percentilen*percentilemarker)  {
              percentilesA.push_back(qAhist->GetBinCenter(i));
              percentilen++;
          }
        if (i == qAhist->GetNbinsX())  {
            allPercentA.push_back(percentilesA);
            percentilesA.resize(0);
            break;}             
        }
      }
  TH1D *perc20 = new TH1D("perc20", "perc20", centdiff, centl, centr);
  TH1D *perc20C = new TH1D("perc20C", "perc20C", centdiff, centl, centr);
  TH1D *perc80 = new TH1D("perc80", "perc80", centdiff, centl, centr);
  TH1D *perc80C = new TH1D("perc8C0", "perc80C", centdiff, centl, centr);
  for (int i = 1; i <= centdiff; i++) {
      perc20->SetBinContent(i, allPercentA[i-1][3]);
      perc20C->SetBinContent(i, allPercent[i-1][6]);
      perc80->SetBinContent(i, allPercentA[i-1][6]);
      perc80C->SetBinContent(i, allPercent[i-1][6]);
  }
  perc20->SetLineColor(kRed);
    perc20->SetLineWidth(2);
  perc20C->SetLineColor(kRed);
    perc20C->SetLineWidth(2);
  perc80->SetLineColor(kRed);
    perc80->SetLineWidth(2);
  perc80C->SetLineColor(kRed);
    perc80C->SetLineWidth(2);


    //Determine Extractor Bin Scaling
    const int extraBins = 9;
    double extraFrac[extraBins] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double extraEdge[extraBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0};
    
    //set branch names
    stringstream Tq2Det, TepDet;
    Tq2Det << "Event_Q2Vector" << q2Det;
    TepDet << "Jet_EPangle" << epDet;

    //Retrieve Tree and Populate Response
    TTree *T = nullptr;
    if (R == 2) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
    if (R == 4) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_Jet", T);
    float jetpT, centra, qvec, ep;
    T->SetBranchAddress("Jet_Pt", &jetpT);
    T->SetBranchAddress("Event_Centrality", &centra);
    T->SetBranchAddress(Tq2Det.str().c_str(), &qvec);    
    T->SetBranchAddress(TepDet.str().c_str(), &ep);    
    int entries = T->GetEntries();
    for (int i = 0; i < entries; i++)   {
       T->GetEntry(i);     
       if (centra < centl || centra > centr)                                    continue;       
       if (jetpT < truncation || jetpT > 120.0)                                 continue;
       if (ese == 0)           {if (qvec > allPercent[int(centra-centl)][2])    continue;}  //20th percentile
       if (ese == 1)           {if (qvec < allPercent[int(centra-centl)][6])    continue;}  //80th percentile
       if (ese == 3)          
           {if (qvec < allPercent[int(centra-centl)][0])    continue;  //80th percentile
           if (qvec > allPercent[int(centra-centl)][3])    continue;}  //80th percentile
       if (epcut == 0)         {if (abs(cos(ep)) > sqrt(2)/2.0)                 continue;}  //out-of-plane
       if (epcut == 1)         {if (abs(cos(ep)) < sqrt(2)/2.0)                 continue;}  //in-plane
       int scaleindex = -1;
         if (jetpT >= extraEdge[0] && jetpT < extraEdge[1])  scaleindex = 0;
         if (jetpT >= extraEdge[1] && jetpT < extraEdge[2])  scaleindex = 1;
         if (jetpT >= extraEdge[2] && jetpT < extraEdge[3])  scaleindex = 2;
         if (jetpT >= extraEdge[3] && jetpT < extraEdge[4])  scaleindex = 3;
         if (jetpT >= extraEdge[4] && jetpT < extraEdge[5])  scaleindex = 4;
         if (jetpT >= extraEdge[5] && jetpT < extraEdge[6])  scaleindex = 5;
         if (jetpT >= extraEdge[6] && jetpT < extraEdge[7])  scaleindex = 6;
         if (jetpT >= extraEdge[7] && jetpT < extraEdge[8])  scaleindex = 7;
         if (jetpT >= extraEdge[8] && jetpT < extraEdge[9])  scaleindex = 8;
           if (scaleindex == 0)   continue;      
       specraw1->Fill(jetpT, 1.0/extraFrac[scaleindex]);
       specraw2->Fill(jetpT, abs(cos(ep)), 1.0/extraFrac[scaleindex]);
       specrawep1->Fill(abs(cos(ep)), 1.0/extraFrac[scaleindex]);
    }

  
  TH2F *hcovariance(0);
  hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);

 
   //h1smeared->Sumw2();
  


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

    //set branch names
    stringstream MCq2Det, MCepDet, MCepDetMatch;
    MCq2Det << "Event_Q2Vector" << q2Det;
    MCepDet << "Jet_EPangle" << epDet;
    MCepDetMatch << "Jet_MC_MatchedPartLevelJet_EPangle" << epDet;

    // set the branches of the mc tree
    mc->SetBranchAddress("Jet_Pt", &ptJet); 
    mc->SetBranchAddress("Jet_Eta", &eta);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &ptPar);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptdet);
    mc->SetBranchAddress("Event_Centrality", &centrality);
    mc->SetBranchAddress("Jet_Area", &area);
    mc->SetBranchAddress(MCq2Det.str().c_str(), &q2);
    mc->SetBranchAddress(MCepDet.str().c_str(), &epangle);
    mc->SetBranchAddress(MCepDetMatch.str().c_str(), &epPar);
    
  
  int countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
      mc->GetEntry(iEntry);
      //****** Get the scaling factor
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      if(ptJet >= -20. && ptJet < 10.)         EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)     EBscale = 1.0/0.01;
      else if(ptJet >= 20. && ptJet < 30.)     EBscale = 1.0/0.04;
      else if(ptJet >= 30. && ptJet < 40.)     EBscale = 1.0/0.15;
      else if(ptJet >= 40. && ptJet < 60.)     EBscale = 1.0/0.25;
      else if(ptJet >= 60. && ptJet < 80.)     EBscale = 1.0/0.35;
      else if(ptJet >= 80. && ptJet < 100.)    EBscale = 1.0/0.35;
      else if(ptJet >= 100. && ptJet < 140.)   EBscale = 1.0/0.15;
      else if(ptJet >= 140. && ptJet < 200.)   EBscale = 1.0/0.05;
      scalefactor*=EBscale;
      // *********

      if (centrality < centl || centrality > centr)  continue;
      if (epcut == 0)    {if (abs(cos(epangle)) > sqrt(2)/2.0) continue;}  //out-plane
      if (epcut == 1)    {if (abs(cos(epangle)) < sqrt(2)/2.0) continue;}  //in-plane
      if (ese == 0)           {if (q2 > allPercent[int(centrality-centl)][2])    continue;}  //20th percentile
      if (ese == 1)           {if (q2 < allPercent[int(centrality-centl)][6])    continue;}  //80th percentile
      if (ese == 3)  
           {if (q2 < allPercent[int(centrality-centl)][0])    continue;  //80th percentile
           if (q2 > allPercent[int(centrality-centl)][3])    continue;}  //80th percentile

      preeta++;      
      if (abs(eta) > (0.9-(double)R*0.1))  continue;  //TPC fiducial cut
      posteta++;
      if(ptPar > 200.) continue;    
  
      
      hKinPre->Fill(ptPar, scalefactor);
      h2KinPre->Fill(ptPar, abs(cos(epPar)), scalefactor);
      if (ptJet < truncation && ptPar > 50.0)  {
          if (pthardbin-1 <= 5)  {
          cout << "BUMP ALERT\n";
          }} 
      if (ptJet < truncation || ptJet > 120) continue;
      hKinPos->Fill(ptPar, scalefactor);
      h2KinPos->Fill(ptPar, abs(cos(epPar)), scalefactor);
      if (ptPar < 10.0 || ptJet < 10.0 )   continue; 

	  h1true->Fill(ptPar, scalefactor);
	  h2true->Fill(ptPar, abs(cos(epPar)), scalefactor);
	  h1raw->Fill(ptJet, scalefactor);
	  h2raw->Fill(ptJet, abs(cos(epangle)), scalefactor);
          h1rawep->Fill(abs(cos(epangle)), scalefactor);
	  response1D.Fill(ptJet, ptPar, scalefactor);
	  response2D.Fill(ptJet, abs(cos(epangle)), ptPar, abs(cos(epPar)), scalefactor);
    }}


    long double avgevents = (long double)totevents/20.0;
    h1raw->Scale(avgevents);
    h1true->Scale(avgevents);
    //////////efficiencies done////////////////////////////////////



  //==================================================================================
  //==========================   Setup the Outputs   =================================   
  //==================================================================================

    stringstream ww;
    if (ese == 0)  ww << "UnfoldingData_2D_non_R0" << R << "_q" << q2Det << "_ep" << epDet << "_" 
                      << Form("%.0f", centl) << Form("%.0f", centr) << "_lo30_" << Form("%.0f", truncation) << "_Dec1.root"; 
    if (ese == 1)  ww << "UnfoldingData_2D_non_R0" << R << "_q" << q2Det << "_ep" << epDet << "_" 
                      << Form("%.0f", centl) << Form("%.0f", centr) << "_hi30_" << Form("%.0f", truncation) << "_Dec1.root"; 
    if (ese == 2)  ww << "UnfoldingData_2D_non_R0" << R << "_q" << q2Det << "_ep" << epDet << "_"
                      << Form("%.0f", centl) << Form("%.0f", centr) << "_allq_" << Form("%.0f", truncation) << "_Dec1.root"; 
    if (ese == 3)  ww << "UnfoldingData_2D_non_R0" << R << "_q" << q2Det << "_ep" << epDet << "_" 
                      << Form("%.0f", centl) << Form("%.0f", centr) << "_1040_" << Form("%.0f", truncation) << "_Dec1.root"; 
    TFile *fout=new TFile (ww.str().c_str(),"RECREATE");
    fout->cd();
    h1raw->SetName("raw");
    h1raw->Write();
    h2raw->SetName("raw2");
    h2raw->Write();
    h1rawep->SetName("rawep1");
    h1rawep->Write();
    h1true->SetName("true1");
    h1true->Write();
    h2true->SetName("true2");
    h2true->Write();
    specraw1->SetName("specraw1");
    specraw1->Write();
    specraw2->SetName("specraw2");
    specraw2->Write();
    specrawep1->SetName("specrawep1");
    specrawep1->Write();
    hKinPre->Write();
    hKinPos->Write();
    h2KinPre->Write();
    h2KinPos->Write();
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
      cout<< "Iteration" << iter << endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&response2D, specraw2, iter, false);    // OR
      TH2D* hunf = (TH2D*)unfold.Hreco(errorTreatment);
      //FOLD BACK
      TH1* hfold = response2D.ApplyToTruth(hunf, "");

      TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
      htempUnf->SetName(Form("Bayesian_Unfolded_%d",iter));
      
      TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
      htempFold->SetName(Form("Bayesian_Folded_%d",iter));        

      htempUnf->Write();
      htempFold->Write();
	  
    }
	  

  //==================================================================================
  //========================== Draw QA Plots =========================================   
  //==================================================================================
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
          //pad1->SetLogz();
       qAhist2D->GetXaxis()->SetRangeUser(centl, centr);
       qAhist2D->GetXaxis()->SetTitle("Centrality (%)");
       qAhist2D->GetYaxis()->SetTitle("q_{2}^{V01}");
       qAhist2D->Draw("same colz");
       perc20->Draw("same");
       perc80->Draw("same");

  TCanvas *c2 = new TCanvas("Q2 vs CentC", "Q2 vs CentC", 500, 500);
  c2->cd();
       gStyle->SetOptStat(0);
       //gStyle->SetPadTickX(1);
       //gStyle->SetPadTickY(1);
       TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
          pad2->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad2->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad2->SetRightMargin(0.15);
          pad2->Draw();
          pad2->cd();
          //pad1->SetLogz();
       q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
       q2hist2D->GetXaxis()->SetTitle("Centrality (%)");
       q2hist2D->GetYaxis()->SetTitle("q_{2}^{V02}");
       q2hist2D->Draw("same colz");
       perc20C->Draw("same");
       perc80C->Draw("same");


}
#ifndef __CINT__
//int main () { RooSimplepTPbPb_data_Caitie(); return 0; }  // Main program when run stand-alone
#endif
