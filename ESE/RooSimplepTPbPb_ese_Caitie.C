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

TH2D* CorrelationHistShape (const TMatrixD& cov,const char* name, const char* title, int na, int nb, int kbin)
{ 
   TH2D* h= new TH2D (name, title, nb, 0, nb, nb, 0, nb); 
     	  for(int l=0;l<nb;l++){
             for(int n=0;n<nb;n++){
                int index1=kbin+na*l;
                int index2=kbin+na*n;
   	        Double_t Vv=cov(index1,index1)*cov(index2,index2);
   		if (Vv>0.0) h->SetBinContent(l+1,n+1,cov(index1,index2)/sqrt(Vv));
  	      }
           }
  return h;
}


TH2D* CorrelationHistPt (const TMatrixD& cov,const char* name, const char* title, int na, int nb, int kbin)
{
   TH2D* h= new TH2D (name, title, na, 0, na, na, 0, na);
     	  for(int l=0;l<na;l++){
             for(int n=0;n<na;n++){
                int index1=l+na*kbin;
                int index2=n+na*kbin;
   	        Double_t Vv=cov(index1,index1)*cov(index2,index2);
  		if (Vv>0.0)h->SetBinContent(l+1,n+1,cov(index1,index2)/sqrt(Vv));
             }
           }
  return h;
}


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

void RooSimplepTPbPb_ese_Caitie(double centl, double centr, int ese, int R = 2, TString cFiles2="files12.txt")
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
  
  //###########################
  //##   Get Raw Spectra     ##
  //###########################
   //Raw Spectra from Trains
   //string rootfilename = "../Spectra/AnalysisResults7283.root";   //semi-central
    string rootfilename = "AnalysisResults7665_2.root";
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
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_V2_TPCEta0dot8");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrFullV0");
      //Fill q2 hist
      cout << "q2hist2D->GetEntries(): " << q2hist2D->GetEntries() <<"\n";
      TH1F *q2all = (TH1F*)q2hist2D->ProjectionY();
      q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
      TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();
      q2hist->Rebin(100);
      cout << "q2hist->GetEntries(): " << q2hist->GetEntries() <<"\n";
      long int percentileticker = 0;                              //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                   //number percentile we're on
      long double percentilemarker = (long double)q2hist->GetEntries()/10.0;        //number of entries that compose 10% of sample
      cout << "Percentile Entries: " << percentilemarker <<"\n";
      vector<long double> percentiles(0);                         //vector that stores 10th, 20th etc percentile of q2
        for (int i = 0; i <= 100; i++)   {
          percentileticker += q2hist->GetBinContent(i);         //increment the percentileticker
          if (percentileticker >= percentilen*percentilemarker)  {
              percentiles.push_back(q2hist->GetBinCenter(i));
              percentilen++;
          }
        if (percentiles.size() >= 10)   break; 
      }
      cout << "20th percentile: " << percentiles[1] <<"\n";
      cout << "80th percentile: " << percentiles[7] <<"\n";


    //Determine Extractor Bin Scaling
    const int extraBins = 8;
    double extraFrac[extraBins] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double extraEdge[extraBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 200.0};
    
    //const int pbins = 18;
    //double pTedgepp[pbins + 1]  = {25, 30, 40, 50, 60, 70, 85, 100};
    TH1D* specraw = new TH1D("specraw", "specraw", kBinsMeasured.size()-1, kBinsMeasured.data());
    int entries = T->GetEntries();
    for (int i = 0; i < entries; i++)   {
       T->GetEntry(i);     
       if (centra < centl || centra > centr)    continue;       
       if (ese == 0)  {if (qvec > percentiles[1]) continue;} //20th percentile
       //if (ese == 0)  {if (qvec > 3.55) continue;} //40th percentile
       //if (ese == 1)  {if (qvec < 4.65) continue;} //60th percentile
       if (ese == 1)  {if (qvec < percentiles[7]) continue;} //80th percentile
       if (jetpT < truncation || jetpT > 120.0)       continue;
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
    }


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
      else if(ptJet >= 10. && ptJet < 20.)    EBscale = 1.0/0.07;
      else if(ptJet >= 20. && ptJet < 30.)    EBscale = 1.0/0.07;
      else if(ptJet >= 30. && ptJet < 40.)    EBscale = 1.0/0.09;
      else if(ptJet >= 40. && ptJet < 60.)    EBscale = 1.0/0.09;
      else if(ptJet >= 60. && ptJet < 80.)    EBscale = 1.0/0.10;
      else if(ptJet >= 80. && ptJet < 100.)   EBscale = 1.0/0.10;
      else if(ptJet >= 100. && ptJet < 200.)  EBscale = 1.0/0.10;
      scalefactor*=EBscale;
      // *********

      //if (abs(cos(epangle)) < sqrt(2)/2.0) continue;  //in-plane
      //if (abs(cos(epangle)) > sqrt(2)/2.0) continue;  //out-plane
      if (centrality < centl || centrality > centr)  continue;
      if (ese == 0)  {if (q2 > percentiles[1]) continue;} //20th percentile
      //if (ese == 0)  {if (q2 > 3.55) continue;} //40th percentile
      //if (ese == 1)  {if (q2 < 4.65) continue;} //60th percentile
      if (ese == 1)  {if (q2 < percentiles[7]) continue;} //80th percentile

      preeta++;      
      if (abs(eta) > (0.9-(double)R*0.1))  continue;  //TPC fiducial cut
      posteta++;
      if(ptPar > 200.) continue;    
  
      
      hKinPre->Fill(ptPar, scalefactor);
      if (ptJet < truncation || ptJet > 120) continue;
      hKinPos->Fill(ptPar, scalefactor);
      if (ptPar < 10.0 || ptJet < 10.0 /*||  ptdet < 10.0*/)   continue; 
      //if (ptJet > ptPar*((double)hvt/10.0))   continue;
      tot++;

	  h1true->Fill(ptPar, scalefactor);
	  h1raw->Fill(ptJet, scalefactor);
	  responsesplit.Fill(ptJet, ptPar, scalefactor);
    }}


    
    TH1F *htrueptd = (TH1F*) h1fulleff->Clone("trueptd");
    TH1F *htruept  = (TH1F*) h1fulleff->Clone("truept"); 

    cout << "Tot: " << tot <<"\n";
    cout << "LargepT: " << largepT <<"\n"; 

    long double avgevents = (long double)totevents/20.0;
    h1raw->Scale(avgevents);
    h1true->Scale(avgevents);
    //////////efficiencies done////////////////////////////////////


    stringstream ww;
    if (ese == 0)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_lo20_Dec2.root"; 
    if (ese == 1)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_hi20_Dec2.root"; 
    if (ese == 2)  ww << "UnfoldingData_non_R02_" << Form("%.0f", centl) << Form("%.0f", centr) << "_allq_Dec2.root"; 
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
	  
      ///HERE I GET THE COVARIANCE MATRIX/////
      /*
      if(iter==8){
	TMatrixD covmat= unfold.Ereco((RooUnfold::ErrorTreatment)RooUnfold::kCovariance);
	for(int k=0;k<h2true->GetNbinsX();k++){
	  TH2D *hCorr= (TH2D*) CorrelationHistShape(covmat, Form("corr%d",k), "Covariance matrix",h2true->GetNbinsX(),h2true->GetNbinsY(),k);
	  TH2D *covshape=(TH2D*)hCorr->Clone("covshape");      
	  covshape->SetName(Form("pearsonmatrix_iter%d_binshape%d",iter,k));
	  covshape->SetDrawOption("colz");
	  covshape->Write();
	}
	
	for(int k=0;k<h2true->GetNbinsY();k++){
	  TH2D *hCorr= (TH2D*) CorrelationHistPt(covmat, Form("corr%d",k), "Covariance matrix",h2true->GetNbinsX(),h2true->GetNbinsY(),k);
	  TH2D *covpt=(TH2D*)hCorr->Clone("covpt");      
	  covpt->SetName(Form("pearsonmatrix_iter%d_binpt%d",iter,k));
	  covpt->SetDrawOption("colz");
	  covpt->Write();
	  }
	  }*/
    }
	  
}
#ifndef __CINT__
//int main () { RooSimplepTPbPb_data_Caitie(); return 0; }  // Main program when run stand-alone
#endif
