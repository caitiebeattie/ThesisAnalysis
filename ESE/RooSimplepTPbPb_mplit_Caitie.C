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
#include <ctime>
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
//#include "RooUnfoldTestHarness2D.h"

#endif

//==============================================================================
// Global definitions
//==============================================================================

const Double_t cutdummy= -99999.0;

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


   



void Normalize2D(TH2* h)
{
   int nbinsYtmp = h->GetNbinsY();
   const int nbinsY = nbinsYtmp;
   Double_t norm[nbinsY];
   for(int biny=1; biny<=nbinsY; biny++)  {
       norm[biny-1] = 0;
       for(int binx=1; binx<=h->GetNbinsX(); binx++) {
         norm[biny-1] += h->GetBinContent(binx,biny);
       }
     }

   for(int biny=1; biny<=nbinsY; biny++)
     {
       for(int binx=1; binx<=h->GetNbinsX(); binx++)
     {
       if(norm[biny-1]==0)  continue;
       else
         {
  h->SetBinContent(binx,biny,h->GetBinContent(binx,biny)/norm[biny-1]);
  h->SetBinError(binx,biny,h->GetBinError(binx,biny)/norm[biny-1]);
         }
     }
     }
}



 vector<const char*> filenames3 = {"TrainOutputs/old/AnalysisResults6846.root", "TrainOutputs/old/AnalysisResults6847.root", "TrainOutputs/old/AnalysisResults6848.root", "TrainOutputs/old/AnalysisResults6849.root", "TrainOutputs/old/AnalysisResults6850.root", "TrainOutputs/old/AnalysisResults6851.root", "TrainOutputs/old/AnalysisResults6852.root", "TrainOutputs/old/AnalysisResults6853.root", "TrainOutputs/old/AnalysisResults6854.root", "TrainOutputs/old/AnalysisResults6855.root", "TrainOutputs/old/AnalysisResults6856.root", "TrainOutputs/old/AnalysisResults6857.root", "TrainOutputs/old/AnalysisResults6858.root", "TrainOutputs/old/AnalysisResults6859.root", "TrainOutputs/old/AnalysisResults6860.root", "TrainOutputs/old/AnalysisResults6861.root", "TrainOutputs/old/AnalysisResults6862.root", "TrainOutputs/old/AnalysisResults6863.root", "TrainOutputs/old/AnalysisResults6864.root", "TrainOutputs/old/AnalysisResults6865.root"};



/*TH2D* CorrelationHist (const TMatrixD& cov,const char* name, const char* title,
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
*/
//=============================================================================
// Example Unfolding
//==============================================================================

void RooSimplepTPbPb_mplit_Caitie(double centl, double centr, int opt, string tag, TString cFiles2="files8.txt")
{
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  int difference=1;
  int Ppol=0;
  cout << "==================================== pick up the responsesplit matrix for background==========================" << endl;
  ///////////////////parameter setting
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;   
  
  TRandom3* rand = new TRandom3();

   //***************************************************

  TH1F* manrem = new TH1F("manrem", "manrem", 90, 30.0, 120.0);

  std::vector<Double_t> kBinsUnfolded = {10, 15, 20, 25, 30, 40, 50, 60, 70, 85, 100, 120, 140, 190/*, 250*/};
  std::vector<Double_t> kBinsMeasured = {30, 40, 50, 60, 70, 85, 100, 120};

  //raw jets used for spectra (0.3)
  TH1F *h1rawdata = new TH1F("rawdata", "rawdata", kBinsMeasured.size()-1, kBinsMeasured.data());
  //raw jets used for response (0.7)
  TH1F *h1rawresponse = h1rawresponse = new TH1F("rawresponse", "rawresponse",kBinsMeasured.size()-1, kBinsMeasured.data());
  //true jets corresponding to spectra
  TH1F *h1true = new TH1F("true", "true", kBinsUnfolded.size()-1, kBinsUnfolded.data());

  //kinematic efficiency
  TH1F *h1kinPre = new TH1F("kinPre", "kinPre", kBinsUnfolded.size()-1, kBinsUnfolded.data()); 
  TH1F *h1kinPos = new TH1F("kinPos", "kinPos", kBinsUnfolded.size()-1, kBinsUnfolded.data());
  //reconstruction efficiency
  //TH1F *recPre = new TH1F("recPre", "recPre", );
  //TH1F *recPos = new TH1F("recPos", "recPos");

  int tot = 0;
  int small = 0;
  int large = 0;
      int largepT = 0;
     
  long int preeta = 0;
  long int posteta = 0;
  long int prematch = 0;
  long int postmatch = 0;

  TH2F *hcovariance(0);
  hcovariance=new TH2F("covariance","covariance",10,0.,1.,10,0,1.);

   h1rawresponse->Sumw2();
   h1true->Sumw2();
   h1rawdata->Sumw2();
   h1kinPre->Sumw2();
    
   TH2D *fullresponse = new TH2D("fullresponse", "fullresponse", 140, 10, 150, 140, 10, 150);

   //branches in the tree that you need in this analysis
   float ptJet,ptPar, ptdet, centrality, matchradius, eta, area, bkgd;
  
   int nEv=0;; 
   //so mcr is correctly normalized to one, not the responsesplit.       
   cout<<"cucu"<<endl;

cout << "a\n"; 
   ifstream infile2;
   infile2.open(cFiles2.Data());
   char filename2[300];
   //split test
   RooUnfoldResponse responsesplit;
   responsesplit.Setup(h1rawresponse, h1true);
   //responsesplit.Setup(h1rawresponse, h1rawdata);

cout <<"b\n";

    int totevents = 0;
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
      //pthardbin = htrials->GetMean();
      double n_event = hnevent->GetBinContent(1);
      double xsect = hcross->GetBinContent(htrials->GetMean() + 1);
      double xscale = xsect*hcross->GetEntries();
      double trials = htrials->GetBinContent(htrials->GetMean() + 1);
      double pTHardscalefactor=xscale/(n_event*trials);


    std::cout << "pT Hard Bin: " << htrials->GetMean() <<"\n";
    std::cout << "	trials: "<< trials << ", n_event:  " << n_event << ", Scale Factor:  " << pTHardscalefactor << std::endl;
    TFile *input2=TFile::Open(filename2);
    // get the mc tree
    TTree *mc=(TTree*)input2->Get("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet"); 
    int nEv=mc->GetEntries(); 

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
    mc->SetBranchAddress("Event_BackgroundDensity", &bkgd); 
 
    //Find Bins to Suppress
    for (int i = 0; i < nEv; i++) {
       mc->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       if(ptPar > 200.) continue;    
       if (ptJet < 10.0 || ptPar < 10.0 )  continue;
       manrem->Fill(ptJet);
     }

    //Store Bins to Suppress
    vector<int> manset(0);
    for (int b = 0; b < 95; b++) {
        if (opt == 1 || opt == 3)   {
           if (manrem->GetBinContent(b+1) != 1)    continue;}
        if (opt == 2 || opt == 4)   { 
           if (manrem->GetBinContent(b+1) == 0)    continue;}
        if (opt == 1 || opt == 2)   {
           if (manrem->GetBinContent(b) == 0 && manrem->GetBinContent(b+2) == 0)  manset.push_back(b+1);}
        if (opt == 3 || opt == 4)   {
        if (manrem->GetBinContent(b) == 0 || manrem->GetBinContent(b+2) == 0)  manset.push_back(b+1);}
    }

   
  
  //Fill Response
  int countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
        mc->GetEntry(iEntry);
          if (abs(eta) > 0.5)   continue;
          if(ptPar > 200.) continue;    
          if (centrality < centl || centrality > centr)  continue;
          if (ptPar < 10.0 /*|| ptdet < 10.0*/ || ptJet < 10.0)   continue; 
       //skip entries that don't meet outlier criteria
       bool outlier = false;
       for (int i = 0; i < manset.size(); i++)   {
           if (ptJet > manrem->GetXaxis()->GetBinLowEdge(manset[i]) && ptJet < manrem->GetXaxis()->GetBinLowEdge(manset[i]+1))  {
              outlier = true;
              break;
           }   
       }
       if (outlier == true)  continue;
      //****** Get the scaling factor
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      // put in if/else statements on the ptJet
      if(ptJet >= -20. && ptJet < 10.)        EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)    EBscale = 1.0/0.03;
      else if(ptJet >= 20. && ptJet < 30.)    EBscale = 1.0/0.03;
      else if(ptJet >= 30. && ptJet < 40.)    EBscale = 1.0/0.05;
      else if(ptJet >= 40. && ptJet < 60.)    EBscale = 1.0/0.05;
      else if(ptJet >= 60. && ptJet < 80.)    EBscale = 1.0/0.09;
      else if(ptJet >= 80. && ptJet < 100.)   EBscale = 1.0/0.09;
      else if(ptJet >= 100. && ptJet < 200.)  EBscale = 1.0/0.10;
      scalefactor*=EBscale;
      // *********
      prematch++;
      postmatch++;
      //if (ptJet > 3.0*ptPar)  continue;
      h1kinPre->Fill(ptPar, scalefactor);  
      if (ptJet < 30 || ptJet > 120) continue;
      h1kinPos->Fill(ptPar, scalefactor);

      tot++;
      //Fill Split Objects
      double split = rand->Rndm();
      if (split < 0.7)	{
	  h1rawresponse->Fill(ptJet, scalefactor);
	  //this is the half split to be the responsesplit 
	  responsesplit.Fill(ptJet, ptPar, scalefactor);          //responsesplit.Fill(ptJet, ptJet, scalefactor);
          fullresponse->Fill(ptJet, ptPar, scalefactor);
	}
      else	{
	  //this is the psuedo data!
          h1rawdata->Fill(ptJet, scalefactor);
	  //this is the generator level distribution for the pseudo data or our answer :)
	  h1true->Fill(ptPar, scalefactor);
	}
    }
   cout << "number of bins to be suppressed: " << manset.size() <<"\n";
   manset.resize(0);

   }


    //////////efficiencies done////////////////////////////////////
 
    
    stringstream ss;
    ss << "UnfoldingSplit_AreaBased_R04_non_" << Form("%.0f%.0f", centl, centr) << "_test_manu" << opt << "_"<< tag << ".root";

    TFile *fout = new TFile (ss.str().c_str(),"RECREATE");
    ss.str("");
    fout->cd();
    h1rawdata->SetName("rawdata");
    h1rawdata->Write();
    //h1rawfull->Write();
    h1rawresponse->SetName("rawresponse");
    h1rawresponse->Write();
    h1true->SetName("true");
    h1true->Write();
   //h1truefull->Write();
    responsesplit.Write();
    fullresponse->Write(); 
   //projection of the responsesplit onto the truth axis
    //TH1D* htruth = (TH1D*)responsesplit.Htruth();
    //htruth->SetName("htruth");
    //htruth->Write();

   cout << "Tot: " << tot <<"\n";


    for(int jar = 1; jar < 16; jar++) {
      int iter = jar;
      cout<<"iteration"<<iter<<endl;
      cout<<"==============Unfold h1====================="<<endl;

      RooUnfoldBayes unfold(&responsesplit, h1rawdata, iter, false);    // OR
      std::cout << "line324\n";
      TH1D* hunf = (TH1D*)unfold.Hreco(errorTreatment);
      //FOLD BACK
      std::cout << "line327\n";
      TH1* hfold = responsesplit.ApplyToTruth(hunf, "");

      TH2D *htempUnf=(TH2D*)hunf->Clone("htempUnf");          
      htempUnf->SetName(Form("Bayesian_Unfolded_%d",iter));
      
      TH2D *htempFold=(TH2D*)hfold->Clone("htempFold");          
      htempFold->SetName(Form("Bayesian_Folded_%d",iter));        

      std::cout << "line340\n";
      htempUnf->Write();
      std::cout << "line342\n";
      htempFold->Write();
	  
    }

cout << ss.str().c_str() <<"\n";
	  
}
#ifndef __CINT__
//int main () { RooSimplepTPbPb_split_Caitie(TString tag); return 0; }  // Main program when run stand-alone
#endif
