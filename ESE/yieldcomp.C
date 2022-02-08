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
void updateR (double &phiscale, double &etascale, double jetR) {
    phiscale = (2*pi)/(1.92 - 2*jetR);
    etascale = (1.4 - 2*jetR);
}

bool cent = false;
bool semi = false;

int centcolor = kViolet+5;
int semicolor = kCyan-4;

vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};



//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void yieldcomp (int iteration, string lofilename, string hifilename, int planestatus = 0) {

  const char* lofile = lofilename.c_str();
  const char* hifile = hifilename.c_str();
  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

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
      semi = true;
      centl = 30.0;
      centr = 50.0;

      //Get File
      TFile *specfile = new TFile("AnalysisResults7665_2.root");
   
      //Get N_event
      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral()/5.0;
      cout << "N_event: " << N_event <<"\n";

      //Get qn percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)specfile->Get("PWGJE_QnVectorTender");
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

      //Get Statistical Uncertainty
      TTree *T = nullptr;
      specfile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
      float pTjet, centra, q2, epangle;
      T->SetBranchAddress("Jet_Pt", &pTjet);
      T->SetBranchAddress("Event_Centrality", &centra);
      T->SetBranchAddress("Event_Q2Vector", &q2);
      T->SetBranchAddress("Jet_EPangle", &epangle);
      vector<double> stat_lo(7);
      vector<double> stat_hi(7);
      for (int i = 0; i < 7; i++)   stat_lo[i] = 0;
      for (int i = 0; i < 7; i++)   stat_hi[i] = 0;
      int entries = T->GetEntries();
      //calculate N per bin
      for (int i = 0; i < entries; i++) {
        T->GetEntry(i);
        if (centra < centl || centra > centr)   continue;
        if (q2 > percentiles[1])   continue;
        if (planestatus == 1) {if (abs(cos(epangle)) < sqrt(2)/2.0) continue;};  //in-plane  
        if (planestatus == 2) {if (abs(cos(epangle)) > sqrt(2)/2.0) continue;};  //out-of-plane
        if (pTjet > 25.0 && pTjet < 30.0 )   stat_lo[0]++;
        if (pTjet > 30.0 && pTjet < 40.0 )   stat_lo[1]++;
        if (pTjet > 40.0 && pTjet < 50.0 )   stat_lo[2]++;
        if (pTjet > 50.0 && pTjet < 60.0 )   stat_lo[3]++;
        if (pTjet > 60.0 && pTjet < 70.0 )   stat_lo[4]++;
        if (pTjet > 70.0 && pTjet < 85.0 )   stat_lo[5]++;
        if (pTjet > 85.0 && pTjet < 100.0)   stat_lo[6]++;
      }
      for (int i = 0; i < entries; i++)   {
        T->GetEntry(i);
        if (centra < centl || centra > centr)   continue;
        if (q2 < percentiles[7])   continue;
        if (planestatus == 1) {if (abs(cos(epangle)) < sqrt(2)/2.0) continue;};  //in-plane  
        if (planestatus == 2) {if (abs(cos(epangle)) > sqrt(2)/2.0) continue;};  //out-of-plane
        if (pTjet > 25.0 && pTjet < 30.0 )   stat_hi[0]++;
        if (pTjet > 30.0 && pTjet < 40.0 )   stat_hi[1]++;
        if (pTjet > 40.0 && pTjet < 50.0 )   stat_hi[2]++;
        if (pTjet > 50.0 && pTjet < 60.0 )   stat_hi[3]++;
        if (pTjet > 60.0 && pTjet < 70.0 )   stat_hi[4]++;
        if (pTjet > 70.0 && pTjet < 85.0 )   stat_hi[5]++;
        if (pTjet > 85.0 && pTjet < 100.0)   stat_hi[6]++;
      }
      //convert to fractional error
      for (int i = 0; i < 7; i++)        {
        stat_lo[i] = (1.0/(double)stat_lo[i])*pow(stat_lo[i], 0.5);
        stat_hi[i] = (1.0/(double)stat_hi[i])*pow(stat_hi[i], 0.5);
      }

  //======================================================
  //========= 4. retrieve unfolded spectra from file
  //=======================================================
    TFile *lefilel = new TFile(lofile);
    TFile *lefileh = new TFile(hifile);
      double phifrac, etafrac;
      updateR(phifrac, etafrac, 0.4);

    //Unfolded Hists for Low q2
    vector<TH1F*> unfoldhistsl(0);
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_1"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_2"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_3"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_4"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_5"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_6"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_7"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_8"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_9"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_10"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_11"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_12"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_13"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_14"));
        unfoldhistsl.push_back((TH1F*)lefilel->Get("Bayesian_Unfolded_15"));

    //Unfolded Hists for High q2
    vector<TH1F*> unfoldhistsh(0);
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_1"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_2"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_3"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_4"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_5"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_6"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_7"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_8"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_9"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_10"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_11"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_12"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_13"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_14"));
        unfoldhistsh.push_back((TH1F*)lefileh->Get("Bayesian_Unfolded_15"));




  //=======================================================
  //========= 5. retrieve efficiencies from file ==========
  //=======================================================
        //Get Kinematic Efficiency
        TH1F *kinpre = (TH1F*)lefilel->Get("KinPre");
        TH1F *kinpos = (TH1F*)lefilel->Get("KinPos");
             TH1F *effic = (TH1F*)kinpos->Clone("effic");
             effic->Divide(kinpre);

       //Get Reconstruction Efficiency
       TFile *recfile = new TFile("../Unfolding/ChargedJetRecEfficiencies.root");
       TH1D  *receffic = (TH1D*)recfile->Get("RecEff_R040_5GeV");


  //=============================================
  //============= Rebin PbPb to match pp ========
  //=============================================
      vector<double> pTedgerebin = {25, 30, 40, 50, 60, 70, 85, 100};
      TH1F* efficiency = (TH1F*)effic->Rebin(pTedgerebin.size()-1, "efficiency", pTedgerebin.data());
      TH1F* recefficiency = (TH1F*)receffic->Rebin(pTedgerebin.size()-1, "recefficiency", pTedgerebin.data());
          recefficiency->Scale(1.0, "width");
      TH1F* lorebin = new TH1F("lorebin", "", pTedgerebin.size()-1, pTedgerebin.data());
          for (int i = 0; i < pTedgerebin.size()-1; i++)   {
                int binnum = unfoldhistsl[iteration]->GetXaxis()->FindBin(pTedgerebin[i]+0.01);
                double weight = unfoldhistsl[iteration]->GetBinContent(binnum);
                lorebin->Fill(pTedgerebin[i]+0.01, weight);
          }
      TH1F* hirebin = new TH1F("hirebin", "", pTedgerebin.size()-1, pTedgerebin.data());
          for (int i = 0; i < pTedgerebin.size()-1; i++)   {
                int binnum = unfoldhistsh[iteration]->GetXaxis()->FindBin(pTedgerebin[i]+0.01);
                double weight = unfoldhistsh[iteration]->GetBinContent(binnum);
                hirebin->Fill(pTedgerebin[i]+0.01, weight);
          }
  


//============================================================
  //========= 6. apply scaling for bin width + N_events ========
  //============================================================
        //Unfolded Spectra 
        for (int i = 0; i < unfoldhistsl.size(); i++)    {
        unfoldhistsl[i]->Scale(1.0, "width");
        //unfoldhistsl[i]->Scale(phifrac/etafrac);
        unfoldhistsl[i]->Scale(1.0/(double)N_event);}


        for (int i = 0; i < unfoldhistsh.size(); i++)    {
        unfoldhistsh[i]->Scale(1.0, "width");
        //unfoldhistsl[i]->Scale(phifrac/etafrac);
        unfoldhistsh[i]->Scale(1.0/(double)N_event);}






 
  //=====================================================
  //============= Calculate the Yield Ratio =============
  //=====================================================
 
      TH1F* yieldratio =(TH1F*)hirebin->Clone("yieldratio");
      yieldratio->Divide(lorebin);

      double lowedge;
      if (cent == true)  lowedge = 60.0;
      if (semi == true)  lowedge = 30.0;
      int startbin = yieldratio->FindBin(lowedge+0.01);
      yieldratio->GetXaxis()->SetRangeUser(lowedge, 100.0);

      
      yieldratio->SetLineColorAlpha(kBlack, 0.0);
      yieldratio->SetMarkerColor(kBlack);
      yieldratio->SetMarkerStyle(20);
      yieldratio->SetMarkerSize(0.6);


  //===============================================================
  //============= Calculate Statistical Uncertainties =============
  //===============================================================
      
      //sum low and high q2 statistical errors in quadrature
      vector<double> stat_comb(7); //stat_lo/hi[0] -> 25-30 GeV
      for (int i = 0; i < 7; i++)   stat_comb[i] = pow(pow(stat_lo[i],2) + pow(stat_hi[i],2), 0.5);
      //multiply fractional errors by yield ratio
      for (int i = 0; i < 7; i++)   stat_comb[i] = stat_comb[i]*yieldratio->GetBinContent(i+yieldratio->FindBin(25.01));


      //create TGraphErrors for statistical errors
      double bincenters010[3] = {65.0, 77.5, 92.5};
      double bincenters3050[6] = {35.0, 45.0, 55.0, 65.0, 77.5, 92.5};
      double binvalues010[3] = {yieldratio->GetBinContent(startbin), yieldratio->GetBinContent(startbin+1), yieldratio->GetBinContent(startbin+2)};
      double binvalues3050[6] = {yieldratio->GetBinContent(startbin), yieldratio->GetBinContent(startbin+1), yieldratio->GetBinContent(startbin+2), yieldratio->GetBinContent(startbin+3), yieldratio->GetBinContent(startbin+4), yieldratio->GetBinContent(startbin+5)};
      double nul010[3] = {5.0, 7.5, 7.5};
      double nul3050[6] = {5.0, 5.0, 5.0, 5.0, 7.5, 7.5};
      //convert stat_comb to correct type
      double stat_comb010[3] = {stat_comb[2], stat_comb[3], stat_comb[4]};
      double stat_comb3050[6] = {stat_comb[1], stat_comb[2], stat_comb[3], stat_comb[4], stat_comb[5], stat_comb[6]};
      //input to TGraphErrors
      TGraphErrors *stat;
      if (cent == true) stat = new TGraphErrors(3, bincenters010, binvalues010, nul010, stat_comb010);
      if (semi == true) stat = new TGraphErrors(6, bincenters3050, binvalues3050, nul3050, stat_comb3050);
           stat->SetMarkerStyle(20);
           stat->SetMarkerSize(0.7);
           stat->SetMarkerColor(kBlack);
           stat->SetLineColor(kBlack);


  //===============================================================
  //============= Calculate Systematic Uncertainties ==============
  //===============================================================
  /*    //get systematic errors from file
      TFile *errorfile;
      if (cent == true) errorfile = new TFile("../Systematics/Systematics010.root");
      if (semi == true) errorfile = new TFile("../Systematics/Systematics3050.root");
      TH1F *errors = (TH1F*)errorfile->Get("h_err"); 
      int startbinerror = errors->FindBin(40.01);
      double syst_PbPb[5] = {errors->GetBinContent(startbinerror), errors->GetBinContent(startbinerror+1), errors->GetBinContent(startbinerror+2), errors->GetBinContent(startbinerror+3), errors->GetBinContent(startbinerror+4)};

      //convert pp errors from absolute to fractional
      double syst_pp[5] = {jetspec_R040_ltb5_syst[0]/pprebin->GetBinContent(0), jetspec_R040_ltb5_syst[1]/pprebin->GetBinContent(1), jetspec_R040_ltb5_syst[2]/pprebin->GetBinContent(2), jetspec_R040_ltb5_syst[3]/pprebin->GetBinContent(3), jetspec_R040_ltb5_syst[4]/pprebin->GetBinContent(4)};

      //combine PbPb and pp errors in quadrature
  //===============================================================
      //get systematic errors from file
      TFile *errorfile;
      if (cent == true) errorfile = new TFile("../Systematics/Systematics010.root");
      if (semi == true) errorfile = new TFile("../Systematics/Systematics3050.root");
      TH1F *errors = (TH1F*)errorfile->Get("h_err"); 
      int startbinerror = errors->FindBin(40.01);
      double syst_PbPb[5] = {errors->GetBinContent(startbinerror), errors->GetBinContent(startbinerror+1), errors->GetBinContent(startbinerror+2), errors->GetBinContent(startbinerror+3), errors->GetBinContent(startbinerror+4)};

      //convert pp errors from absolute to fractional
      double syst_pp[5] = {jetspec_R040_ltb5_syst[0]/pprebin->GetBinContent(0), jetspec_R040_ltb5_syst[1]/pprebin->GetBinContent(1), jetspec_R040_ltb5_syst[2]/pprebin->GetBinContent(2), jetspec_R040_ltb5_syst[3]/pprebin->GetBinContent(3), jetspec_R040_ltb5_syst[4]/pprebin->GetBinContent(4)};

      //combine PbPb and pp errors in quadrature
      vector<double> systerrors(5);
      for (int i = 0; i < 5; i++)  {
          systerrors[i] = pow(pow(syst_PbPb[i], 2) + pow(syst_pp[i], 2), 0.5);
          cout << "systerrors: " << systerrors[i] <<"\n";
          cout << "	pp: " << syst_pp[i] <<"\n";
          cout << "	PbPb: " << syst_PbPb[i] <<"\n";
          }

      //multiply fractional errors by RAA values   
      double sys010[3] = {Raa->GetBinContent(startbin)*systerrors[2], Raa->GetBinContent(startbin+1)*systerrors[3], Raa->GetBinContent(startbin+2)*systerrors[4]};
      double sys3050[5] = {Raa->GetBinContent(startbin)*systerrors[0], Raa->GetBinContent(startbin+1)*systerrors[1], Raa->GetBinContent(startbin+2)*systerrors[2], Raa->GetBinContent(startbin+3)*systerrors[3], Raa->GetBinContent(startbin+4)*systerrors[4]};

      //ctraye TGraph Errors for systematic errors
      auto syst010 = new TGraphErrors(3, bincenters, binvalues, nul, sys010);
      auto syst3050 = new TGraphErrors(5, bincenters, binvalues, nul, sys3050);
      syst010->SetLineColor(kBlack);
      syst010->SetFillStyle(0);
      syst3050->SetLineColor(kBlack);
      syst3050->SetFillStyle(0);   
 

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
   auto gen = new TLegend(0.2, 0.7, 0.35, 0.85);
      gen->SetTextSize(0.035);
      gen->SetBorderSize(0);
        if (cent == true)    gen->AddEntry((TObject*)0, Form("ALICE Pb#font[122]{-}Pb, 5.02 TeV, %.0f#font[122]{-}%.0f%%", centl, centr), "");
        if (semi == true)    gen->AddEntry((TObject*)0, Form("ALICE Pb#font[122]{-}Pb, 5.02 TeV, %.0f#font[122]{-}%.0f%%", centl, centr), "");
        gen->AddEntry((TObject*)0, "Charged Jets, anti-#it{k}_{T}, R = 0.2, |#it{#eta}_{jet}| < 0.7", "");
        gen->AddEntry((TObject*)0, "Work in Progress", "");
   auto datadetails = new TLegend(0.55, 0.25, 0.75, 0.4);
      datadetails->SetTextSize(0.035);
      datadetails->SetBorderSize(0);
        datadetails->AddEntry((TObject*)0, "LHC18q_pass1", "");
        datadetails->AddEntry((TObject*)0, "20\% high/low q_{2}", "");

cout <<"line 454\n";

  //====================================
  //========== draw histograms =========
  //====================================
  //Draw Yield Ratio
  TCanvas *c1 = new TCanvas("Yield Ratio", lofile, 600, 500);
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
             yieldratio->SetTitle("");
             yieldratio->SetMinimum(0.01);
             yieldratio->SetMaximum(1.99);
             yieldratio->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             yieldratio->GetXaxis()->SetTitleSize(0.04);
             yieldratio->GetYaxis()->SetTitleSize(0.04);
             yieldratio->SetYTitle("high q_{2}/low q_{2}");
         yieldratio->Draw("same");
         //RaaRaw->Draw("hist same");
      //stat->Draw("e same");
     stat->Draw("ez same");
     //caitie->Draw("ezp 2 same");
    //syst3050->Draw("ezp 2 same");
    //stat->Draw("e 3 same");
     //stat->Draw("ezp same");
   //Raa->Draw("same");
   gen->Draw("same");
   datadetails->Draw("same");
   mline->Draw("same");


  TCanvas *c2 = new TCanvas("Raw Style", "raw style", 600, 500);
  c2->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad2a = new TPad("pad2a", "", 0.0, 0.4, 1.0, 1.0);
          pad2a->SetBottomMargin(0.0); // Upper and lower plot are joined
          pad2a->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad2a->Draw();
       TPad *pad2b = new TPad("pad2b", "", 0.0, 0.05, 1.0, 0.4);
          pad2b->SetBottomMargin(0.2); // Upper and lower plot are joined
          pad2b->SetTopMargin(0.0); // Upper and lower plot are joined
          pad2b->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad2b->Draw();
          pad2a->cd();               // pad1 becomes the current pad   
         //pad1a->SetLogx();
          pad2a->SetLogy();
       setcolor(lorebin, kBlue);
       setcolor(hirebin, kRed);
       hirebin->SetMaximum(hirebin->GetBinContent(2)*10);
       
       hirebin->SetTitle("");
       //hirebin->GetXaxis()->SetRangeUser(lbound, rbound);

       hirebin->Scale(1.0, "width");
       lorebin->Scale(1.0, "width");

       hirebin->Draw("hist same");
       lorebin->Draw("hist same");
 
       //histkey->Draw("same");

       c2->cd();
       pad2b->cd();
       yieldratio->SetTitle("");
       yieldratio->SetXTitle("p_{T}");
       yieldratio->GetXaxis()->SetLabelSize(0.05);
       yieldratio->GetYaxis()->SetLabelSize(0.05);
       yieldratio->GetXaxis()->SetTitleSize(0.06);
       yieldratio->SetMinimum(0.01);
       yieldratio->SetMaximum(1.99); 
       yieldratio->Draw("hist same");
       mline->Draw("same"); 


 }


void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}
 
