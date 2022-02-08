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
double rbound = 140.0;

void setcolor(TH1F* h, int kcolor);
void updateR (double &phiscale, double &etascale, double jetR) {
    phiscale = (2*pi)/(1.92 - 2*jetR);
    etascale = (1.4 - 2*jetR);
}

bool cent = false;
bool semi = false;

int centcolor = kViolet+5;
int semicolor = kCyan-4;

std::vector<double> kBinsUnfolded = {10.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 140.0};

//TH1F* q2hist = new TH1F("q2hist", "q2hist", 100, 0, 20);
TH1F* lowq2 = new TH1F("lowq2", "lowq2", kBinsUnfolded.size()-1, kBinsUnfolded.data());
TH1F* highq2 = new TH1F("highq2", "highq2", kBinsUnfolded.size()-1, kBinsUnfolded.data());

//TH1F* lowqcent = new TH1F("lowqcent", "lowqcent", 10, 30.0, 50.0);
//TH1F* highqcent = new TH1F("highqcent", "highqcent", 10, 30.0, 50.0);

//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void rawcomp (string rootfilename) {

  const char* lerootfile = rootfilename.c_str();
  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

   //==================================================
   //======== 1. get file tag for legend 
   //==================================================
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = "_";
   for (int i = 0; i < rootfilename.length(); i++)    {
       const char *a = &rootfilename[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
       else ss << rootfilename[i];
       }

  
   //========================================================
   //====== 2. establish centrality from unfolding filename
   //========================================================  
    string tags;
       vector<string> words(0);
       while (ss >> tags)   words.push_back(tags);
       stringstream s2;
    const char* ccent = "010";
    const char* scent = "3050";
      for (int i = 2; i < words.size() - 1; i++) {
        s2 << words[i] << " ";
        const char *newword = words[i].c_str();
        if (strncmp(newword, ccent, 3) == 0)   cent = true;
        if (strncmp(newword, scent, 3) == 0)   semi = true;
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

  //===========================================================
  //========= 3. retrieve info from data file
  //===========================================================
      //GetFile
      TFile *lefile = new TFile(lerootfile);
   
      centl = 30.0;
      centr = 50.0;
      //Get N_event
      TDirectoryFile* tdf = (TDirectoryFile*)lefile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral();
      cout << "N_event: " << N_event <<"\n";

      //Calculate q2 percentiles
      TTree *T = nullptr;
      lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
      TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_V2_TPCEta0dot8");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrFullV0");
      float pTjet, centra, q2;
      T->SetBranchAddress("Jet_Pt", &pTjet);
      T->SetBranchAddress("Event_Centrality", &centra);
      T->SetBranchAddress("Event_Q2Vector", &q2);
      //Fill q2 hist
      int entries = T->GetEntries();
      cout << "q2hist2D->GetEntries(): " << q2hist2D->GetEntries() <<"\n";
      TH1F *q2all = (TH1F*)q2hist2D->ProjectionY();
      q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
      TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();
      q2hist->Rebin(100);
      //for (int i = 0; i < entries; i++) {
      //  T->GetEntry(i);
      //  if (centra < centl || centra > centr)   continue;
      //  q2hist->Fill(q2);
      //}
      cout << "q2hist->GetEntries(): " << q2hist->GetEntries() <<"\n";
      long int percentileticker = 0;                              //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                   //number percentile we're on
      cout << "entries (#jets in tree): " << entries <<"\n";
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

      for (int i = 0; i < percentiles.size(); i++)   cout << "Percentile: " << i << ", q2: " << percentiles[i] << "\n";
      TH2F* loq = (TH2F*)q2hist2D->Clone("loq");
      TH2F* hiq = (TH2F*)q2hist2D->Clone("hiq");
      loq->GetYaxis()->SetRangeUser(0.0,percentiles[3]);
      hiq->GetYaxis()->SetRangeUser(percentiles[5], 15.0);
      TH1F* lowqcent = (TH1F*)loq->ProjectionX();
      TH1F* highqcent = (TH1F*)hiq->ProjectionX();

  //===========================================================
  //========= 3. plot initial spectra
  //===========================================================
     long int low = 0;
     long int high = 0;

     for (int i = 0; i < entries; i++)   {
        T->GetEntry(i);
        if (centra < centl || centra > centr)   continue; 
        if (q2 < percentiles[0]) {
                lowq2->Fill(pTjet);
                }
        if (q2 > percentiles[8]) {
                highq2->Fill(pTjet); 
                }
     }

      cout << "low: " << low <<"\n";
      cout << "high: " << high <<"\n";

      lowq2->Scale(1.0, "width");
      highq2->Scale(1.0, "width");

      TH1F* rat = (TH1F*)highq2->Clone("rat");
      rat->Divide(lowq2);
      rat->GetXaxis()->SetRangeUser(lbound, rbound);
/*
  //============================================================
  //========= 6. apply scaling for bin width + N_events ========
  //============================================================
        //Unfolded Spectra 
        for (int i = 0; i < unfoldhists.size(); i++)    {
        unfoldhists[i]->Scale(1.0, "width");
        //unfoldhists[i]->Scale(phifrac/etafrac);
        unfoldhists[i]->Scale(1.0/(double)N_event);}

        //Refolded Spectrum
        TH1F* refolded = (TH1F*)lefile->Get("Bayesian_Folded_6");
        refolded->Scale(1.0, "width");
        refolded->Scale(1.0/(double)N_event);

        //Raw from Data
        TH1F* raw = (TH1F*)lefile->Get("specraw");
        raw->Scale(1.0, "width");
        //raw->Scale(phifrac/etafrac);
        raw->Scale(1.0/(double)N_event);

        //Raw from Response
        TH1F* respraw = (TH1F*)lefile->Get("raw");
        respraw->Scale(1.0, "width");
        //respraw->Scale(phifrac/etafrac);
        respraw->Scale(1.0/(double)N_event);

        //Truth from Response
        TH1F* truth = (TH1F*)lefile->Get("true");
        truth->Scale(1.0, "width");
        //truth->Scale(phifrac/etafrac);
        truth->Scale(1.0/(double)N_event);



 
TH1F* closure = (TH1F*)refolded->Clone("closure");
   closure->Divide(raw);





 
  //=============================================
  //============= Calculate the RAA =============
  //=============================================
      //divide PbPb by pp
      TH1F* RaaLo = (TH1F*)lowq2->Clone("RaaLo");
      RaaLo->Divide(pprebin);

      //apply efficiencies
      RaaLo->Divide(efficiency);
      RaaLo->Divide(recefficiency);

      //select and apply TAA
      double TAA010 = 23.3;
      double TAA3050 = 3.9;
      double TAA8090 = 0.084;
      double TAA = 1.0;
      if (cent == true)   TAA = TAA010;
      if (semi == true)   TAA = TAA3050;
      Raa->Scale(1.0/TAA);

      //set aesthetics
      Raa->SetLineColorAlpha(kWhite, 0.0);
      Raa->SetMarkerStyle(20);
      Raa->SetMarkerColorAlpha(kBlack, 0.0);
      Raa->SetMarkerSize(0.7);
      //lower edge based on unfolding stability/kinematic efficiency 
      double lowedge;
      if (cent == true)  lowedge = 60.0;
      if (semi == true)  lowedge = 40.0;
      int startbin = Raa->FindBin(lowedge+0.01);
      Raa->GetXaxis()->SetRangeUser(lowedge, 100.0);

      //calculate RAA pre-unfolding
      TH1F* RaaRaw = (TH1F*)rawrebin->Clone("RaaRaw");
      RaaRaw->Divide(pprebin);
      RaaRaw->Scale(1.0/TAA);
      RaaRaw->GetXaxis()->SetRangeUser(20.0, 100.0);
      RaaRaw->SetLineColorAlpha(kWhite, 0.0);
      RaaRaw->SetMarkerColorAlpha(kWhite, 0.0);

*/


  //====================================
  //========== create lines ============
  //====================================
  TLine *line = new TLine(lbound, 1, rbound, 1);
  line->SetLineStyle(2);
  TLine *sline = new TLine(40.0, 1, 100.0, 1);
  sline->SetLineStyle(2);
  TLine *mline = new TLine(15.0, 1, 100.0, 1);
  mline->SetLineStyle(2);

  TLine *loline = new TLine(percentiles[0], 1.0, percentiles[0], pow(10,7));
  loline->SetLineStyle(2);
  TLine *hiline = new TLine(percentiles[8], 1.0, percentiles[8], pow(10,7));
  hiline->SetLineStyle(2);




  //======================================
  //========== create legends ============
  //======================================
   auto fileleg = new TLegend(0.6, 0.8, 0.85, 0.85);
     fileleg->SetTextSize(0.055);
     fileleg->SetBorderSize(0);
     fileleg->AddEntry((TObject*)0, s2.str().c_str(), "");
   auto histkey = new TLegend(0.6, 0.675, 0.85, 0.85);
     histkey->SetTextSize(0.04);
     histkey->SetBorderSize(0); 
     histkey->AddEntry((TObject*)0, "LHC18q pass1, 40-50%", "");
     histkey->AddEntry(lowq2, "smallest 10\% q_{2}", "pl");
     histkey->AddEntry(highq2, "largest 10\% q_{2}", "pl");
/*

  //====================================
  //========== TAA uncertainty =========
  //====================================
  double systTAA = 0.0165;
  auto TAAsyst = new TGraphErrors(4, hcenters, hvalues010, hnul, herr010);
     TAAsyst->SetFillColor(kGray+2);
     TAAsyst->SetLineColor(kGray+2);
  TBox* uncertaintyTAA = new TBox(96, 1-systTAA , 100, 1+systTAA);
    uncertaintyTAA->SetFillStyle(1001);
    uncertaintyTAA->SetFillColor(kGray+2);
   auto TAAunc = new TLegend(0.15, 0.425, 0.35, 0.475);
     TAAunc->SetTextSize(0.03);
     TAAunc->SetBorderSize(0);
     TAAunc->AddEntry(TAAsyst, "T_{AA} normalization uncertainty", "plf"); 
 

*/

  //====================================
  //========== draw histograms =========
  //====================================
  //Draw RAA
  TCanvas *c1 = new TCanvas("Raw Spectra", lerootfile, 600, 500);
  c1->cd();
       gStyle->SetOptStat(0);
       gStyle->SetPadTickX(1);
       gStyle->SetPadTickY(1);
       TPad *pad1a = new TPad("pad1a", "", 0.0, 0.4, 1.0, 1.0);
          pad1a->SetBottomMargin(0.0); // Upper and lower plot are joined
          pad1a->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad1a->Draw();
       TPad *pad1b = new TPad("pad1b", "", 0.0, 0.05, 1.0, 0.4);
          pad1b->SetBottomMargin(0.2); // Upper and lower plot are joined
          pad1b->SetTopMargin(0.0); // Upper and lower plot are joined
          pad1b->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad1b->Draw();
          pad1a->cd();               // pad1 becomes the current pad   
         //pad1a->SetLogx();
          pad1a->SetLogy();
       setcolor(lowq2, kBlue);
       setcolor(highq2, kRed);
       highq2->SetMaximum(highq2->GetBinContent(2)*10);
       
       highq2->SetTitle("");
       highq2->GetXaxis()->SetRangeUser(lbound, rbound);

       highq2->Draw("same e");
       lowq2->Draw("same e");
 
       histkey->Draw("same");

       c1->cd();
       pad1b->cd();
       rat->SetTitle("");
       rat->SetXTitle("p_{T}");
       rat->GetXaxis()->SetLabelSize(0.08);
       rat->GetYaxis()->SetLabelSize(0.08);
       rat->GetXaxis()->SetTitleSize(0.09);
       rat->SetMinimum(0.01);
       rat->SetMaximum(3.0); 
       rat->Draw("same");
       line->Draw("same"); 

  /*      //unfoldhists[0]->SetMaximum(pow(10,-1));
         //unfoldhists[0]->SetTitle("");
         //for (int i = 0; i < 15; i++)         unfoldhists[i]->Draw("same");
             RaaRaw->SetTitle("");
             RaaRaw->SetMinimum(0.01);
             RaaRaw->SetMaximum(1.99);
             RaaRaw->SetXTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
             RaaRaw->GetXaxis()->SetTitleSize(0.05);
             RaaRaw->GetYaxis()->SetTitleSize(0.05);
             RaaRaw->SetYTitle("#it{R}_{AA}");
         RaaRaw->Draw("same");
         Raa->Draw("same");
         //RaaRaw->Draw("hist same");
     stat->Draw("ez same");
     caitie->Draw("ezp 2 same");
    //syst3050->Draw("ezp 2 same");
    //stat->Draw("e 3 same");
     //if (cent == true)     hsyst010->Draw("ezp 2 same");
     //if (semi == true)     hsyst3050->Draw("ezp 2 same");
     //if (cent == true)     mlsyst010->Draw("ezp 2 same");
     //if (semi == true)     mlsyst3050->Draw("ezp 2 same");
     //stat->Draw("ezp same");
   //Raa->Draw("same");
   if (ESE == true)  {
     linFit_loL->Draw("same");
     linFit_loL2->Draw("same");
     linFit_hiL->Draw("same");
     linFit_hiL2->Draw("same");}
   
   uncertaintyTAA->Draw("same");
   TAAunc->Draw("same");
   histkey->Draw("same");
   gen->Draw("same");
   mline->Draw("same");
  */
  //Draw q2 Hist
  TCanvas *c2 = new TCanvas("q2 Dist", "", 500, 500);
   c2->cd();
   gStyle->SetOptStat(0);
   TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
     pad2->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad2->SetLogy();
     pad2->Draw();
   pad2->cd();               // pad1 becomes the current pad
   q2hist->SetXTitle("q_{2}");
   q2hist->SetTitle("");
   q2hist->SetMinimum(1.0);
   q2hist->SetMaximum(pow(10,7));

    setcolor(q2hist, kRed);
    setcolor(q2all, kBlue);

   q2hist->Draw("same");
   q2all->Draw("same");
   loline->Draw("same");
   hiline->Draw("same");


  TCanvas *c3 = new TCanvas("LHC18q_pass1", lerootfile, 500, 500);
   c3->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad3 = new TPad("pad3", "", 0.0, 0.05, 1.0, 1.0);
     pad3->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad3->SetRightMargin(0.15);
     pad3->Draw();
   pad3->cd();               // pad1 becomes the current pad   
   q2hist2D->SetXTitle("centrality");
   q2hist2D->SetYTitle("q_{2}");
   q2hist2D->Draw("same colz");

  TCanvas *c4 = new TCanvas("Centrality Bias", "Centrality Bias", 500, 500);
   c4->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
     pad4->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad4->SetRightMargin(0.15);
     pad4->Draw();
     //pad4->SetLogy();
     pad4->cd();
     setcolor(lowqcent, kBlue);
     setcolor(highqcent, kRed);
     highqcent->SetMinimum(0.01);
     highqcent->Draw("same");
     lowqcent->Draw("same");

   auto lemean = new TLegend(0.45, 0.2, 0.8, 0.35);
     lemean->SetTextSize(0.045);
     lemean->SetBorderSize(0);
     lemean->AddEntry(highqcent, Form("High Mean: %.2f", highqcent->GetMean()), "pl");
     lemean->AddEntry(lowqcent, Form("Low Mean: %.2f", lowqcent->GetMean()), "pl");
     lemean->Draw("same");

   //pad1->SetLogx();
/*

   c4->cd();
   TPad *pad6 = new TPad("pad6", "", 0.5, 0.4, 1.0, 1.0);
        pad6->SetBottomMargin(0.0);
        pad6->SetLeftMargin(0.0);
        pad6->Draw();
      pad6->cd();
      pad6->SetLogy();
         unfoldhists[iteration]->SetMaximum(pow(10, -1));
         unfoldhists[iteration]->SetMinimum(pow(10, -9));
         unfoldhists[iteration]->SetTitle("");
            setcolor(unfoldhists[iteration], kGreen+2);
         unfoldhists[iteration]->GetXaxis()->SetRangeUser(25.0, 120.0);
         unfoldhists[iteration]->Draw("same");
            setcolor(raw, kRed+2);
         raw->Draw("same");
            setcolor(truth, kBlue+2);
         truth->Draw("same");
            setcolor(respraw, kViolet);	
         respraw->Draw("same");
       close->Draw("same");
      c4->cd();

   TPad *pad7 = new TPad("pad7", "", 0.5, 0.05, 1.0, 0.4);
        pad7->SetTopMargin(0.0);
        pad7->SetBottomMargin(0.2);
        pad7->Draw();
        pad7->SetLeftMargin(0.0);
      pad7->cd();
        closure->SetMinimum(0.5);
        closure->SetMaximum(1.5);
      setcolor(closure, kBlack);
      closure->Draw("same");      
      refo->Draw("same");
     closure->SetTitle("");
         closure->SetXTitle("p_{T} [GeV/c]");
         closure->GetXaxis()->SetTitleSize(0.09);
         closure->GetXaxis()->SetLabelSize(0.07);
         closure->GetYaxis()->SetLabelSize(0.06);
       mline->Draw("same");

  //Draw Reconstruction Efficiency
  TCanvas *c9 = new TCanvas("c9", "rec effic rebinned", 500, 500);
  c9->cd();
    TPad *pad9 = new TPad("pad9", "pad9", 0.0, 0.05, 1.0, 1.0);
    pad9->Draw();
    pad9->cd();
      recefficiency->SetMinimum(0.01);
      recefficiency->SetMaximum(0.99);
      recefficiency->Draw("hist");
   auto rec = new TLegend(0.45, 0.3, 0.8, 0.55);
     rec->SetBorderSize(0);
     rec->SetTextSize(0.035);
     rec->AddEntry((TObject*)0, "Reconstruction Efficiency", "");
     rec->Draw("same");
*/


}

void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}

