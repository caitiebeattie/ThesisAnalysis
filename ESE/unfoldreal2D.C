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
double lbound = 30.0;
double rbound = 120.0;

void setcolor(TH1F* h, int kcolor);

vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};


bool cent = false;
bool semi = false; 
bool r020 = false;
bool r040 = false;



//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void unfoldreal2D (int iteration, string rootfilename, bool zoom = false) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);


   //==================================================
   //======== 1. get file tag for legend 
   //==================================================
   
   //Get File tag for legend
   vector<const char*> letters(0);
   stringstream ss;
   const char* delim = "_";
   for (int i = 0; i < rootfilename.length(); i++)    {
       const char *a = &rootfilename[i];
       letters.push_back(a);
       if (strncmp(letters[i], delim, 1) == 0)    ss << " ";
       else ss << rootfilename[i];
       }
    string tags;
    vector<string> words(0);
    while (ss >> tags)  {
       words.push_back(tags);
    }  


   //========================================================
   //======  establish centrality from unfolding filename
   //========================================================  

    stringstream s2;
    const char* ccent = "010";
    const char* scent = "3050";
    const char* r02 = "R02";
    const char* r04 = "R04";
    for (int i = 2; i < words.size() - 1; i++) {
        s2 << words[i] << " ";
        const char *newword = words[i].c_str();
        if (strncmp(newword, ccent, 3) == 0)   cent = true;
        if (strncmp(newword, scent, 3) == 0)   semi = true;
        if (strncmp(newword, r02, 3) == 0)     r020 = true;
        if (strncmp(newword, r04, 3) == 0)     r040 = true;
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
  //=========  retrieve info from data file used to unfold
  //===========================================================
      //Get File
      string specfilename;
      //if (cent == true)     specfilename = "../Spectra/AnalysisResults7500.root";
      if (semi == true)     specfilename = "AnalysisResults7812.root";
      const char *lespecfile = specfilename.c_str(); 
      TFile *specfile = new TFile(lespecfile);
   
      //Get N_event
      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral();
      cout << "N_event: " << N_event <<"\n";



  //======================================================
  //=========  retrieve unfolded spectra from file
  //=======================================================
 
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);

    vector<TH1F*> unfoldhists(0);
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_1"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_2"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_3"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_4"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_5"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_6"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_7"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_8"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_9"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_10"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_11"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_12"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_13"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_14"));
    unfoldhists.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_15"));
 

   //Convert 2D hists into 1D projections
    vector<TH1F*> unfoldHistsIP(0); 
    vector<TH1F*> unfoldHistsOP(0); 
    vector<TH1F*> unfoldHists1D(0); 
    vector<TH1F*> unfoldHistsEP(0); 
    for (int i = 0; i < 15; i++) {
       stringstream ssip, ssop, ss1D, ssEP;
       ssip << "tempHistIP" << i;
       ssop << "tempHistOP" << i;
       ss1D << "tempHist1D" << i;
       ssEP << "tempHistEP" << i;
       TH2F* tempHistIP = (TH2F*)unfoldhists[i]->Clone(ssip.str().c_str());
       TH2F* tempHistOP = (TH2F*)unfoldhists[i]->Clone(ssop.str().c_str());
       TH2F* tempHist1D = (TH2F*)unfoldhists[i]->Clone(ss1D.str().c_str());
       TH2F* tempHistEP = (TH2F*)unfoldhists[i]->Clone(ssEP.str().c_str());
       tempHistIP->GetYaxis()->SetRangeUser(0.0, sqrt(2)/2.0);
       tempHistOP->GetYaxis()->SetRangeUser(sqrt(2)/2.0, 1.0);
       tempHistEP->GetXaxis()->SetRangeUser(40.0, 100.0);
       unfoldHistsIP.push_back((TH1F*)tempHistIP->ProjectionX());
       unfoldHistsOP.push_back((TH1F*)tempHistOP->ProjectionX());
       unfoldHists1D.push_back((TH1F*)tempHist1D->ProjectionX());
       unfoldHistsEP.push_back((TH1F*)tempHistEP->ProjectionY());
    }
    for (int i = 0; i < unfoldhists.size(); i++)  {
        unfoldHistsIP[i]->Scale(1.0, "width");
        unfoldHistsOP[i]->Scale(1.0, "width");
        unfoldHists1D[i]->Scale(1.0, "width");
        unfoldHists1D[i]->Scale(1.0/(double)N_event);
        unfoldHistsEP[i]->Scale(1.0, "width");
        unfoldHistsEP[i]->Scale(1.0/(double)N_event);

     }

   
        //Refolded Spectrum
        TH2D* refolded2 = (TH2D*)lefile->Get("Bayesian_Folded_6");
        TH1D* refolded = (TH1D*)refolded2->ProjectionX();
        TH1D* refoldedEP = (TH1D*)refolded2->ProjectionY();
        refolded->Scale(1.0, "width");
        refoldedEP->Scale(1.0, "width");
        refolded->Scale(1.0/(double)N_event);
        refoldedEP->Scale(1.0/(double)N_event);


        //Raw from Data
        TH2F* raw2 = (TH2F*)lefile->Get("specraw2");
        TH1F* raw = (TH1F*)raw2->ProjectionX();
        TH1F* rawEP = (TH1F*)raw2->ProjectionY();
        raw->Scale(1.0, "width");
        rawEP->Scale(1.0, "width");
        //raw->Scale(phifrac/etafrac);
        raw->Scale(1.0/(double)N_event);
        rawEP->Scale(1.0/(double)N_event);

        //Raw from Response
        TH2F* respraw2 = (TH2F*)lefile->Get("raw2");
        TH1F* respraw = (TH1F*)respraw2->ProjectionX();
        TH1F* resprawEP = (TH1F*)respraw2->ProjectionY();
        respraw->Scale(1.0, "width");
        resprawEP->Scale(1.0, "width");
        //respraw->Scale(phifrac/etafrac);
        respraw->Scale(1.0/(double)N_event);
        respraw->Scale(1.0/(double)N_event);

        //Truth from Response
        TH2F* truth2 = (TH2F*)lefile->Get("true2");
        TH1F* truth = (TH1F*)truth2->ProjectionX();
        TH1F* truthEP = (TH1F*)truth2->ProjectionY();
        truth->Scale(1.0, "width");
        truthEP->Scale(1.0, "width");
        //truth->Scale(phifrac/etafrac);
        truth->Scale(1.0/(double)N_event);
        truthEP->Scale(1.0/(double)N_event);

  //=======================================================
  //========= 5. retrieve efficiencies from file ==========
  //=======================================================
        //Get Kinematic Efficiency
        //File Low
        TH1F *kinpre1 = (TH1F*)lefile->Get("KinPre");
        TH1F *kinpos1 = (TH1F*)lefile->Get("KinPos");
             TH1F *effic1 = (TH1F*)kinpos1->Clone("effic1");
             effic1->Divide(kinpre1);

       //Get Reconstruction Efficiency
       TFile *recfile = new TFile("../Unfolding/ChargedJetRecEfficiencies.root");
       TH1D  *receffic = (TH1D*)recfile->Get("RecEff_R020_5GeV");
 
  
  //=======================================================
  //========= 6. Get Ratios ===============================
  //=======================================================

    TH1F* closure = (TH1F*)refolded->Clone("closure");
    closure->Divide(raw);

    TH1F* closureEP = (TH1F*)refoldedEP->Clone("closureEP");
    closureEP->Divide(rawEP);


    vector<TH1F*> ratiohists(0);
    vector<TH1F*> ratiohistsEP(0);
    for (int i = 0; i < unfoldHists1D.size(); i++)   {
      stringstream nameep;
      nameep << names[i] << "_ep";
      
      ratiohists.push_back((TH1F*)unfoldHists1D[i]->Clone(names[i]));
      ratiohists[i]->Divide(unfoldHists1D[iteration]);
      ratiohistsEP.push_back((TH1F*)unfoldHistsEP[i]->Clone(nameep.str().c_str()));
      ratiohistsEP[i]->Divide(unfoldHistsEP[iteration]);
     // ratiohists[i]->GetXaxis()->SetRangeUser(lbound, rbound); 
      }

vector<int> colorlist = {kRed-7, kOrange+1, kOrange-2, kYellow-7, kSpring+5, kGreen-6, kGreen-2, kTeal+2, kTeal-5, kCyan-3, kCyan+2, kAzure-9, kAzure-4, kBlue-4, kBlue-6, kViolet+6, kViolet-8, kMagenta-8};
for (int i = 0; i < unfoldhists.size(); i++)    {
     setcolor(unfoldHists1D[i], colorlist[i]);
     setcolor(ratiohists[i], colorlist[i]);
     setcolor(unfoldHistsEP[i], colorlist[i]);
     setcolor(ratiohistsEP[i], colorlist[i]);}
     unfoldHists1D[iteration]->SetLineWidth(3);
  
  TLine *line = new TLine(lbound, 1, rbound, 1);
  line->SetLineStyle(2);
  TLine *mline = new TLine(15.0, 1, 100.0, 1);
  mline->SetLineStyle(2);




  
  //======================================
  //========== create legends ============
  //======================================

    auto fileleg = new TLegend(0.4, 0.2, 0.85, 0.45);
         fileleg->SetTextSize(0.055);
         fileleg->SetBorderSize(0);
         fileleg->AddEntry((TObject*)0, s2.str().c_str(), "");

    auto itleg = new TLegend(0.6, 0.3, 0.85, 0.85);
         itleg->SetTextSize(0.045);
         itleg->SetBorderSize(0); 
         for (int i = 0; i < unfoldHists1D.size(); i++)    itleg->AddEntry(unfoldHists1D[i], Form("%d", i+1), "pl");

    auto stab = new TLegend(0.15, 0.65, 0.4, 0.85);
         stab->SetTextSize(0.05);
         stab->SetBorderSize(0); 
         stab->AddEntry((TObject*)0, "Stability", "");
         if (cent == true)  stab->AddEntry((TObject*)0, "R = 0.2, 0-10%", "");
         if (semi == true)  stab->AddEntry((TObject*)0, "R = 0.2, 30-50%", "");

    auto split = new TLegend(0.15, 0.65, 0.4, 0.85);
         split->SetTextSize(0.05);
         split->SetBorderSize(0); 
         if (cent == true)  split->AddEntry((TObject*)0, "R = 0.2, 0-10%", "");
         if (semi == true)  split->AddEntry((TObject*)0, "R = 0.2, 30-50%", "");

   auto *close = new TLegend(0.5, 0.6, 0.85, 0.85);
        close->SetBorderSize(0);
        close->SetTextSize(0.04);
        close->AddEntry(unfoldHists1D[iteration], "unfolded", "pl");
        close->AddEntry(raw, "raw", "pl");
        close->AddEntry(truth, "truth (from embedding)", "pl");
        close->AddEntry(respraw, "hybrid (from embedding)", "pl");

   auto *refo = new TLegend(0.6, 0.75, 0.85, 0.85);
        refo->SetBorderSize(0);
        refo->SetTextSize(0.055);
        refo->AddEntry((TObject*)0, "Refolded/Raw", "");


  //==================================================
  //========== determine viewing properties  =========
  //==================================================

    float ymin, ymax;
    if (zoom == false)  {
        ymin = 0.5;
        ymax = 1.5;}
    if (zoom == true)   {
      ymin = 0.8;
      ymax = 1.2;}



  //====================================
  //========== draw histograms =========
  //====================================


  TCanvas *c1= new TCanvas("Stability Test", lerootfile, 1200, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1a = new TPad("pad1a", "", 0.0, 0.4, 0.5, 1.0);
   TPad *pad1b = new TPad("pad1b", "", 0.0, 0.05, 0.5, 0.4);
     pad1a->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad1a->SetRightMargin(0.0);
     pad1b->SetTopMargin(0.0);
     pad1b->SetRightMargin(0.0);
     pad1b->SetBottomMargin(0.2);
     pad1a->Draw();
     pad1b->Draw();
   pad1a->cd();               // pad1 becomes the current pad   
     //pad1->SetLogx();
     pad1a->SetLogy();
     double topmax = unfoldHists1D[0]->GetBinContent(3)*pow(10,3);
     double topmin = unfoldHists1D[0]->GetBinContent(9)*pow(10,-2);
     unfoldHists1D[0]->SetMaximum(topmax);
     unfoldHists1D[0]->SetMinimum(topmin);
     unfoldHists1D[0]->SetTitle("");
     unfoldHists1D[0]->GetXaxis()->SetRangeUser(30.0, 120.0);
     for (int i = 0; i < unfoldHists1D.size(); i++)         unfoldHists1D[i]->Draw("same");
     itleg->Draw("same");
     c1->cd();
   pad1b->cd();
     ratiohists[0]->GetXaxis()->SetRangeUser(30.0, 119.9);
     ratiohists[0]->GetXaxis()->SetLabelSize(0.07);
     ratiohists[0]->GetYaxis()->SetLabelSize(0.06);
     ratiohists[0]->SetTitle("");
     ratiohists[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
     ratiohists[0]->GetXaxis()->SetTitleSize(0.09);
     ratiohists[0]->SetMinimum(0.01);
     ratiohists[0]->SetMaximum(1.99);
     for (int i = 0; i < ratiohists.size(); i++)         ratiohists[i]->Draw("same");
     line->Draw("same");  
     fileleg->Draw("same");
     c1->cd();
   TPad *pad1c = new TPad("pad1c", "", 0.5, 0.4, 1.0, 1.0);
        pad1c->SetBottomMargin(0.0);
        pad1c->SetLeftMargin(0.0);
        pad1c->Draw();
        pad1c->cd();
        pad1c->SetLogy();
         unfoldHists1D[iteration]->SetMaximum(topmax);
         unfoldHists1D[iteration]->SetMinimum(topmin);
         unfoldHists1D[iteration]->SetTitle("");
            setcolor(unfoldHists1D[iteration], kGreen+2);
         unfoldHists1D[iteration]->GetXaxis()->SetRangeUser(30.0, 120.0);
         unfoldHists1D[iteration]->Draw("same");
            setcolor(raw, kRed+2);
         raw->Draw("same");
            setcolor(truth, kBlue+2);
         truth->Draw("same");
            setcolor(respraw, kViolet);	
         respraw->Draw("same");
       close->Draw("same");
      c1->cd();
   TPad *pad1d = new TPad("pad1d", "", 0.5, 0.05, 1.0, 0.4);
        pad1d->SetTopMargin(0.0);
        pad1d->SetBottomMargin(0.2);
        pad1d->Draw();
        pad1d->SetLeftMargin(0.0);
      pad1d->cd();
        closure->SetMinimum(0.01);
        closure->SetMaximum(1.99);
      setcolor(closure, kBlack);
      closure->GetXaxis()->SetRangeUser(30.01, 120.0);
      closure->Draw("same");      
      refo->Draw("same");
     closure->SetTitle("");
         closure->SetXTitle("#it{p}_{T} (GeV/#it{c})");
         closure->GetXaxis()->SetTitleSize(0.09);
         closure->GetXaxis()->SetLabelSize(0.07);
         closure->GetYaxis()->SetLabelSize(0.06);
       line->Draw("same");




  TCanvas *c2= new TCanvas("Stability EP", "Stability EP", 1200, 500);
   c2->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad2a = new TPad("pad2a", "", 0.0, 0.4, 0.5, 1.0);
   TPad *pad2b = new TPad("pad2b", "", 0.0, 0.05, 0.5, 0.4);
     pad2a->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad2a->SetRightMargin(0.0);
     pad2b->SetTopMargin(0.0);
     pad2b->SetRightMargin(0.0);
     pad2b->SetBottomMargin(0.2);
     pad2a->Draw();
     pad2b->Draw();
   pad2a->cd();               // pad1 becomes the current pad   
   //pad1->SetLogx();
   pad2a->SetLogy();
   //topmax = unfoldHistsEP[0]->GetBinContent(3)*pow(10,3);
   //topmin = unfoldHistsEP[0]->GetBinContent(9)*pow(10,-2);
   unfoldHistsEP[0]->SetMaximum(topmax);
   unfoldHistsEP[0]->SetMinimum(topmin);
   unfoldHistsEP[0]->SetTitle("");
   //unfoldHistsEP[0]->GetXaxis()->SetRangeUser(30.0, 120.0);
   for (int i = 0; i < unfoldHistsEP.size(); i++)         unfoldHistsEP[i]->Draw("same");
   itleg->Draw("same");
   pad2b->cd();
   ratiohistsEP[0]->GetXaxis()->SetLabelSize(0.07);
   ratiohistsEP[0]->GetYaxis()->SetLabelSize(0.06);
   ratiohistsEP[0]->SetTitle("");
   ratiohistsEP[0]->SetXTitle("|#Delta#varphi|");
   ratiohistsEP[0]->GetXaxis()->SetTitleSize(0.09);
   ratiohistsEP[0]->SetMinimum(0.01);
   ratiohistsEP[0]->SetMaximum(1.99);
   for (int i = 0; i < ratiohists.size(); i++)         ratiohistsEP[i]->Draw("same");
   line->Draw("same");  
   fileleg->Draw("same");
   c2->cd();
   TPad *pad2c = new TPad("pad2c", "", 0.5, 0.4, 1.0, 1.0);
        pad2c->SetBottomMargin(0.0);
        pad2c->SetLeftMargin(0.0);
        pad2c->Draw();
      pad2c->cd();
      pad2c->SetLogy();
         unfoldHistsEP[iteration]->SetMaximum(topmax);
         unfoldHistsEP[iteration]->SetMinimum(topmin);
         unfoldHistsEP[iteration]->SetTitle("");
            setcolor(unfoldHistsEP[iteration], kGreen+2);
         unfoldHistsEP[iteration]->Draw("same");
            setcolor(rawEP, kRed+2);
         rawEP->Draw("same");
            setcolor(truthEP, kBlue+2);
         truthEP->Draw("same");
            setcolor(resprawEP, kViolet);	
         resprawEP->Draw("same");
       close->Draw("same");
      c2->cd();
   TPad *pad2d = new TPad("pad2d", "", 0.5, 0.05, 1.0, 0.4);
        pad2d->SetTopMargin(0.0);
        pad2d->SetBottomMargin(0.2);
        pad2d->Draw();
        pad2d->SetLeftMargin(0.0);
      pad2d->cd();
        closureEP->SetMinimum(0.01);
        closureEP->SetMaximum(1.99);
      setcolor(closureEP, kBlack);
      closureEP->Draw("same");      
      refo->Draw("same");
     closureEP->SetTitle("");
         closureEP->SetXTitle("|#Delta#varphi|");
         closureEP->GetXaxis()->SetTitleSize(0.09);
         closureEP->GetXaxis()->SetLabelSize(0.07);
         closureEP->GetYaxis()->SetLabelSize(0.06);
       line->Draw("same");
  


  TCanvas *c3= new TCanvas("Reconstruction Efficiency", "Reconstruction Efficiency", 500, 500);
   c3->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad3 = new TPad("pad3", "", 0.0, 0.05, 1.0, 1.0);
     pad3->SetBottomMargin(0.15); // Upper and lower plot are joined
     pad3->Draw();
     pad3->Draw();
   pad3->cd();               // pad1 becomes the current pad   
   receffic->SetTitle("");
   receffic->SetMinimum(0.01);
   receffic->SetMaximum(0.99);
   receffic->GetXaxis()->SetRangeUser(10.0, 120.0);
   receffic->Draw("same");

  TCanvas *c4= new TCanvas("Kinematic Efficiency", "Kinematic Efficiency", 500, 500);
   c4->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
     pad4->SetBottomMargin(0.15); // Upper and lower plot are joined
     pad4->Draw();
     pad4->Draw();
   pad4->cd();               // pad1 becomes the current pad   
   effic1->SetTitle("");
   effic1->SetMinimum(0.01);
   effic1->SetMaximum(0.99);
   effic1->GetXaxis()->SetRangeUser(10.0, 120.0);
   effic1->SetXTitle("p_{T} (GeV/#it{c})");
   effic1->Draw("same");



/*

  TCanvas *c1 = new TCanvas("Stability Test", lerootfile, 800, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad1 = new TPad("pad1", "", 0.0, 0.4, 0.5, 1.0);
   TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 0.5, 0.4);
     pad1->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad1->SetRightMargin(0.0);
     pad2->SetTopMargin(0.0);
     pad2->SetRightMargin(0.0);
     pad2->SetBottomMargin(0.2);
     pad1->Draw();
     pad2->Draw();
   pad1->cd();               // pad1 becomes the current pad   
   //pad1->SetLogx();
   pad1->SetLogy();
   double topmax = unfoldHists1D[0]->GetBinContent(3)*pow(10,3);
   double topmin = unfoldHists1D[0]->GetBinContent(9)*pow(10,-2);
   unfoldHists1D[0]->SetMaximum(topmax);
   unfoldHists1D[0]->SetMinimum(topmin);
   //unfoldHists1D[0]->SetMinimum(pow(10,-9));
   unfoldHists1D[0]->SetTitle("");
   for (int i = 0; i < 15; i++)         unfoldHists1D[i]->Draw("same");
   itleg->Draw("same");
   stab->Draw("same");
   pad2->cd();
   ratiohists[0]->GetXaxis()->SetLabelSize(0.05);
   ratiohists[0]->GetYaxis()->SetLabelSize(0.05);
   ratiohists[0]->SetTitle("");
   ratiohists[0]->SetXTitle("#it{p}_{T} (GeV/#it{c})");
   ratiohists[0]->GetXaxis()->SetTitleSize(0.09);
   ratiohists[0]->SetMinimum(ymin);
   ratiohists[0]->SetMaximum(ymax);
   for (int i = 0; i < unfoldhists.size(); i++)         ratiohists[i]->Draw("same");
   line->Draw("same");  

  //TCanvas *c3 = new TCanvas("Split Test", "Split Test", 500, 500);
   //c3->cd();
   gStyle->SetOptStat(0);
   c1->cd();
   TPad *pad3 = new TPad("pad3", "", 0.5, 0.4, 1.0, 1.0);
   TPad *pad4 = new TPad("pad4", "", 0.5, 0.05, 1.0, 0.4);
   TLegend *leg = new TLegend(0.6, 0.6, 0.85, 0.85);
     setcolor(unfoldHists1D[iteration], kBlue);
     setcolor(truth, kGreen);
     setcolor(raw, kRed);
     setcolor(trivial, kBlack);
   leg->AddEntry(unfoldHists1D[iteration], "unfolded", "pl");
   leg->AddEntry(truth, "truth", "pl");
   leg->AddEntry(raw, "raw", "pl");
   leg->SetBorderSize(0);
     pad3->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad4->SetTopMargin(0.0);
     pad4->SetBottomMargin(0.2);
     pad3->SetLeftMargin(0.0);
     pad4->SetLeftMargin(0.0);
     pad3->Draw();
     pad4->Draw();
   pad3->cd();               // pad1 becomes the current pad   
     //pad1->SetLogx();
     pad3->SetLogy();
     unfoldHists1D[iteration]->GetXaxis()->SetRangeUser(lbound, rbound);
     //unfoldHists1D[iteration]->SetMaximum(pow(10,-1));
     //unfoldHists1D[iteration]->SetMinimum(pow(10,-9));
     unfoldHists1D[iteration]->SetTitle("");
     unfoldHists1D[iteration]->SetMaximum(topmax);
     unfoldHists1D[iteration]->SetMinimum(topmin);
     unfoldHists1D[iteration]->Draw("same");      
     truth->GetXaxis()->SetRangeUser(lbound, rbound);
     truth->Draw("same");
     raw->GetXaxis()->SetRangeUser(lbound, rbound);
     raw->Draw("same");
     leg->Draw("same"); 
     split->Draw("same"); 
   pad4->cd();
     trivial->SetTitle("");
     trivial->GetXaxis()->SetRangeUser(lbound, rbound);
     trivial->SetMinimum(ymin);
     trivial->SetMaximum(ymax);
     trivial->Draw("same");
     trivial->SetXTitle("#it{p}_{T} (GeV/#it{c})");
     trivial->GetXaxis()->SetTitleSize(0.09);
     trivial->GetXaxis()->SetLabelSize(0.05);
     trivial->GetYaxis()->SetLabelSize(0.05);
     line->Draw("same");
     fileleg->Draw("same");
*/

}


void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}

