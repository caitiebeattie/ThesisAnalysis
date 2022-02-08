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

float q2cutlo = 2.45;
float q2cuthi = 5.85;

vector<const char*> names = {"i1", "i2", "i3", "i4", "i5", "i6", "i7", "i8", "i9", "i10", "i11", "i12","i13", "i14", "i15"};



//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void unfoldreal (int iteration, string rootfilename, int qsel, int R = 4) {

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

   semi = true;
   centl = 40.0;
   centr = 45.0;
  //==============================================================
  //========= 3. retrieve info from data file used to unfold
  //==============================================================
      //Get File
      string specfilename;
      if (cent == true)     specfilename = "AnalysisResults7696_2.root";
      if (semi == true)     specfilename = "AnalysisResults7665_2.root";
      const char *lespecfile = specfilename.c_str(); 
      TFile *specfile = new TFile(lespecfile);
   
      //Get N_event
      TDirectoryFile* tdf = (TDirectoryFile*)specfile->Get("ChargedJetsHadronCF");
      AliEmcalList *tlist;
      if (R == 2) tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      if (R == 4) tlist = (AliEmcalList*)tdf->Get("AliAnalysisTaskJetExtractor_Jet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_Jet_histos"); 
      TH1F *centhist = (TH1F*)tlist->FindObject("fHistCentrality");
      centhist->GetXaxis()->SetRangeUser(centl, centr);
      double N_event = centhist->Integral()/5.0;
      cout << "N_event: " << N_event <<"\n";


      //Get Statistical Uncertainty
      TTree *T = nullptr;
      if (R == 2) specfile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
      if (R == 4) specfile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_Jet", T);
      float pTjet, centra, q2;
      T->SetBranchAddress("Jet_Pt", &pTjet);
      T->SetBranchAddress("Event_Centrality", &centra);
      T->SetBranchAddress("Event_Q2Vector", &q2);
      vector<double> stat_PbPb(7);
      for (int i = 0; i < 7; i++)   stat_PbPb[i] = 0;
      int entries = T->GetEntries();
      for (int i = 0; i < entries; i++) {
        T->GetEntry(i);
        if (centra < centl || centra > centr)   continue;
        if (qsel == 0)      {if (q2 < q2cutlo)   continue;};
        if (qsel == 1)      {if (q2 > q2cuthi)   continue;};
        if (pTjet > 25.0 && pTjet < 30.0 )   stat_PbPb[0]++;
        if (pTjet > 30.0 && pTjet < 40.0 )   stat_PbPb[1]++;
        if (pTjet > 40.0 && pTjet < 50.0 )   stat_PbPb[2]++;
        if (pTjet > 50.0 && pTjet < 60.0 )   stat_PbPb[3]++;
        if (pTjet > 60.0 && pTjet < 70.0 )   stat_PbPb[4]++;
        if (pTjet > 70.0 && pTjet < 85.0 )   stat_PbPb[5]++;
        if (pTjet > 85.0 && pTjet < 100.0)   stat_PbPb[6]++;
      }
      for (int i = 0; i < 7; i++)        stat_PbPb[i] = (1.0/(double)stat_PbPb[i])*pow(stat_PbPb[i], 0.5);
     


  //======================================================
  //========= 4. retrieve unfolded spectra from file
  //=======================================================
    TFile *lefile = new TFile(lerootfile);
    cout << "Unfolding obtained from " << rootfilename <<"\n";
      double phifrac, etafrac;
      updateR(phifrac, etafrac, 0.4);
    vector<TH1F*> unfoldhistsl(0);
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_1"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_2"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_3"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_4"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_5"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_6"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_7"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_8"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_9"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_10"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_11"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_12"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_13"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_14"));
        unfoldhistsl.push_back((TH1F*)lefile->Get("Bayesian_Unfolded_15"));

             //divide unfolded hists by chosen iteration
             vector<TH1F*> ratiohists(0);
             for (int i = 0; i < unfoldhistsl.size(); i++)   {
                 ratiohists.push_back((TH1F*)unfoldhistsl[i]->Clone(names[i]));
                 ratiohists[i]->Divide(unfoldhistsl[iteration]);
             // ratiohists[i]->GetXaxis()->SetRangeUser(lbound, rbound); 
             }

             //set aesthetics for ratio hists
             vector<int> colorlist = {kRed-7, kOrange+1, kOrange-2, kYellow-7, kSpring+5, kGreen-6, kGreen-2, kTeal+2, kTeal-5, kCyan-3, kCyan+2, kAzure-9, kAzure-4, kBlue-4, kBlue-6, kViolet+6, kViolet-8, kMagenta-8};
             //vector<int> colorlist = {kRed-7, kOrange-2, kSpring+5, kGreen-2, kTeal-5, kCyan+2, kAzure-4, kBlue-6, kViolet-8, kMagenta-8};
             for (int i = 0; i < unfoldhistsl.size(); i++)    {
                 setcolor(unfoldhistsl[i], colorlist[i]);
                 setcolor(ratiohists[i], colorlist[i]);
             }
             unfoldhistsl[iteration]->SetLineWidth(2);


  //=======================================================
  //========= 5. retrieve efficiencies from file ==========
  //=======================================================
        //Get Kinematic Efficiency
        TH1F *kinpre = (TH1F*)lefile->Get("KinPre");
        TH1F *kinpos = (TH1F*)lefile->Get("KinPos");
             TH1F *effic = (TH1F*)kinpos->Clone("effic");
             effic->Divide(kinpre);

       //Get Reconstruction Efficiency
       TFile *recfile = new TFile("../Unfolding/ChargedJetRecEfficiencies.root");
       TH1D  *receffic = (TH1D*)recfile->Get("RecEff_R040_5GeV");


  //============================================================
  //========= 6. apply scaling for bin width + N_events ========
  //============================================================
        //Unfolded Spectra 
        for (int i = 0; i < unfoldhistsl.size(); i++)    {
        unfoldhistsl[i]->Scale(1.0, "width");
        //unfoldhistsl[i]->Scale(phifrac/etafrac);
        unfoldhistsl[i]->Scale(1.0/(double)N_event);}

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



  //==========================================================
  //============= published pp spectra and errors
  //==========================================================
      double pTedgepp[19]  = {5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 60, 70, 85, 100};

      //pp spectra
      double jetspec_R040_ltb5[19] = {0.0381711, 0.034356, 0.0271665, 0.0216575, 0.0171139, 0.0119927, 0.00739304, 0.00470566, 0.0030287, 0.00197544, 0.00103004, 0.000412768, 0.000139049, 3.81726*pow(10, -5), 1.33766*pow(10, -5), 5.4954*pow(10, -6), 2.18345*pow(10, -6), 8.16514*pow(10, -7)};
           TH1F* refpp= new TH1F("refpp_R040_ltb5", "", 18, pTedgepp);
           for (int i = 0; i < 18; i++) { refpp->Fill(pTedgepp[i] + 0.01, jetspec_R040_ltb5[i]); }

      //pp statistical errors (absolute)
      double jetspec_R040_ltb5_stat[/*19*/5] = {/*0.000160736, 0.000128482, 0.000120857, 0.000106568, 8.90741*pow(10,-5), 5.94849*pow(10,-5), 4.19683*pow(10,-5), 3.10298*pow(10,-5), 2.25445*pow(10,-5), 1.54392*pow(10,-5), 8.88172*pow(10,-6), 4.11109*pow(10,-6), 1.84098*pow(10,-6), */7.91806*pow(10,-7), 4.85244*pow(10,-7), 3.2539*pow(10,-7), 1.8771*pow(10,-7), 9.15283*pow(10,-8)};

      //pp systematic errors (absolute)
      double jetspec_R040_ltb5_syst[/*19*/] = {/*0.00179404, 0.00164909, 0.00119533, 0.000909613, 0.0007359, 0.000479709, 0.000325294, 0.000216461, 0.000148406, 0.000100747, 5.76822*pow(10,-5), 2.51789*pow(10,-5), 9.59436*pow(10,-6), */2.97746*pow(10,-6), 1.17714*pow(10,-6), 5.33054*pow(10,-7), 2.37996*pow(10,-7), 9.96147*pow(10,-8)};





  raw->Sumw2();
  unfoldhistsl[iteration]->Sumw2();

  //=============================================
  //============= Rebin PbPb to match pp ========
  //=============================================
      vector<double> pTedgerebin = {25, 30, 40, 50, 60, 70, 85, 100};
      TH1F* efficiency = (TH1F*)effic->Rebin(pTedgerebin.size()-1, "efficiency", pTedgerebin.data());
      TH1F* recefficiency = (TH1F*)receffic->Rebin(pTedgerebin.size()-1, "recefficiency", pTedgerebin.data());
          recefficiency->Scale(1.0, "width");
      TH1F* pprebin =(TH1F*)refpp->Rebin(pTedgerebin.size()-1, "pprebin", pTedgerebin.data()); 
      TH1F* unfrebin = new TH1F("unfrebin", "", pTedgerebin.size()-1, pTedgerebin.data());
          for (int i = 0; i < pTedgerebin.size()-1; i++)   {
                int binnum = unfoldhistsl[iteration]->GetXaxis()->FindBin(pTedgerebin[i]+0.01);
                double weight = unfoldhistsl[iteration]->GetBinContent(binnum);
                unfrebin->Fill(pTedgerebin[i]+0.01, weight);
          }
      TH1F* rawrebin = new TH1F("rawrebin", "", pTedgerebin.size()-1, pTedgerebin.data());
          for (int i = 0; i < pTedgerebin.size()-1; i++)   { 
                int binnum = raw->GetXaxis()->FindBin(pTedgerebin[i]+0.01);
                double weight = raw->GetBinContent(binnum);
                rawrebin->Fill(pTedgerebin[i]+0.01, weight);
          }
 

 
  //=============================================
  //============= Calculate the RAA =============
  //=============================================
      //divide PbPb by pp
      TH1F* Raa = (TH1F*)unfrebin->Clone("Raa");
      Raa->Divide(pprebin);

      //apply efficiencies
      Raa->Divide(efficiency);
      Raa->Divide(recefficiency);

      //select and apply TAA
      double TAA010 = 23.3;
      double TAA3050 = 3.9;
      double TAA8090 = 0.084;
      double TAA = 1.0;
      if (cent == true)   TAA = TAA010;
      if (semi == true)   TAA = TAA3050;
      Raa->Scale(1.0/TAA);

      //set aesthetics
      //Raa->SetLineColorAlpha(kWhite, 0.0);
      Raa->SetLineColor(kBlack);
      Raa->SetMarkerStyle(20);
      //Raa->SetMarkerColorAlpha(kBlack, 0.0);
      Raa->SetMarkerColor(kBlack);
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

  //================================================
  //============= Calculate ESE Bounds =============
  //================================================

      TH1F* Raa_lo_L = new TH1F("Raa_lo_L", "", pTedgerebin.size()-1, pTedgerebin.data());
      TH1F* Raa_lo_L2 = new TH1F("Raa_lo_L2", "", pTedgerebin.size()-1, pTedgerebin.data());
      TH1F* Raa_hi_L = new TH1F("Raa_hi_L", "", pTedgerebin.size()-1, pTedgerebin.data());
      TH1F* Raa_hi_L2 = new TH1F("Raa_hi_L2", "", pTedgerebin.size()-1, pTedgerebin.data());
      
      int ESEend;
      if (cent == true) ESEend = startbin+3;
      if (semi == true) ESEend = startbin+5;

      for (int i = startbin; i < ESEend; i++) {
        //calculate E-loss per bin
        double dE = 1 - pow(Raa->GetBinContent(i), 0.25); 
        cout << "Low Edge: " << Raa->GetXaxis()->GetBinLowEdge(i) << ", High Edge: " << Raa->GetXaxis()->GetBinLowEdge(i+1) << ", dE: " << dE <<"\n";

        //calculate R_AA+
        Raa_lo_L->SetBinContent(i, pow(1-1.16*dE, 4)); // ~L
        Raa_lo_L2->SetBinContent(i, pow(1-1.34*dE, 4)); // ~L^2

        //calculate R_AA-
        Raa_hi_L->SetBinContent(i, pow(1-0.83*dE, 4)); // ~L
        Raa_hi_L2->SetBinContent(i, pow(1-0.69*dE, 4)); // ~L^2
      }
      Raa_lo_L->SetLineColor(kMagenta+1);
      Raa_hi_L->SetLineColor(kMagenta+1);
      Raa_lo_L2->SetLineColor(kOrange-3);
      Raa_hi_L2->SetLineColor(kOrange-3);

      Raa_lo_L->SetMarkerColor(kMagenta+1);
      Raa_hi_L->SetMarkerColor(kMagenta+1);
      Raa_lo_L2->SetMarkerColor(kOrange-3);
      Raa_hi_L2->SetMarkerColor(kOrange-3);

      
      TF1* linFit_loL = new TF1("linFitloL", "pol2", 40, 100);
        linFit_loL->SetRange(lowedge, 100);
        Raa_lo_L->Fit(linFit_loL, "N R I", "");
        linFit_loL->SetLineColorAlpha(kMagenta+1, 0.3);
        linFit_loL->SetLineWidth(3);
      TF1* linFit_loL2 = new TF1("linFitloL2", "pol2", 40, 100);
        linFit_loL2->SetRange(lowedge, 100);
        Raa_lo_L2->Fit(linFit_loL2, "N R I", "");
        linFit_loL2->SetLineColorAlpha(kOrange-3, 0.3);
        linFit_loL2->SetLineWidth(3);
      TF1* linFit_hiL = new TF1("linFithiL", "pol2", 40, 100);
        linFit_hiL->SetRange(lowedge, 100);
        Raa_hi_L->Fit(linFit_hiL, "N R I", "");
        linFit_hiL->SetLineColorAlpha(kMagenta+1, 0.3);
        linFit_hiL->SetLineWidth(3);
      TF1* linFit_hiL2 = new TF1("linFithiL2", "pol2", 40, 100);
        linFit_hiL2->SetRange(lowedge, 100);
        Raa_hi_L2->Fit(linFit_hiL2, "N R I", "");
        linFit_hiL2->SetLineColorAlpha(kOrange-3, 0.3);
        linFit_hiL2->SetLineWidth(3);


  //===============================================================
  //============= Calculate Statistical Uncertainties =============
  //===============================================================
      //convert pp errors from absolute to fractional
      double stat_pp[5] = {jetspec_R040_ltb5_stat[0]/pprebin->GetBinContent(0), jetspec_R040_ltb5_stat[1]/pprebin->GetBinContent(1), jetspec_R040_ltb5_stat[2]/pprebin->GetBinContent(2), jetspec_R040_ltb5_stat[3]/pprebin->GetBinContent(3), jetspec_R040_ltb5_stat[4]/pprebin->GetBinContent(4)};
      
      //sum PbPb and pp statistical errors in quadrature
      vector<double> stat_comb(5);
      for (int i = 0; i < 5; i++)   stat_comb[i] = pow(pow(stat_PbPb[i],2) + pow(stat_pp[i],2), 0.5);
      //multiply fractional errors by RAA values
      for (int i = 0; i < 5; i++)   stat_comb[i] = stat_comb[i]*Raa->GetBinContent(i+startbin);

      //create TGraphErrors for statistical errors
      /*int nbins;
      if (cent == true) nbins = 3;
      if (semi == true) nbins = 5;*/
      double bincenters010[3] = {65.0, 77.5, 92.5};
      double bincenters3050[5] = {45.0, 55.0, 65.0, 77.5, 92.5};
      double binvalues010[3] = {Raa->GetBinContent(startbin), Raa->GetBinContent(startbin+1), Raa->GetBinContent(startbin+2)};
      double binvalues3050[5] = {Raa->GetBinContent(startbin), Raa->GetBinContent(startbin+1), Raa->GetBinContent(startbin+2), Raa->GetBinContent(startbin+3), Raa->GetBinContent(startbin+4)};
      double nul010[3] = {5.0, 7.5, 7.5};
      double nul3050[5] = {5.0, 5.0, 5.0, 7.5, 7.5};
      //convert stat_comb to correct type
      double stat_comb010[3] = {stat_comb[2], stat_comb[3], stat_comb[4]};
      double stat_comb3050[5] = {stat_comb[1], stat_comb[2], stat_comb[3], stat_comb[4], stat_comb[5]};
      //input to TGraphErrors
      TGraphErrors *stat;
      if (cent == true) stat = new TGraphErrors(3, bincenters010, binvalues010, nul010, stat_comb010);
      if (semi == true) stat = new TGraphErrors(5, bincenters3050, binvalues3050, nul3050, stat_comb3050);
           stat->SetMarkerStyle(20);
           stat->SetMarkerSize(0.7);
           if (cent == true) {
              stat->SetMarkerColor(centcolor);
              stat->SetLineColor(centcolor);}
           if (semi == true) {
              stat->SetMarkerColor(kCyan+3);
              stat->SetLineColor(kCyan+3);
           }


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


  //===================================================
  //============= reference values from Hannah
  //===================================================
      //Hannah Area Based RAA
      double hannahedge[5] = {40, 60, 70, 85, 100};
      double hannah010[4]  = {0.306683, 0.300064, 0.324841, 0.368015};
      double hannah3050[4] = {0.536177, 0.602989, 0.592376, 0.596467};
      TH1F* hannahRaa = new TH1F("hannahraa", "", 4, hannahedge);
      if (cent == true)  {
         for (int i = 0; i < 4; i++)   hannahRaa->Fill(hannahedge[i]+0.01, hannah010[i]);}
      if (semi == true)  {
         for (int i = 0; i < 4; i++)   hannahRaa->Fill(hannahedge[i]+0.01, hannah3050[i]);}

      //Hannah Machine Learning RAA
      double mledge[8] = {20, 30, 40, 50, 60, 70, 85, 100};
      double ml010[7] = {0.167145, 0.195948, 0.232806, 0.267605, 0.295418, 0.333034, 0.307735};
      double ml3050[7] = {0.396661, 0.441344, 0.518221, 0.54668, 0.533559, 0.608539, 0.644078};
      TH1F* mlRaa = new TH1F("mlRaa", "", 7, mledge);
      if (cent == true)  {
         for (int i = 0; i < 4; i++)   mlRaa->Fill(mledge[i]+0.01, ml010[i]);}
      if (semi == true)  {
         for (int i = 0; i < 4; i++)   mlRaa->Fill(mledge[i]+0.01, ml3050[i]);}




double hcenters[4] = {50, 65, 77.5, 92.5};
double hvalues010[4] = {0.306683, 0.300064, 0.324841, 0.368015};
double hvalues3050[4] = {0.536177, 0.602989, 0.592376, 0.596467};
double hnul[4] = {10, 5, 7.5, 7.5};
double herr010[4]  = {0.0322854, 0.0327882, 0.0300957, 0.0361132};
double herr3050[4] = {0.0373936, 0.0527642, 0.0548663, 0.0592145};
//double hvalues[4];
//double herr[4];
    auto hsyst010 = new TGraphErrors(4, hcenters, hvalues010, hnul, herr010);
    auto hsyst3050 = new TGraphErrors(4, hcenters, hvalues3050, hnul, herr3050);
     hsyst010->SetFillColorAlpha(kCyan, 0.2); 
     hsyst010->SetLineColor(kAzure);
     hsyst010->SetMarkerSize(0.6);
     hsyst010->SetMarkerStyle(20);
     hsyst010->SetMarkerColor(kAzure);
     hsyst3050->SetFillColorAlpha(kCyan, 0.2); 
     hsyst3050->SetLineColor(kAzure);
     hsyst3050->SetMarkerSize(0.6);
     hsyst3050->SetMarkerStyle(20);
     hsyst3050->SetMarkerColor(kAzure);
double mlcenters[6] = {35, 45, 55, 65, 77.5, 92.5};
double mlvalues010[6] = {0.195948, 0.232806, 0.267605, 0.295418,  0.333034 ,  0.307735};
double mlvalues3050[6] = {0.441344, 0.518221, 0.54668, 0.533559, 0.608539, 0.644078};
double mlnul[6] =  {5, 5, 5, 5, 7.5, 7.5};
double mlerr010[6] = {0.0284162, 0.0241808, 0.0235063, 0.0380027, 0.036138, 0.0532162};
double mlerr3050[6] = {0.0549399, 0.0743638, 0.103646, 0.090969, 0.100182, 0.134028};
   auto mlsyst010 = new TGraphErrors(6, mlcenters, mlvalues010, mlnul, mlerr010);
   auto mlsyst3050 = new TGraphErrors(6, mlcenters, mlvalues3050, mlnul, mlerr3050);
     mlsyst010->SetFillColorAlpha(kMagenta, 0.2);
     mlsyst010->SetLineColor(kMagenta+1);
     mlsyst010->SetMarkerSize(0.6);
     mlsyst010->SetMarkerStyle(21);
     mlsyst010->SetMarkerColor(kMagenta+1);
     mlsyst3050->SetFillColorAlpha(kAzure, 0.2);
     mlsyst3050->SetLineColor(kAzure);
     mlsyst3050->SetMarkerSize(0.6);
     mlsyst3050->SetMarkerStyle(21);
     mlsyst3050->SetMarkerColor(kAzure);


//TH1F* trivial  = (TH1F*)unfoldhistsl[iteration]->Clone("trivial");
//trivial->Divide(truth);
//trivial->Scale(1.0, "width");

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
   auto gen = new TLegend(0.15, 0.7, 0.3, 0.85);
      gen->SetTextSize(0.035);
      gen->SetBorderSize(0);
        if (cent == true)    gen->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb, 5.02 TeV, 0#font[122]{-}10%", "");
        if (semi == true)    gen->AddEntry((TObject*)0, "ALICE Pb#font[122]{-}Pb, 5.02 TeV, 40#font[122]{-}45%", "");
        gen->AddEntry((TObject*)0, Form("Charged Jets, anti-#it{k}_{T}, |#it{#eta}_{jet}| < %.1f", 0.9-(double)R*0.1), "");
        gen->AddEntry((TObject*)0, "Work in Progress", "");
     //if (ESE == false)    histkey->AddEntry(caitie, "#it{R} = 0.4, #it{p}_{T}^{leading track} > 5 GeV/#it{c}", "plf");
    // histkey->AddEntry(caitie, "ALICE Work in Progress", "plf");
     //histkey->AddEntry(mlsyst3050, "#it{R} = 0.4, ML-Bkgd. Sub. (2015)", "plf");
     //if (cent == true)   histkey->AddEntry(hsyst010, "Hannah Area-Based, p_{T}^{leading track} > 7 GeV (2015)", "pl");
     //if (semi == true)   histkey->AddEntry(hsyst3050, "Hannah Area-Based, p_{T}^{leading track} > 7 GeV (2015)", "pl");
     //if (cent == true)   histkey->AddEntry(mlsyst010, "Hannah ML (2015)", "pl");
     //if (semi == true)   histkey->AddEntry(mlsyst3050, "Hannah ML (2015)", "pl");
   auto *close = new TLegend(0.5, 0.6, 0.85, 0.85);
     close->SetBorderSize(0);
     close->SetTextSize(0.04);
     close->AddEntry(unfoldhistsl[iteration], "unfolded", "pl");
     close->AddEntry(raw, "raw", "pl");
     close->AddEntry(truth, "truth (from embedding)", "pl");
     close->AddEntry(respraw, "hybrid (from embedding)", "pl");
   auto *refo = new TLegend(0.6, 0.75, 0.85, 0.85);
     refo->SetBorderSize(0);
     refo->SetTextSize(0.055);
     refo->AddEntry((TObject*)0, "Refolded/Raw", "");


  //====================================
  //========== TAA uncertainty =========
  //====================================
  /*double systTAA = 0.0165;
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
  TCanvas *c1 = new TCanvas("RAA", lerootfile, 600, 500);
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
     //caitie->Draw("ezp 2 same");
    //syst3050->Draw("ezp 2 same");
    //stat->Draw("e 3 same");
     //if (cent == true)     hsyst010->Draw("ezp 2 same");
     //if (semi == true)     hsyst3050->Draw("ezp 2 same");
     //if (cent == true)     mlsyst010->Draw("ezp 2 same");
     //if (semi == true)     mlsyst3050->Draw("ezp 2 same");
     //stat->Draw("ezp same");
   //Raa->Draw("same");
   
   //uncertaintyTAA->Draw("same");
   //TAAunc->Draw("same");
   gen->Draw("same");
   mline->Draw("same");
  
  //Draw Kinematic Efficiency
  TCanvas *c2 = new TCanvas("Kinematic Efficiency", "", 500, 500);
   c2->cd();
   gStyle->SetOptStat(0);
   TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
     pad2->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad2->Draw();
   pad2->cd();               // pad1 becomes the current pad
   efficiency->SetMinimum(0.01);
   efficiency->SetMaximum(.99);   
   efficiency->SetXTitle("p_{T} [GeV/c]");
   efficiency->SetTitle("");
   efficiency->Draw("hist same");
   auto kin = new TLegend(0.45, 0.3, 0.8, 0.55);
     kin->SetBorderSize(0);
     kin->SetTextSize(0.035);
     kin->AddEntry((TObject*)0, "Kinematic Efficiency", "");
     kin->Draw("same");


  //Draw Stability/Refolding Tests
  TCanvas *c4= new TCanvas("Stability Test", lerootfile, 1000, 700);
   c4->cd();
   gStyle->SetOptStat(0);
   gStyle->SetPadTickY(1);
   TPad *pad4 = new TPad("pad4", "", 0.0, 0.4, 0.5, 1.0);
   TPad *pad5 = new TPad("pad5", "", 0.0, 0.05, 0.5, 0.4);
     pad4->SetBottomMargin(0.0); // Upper and lower plot are joined
     pad4->SetRightMargin(0.0);
     pad5->SetTopMargin(0.0);
     pad5->SetRightMargin(0.0);
     pad5->SetBottomMargin(0.2);
     pad4->Draw();
     pad5->Draw();
   pad4->cd();               // pad1 becomes the current pad   
   //pad1->SetLogx();
   pad4->SetLogy();
   unfoldhistsl[0]->SetMaximum(pow(10, -1));
   unfoldhistsl[0]->SetMinimum(pow(10, -9));
   unfoldhistsl[0]->SetTitle("");
   for (int i = 0; i < unfoldhistsl.size(); i++)         unfoldhistsl[i]->Draw("same");
   pad5->cd();
   ratiohists[0]->GetXaxis()->SetLabelSize(0.07);
   ratiohists[0]->GetYaxis()->SetLabelSize(0.06);
   ratiohists[0]->SetTitle("");
   ratiohists[0]->SetXTitle("p_{T} [GeV/c]");
   ratiohists[0]->GetXaxis()->SetTitleSize(0.09);
   ratiohists[0]->SetMinimum(0.5);
   ratiohists[0]->SetMaximum(1.5);
   for (int i = 0; i < ratiohists.size(); i++)         ratiohists[i]->Draw("same");
   line->Draw("same");  
     dcut->Draw("same");
     dcut2->Draw("same");
     dcut3->Draw("same");
     ucut->Draw("same");
     ucut2->Draw("same");
     ucut3->Draw("same");
   fileleg->Draw("same");


   c4->cd();
   TPad *pad6 = new TPad("pad6", "", 0.5, 0.4, 1.0, 1.0);
        pad6->SetBottomMargin(0.0);
        pad6->SetLeftMargin(0.0);
        pad6->Draw();
      pad6->cd();
      pad6->SetLogy();
         unfoldhistsl[iteration]->SetMaximum(pow(10, -1));
         unfoldhistsl[iteration]->SetMinimum(pow(10, -9));
         unfoldhistsl[iteration]->SetTitle("");
            setcolor(unfoldhistsl[iteration], kGreen+2);
         unfoldhistsl[iteration]->GetXaxis()->SetRangeUser(25.0, 120.0);
         unfoldhistsl[iteration]->Draw("same");
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



}

void setcolor(TH1F* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(1);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.6);
}

