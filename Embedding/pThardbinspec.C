
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;

void getq2(const char *anfile);
void getbin(const char *analysisfile, double centl, double centr, int jetR, vector<vector<long double>> q2vec);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);

int lbin = 10;
int rbin = 200;
int binres = 95;
  
vector<vector<long double>> q2percentiles(0);
    
TH1D *binspec = new TH1D("bin", "", binres, lbin, rbin);
TH1D *genspec = new TH1D("gen", "", binres, lbin, rbin);
TH1D *detspec = new TH1D("det", "", binres, lbin, rbin);

TH1D *binsum1 = new TH1D("bin1sum", "", binres, lbin, rbin);
TH1D *gensum1 = new TH1D("gen1sum", "", binres, lbin, rbin);
TH1D *detsum1 = new TH1D("det1sum", "", binres, lbin, rbin);
TH1D *binsum2 = new TH1D("bin2sum", "", binres, lbin, rbin);
TH1D *gensum2 = new TH1D("gen2sum", "", binres, lbin, rbin);
TH1D *detsum2= new TH1D("det2sum", "", binres, lbin, rbin);

long int totevents = 0;

vector<const char*> filenames4 = {"TrainOutputs/AnalysisResults7158.root", "TrainOutputs/AnalysisResults7121.root", "TrainOutputs/AnalysisResults7122.root", "TrainOutputs/AnalysisResults7123.root", "TrainOutputs/AnalysisResults7124.root", "TrainOutputs/AnalysisResults7126.root", "TrainOutputs/AnalysisResults7127.root", "TrainOutputs/AnalysisResults7128.root", "TrainOutputs/AnalysisResults7129.root", "TrainOutputs/AnalysisResults7131.root", "TrainOutputs/AnalysisResults7132.root", "TrainOutputs/AnalysisResults7133.root", "TrainOutputs/AnalysisResults7134.root", "TrainOutputs/AnalysisResults7135.root", "TrainOutputs/AnalysisResults7136.root", "TrainOutputs/AnalysisResults7159.root", "TrainOutputs/AnalysisResults7160.root"};


//R = 0.4 ESE, 30-50%, pass1q
//vector<const char*> filenames3 = {"../ESE/Trains/AnalysisResults7637_2.root", "../ESE/Trains/AnalysisResults7638_2.root", "../ESE/Trains/AnalysisResults7639_2.root", "../ESE/Trains/AnalysisResults7640_2.root", "../ESE/Trains/AnalysisResults7641_2.root", "../ESE/Trains/AnalysisResults7642_2.root", "../ESE/Trains/AnalysisResults7643_2.root", "../ESE/Trains/AnalysisResults7644_2.root", "../ESE/Trains/AnalysisResults7645_2.root", "../ESE/Trains/AnalysisResults7646_2.root", "../ESE/Trains/AnalysisResults7647_2.root", "../ESE/Trains/AnalysisResults7648_2.root", "../ESE/Trains/AnalysisResults7650_2.root", "../ESE/Trains/AnalysisResults7651_2.root", "../ESE/Trains/AnalysisResults7652_2.root", "../ESE/Trains/AnalysisResults7653_2.root", "../ESE/Trains/AnalysisResults7654_2.root", "../ESE/Trains/AnalysisResults7655_2.root", "../ESE/Trains/AnalysisResults7656_2.root"};

//R = 0.2 ESE, 30-50%, pass1q
/*vector<const char*> filenames3 = {"../ESE/Trains/AnalysisResults7671_2.root", "../ESE/Trains/AnalysisResults7694_2.root", 
       "../ESE/Trains/AnalysisResults7673_2.root", "../ESE/Trains/AnalysisResults7674_2.root", "../ESE/Trains/AnalysisResults7675_2.root",
       "../ESE/Trains/AnalysisResults7676_2.root", "../ESE/Trains/AnalysisResults7677_2.root", "../ESE/Trains/AnalysisResults7678_2.root",
       "../ESE/Trains/AnalysisResults7695_2.root", "../ESE/Trains/AnalysisResults7680_2.root", "../ESE/Trains/AnalysisResults7681_2.root",
       "../ESE/Trains/AnalysisResults7682_2.root", "../ESE/Trains/AnalysisResults7683_2.root", "../ESE/Trains/AnalysisResults7685_2.root", 
       "../ESE/Trains/AnalysisResults7686_2.root", "../ESE/Trains/AnalysisResults7687_2.root", "../ESE/Trains/AnalysisResults7690_2.root",
       "../ESE/Trains/AnalysisResults7691_2.root", "../ESE/Trains/AnalysisResults7692_2.root", "../ESE/Trains/AnalysisResults7693_2.root"};*/

//R = 0.2, ESE, 30-50%, pass3q
vector<const char*> filenames3 = {"../ESE/LHC18q/AnalysisResults7813.root", "../ESE/LHC18q/AnalysisResults7814.root", 
       "../ESE/LHC18q/AnalysisResults7815.root", "../ESE/LHC18q/AnalysisResults7816.root", "../ESE/LHC18q/AnalysisResults7835.root",
       "../ESE/LHC18q/AnalysisResults7818.root", "../ESE/LHC18q/AnalysisResults7819.root", "../ESE/LHC18q/AnalysisResults7820.root",
       "../ESE/LHC18q/AnalysisResults7821.root", "../ESE/LHC18q/AnalysisResults7822.root", "../ESE/LHC18q/AnalysisResults7823.root", 
       "../ESE/LHC18q/AnalysisResults7824.root", "../ESE/LHC18q/AnalysisResults7825.root", "../ESE/LHC18q/AnalysisResults7826.root",
       "../ESE/LHC18q/AnalysisResults7827.root", "../ESE/LHC18q/AnalysisResults7828.root", "../ESE/LHC18q/AnalysisResults7829.root",
       "../ESE/LHC18q/AnalysisResults7830.root", "../ESE/LHC18q/AnalysisResults7831.root", "../ESE/LHC18q/AnalysisResults7832.root"};

//R = 0.2 ESE, 0-10%, pass1q
//vector<const char*> filenames3 = {"../ESE/Trains/AnalysisResults7697_2.root", "../ESE/Trains/AnalysisResults7698_2.root", "../ESE/Trains/AnalysisResults7699_2.root", "../ESE/Trains/AnalysisResults7700_2.root", "../ESE/Trains/AnalysisResults7701_2.root", "../ESE/Trains/AnalysisResults7702_2.root", "../ESE/Trains/AnalysisResults7703_2.root", "../ESE/Trains/AnalysisResults7704_2.root", "../ESE/Trains/AnalysisResults7705_2.root", "../ESE/Trains/AnalysisResults7706_2.root", "../ESE/Trains/AnalysisResults7707_2.root", "../ESE/Trains/AnalysisResults7708_2.root", "../ESE/Trains/AnalysisResults7709_2.root", "../ESE/Trains/AnalysisResults7710_2.root", "../ESE/Trains/AnalysisResults7711_2.root", "../ESE/Trains/AnalysisResults7712_2.root", "../ESE/Trains/AnalysisResults7713_2.root", "../ESE/Trains/AnalysisResults7714_2.root", "../ESE/Trains/AnalysisResults7715_2.root", "../ESE/Trains/AnalysisResults7716_2.root"};

int j = 0;
int start = 0;

  vector<TH1D*> bin(0);
  vector<TH1D*> gen(0);
  vector<TH1D*> det(0);

  vector<const char*> bins = {"bin1", "bin2", "bin3", "bin4", "bin5", "bin6", "bin7", "bin8", "bin9", "bin10", "bin11", "bin12", "bin13", "bin14", "bin15", "bin16", "bin17", "bin18", "bin19", "bin20"};
  vector<const char*> gens = {"gen1", "gen2", "gen3", "gen4", "gen5", "gen6", "gen7", "gen8", "gen9", "gen10", "gen11", "gen12", "gen13", "gen14", "gen15", "gen16", "gen17", "gen18", "gen19", "gen20"};
  vector<const char*> dets = {"det1", "det2", "det3", "det4", "det5", "det6", "det7", "det8", "det9", "det10", "det11", "det12", "det13", "det14", "det15", "det16", "det17", "det18", "det19", "det20"};




//========================================
//========================================
//=========== main function ==============
//========================================
//========================================

void pThardbinspec (double cl, double cr, int jetR = 2) {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  getq2("../ESE/AnalysisResults7812.root");

  for (int i = 0; i < filenames3.size(); i++)    { 
      getbin(filenames3[i], cl, cr, jetR, q2percentiles);
        bin.push_back((TH1D*)binspec->Clone(bins[i]));
        gen.push_back((TH1D*)genspec->Clone(gens[i]));
        det.push_back((TH1D*)detspec->Clone(dets[i]));
          binsum1->Add(bin[i]);
          gensum1->Add(gen[i]);
          detsum1->Add(det[i]);
      }
  bin.push_back(binsum1);
  gen.push_back(gensum1);
  det.push_back(detsum1);


  cout << "Total events: " << totevents <<"\n";
  double avgevents = (double)totevents/20.0;
  cout << "Average Events: " << avgevents <<"\n";
  
  //Scale by Average Events
  for (int i = 0; i <= filenames3.size(); i++) {
      bin[i]->Scale(avgevents);  
      gen[i]->Scale(avgevents);
      det[i]->Scale(avgevents);
   }
  
  setstyle(bin, gen, det);




  //==========================================================
  //============= Legends  ===================================
  //==========================================================

  //pThard bin legend
  vector<string> newleg(0);
  for (int i = 0; i < filenames3.size(); i++) {
    int index = i+start;
    string nl = Form("bin %d", index);      //bin key
    newleg.push_back(nl);
  }
  vector<const char*> newlegchar(0);
  for (int i = 0; i < newleg.size(); i++) newlegchar.push_back(newleg[i].c_str());

  auto leg1 = new TLegend(0.45, 0.55, 0.65, 0.85);
  leg1->SetTextSize(0.035);
  leg1->SetBorderSize(0);
  auto leg2 = new TLegend(0.65, 0.55, 0.85, 0.85);
  leg2->SetTextSize(0.035);
  leg2->SetBorderSize(0);
  for (int i = 0; i < filenames3.size(); i++) {
      if (i < 10)      leg1->AddEntry(bin[i], newlegchar[i], "pl");
      if (i >= 10)     leg2->AddEntry(bin[i], newlegchar[i], "pl");
  }

  auto info = new TLegend(0.15, 0.7, 0.35, 0.85);
    info->SetTextSize(0.035);
    info->SetBorderSize(0);
    info->AddEntry((TObject*)0, "ALICE, Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV, 30#font[122]{-}50%", "");
    info->AddEntry((TObject*)0, "Embedding Spectra, #it{p}_{T, true} > 10 GeV/#it{c}", "");


  auto datarun = new TLegend(0.15, 0.75, 0.35, 0.85);
  datarun->SetTextSize(0.035);
  datarun->SetBorderSize(0);
  datarun->AddEntry((TObject*)0, "Charged Jets, R = 0.2, anti-#it{k}_{T}, |#it{#eta}_{jet}| < 0.7", "");
  datarun->AddEntry((TObject*)0, "30\% largest #it{q}_{2}", ""); 


  TLine *line = new TLine(10.0, pow(10,4), 10.0, pow(10,16));
  line->SetLineStyle(2);

  TCanvas *c1 = new TCanvas("c1", "new", 1650, 550);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad11 = new TPad("pad11", "pad11", 0.0, 0.05, 0.33, 1.0);
     pad11->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad11->Draw();             // Draw the upper pad: pad1
     pad11->cd();               // pad1 becomes the current pad
     pad11->SetLogy();
     gensum1->SetXTitle("#it{p}_{T}^{truth} (GeV/#it{c})");
     double ymax = gensum1->GetBinContent(2);
     double ymin = gensum1->GetBinContent(40);
       cout << "ymax: " << ymax <<"\n";
       cout << "ymin: " << ymin <<"\n";
     gensum1->SetMinimum(pow(10,-6)*ymin);
     gensum1->SetMaximum(pow(10,5)*ymax);
     gensum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   gen[i]->Draw("same");  
       info->Draw("same");
     c1->cd();
   TPad *pad12 = new TPad("pad12", "pad12", 0.33, 0.05, 0.67, 1.0);
     pad12->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad12->Draw();             // Draw the upper pad: pad1
     pad12->cd();               // pad1 becomes the current pad
     pad12->SetLogy();
     detsum1->SetXTitle("#it{p}_{T}^{detector} (GeV/#it{c})");
     detsum1->SetMinimum(pow(10,-6)*ymin);
     detsum1->SetMaximum(pow(10,5)*ymax);
     detsum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   det[i]->Draw("same");  
       datarun->Draw("same");
     c1->cd();
   TPad *pad13 = new TPad("pad13", "pad13", 0.67, 0.05, 1.0, 1.0);
     pad13->SetBottomMargin(0.1); // Upper and lower plot are joined
     pad13->Draw();             // Draw the upper pad: pad1
     pad13->cd();               // pad1 becomes the current pad
     pad13->SetLogy();
     binsum1->SetXTitle("#it{p}_{T}^{hybrid} (GeV/#it{c})");
     binsum1->SetMinimum(pow(10,-6)*ymin);
     binsum1->SetMaximum(pow(10,5)*ymax);
     binsum1->Draw("same");
     for (int i = 0; i <= filenames3.size(); i++)   bin[i]->Draw("same");  
       leg1->Draw("same");
       if (filenames3.size() > 10) leg2->Draw("same");
 

  
}






//-----------------------------------------------------------------------------------------------------



void getq2(const char *anfile)  {

      TFile *fin = new TFile(anfile);
      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)fin->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2V0C_EPV0C");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                             //number percentile we're on
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      double centl = 30.0;
      double centr = 50.0;
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
  /*TH1D *perc20 = new TH1D("perc20", "perc20", centdiff, centl, centr);
  TH1D *perc80 = new TH1D("perc80", "perc80", centdiff, centl, centr);
  for (int i = 1; i <= centdiff; i++) {
      perc20->SetBinContent(i, allPercent[i-1][1]);
      perc80->SetBinContent(i, allPercent[i-1][7]);
  }
  perc20->SetLineColor(kRed);
  perc80->SetLineColor(kRed);
*/
  q2percentiles = allPercent;

}




//----------------------------------------------------------------------------------





void getbin(const char *analysisfile, double centl, double centr, int jetR, vector<vector<long double>> q2vec) {

    j++;    

    binspec->Reset();
    genspec->Reset();
    detspec->Reset(); 


    string rootfilename = analysisfile;
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);
 
    AliEmcalList *EmbeddingHelper = (AliEmcalList*)lefile->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");    
    TProfile *pThardscale = (TProfile*)EmbeddingHelper->FindObject("fHistXsection");
    TH1D *pThardbin = (TH1D*)EmbeddingHelper->FindObject("fHistTrials");
    TH1D *eventcount = (TH1D*)EmbeddingHelper->FindObject("fHistEventCount");


    TTree *T = nullptr;
    if (jetR == 2) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_hybJet", T);
    if (jetR == 4) lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_hybJet", T);
    float hybridpT, detpT, partpT, centrality, area, matchradius, q2;
    T->SetBranchAddress("Jet_Pt", &hybridpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &partpT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Jet_Area", &area);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius);
    T->SetBranchAddress("Event_Q2Vector", &q2); 

    int entries = T->GetEntries();

    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    double pThardFrac = xscale/(n_event*trials);

    if (j == 1)  start = pThardbin->GetMean();  


    //Determine Extractor Bin Scaling
    const int extractBins = 9;
    //double extractFrac[extractBins] = {0.0, 0.01, 0.01, 0.02, 0.02, 0.05, 0.05, 0.05};
    //double extractFrac[extractBins] = {0.0, 0.02, 0.03, 0.1, 0.12, 0.12, 0.12, 0.07};  //R = 0.2, 0-10%, ESE, lhc18q
    double extractFrac[extractBins] = {0.0, 0.01, 0.03, 0.15, 0.20, 0.20, 0.20, 0.15, 0.05};  //R = 0.2, 30-50%, ESE, lhc18q
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0};
    

    int fillcount = 0;
    for (int i = 0; i < entries; i++) {
       T->GetEntry(i);
       if (centrality < centl || centrality > centr)  continue;
       if (hybridpT < 10.0 /*|| detpT < 10.0*/ || partpT < 10.0 )  continue;
       //if (q2 > 2.45) continue; //20th percentile
       if (q2 < q2vec[(int)centrality-30][6])  continue;
       //if (q2 < 5.85) continue; //80th percentile
       int scaleindex = -1;
         if (hybridpT >= extractEdge[0] && hybridpT < extractEdge[1])  scaleindex = 0;
         if (hybridpT >= extractEdge[1] && hybridpT < extractEdge[2])  scaleindex = 1;
         if (hybridpT >= extractEdge[2] && hybridpT < extractEdge[3])  scaleindex = 2;
         if (hybridpT >= extractEdge[3] && hybridpT < extractEdge[4])  scaleindex = 3;
         if (hybridpT >= extractEdge[4] && hybridpT < extractEdge[5])  scaleindex = 4;
         if (hybridpT >= extractEdge[5] && hybridpT < extractEdge[6])  scaleindex = 5;
         if (hybridpT >= extractEdge[6] && hybridpT < extractEdge[7])  scaleindex = 6;
         if (hybridpT >= extractEdge[7] && hybridpT < extractEdge[8])  scaleindex = 7;
         if (hybridpT >= extractEdge[8] && hybridpT < extractEdge[9])  scaleindex = 8;
       if (scaleindex == 0)   continue;
       binspec->Fill(hybridpT, pThardFrac/extractFrac[scaleindex]);
       genspec->Fill(partpT, pThardFrac/extractFrac[scaleindex]);
       detspec->Fill(detpT, pThardFrac/extractFrac[scaleindex]);
   }

    cout << "pT hard bin: " << pThardbin->GetMean() <<"\n";
    //cout << "	Scale Factor: " << pThardFrac << "*avg_events" << "\n";    

    totevents += n_event;
    int pTindex = pThardbin->GetMean() - 1;
   
    lefile->Close();

}

void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3) {
  vector<int> colors = {kRed+2, kRed, kOrange+8, kOrange+6, kOrange-4, kYellow, kSpring+10, kSpring, kGreen-3, /*kGreen+2,*/ kTeal+6, kTeal, kCyan-3, kAzure+7, kBlue+1, kBlue-2, kViolet+2, kViolet, kMagenta-4, kPink+6, kPink+10};
  int entryremoval = colors.size() - filenames3.size();
  cout << "Entry removal" << entryremoval <<"\n"; 
  for (int i = 0; i < entryremoval; i++)   colors.pop_back();
  colors.push_back(kBlack);

  for (int i = 0; i < vec1.size(); i++) {
      vec1[i]->SetLineColor(colors[i]);
      vec1[i]->SetMarkerColor(colors[i]);
        vec1[i]->SetLineWidth(2);
        vec1[i]->SetMarkerSize(0.5);
        vec1[i]->SetMarkerStyle(20);
      vec2[i]->SetLineColor(colors[i]);
      vec2[i]->SetMarkerColor(colors[i]);
        vec2[i]->SetLineWidth(2);
        vec2[i]->SetMarkerSize(0.5);
        vec2[i]->SetMarkerStyle(20);
      vec3[i]->SetLineColor(colors[i]);
      vec3[i]->SetMarkerColor(colors[i]);
        vec3[i]->SetLineWidth(2);
        vec3[i]->SetMarkerSize(0.5);
        vec3[i]->SetMarkerStyle(20);
  }



}
