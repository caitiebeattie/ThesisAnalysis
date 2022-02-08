#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TColor.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRatioPlot.h"
#include "TLatex.h"

using namespace std;


    int xres = 9;
    int yres = 100;
    int xl = 10.0;
    int xr = 100.0;


//R = 0.2, 30-50%, ESE, lhc18q, pass1
/*vector<const char*> filenames = {"../ESE/Trains/AnalysisResults7671_2.root", "../ESE/Trains/AnalysisResults7694_2.root",
          "../ESE/Trains/AnalysisResults7673_2.root", "../ESE/Trains/AnalysisResults7674_2.root", "../ESE/Trains/AnalysisResults7675_2.root", 
          "../ESE/Trains/AnalysisResults7676_2.root", "../ESE/Trains/AnalysisResults7677_2.root", "../ESE/Trains/AnalysisResults7678_2.root", 
          "../ESE/Trains/AnalysisResults7695_2.root", "../ESE/Trains/AnalysisResults7680_2.root", "../ESE/Trains/AnalysisResults7681_2.root",
          "../ESE/Trains/AnalysisResults7682_2.root", "../ESE/Trains/AnalysisResults7683_2.root", "../ESE/Trains/AnalysisResults7685_2.root", 
          "../ESE/Trains/AnalysisResults7686_2.root", "../ESE/Trains/AnalysisResults7687_2.root", "../ESE/Trains/AnalysisResults7690_2.root",
          "../ESE/Trains/AnalysisResults7691_2.root", "../ESE/Trains/AnalysisResults7692_2.root"};*/

//R = 0.2, 30-50%, ESE, lhc18q, pass3
vector<const char*> filenames = {"../ESE/LHC18q/AnalysisResults7813.root", "../ESE/LHC18q/AnalysisResults7814.root", 
       "../ESE/LHC18q/AnalysisResults7815.root", "../ESE/LHC18q/AnalysisResults7816.root", "../ESE/LHC18q/AnalysisResults7835.root",
       "../ESE/LHC18q/AnalysisResults7818.root", "../ESE/LHC18q/AnalysisResults7819.root", "../ESE/LHC18q/AnalysisResults7820.root",
       "../ESE/LHC18q/AnalysisResults7821.root", "../ESE/LHC18q/AnalysisResults7822.root", "../ESE/LHC18q/AnalysisResults7823.root", 
       "../ESE/LHC18q/AnalysisResults7824.root", "../ESE/LHC18q/AnalysisResults7825.root", "../ESE/LHC18q/AnalysisResults7826.root",
       "../ESE/LHC18q/AnalysisResults7827.root", "../ESE/LHC18q/AnalysisResults7828.root", "../ESE/LHC18q/AnalysisResults7829.root",
       "../ESE/LHC18q/AnalysisResults7830.root", "../ESE/LHC18q/AnalysisResults7831.root", "../ESE/LHC18q/AnalysisResults7832.root"};

  vector<const char*> hds = {"hd1", "hd2", "hd3", "hd4", "hd5", "hd6", "hd7", "hd8", "hd9", "hd10", "hd11", "hd12", "hd13", "hd14", "hd15", "hd16", "hd17", "hd18", "hd19", "hd20"};
  vector<const char*> dts = {"dt1", "dt2", "dt3", "dt4", "dt5", "dt6", "dt7", "dt8", "dt9", "dt10", "dt11", "dt12", "dt13", "dt14", "dt15", "dt16", "dt17", "dt18", "dt19", "dt20"};
  vector<const char*> hts = {"ht1", "ht2", "ht3", "ht4", "ht5", "ht6", "ht7", "ht8", "ht9", "ht10", "ht11", "ht12", "ht13", "ht14", "ht15", "ht16", "ht17", "ht18", "ht19", "ht20"};
  
  vector<const char*> reshds = {"reshd1", "reshd2", "reshd3", "reshd4", "reshd5", "reshd6", "reshd7", "reshd8", "reshd9", "reshd10", "reshd11", "reshd12", "reshd13", "reshd14", "reshd15", "reshd16", "reshd17", "reshd18", "reshd19", "reshd20"};
  vector<const char*> resdts = {"resdt1", "resdt2", "resdt3", "resdt4", "resdt5", "resdt6", "resdt7", "resdt8", "resdt9", "resdt10", "resdt11", "resdt12", "resdt13", "resdt14", "resdt15", "resdt16", "resdt17", "resdt18", "resdt19", "resdt20"};
  vector<const char*> reshts = {"resht1", "resht2", "resht3", "resht4", "resht5", "resht6", "resht7", "resht8", "resht9", "resht10", "resht11", "resht12", "resht13", "resht14", "resht15", "resht16", "resht17", "resht18", "resht19", "resht20"};
 
    vector<bool> skip(0);
 
    vector<TH1D*> hdvec(0);
    vector<TH1D*> dtvec(0);
    vector<TH1D*> htvec(0);
 
    vector<TH2D*> reshd(0);
    vector<TH2D*> resdt(0);
    vector<TH2D*> resht(0);
 
    //pT^reco/pT^part
    TH1D *jetEnergyHT = new TH1D("jetEnergyHT", "", xres, 0.0, 5.0); 
    TH1D *jetEnergyHD = new TH1D("jetEnergyHD", "", xres, 0.0, 5.0);
    TH1D *jetEnergyDT = new TH1D("jetEnergyDT", "", xres, 0.0, 5.0);

    TH1D *HT = new TH1D("HT", "", xres, 0.0, 5.0); 
    TH1D *HD = new TH1D("HD", "", xres, 0.0, 5.0);
    TH1D *DT = new TH1D("DT", "", xres, 0.0, 5.0);

    //delta pT
    TH1D *pTht = new TH1D("pTht", "", 120, -60.0, 60.0);

    //residuals
    TH2D *residualht = new TH2D("residualht", "", xres, xl, xr, yres, -2.0, 2.0);
    TH2D *residualhd = new TH2D("residualhd", "", xres, xl, xr, yres, -2.0, 2.0);
    TH2D *residualdt = new TH2D("residualdt", "", xres, xl, xr, yres, -2.0, 2.0);


long int totevents = 0;
void gethd(const char *analysisfile);
void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3);
void setcolor(TH1D* h, int kcolor);

void performance () {

  gStyle->SetTitleFontSize(0.1);
  //gStyle->SetOptStat(0);

  
  cout << "line 72\n";
  for (int i = 0; i < filenames.size(); i++)    { 
      gethd(filenames[i]);
        if (skip[i] == true) {
           //htvec.push_back(nullptr);
           //hdvec.push_back(nullptr);
           //dtvec.push_back(nullptr);
           continue;}
        htvec.push_back((TH1D*)jetEnergyHT->Clone(hds[i]));
        hdvec.push_back((TH1D*)jetEnergyHD->Clone(dts[i]));
        dtvec.push_back((TH1D*)jetEnergyDT->Clone(hts[i]));
      }
 
   long double avgevents = (long double)totevents/20.0;
   cout << "	avgevents: " << avgevents <<"\n";

   for (int i = 0; i < htvec.size(); i++) {  
     //if (skip[i] == true)   continue;
       htvec[i]->Scale(avgevents, "width");
       hdvec[i]->Scale(avgevents, "width");
       dtvec[i]->Scale(avgevents, "width");
       if (htvec[i]->Integral() == 0)   continue;
       htvec[i]->Scale(1.0/htvec[i]->Integral());
       hdvec[i]->Scale(1.0/hdvec[i]->Integral());
       dtvec[i]->Scale(1.0/dtvec[i]->Integral());
   }

   HT->Scale(avgevents, "width");
   HD->Scale(avgevents, "width");
   DT->Scale(avgevents, "width");

   HT->Scale(1.0/HT->Integral());
   HD->Scale(1.0/HD->Integral());
   DT->Scale(1.0/DT->Integral());

  residualht->Scale(avgevents);
  residualhd->Scale(avgevents);
  residualdt->Scale(avgevents);

  //Plot projections of residuals by pT
  TH1D *residualslicedt = (TH1D*)residualdt->Clone("residualslicedt");
  TH1D *residualslicehd = (TH1D*)residualhd->Clone("residualslicehd");
  TH1D* JESov = new TH1D("JESov", "", xres, xl, xr);
  TH1D* JERov = new TH1D("JERov", "", xres, xl, xr);
    stringstream ssht;
    for (int i = 1; i <= xres; i++)  {   
      ssht << "slice" << i;
      TH1D* residualsliceht = (TH1D*)residualht->ProjectionY(ssht.str().c_str(), i, i);   
      JESov->SetBinContent(i, residualsliceht->GetMean());
      JESov->SetBinError(i, residualsliceht->GetMeanError());
      JERov->SetBinContent(i, residualsliceht->GetStdDev());
      JERov->SetBinError(i, residualsliceht->GetRMSError());
      residualsliceht->Reset();
      }
    setcolor(JESov, kBlack);
    setcolor(JERov, kBlack);
  TH1D* JESde = new TH1D("JESde", "", xres, xl, xr);
  TH1D* JERde = new TH1D("JERde", "", xres, xl, xr);
    stringstream ssdt;
    for (int i = 1; i <= xres; i++)  {
      ssdt << "slice" << i;
      TH1D* residualslicedt = (TH1D*)residualdt->ProjectionY(ssdt.str().c_str(), i, i);   
      JESde->SetBinContent(i, residualslicedt->GetMean());
      JESde->SetBinError(i, residualslicedt->GetMeanError());
      JERde->SetBinContent(i, residualslicedt->GetStdDev());
      JERde->SetBinError(i, residualslicedt->GetRMSError());
      residualslicedt->Reset();
      }
    setcolor(JESde, kBlue);
    setcolor(JERde, kBlue);
  TH1D* JESba = new TH1D("JESba", "", xres, xl, xr);;
  TH1D* JERba = new TH1D("JERba", "", xres, xl, xr);
    stringstream sshd;
    for (int i = 1; i <= xres; i++)  {
      sshd << "slice" << i;
      TH1D* residualslicehd = (TH1D*)residualhd->ProjectionY(sshd.str().c_str(), i, i);   
      JESba->SetBinContent(i, residualslicehd->GetMean());
      JESba->SetBinError(i, residualslicehd->GetMeanError());
      JERba->SetBinContent(i, residualslicehd->GetStdDev());
      JERba->SetBinError(i, residualslicehd->GetRMSError());
      residualslicehd->Reset();
      }
    setcolor(JESba, kRed);
    setcolor(JERba, kRed);
  cout << "line 151\n";


  setstyle(htvec, hdvec, dtvec);
   
/*
//Differential Scaling
y0->Scale(1.0/y0->Integral(), "width");
y5->Scale(1.0/y5->Integral(), "width");
y7->Scale(1.0/y7->Integral(), "width");
*/
  

  HT->SetLineWidth(3);
  HT->SetMarkerStyle(20);
  HT->SetMarkerSize(0.5);
  HD->SetLineWidth(3);
  HD->SetMarkerStyle(20);
  HD->SetMarkerSize(0.5);
  DT->SetLineWidth(3);
  DT->SetMarkerStyle(20);
  DT->SetMarkerSize(0.5);
  
    //Delta pT
    pTht->SetLineColor(kBlack);
    pTht->SetMarkerColor(kBlack);
    pTht->SetLineWidth(1);
    pTht->SetMarkerSize(0.5);
    pTht->SetMarkerStyle(20);
  

  
  string jesht = Form("JES = %.1f", HT->GetMean());
  const char *JES_HT = jesht.c_str();
  string jeshd = Form("JES = %.1f", HD->GetMean());
  const char *JES_HD= jeshd.c_str();
  string jesdt = Form("JES = %.1f", DT->GetMean());
  const char *JES_DT = jesdt.c_str();
  
  string jerht = Form("JER = %.1f", HT->GetStdDev());
  const char *JER_HT = jerht.c_str();
  string jerhd = Form("JER = %.1f", HD->GetStdDev());
  const char *JER_HD= jerhd.c_str();
  string jerdt = Form("JER = %.1f", DT->GetStdDev());
  const char *JER_DT = jerdt.c_str();
  


  //================================
  //=======  Legends  ==============
  //================================

  auto datarun = new TLegend(0.35, 0.725, 0.55, 0.85);
    datarun->SetTextSize(0.025);
    datarun->SetBorderSize(0);
    datarun->AddEntry((TObject*)0, "LHC18q + LHC18r", "");
    datarun->AddEntry((TObject*)0, "Pb-Pb #sqrt{s_{NN} = 5.02 TeV}", "");
    datarun->AddEntry((TObject*)0, "charged jets, anti-k_{T}, R = 0.4", "");
  
  auto residues = new TLegend(0.5, 0.6, 0.85, 0.85);
    residues->SetTextSize(0.03);
    residues->SetBorderSize(0);
    residues->AddEntry(JESov, "Overall", "pl");
    residues->AddEntry(JESba, "Fluctuations", "pl");
    residues->AddEntry(JESde, "Detector", "pl");

 
  //Key: Overall JES, JER 
  auto legHTS = new TLegend(0.45, 0.75, 0.65, 0.85);
    legHTS->SetTextSize(0.03);
    legHTS->SetBorderSize(0);
    legHTS->AddEntry((TObject*)0, JES_HT, "");
    auto legHTR = new TLegend(.65, 0.75, 0.85, 0.85);
    legHTR->SetTextSize(0.03);
    legHTR->SetBorderSize(0);
    legHTR->AddEntry((TObject*)0, JER_HT, "");
  auto legHDS = new TLegend(0.45, 0.75, 0.65, 0.85);
    legHDS->SetTextSize(0.035);
    legHDS->SetBorderSize(0);
    legHDS->AddEntry((TObject*)0, JES_HD, "");
    auto legHDR = new TLegend(0.65, 0.75, 0.85, 0.85);
    legHDR->SetTextSize(0.035);
    legHDR->SetBorderSize(0);
    legHDR->AddEntry((TObject*)0, JER_HD, "");
  auto legDTS = new TLegend(0.45, 0.75, 0.65, 0.85);
    legDTS->SetTextSize(0.035);
    legDTS->SetBorderSize(0);
    legDTS->AddEntry((TObject*)0, JES_DT, "");
    auto legDTR = new TLegend(0.65, 0.75, 0.85, 0.85);
    legDTR->SetTextSize(0.035);
    legDTR->SetBorderSize(0);
    legDTR->AddEntry((TObject*)0, JER_DT, "");
 

  double yd = 0.15;
  double yu = 0.75;  

  //Key: pT hard bin column
  double xl_pT = 0.55;
  double xr_pT = 0.65;
  auto hardbinHT = new TLegend(xl_pT, yd, xr_pT, yu);
    hardbinHT->SetTextSize(0.03);
    hardbinHT->SetBorderSize(0);
    hardbinHT->AddEntry((TObject*)0, "#bf{Bin}", "");
    for (int i = 0; i < htvec.size(); i++)   hardbinHT->AddEntry(htvec[i], Form("%d", i+1), "pl");
  auto hardbinDT = new TLegend(xl_pT, yd, xr_pT, yu);
    hardbinDT->SetTextSize(0.03);
    hardbinDT->SetBorderSize(0);
    hardbinDT->AddEntry((TObject*)0, "#bf{Bin}", "");
    for (int i = 0; i < dtvec.size(); i++)   hardbinDT->AddEntry(dtvec[i], Form("%d", i+1), "pl"); 
  auto hardbinHD = new TLegend(xl_pT, yd, xr_pT, yu);
    hardbinHD->SetTextSize(0.03);
    hardbinHD->SetBorderSize(0);
    hardbinHD->AddEntry((TObject*)0, "#bf{Bin}", "");
    for (int i = 0; i < hdvec.size(); i++)   hardbinHD->AddEntry(hdvec[i], Form("%d", i+1), "pl"); 
  
  //Key: JES Column
  double xl_jes = xr_pT;
  double xr_jes = 0.75;
  auto legJESht = new TLegend(xl_jes, yd, xr_jes, yu);
    legJESht->SetTextSize(0.03);
    legJESht->SetBorderSize(0); 
    legJESht->AddEntry((TObject*)0, "#bf{Scale}", "");
    for (int i = 0; i < htvec.size(); i++)    legJESht->AddEntry(htvec[i], Form("%.2f", htvec[i]->GetMean()), "");
  auto legJEShd = new TLegend(xl_jes, yd, xr_jes, yu);
    legJEShd->SetTextSize(0.03);
    legJEShd->SetBorderSize(0); 
    legJEShd->AddEntry((TObject*)0, "#bf{Scale}", "");
    for (int i = 0; i < hdvec.size(); i++)     legJEShd->AddEntry(hdvec[i], Form("%.2f", hdvec[i]->GetMean()), "");
  auto legJESdt = new TLegend(xl_jes, yd, xr_jes, yu);
    legJESdt->SetTextSize(0.03);
    legJESdt->SetBorderSize(0); 
    legJESdt->AddEntry((TObject*)0, "#bf{Scale}", "");
    for (int i = 0; i < dtvec.size(); i++)     legJESdt->AddEntry(dtvec[i], Form("%.2f", dtvec[i]->GetMean()), "");
  cout << "line 274\n";
  

  //Key: JER Column
  double xl_jer = xr_jes;
  auto legJERht = new TLegend(xl_jer, yd, 0.85, yu);
    legJERht->SetTextSize(0.03);
    legJERht->SetBorderSize(0); 
    legJERht->AddEntry((TObject*)0, "#bf{Res}", "");
    for (int i = 0; i < htvec.size(); i++)    legJERht->AddEntry(htvec[i], Form("%.2f", htvec[i]->GetStdDev()), "");
  auto legJERhd = new TLegend(xl_jer, yd, 0.85, yu);
    legJERhd->SetTextSize(0.03);
    legJERhd->SetBorderSize(0); 
    legJERhd->AddEntry((TObject*)0, "#bf{Res}", "");
    for (int i = 0; i < hdvec.size(); i++)     legJERhd->AddEntry(hdvec[i], Form("%.2f", hdvec[i]->GetStdDev()), "");
  auto legJERdt = new TLegend(xl_jer, yd, 0.85, yu);
    legJERdt->SetTextSize(0.03);
    legJERdt->SetBorderSize(0); 
    legJERdt->AddEntry((TObject*)0, "#bf{Res}", "");
    for (int i = 0; i < dtvec.size(); i++)     legJERdt->AddEntry(dtvec[i], Form("%.2f", dtvec[i]->GetStdDev()), "");





  //===============================================
  //==============  Draw Objects  =================
  //===============================================

  // y0->SetYTitle("#frac{1}{N_{jet}}#frac{dN_{jet}}{dp_{T}}"); 
  TCanvas *c1 = new TCanvas("c1", "Hybrid <-> Truth", 500, 500);
   c1->cd();
   gStyle->SetOptStat(0);
   TPad *pad1 = new TPad("pad1", "pad1", 0.0, 0.05, 1.0, 1.0);
   pad1->SetBottomMargin(0.1); // Upper and lower plot are joined
   pad1->SetLeftMargin(0.15);
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   //pad1->SetLogy();
   //jetEnergyHT->SetXTitle("p_{T}^{hybrid}");
   HT->SetXTitle("p_{T}^{hybrid}/p_{T}^{truth}");
   HT->SetMaximum(0.5);
   HT->Draw("same");
   legHTS->Draw("same");
   legHTR->Draw("same");
   hardbinHT->Draw("same");
   legJESht->Draw("same");
   legJERht->Draw("same");
   for (int i = 0; i < htvec.size(); i++)   htvec[i]->Draw("same");
   HT->Draw("hist same");
   //datarun->Draw("same");

  TCanvas *c2 = new TCanvas("c2", "Detector <-> Truth", 500, 500);
   c2->cd();
   gStyle->SetOptStat(0);
   TPad *pad2 = new TPad("pad2", "pad2", 0.0, 0.05, 1.0, 1.0);
   pad2->SetBottomMargin(0.1);
   pad2->SetLeftMargin(0.15);
   pad2->Draw();
   pad2->cd();
   //pad2->SetLogy();
   DT->SetXTitle("p_{T}^{detector}/p_{T}^{truth}");
   DT->SetMaximum(0.5);
   DT->Draw("same");
   legDTS->Draw("same");
   legDTR->Draw("same");
   hardbinDT->Draw("same");
   legJESdt->Draw("same");
   legJERdt->Draw("same");
   for (int i = 0; i < dtvec.size(); i++)   hdvec[i]->Draw("same");       
   DT->Draw("hist same");
   //datarun->Draw("same");

  TCanvas *c3 = new TCanvas("c3", "Hybrid <-> Detector", 500, 500);
   c3->cd();
   gStyle->SetOptStat(0);
   TPad *pad3 = new TPad("pad3", "pad3", 0.0, 0.05, 1.0, 1.0);
   pad3->SetBottomMargin(0.1);
   pad3->SetLeftMargin(0.15);
   pad3->Draw();
   pad3->cd();
   //pad3->SetLogy();
   HD->SetXTitle("p_{T}^{hybrid}/p_{T}^{detector}");
   //jetEnergyHD->SetYTitle("p_{T}^{htector}");
   HD->SetMaximum(0.5);
   HD->Draw("same");
   legHDS->Draw("same");
   legHDR->Draw("same");
   hardbinHD->Draw("same");
   legJEShd->Draw("same");
   legJERhd->Draw("same");
   for (int i = 0; i < dtvec.size(); i++)    hdvec[i]->Draw("same");
   HD->Draw("hist same");
   //datarun->Draw("same");

  TCanvas *c4 = new TCanvas("c4", "Residual: Hybrid <-> Truth", 1200, 500);
   c4->cd();
   gStyle->SetOptStat(0);
   TPad *pad4 = new TPad("pad4", "pad4", 0.0, 0.05, 0.33, 1.0);
   pad4->SetBottomMargin(0.1);
   pad4->SetLeftMargin(0.15);
   pad4->SetRightMargin(0.15);
   pad4->Draw();
   pad4->cd();
   pad4->SetLogz();
   residualht->SetXTitle("p_{T}^{truth}");
   residualht->SetYTitle("#frac{p_{T}^{hyrbid}-p_{T}^{truth}}{p_{T}^{truth}}");
   residualht->Draw("colz");
     JESov->Draw("same");
     c4->cd();
   TPad *pad5 = new TPad("pad5", "pad5", 0.33, 0.05, 0.67, 1.0);
   pad5->SetBottomMargin(0.1);
   pad5->SetLeftMargin(0.15);
   pad5->SetRightMargin(0.15);
   pad5->Draw();
   pad5->cd();
   pad5->SetLogz();
   residualhd->SetXTitle("p_{T}^{detector}");
   residualhd->SetYTitle("#frac{p_{T}^{hyrbid}-p_{T}^{detector}}{p_{T}^{detector}}");
   residualhd->Draw("colz");
   JESba->Draw("same");
     c4->cd();
   TPad *pad6 = new TPad("pad6", "pad6", 0.67, 0.05, 1.0, 1.0);
   pad6->SetBottomMargin(0.1);
   pad6->SetLeftMargin(0.15);
   pad6->SetRightMargin(0.15);
   pad6->Draw();
   pad6->cd();
   pad6->SetLogz();
   residualdt->SetXTitle("p_{T}^{truth}");
   residualdt->SetYTitle("#frac{p_{T}^{detector}-p_{T}^{truth}}{p_{T}^{truth}}");
   residualdt->Draw("colz");
     JESde->Draw("same");


  TLine *line = new TLine(10.0, 0.0, 100.0, 0.0);
  line->SetLineStyle(2);
  TCanvas *c7 = new TCanvas("c7", "Jet Energy Scale", 500, 500);
   c7->cd();
   gStyle->SetOptStat(0);
   TPad *pad7 = new TPad("pad7", "pad7", 0.0, 0.05, 1.0, 1.0);
   pad7->SetBottomMargin(0.1);
   pad7->SetLeftMargin(0.15);
   pad7->SetTicky();
   pad7->Draw();
   pad7->cd();
   JESov->SetXTitle("p_{T}^{true}");
   JESov->SetYTitle("JES");
   JESov->SetMaximum(1.0);
   JESov->SetMinimum(-1.0);
   JESov->Draw("e");
   JESde->Draw("e same");
   JESba->Draw("e same");
   line->Draw("same");
     residues->Draw("same"); 


 
  //Plot Jet Energy Resolution
  TCanvas *c8 = new TCanvas("c8", "Jet Energy Resolution", 500, 500);
   c8->cd();
   gStyle->SetOptStat(0);
   TPad *pad8 = new TPad("pad8", "pad8", 0.0, 0.05, 1.0, 1.0);
   pad8->SetBottomMargin(0.1);
   pad8->SetLeftMargin(0.15);
   pad8->SetTicky();
   pad8->Draw();
   pad8->cd();
   JERov->SetXTitle("p_{T}^{true}");
   JERov->SetYTitle("JER");
   JERov->SetMaximum(0.5);
   JERov->SetMinimum(0.0);
   JERov->Draw("e");
   JERde->Draw("e same");
   JERba->Draw("e same");
     residues->Draw("same");

 //Plot Delta pT
 TCanvas *c9 = new TCanvas("c9", "#Delta pT", 500, 500);
   c9->cd();
   pTht->Draw();
     pTht->SetXTitle("p_{T}^{hybrid} - p_{T}^{true}");
   auto dpt = new TLegend(0.62, 0.65, 0.85, 0.85);
     dpt->SetBorderSize(0);
     dpt->SetTextSize(0.04);
     dpt->AddEntry((TObject*)0, Form("#bf{Mean}: %.2f", pTht->GetMean()), "");
     dpt->AddEntry((TObject*)0, Form("#bf{Width}: %.2f", pTht->GetStdDev()), "");
     dpt->Draw("same");
   auto specs = new TLegend(0.15, 0.65, 0.35, 0.85); 
     specs->SetBorderSize(0);
     specs->SetTextSize(0.0275);
     specs->AddEntry((TObject*)0, "40 GeV < p_{T}^{det} < 120 GeV", "");
     //specs->Draw("same");

cout << "Std Dev Delta pT: " << pTht->GetStdDev() <<"\n";
cout << "Mean: " << pTht->GetMean() <<"\n";

}
  
void gethd(const char *analysisfile) {

    jetEnergyHT->Reset();
    jetEnergyHD->Reset();
    jetEnergyDT->Reset();


    string rootfilename = analysisfile;
    const char *lerootfile = rootfilename.c_str(); 
    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);

    AliEmcalList *EmbeddingHelper = (AliEmcalList*)lefile->Get("AliAnalysisTaskEmcalEmbeddingHelper_histos");    
    TProfile *pThardscale = (TProfile*)EmbeddingHelper->FindObject("fHistXsection");
    TH1D *pThardbin = (TH1D*)EmbeddingHelper->FindObject("fHistTrials");
    TH1D *eventcount = (TH1D*)EmbeddingHelper->FindObject("fHistEventCount");
    //Determine pT hard scaling
    double xsect = pThardscale->GetBinContent(pThardbin->GetMean() + 1);
    double xscale = xsect*pThardscale->GetEntries();
    double trials = pThardbin->GetBinContent(pThardbin->GetMean() + 1);
 
    double n_event = eventcount->GetBinContent(1);
    totevents += n_event;
    double pThardFrac = xscale/(n_event*trials);
     


    string qnfilename = "../ESE/AnalysisResults7812.root";             // R=0.2, 30-50%, pass1
    const char *leqnfile = qnfilename.c_str(); 

    //Load File/Tree into system
    TFile *qnfile = new TFile(leqnfile);

      //Calculate q2 percentiles
      TDirectoryFile *qndirfile = (TDirectoryFile*)qnfile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qndirfile->Get("coutputQnVectorTender_Q2V0C_EPV0C");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                             //number percentile we're on
      double centl = 30.0;
      double centr = 50.0;
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

    float hybpT, detpT, trupT;
    float centrality, q2, epangle, epPar;

    TTree* T = nullptr;
    lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_hybJet",T);
    T->SetBranchAddress("Jet_Pt", &hybpT);
    T->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &detpT);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &trupT);
    T->SetBranchAddress("Event_Centrality", &centrality);
    T->SetBranchAddress("Event_Q2Vector", &q2);
    T->SetBranchAddress("Jet_EPangle", &epangle);
    T->SetBranchAddress("Jet_MC_MatchedPartLevelJet_EPangle", &epPar);

    int entries	= T->GetEntries();
    int passnum = 0;
    //cout << "Entries: " << entries <<"\n";
    
    const int extractBins = 9;
    double extractFrac[extractBins] = {0.0, 0.01, 0.03, 0.15, 0.20, 0.20, 0.20, 0.15, 0.05};
    double extractEdge[extractBins + 1] = {-20.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 140.0, 200.0};
    for (int i = 0; i < entries; i++) {
        T->GetEntry(i);
           int scaleindex = -1;
           if (hybpT >= extractEdge[0] && hybpT < extractEdge[1])  scaleindex = 0;
           if (hybpT >= extractEdge[1] && hybpT < extractEdge[2])  scaleindex = 1;
           if (hybpT >= extractEdge[2] && hybpT < extractEdge[3])  scaleindex = 2;
           if (hybpT >= extractEdge[3] && hybpT < extractEdge[4])  scaleindex = 3;
           if (hybpT >= extractEdge[4] && hybpT < extractEdge[5])  scaleindex = 4;
           if (hybpT >= extractEdge[5] && hybpT < extractEdge[6])  scaleindex = 5;
           if (hybpT >= extractEdge[6] && hybpT < extractEdge[7])  scaleindex = 6;
           if (hybpT >= extractEdge[7] && hybpT < extractEdge[8])  scaleindex = 7;
           if (hybpT >= extractEdge[8] && hybpT < extractEdge[9])  scaleindex = 8;
        if (scaleindex == 0)   continue;
        if (hybpT < 10.0 || trupT < 10.0 /*|| detpT < 10.0*/)   continue;
        if (centrality < 30.0 || centrality > 50.0)       continue;
//                   if (q2 > allPercent[int(centrality-centl)][2])    continue;  //20th percentile
//                   if (q2 < allPercent[int(centrality-centl)][6])    continue;  //80th percentile
        //if (abs(cos(epangle)) < sqrt(3)/2.0) continue;  //in-plane
        //if ( detpT < 40.0 || detpT > 120.0)   continue;
        //if (trupT < 40.0 || trupT > 60.0)   continue;
          passnum++;
        jetEnergyHT->Fill(hybpT/trupT, pThardFrac/extractFrac[scaleindex]);
        jetEnergyHD->Fill(hybpT/detpT, pThardFrac/(extractFrac[scaleindex])); 
        jetEnergyDT->Fill(detpT/trupT, pThardFrac/(extractFrac[scaleindex]));
 
        pTht->Fill(hybpT-trupT, pThardFrac/extractFrac[scaleindex]);

        residualht->Fill(trupT, (hybpT-trupT)/trupT, pThardFrac/(extractFrac[scaleindex]));
        residualhd->Fill(detpT, (hybpT-detpT)/detpT, pThardFrac/(extractFrac[scaleindex]));
        residualdt->Fill(trupT, (detpT-trupT)/trupT, pThardFrac/extractFrac[scaleindex]);
   }
  
  HT->Add(jetEnergyHT);
  HD->Add(jetEnergyHD);
  DT->Add(jetEnergyDT);

  if (passnum == 0)  {
     skip.push_back(true);
   }

  if (passnum != 0)    {
    skip.push_back(false);
   
  jetEnergyHT->Scale(1.0/jetEnergyHT->Integral());
  jetEnergyHD->Scale(1.0/jetEnergyHD->Integral());
  jetEnergyDT->Scale(1.0/jetEnergyDT->Integral());
  }
  double frac = (double)passnum/(double)entries;
  cout << "pT hard hd: " << pThardbin->GetMean() <<"\n";
  cout << "	Pass Frac: " << Form("%.5f", frac) <<"\n";

  lefile->Close();

}


void setcolor(TH1D* h, int kcolor) {
  h->SetLineColor(kcolor);
  h->SetMarkerColor(kcolor);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.7);
}


void setstyle(vector<TH1D*> vec1, vector<TH1D*> vec2, vector<TH1D*> vec3) {
  vector<int> colors = {kRed+2, kRed, kOrange+10, kOrange-4, kYellow, kSpring+10, kSpring, kGreen-3, kGreen+2, kTeal+6, kTeal, kCyan-3, kAzure+7, kBlue+1, kBlue-2, kViolet+2, kViolet, kMagenta-4, kPink+6, kPink+10, kBlack};
  for (int i = 0; i < vec1.size(); i++) {
      if (skip[i] == 0)   continue;
      vec1[i]->SetLineColor(colors[i]);
      vec1[i]->SetMarkerColor(colors[i]);
        vec1[i]->SetLineWidth(1);
        vec1[i]->SetMarkerSize(0.5);
        vec1[i]->SetMarkerStyle(20);
      vec2[i]->SetLineColor(colors[i]);
      vec2[i]->SetMarkerColor(colors[i]);
        vec2[i]->SetLineWidth(1);
        vec2[i]->SetMarkerSize(0.5);
        vec2[i]->SetMarkerStyle(20);
      vec3[i]->SetLineColor(colors[i]);
      vec3[i]->SetMarkerColor(colors[i]);
        vec3[i]->SetLineWidth(1);
        vec3[i]->SetMarkerSize(0.5);
        vec3[i]->SetMarkerStyle(20);
  }
}

