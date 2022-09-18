#include "TRandom.h"
#include "TRandom3.h"
  
void QAplots (int lql = 0, int lqr = 3, int hql = 7, int hqr = 10, int R = 4)  {

   //gStyle->SetPalette(kCool);  
   TRandom3* rand = new TRandom3();   

   double centl = 30.0;
   double centr = 50.0;

   TString cFiles2="files21.txt";

  //============================================
  //==========   Get Raw Spectra    ============
  //============================================
   //Raw Spectra from Trains
    string rootfilename;
    if (R == 2) rootfilename = "AnalysisResults8144.root";
    if (R == 4) rootfilename = "AnalysisResults8294.root";
    const char *lerootfile = rootfilename.c_str(); 

    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);
    TTree *T = nullptr;
    if (R == 2)  lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
    if (R == 4)  lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR040_tracks_pT0150_pt_scheme_Rho_Jet", T);
    float jetpT, centra, qvec, ep, N;
    T->SetBranchAddress("Jet_Pt", &jetpT);
    T->SetBranchAddress("Event_Centrality", &centra);
    T->SetBranchAddress("Event_Q2VectorV0M", &qvec);    
    T->SetBranchAddress("Jet_EPangleV0M", &ep);    
    T->SetBranchAddress("Jet_Eta", &N);

    //blank histogram for plot scaling
    TH2F* rangeup = new TH2F("rangeup", "rangeup", 20, 30., 50., 20000, 0.0, 20.0);




  //=====================================================
  //==========   Calculate q2 Percentiles    ============
  //=====================================================

     //Get Hists from File
     TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2_EP");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrFullV0");
      //plots of cos(2(w1-w2)) for different q2 measurements
      TH3F *resoACqA = (TH3F*)qnlist->FindObject("fHistResolution_epV0AV0C_qV0A");
      TH3F *resoCTqA = (TH3F*)qnlist->FindObject("fHistResolution_epV0CTPC_qV0A");
      TH3F *resoATqA = (TH3F*)qnlist->FindObject("fHistResolution_epV0ATPC_qV0A");
      TH3F *resoACqC = (TH3F*)qnlist->FindObject("fHistResolution_epV0AV0C_qV0C");
      TH3F *resoCTqC = (TH3F*)qnlist->FindObject("fHistResolution_epV0CTPC_qV0C");
      TH3F *resoATqC = (TH3F*)qnlist->FindObject("fHistResolution_epV0AV0C_qV0C");
      TH3F *resoTTqM = (TH3F*)qnlist->FindObject("fHistResolution_epTPCpTPCn_qV0M");
      TH3F *resoMTpqM = (TH3F*)qnlist->FindObject("fHistResolution_epV0MTPCp_qV0M");
      TH3F *resoMTnqM = (TH3F*)qnlist->FindObject("fHistResolution_epV0MTPCn_qV0M");

      //plots to check flatness of EP
      TH1F *EPcalibV0M = (TH1F*)qnlist->FindObject("fHistEventPlaneFullV0");
      TH1F *EPcalibV0A = (TH1F*)qnlist->FindObject("fHistEventPlaneV0A");
      TH1F *EPcalibV0C = (TH1F*)qnlist->FindObject("fHistEventPlaneV0C");
      //normalize plots
      EPcalibV0M->Scale(1.0/EPcalibV0M->Integral());
      EPcalibV0A->Scale(1.0/EPcalibV0A->Integral());
      EPcalibV0C->Scale(1.0/EPcalibV0C->Integral());


      //initialize variables used to calculate percentiles
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      double centdiff = (centr-centl);
      int leftc, rightc;
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < (double)centdiff; j++) {
        //reset percentile ticker for each centrality
        long int percentileticker = 0;                                //number of entries, starting from lowest q2 and counting upward
        int percentilen = 1;                                          //number percentile we're on
        //update centrality
        leftc = (int)centl + j; 
        rightc = (int)centl + j + 1;    
        //cout << "Centrality " << leftc << "-" << rightc <<"\n"; 
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);       //restrict 2D hist to current 1% centrality bin
        //project slice onto q2 axis
        TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();           
        TH1F *q2copy = (TH1F*)q2hist2D->ProjectionY();
        //cout << "	Events in centrality bin: " << q2hist->GetEntries() <<"\n";
        long double percentilemarker = (long double)q2hist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        //increment along q2 axis until percentile is surpassed
        for (int i = 1; i <= q2hist->GetNbinsX(); i++)   {  
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

      //centrality independent percentiles for q2 plot
      vector<double> centind(0);
      //project onto q2 axis
      q2hist2D->GetXaxis()->SetRangeUser(30.0, 50.0);
      TH1F* q2centind = (TH1F*)q2hist2D->ProjectionY();
      int pticker = 0;
      int pn = 1;
      double pmarker = (double)q2centind->GetEntries()/10.0; 
      vector<TLine*> lines(0);
      for (int i = 1; i <= q2centind->GetNbinsX(); i++) {
          pticker += q2centind->GetBinContent(i);
          if (pticker >= pn*pmarker) {
             centind.push_back(q2centind->GetBinCenter(i));
             cout << "	Percentile: " << pn << ", q2: " << q2centind->GetBinCenter(i) <<"\n";
             lines.push_back(new TLine(centind[pn-1], 0.01, centind[pn-1], 0.99));
             lines[pn-1]->SetLineStyle(2);
             pn++;
          }
      }


   //centrality independent q2 percentiles test
   vector<double> cindq2test(0);
    for (int j = 0; j < 10; j++)  {
       double ctot = 0;
       for (int i = 1; i <= centdiff; i++) {
           ctot += allPercent[i-1][j];}
       double avg = ctot/20.;
       cindq2test.push_back(avg);  }
        for (int k = 0; k < 10; k++)  cout << "PercentileTest " << k+1 << ": " << cindq2test[k] << "\n";

     TCanvas *ci = new TCanvas("ci", "ci", 500, 500);
     ci->cd();
       q2centind->Draw("same");



  //=====================================================
  //==========   Calculate EP Resolution    =============
  //=====================================================

      //calculate R2
      // res1 = V0C - V0A
      // res2 = V0C - TPC
      // res3 = V0A - TPC
      // res4 = TPCp - TPCn
      // res5 = V0M - TPCp
      // res6 = V0M - TPCn

     TH1F* hAClow = new TH1F("hAClow", "hAClow", 200, -1., 1.);
        TH1F* hAChigh = new TH1F("hAChigh", "hAChigh", 200, -1., 1.); 
        TH1F* hCTlow = new TH1F("hCTlow", "hCTlow", 200, -1., 1.);
        TH1F* hCThigh = new TH1F("hCThigh", "hCThigh", 200, -1., 1.);
        TH1F* hATlow = new TH1F("hATlow", "hATlow", 200, -1., 1.);
        TH1F* hAThigh = new TH1F("hAThigh", "hAThigh", 200, -1., 1.);
      TH1F* hTTlow = new TH1F("hTTlow", "hTTlow", 200, -1., 1.);
        TH1F* hTThigh = new TH1F("hTThigh", "hTThigh", 200, -1., 1.);
        TH1F* hMTplow = new TH1F("hMTplow", "hMTplow", 200, -1., 1.);
        TH1F* hMTphigh = new TH1F("hMTphigh", "hMTphigh", 200, -1., 1.);
        TH1F* hMTnlow = new TH1F("hMTnlow", "hMTnlow", 200, -1., 1.);
        TH1F* hMTnhigh = new TH1F("hMTnhigh", "hMTnhigh", 200, -1., 1.);

      //iterate through centrality bins
      for (int i = 31; i <= 50; i++)   {
           //restrict centrality to 1% bins
           resoACqA->GetZaxis()->SetRange(i, i+1);
           resoCTqA->GetZaxis()->SetRange(i, i+1); 
           resoATqA->GetZaxis()->SetRange(i, i+1);
           resoTTqM->GetZaxis()->SetRange(i, i+1);
           resoMTpqM->GetZaxis()->SetRange(i, i+1);
           resoMTnqM->GetZaxis()->SetRange(i, i+1);
           //project to q (y) vs resolution(cos2(w1-w2)) (x)
           TH2F* tempResvQ1 = (TH2F*)resoACqA->Project3D("yx");
           TH2F* tempResvQ2 = (TH2F*)resoCTqA->Project3D("yx");
           TH2F* tempResvQ3 = (TH2F*)resoATqA->Project3D("yx");
           TH2F* tempResvQ4 = (TH2F*)resoTTqM->Project3D("yx");
           TH2F* tempResvQ5 = (TH2F*)resoMTpqM->Project3D("yx");
           TH2F* tempResvQ6 = (TH2F*)resoMTnqM->Project3D("yx");
           //project onto resolution (x) axis for restricted high and low q2
           double q2value_lql;                              //low percentile left
             if (lql == 0)  q2value_lql = 0.0;
             if (lql != 0)  q2value_lql = allPercent[i-31][lql-1]; 
           double q2value_lqr = allPercent[i-31][lqr-1];    //low percentile right
           double q2value_hql = allPercent[i-31][hql-1];    //high percentile left
           double q2value_hqr;                              //high percentile right
             if (hqr != 10)  q2value_hqr = allPercent[i-31][hqr-1]; 
             if (hqr == 10)  q2value_hqr = 20.0;
           //projections where q determined using V0A
           TH1F* tempL1 = (TH1F*)tempResvQ1->ProjectionX("tempL1", tempResvQ1->GetYaxis()->FindBin(q2value_lql),  //
                                                                   tempResvQ1->GetYaxis()->FindBin(q2value_lqr));
             TH1F* tempH1 = (TH1F*)tempResvQ1->ProjectionX("tempH1", tempResvQ1->GetYaxis()->FindBin(q2value_hql), -1);
             TH1F* tempL2 = (TH1F*)tempResvQ2->ProjectionX("tempL2", tempResvQ1->GetYaxis()->FindBin(q2value_lql),
                                                                     tempResvQ2->GetYaxis()->FindBin(q2value_lqr));
             TH1F* tempH2 = (TH1F*)tempResvQ2->ProjectionX("tempH2", tempResvQ2->GetYaxis()->FindBin(q2value_hql), -1);
             TH1F* tempL3 = (TH1F*)tempResvQ3->ProjectionX("tempL3", tempResvQ1->GetYaxis()->FindBin(q2value_lql), 
                                                                     tempResvQ3->GetYaxis()->FindBin(q2value_lqr ));
             TH1F* tempH3 = (TH1F*)tempResvQ3->ProjectionX("tempH3", tempResvQ3->GetYaxis()->FindBin(q2value_hql), -1);
           //projections where q determined using V0M
           TH1F* tempL4 = (TH1F*)tempResvQ4->ProjectionX("tempL4", tempResvQ4->GetYaxis()->FindBin(q2value_lql),
                                                                     tempResvQ4->GetYaxis()->FindBin(q2value_lqr));
             TH1F* tempH4 = (TH1F*)tempResvQ4->ProjectionX("tempH4", tempResvQ4->GetYaxis()->FindBin(q2value_hql), -1);
             TH1F* tempL5 = (TH1F*)tempResvQ5->ProjectionX("tempL5", tempResvQ5->GetYaxis()->FindBin(q2value_lql), 
                                                                     tempResvQ5->GetYaxis()->FindBin(q2value_lqr));
             TH1F* tempH5 = (TH1F*)tempResvQ5->ProjectionX("tempH5", tempResvQ5->GetYaxis()->FindBin(q2value_hql), -1);
             TH1F* tempL6 = (TH1F*)tempResvQ6->ProjectionX("tempL6", tempResvQ6->GetYaxis()->FindBin(q2value_lql),
                                                                     tempResvQ6->GetYaxis()->FindBin(q2value_lqr));
             TH1F* tempH6 = (TH1F*)tempResvQ6->ProjectionX("tempH6", tempResvQ6->GetYaxis()->FindBin(q2value_hql), -1);
           //add temporary hists to final result
           hAClow->Add(tempL1);
             hAChigh->Add(tempH1);
             hCTlow->Add(tempL2);
             hCThigh->Add(tempH2);
             hATlow->Add(tempL3);
             hAThigh->Add(tempH3);
           hTTlow->Add(tempL4);
             hTThigh->Add(tempH4);
             hMTplow->Add(tempL5);
             hMTphigh->Add(tempH5);
             hMTnlow->Add(tempL6);
             hMTnhigh->Add(tempH6);
   }

   //Get Averages of Distributions
   double AClow = hAClow->GetMean();
     double AChigh = hAChigh->GetMean();
     double CTlow = hCTlow->GetMean();
     double CThigh = hCThigh->GetMean();
     double ATlow = hATlow->GetMean();
     double AThigh = hAThigh->GetMean();
   double TTlow = hTTlow->GetMean();
     double TThigh = hTThigh->GetMean();
     double MTplow = hMTplow->GetMean();
     double MTphigh = hMTphigh->GetMean();
     double MTnlow = hMTnlow->GetMean();
     double MTnhigh = hMTnhigh->GetMean();
   //get R2 value (R2^A = sqrt(AC*AB/BC))
   double R2Alow = sqrt(AClow*ATlow/CTlow);
   double R2Ahigh = sqrt(AChigh*AThigh/CThigh);
   double R2Clow = sqrt(AClow*CTlow/ATlow);
   double R2Chigh = sqrt(AChigh*CThigh/AThigh);
   double R2Mlow = sqrt(MTplow*MTnlow/TTlow);
   double R2Mhigh = sqrt(MTphigh*MTnhigh/TThigh);

   cout << "R2Mlow: " << R2Mlow <<"\n";
   cout << "R2Mhigh: " << R2Mhigh <<"\n";





  //==========================================================================
  //==========   Consider Resolution Plots for Other EP Combos   =============
  //==========================================================================


   //Get full 2D plots
   resoACqA->GetZaxis()->SetRange(31, 50);
   resoCTqA->GetZaxis()->SetRange(31, 50); 
   resoATqA->GetZaxis()->SetRange(31, 50);
   resoACqC->GetZaxis()->SetRange(31, 50);
   resoCTqC->GetZaxis()->SetRange(31, 50); 
   resoATqC->GetZaxis()->SetRange(31, 50);
   resoTTqM->GetZaxis()->SetRange(31, 50);
   resoMTpqM->GetZaxis()->SetRange(31, 50);
   resoMTnqM->GetZaxis()->SetRange(31, 50);
   vector<double> resMean1A (0);   //V0C-V0A
   vector<double> resMean2A (0);   //V0C-TPC
   vector<double> resMean3A (0);   //V0A-TPC
   vector<double> resMean1C (0);   //V0C-V0A
   vector<double> resMean2C (0);   //V0C-TPC
   vector<double> resMean3C (0);   //V0A-TPC
   vector<double> resMean1M (0);   //TPC-TPC
   vector<double> resMean2M (0);   //V0M-TPCp
   vector<double> resMean3M (0);   //V0M-TPCn
      //project to reolution(y) vs q(x)
      TH2F* tQ1vResA = (TH2F*)resoACqA->Project3D("xy");
      TH2F* tQ2vResA = (TH2F*)resoCTqA->Project3D("xy");
      TH2F* tQ3vResA = (TH2F*)resoATqA->Project3D("xy");
      TH2F* tQ1vResC = (TH2F*)resoACqC->Project3D("xy");
      TH2F* tQ2vResC = (TH2F*)resoCTqC->Project3D("xy");
      TH2F* tQ3vResC = (TH2F*)resoATqC->Project3D("xy");
      TH2F* tQ1vResM = (TH2F*)resoTTqM->Project3D("xy");
      TH2F* tQ2vResM = (TH2F*)resoMTpqM->Project3D("xy");
      TH2F* tQ3vResM = (TH2F*)resoMTnqM->Project3D("xy");
      //project to resolution hists to get mean values per q2
      for (int i = 0; i < tQ1vResA->GetNbinsX(); i++) {
          tQ1vResA->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ2vResA->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ3vResA->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ1vResC->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ2vResC->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ3vResC->GetXaxis()->SetRange(i, i+1);  //restrict q value before projecting
          tQ1vResM->GetXaxis()->SetRange(i, i+1);
          tQ2vResM->GetXaxis()->SetRange(i, i+1);
          tQ3vResM->GetXaxis()->SetRange(i, i+1);
          TH1F* t1A = (TH1F*)tQ1vResA->ProjectionY("t1A");
          TH1F* t2A = (TH1F*)tQ2vResA->ProjectionY("t2A");
          TH1F* t3A = (TH1F*)tQ3vResA->ProjectionY("t3A");
          TH1F* t1C = (TH1F*)tQ1vResC->ProjectionY("t1C");
          TH1F* t2C = (TH1F*)tQ2vResC->ProjectionY("t2C");
          TH1F* t3C = (TH1F*)tQ3vResC->ProjectionY("t3C");
          TH1F* t1M = (TH1F*)tQ1vResM->ProjectionY("t1M");
          TH1F* t2M = (TH1F*)tQ2vResM->ProjectionY("t2M");
          TH1F* t3M = (TH1F*)tQ3vResM->ProjectionY("t3M");
          //cout << resMean1[i] <<"\n";
          resMean1A.push_back(t1A->GetMean());
          resMean2A.push_back(t2A->GetMean());
          resMean3A.push_back(t3A->GetMean());
          resMean1C.push_back(t1C->GetMean());
          resMean2C.push_back(t2C->GetMean());
          resMean3C.push_back(t3C->GetMean());
          resMean1M.push_back(t1M->GetMean());
          resMean2M.push_back(t2M->GetMean());
          resMean3M.push_back(t3M->GetMean());
      }
    //fill hists with calculated resolution
    TH1F* qAepC = new TH1F("qAepC", "qAepC", tQ1vResA->GetNbinsX(), 0, 20);
    TH1F* qAepA = new TH1F("qAepA", "qAepA", tQ1vResA->GetNbinsX(), 0, 20);
    TH1F* qCepC = new TH1F("qCepC", "qCepC", tQ1vResC->GetNbinsX(), 0, 20);
    TH1F* qCepA = new TH1F("qCepA", "qCepA", tQ1vResC->GetNbinsX(), 0, 20);
    TH1F* qMepM = new TH1F("qMepM", "qMepM", tQ1vResM->GetNbinsX(), 0, 20);
    //avoid dividing by 0
    for (int i = 0; i < 130/*tQ1vRes->GetNbinsX()*/; i++) {
      //qAepC->SetBinContent(i, 0.5);
      //qAepA->SetBinContent(i, 0.5);
        if (resMean3A[i] != 0) qAepC->SetBinContent(i, sqrt(resMean1A[i]*resMean2A[i]/resMean3A[i]));
        if (resMean3A[i] == 0) qAepC->SetBinContent(i, 0);
        if (resMean2A[i] != 0) qAepA->SetBinContent(i, sqrt(resMean1A[i]*resMean3A[i]/resMean2A[i]));
        if (resMean2A[i] == 0) qAepA->SetBinContent(i, 0);
        if (resMean3C[i] != 0) qCepC->SetBinContent(i, sqrt(resMean1C[i]*resMean2C[i]/resMean3C[i]));
        if (resMean3C[i] == 0) qCepC->SetBinContent(i, 0);
        if (resMean2C[i] != 0) qCepA->SetBinContent(i, sqrt(resMean1C[i]*resMean3C[i]/resMean2C[i]));
        if (resMean2C[i] == 0) qCepA->SetBinContent(i, 0);
        if (resMean1M[i] != 0) qMepM->SetBinContent(i, sqrt(resMean2M[i]*resMean3M[i]/resMean1M[i]));
        if (resMean1M[i] == 0) qMepM->SetBinContent(i, 0);
    }   

      int kSunPink, kSunPurple, kSunOrange, kSunBlue;
      kSunPink = TColor::GetColor("#fc03ca");
      kSunPurple = TColor::GetColor("#7100b8");
      kSunOrange = TColor::GetColor("#2597fa");
      kSunBlue = TColor::GetColor("#251bb3");


   qAepC->SetLineColor(kSunPink);
      qAepC->SetLineWidth(3);
      qAepC->SetMarkerStyle(20);
      qAepC->SetMarkerSize(0.6);
      qAepC->SetMarkerColor(kSunPink);
   qAepA->SetLineColor(kSunPurple);
      qAepA->SetLineWidth(3);
      qAepA->SetMarkerStyle(20);
      qAepA->SetMarkerSize(0.6);
      qAepA->SetMarkerColor(kSunPurple);
   qCepC->SetLineColor(kSunOrange);
      qCepC->SetLineWidth(3);
      qCepC->SetMarkerStyle(20);
      qCepC->SetMarkerSize(0.6);
      qCepC->SetMarkerColor(kSunOrange);
   qCepA->SetLineColor(kSunBlue);
      qCepA->SetLineWidth(3);
      qCepA->SetMarkerStyle(20);
      qCepA->SetMarkerSize(0.6);
      qCepA->SetMarkerColor(kSunBlue);
  
    qMepM->SetLineColor(kAzure+2);
      qMepM->SetLineWidth(3);
      qMepM->SetMarkerStyle(20);
      qMepM->SetMarkerSize(0.6);
      qMepM->SetMarkerColor(kAzure+2);



   /* cout << "R2Alow: " << R2Alow <<"\n";
   cout << "R2Ahigh: " << R2Ahigh <<"\n";
   cout << "R2Clow: " << R2Clow <<"\n";
   cout << "R2Chigh: " << R2Chigh <<"\n";
   */
   vector<int> colors = {};


  TH1D *perc10 = new TH1D("perc10", "perc10", centdiff, centl, centr);
  TH1D *perc40 = new TH1D("perc40", "perc40", centdiff, centl, centr);
  TH1D *perc70 = new TH1D("perc70", "perc70", centdiff, centl, centr);
  for (int i = 1; i <= centdiff; i++) {
      perc10->SetBinContent(i, allPercent[i-1][0]);
      perc40->SetBinContent(i, allPercent[i-1][3]);
      perc70->SetBinContent(i, allPercent[i-1][6]);
  }
  perc10->SetLineColor(kRed);
    perc10->SetLineWidth(2);
  perc40->SetLineColor(kRed);
    perc40->SetLineWidth(2);
  perc70->SetLineColor(kRed);
    perc70->SetLineWidth(2);




  TLine *line = new TLine(0.0, 1, 1.0, 1);
  line->SetLineStyle(2);


  //--------------------------------------------------------------------------
/*
            cout << "q2hist2D->GetNbinsX(): " << q2hist2D->GetNbinsX() <<"\n";
            cout << "q2hist2D->GetNbinsY(): " << q2hist2D->GetNbinsY() <<"\n";
        for (int i = 1; i < q2hist2D->GetNbinsX(); i++)  {     //x bins
            for (int j = 1; j < q2hist2D->GetNbinsY(); j++)  { //y bins
                rangeup->SetBinContent(i, j, q2hist2D->GetBinContent(i, j));

                }
            }  
*/

  TLegend* perfleg = new TLegend(0.15, 0.7, 0.5, 0.85);
           perfleg->SetBorderSize(0);
           perfleg->SetTextSize(0.035);
           perfleg->AddEntry((TObject*)0, "Work in Progress, Pb#font[122]{-}Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
           perfleg->AddEntry(perc40, "10^{th}, 40^{th}, and 70^{th} #it{q}_{2}^{V0M} percentiles", "l");
  TLegend* q2perms = new TLegend(0.65, 0.15, 0.8, 0.4);
           q2perms->SetBorderSize(0);
           q2perms->SetTextSize(0.035);
           //q2perms->AddEntry(qCepA, "#it{R}_{2}^{V0A} (#it{q}_{2}^{V0C})", "pl");
           //q2perms->AddEntry(qAepA, "#it{R}_{2}^{V0A} (#it{q}_{2}^{V0A})", "pl");
           //q2perms->AddEntry(qCepC, "#it{R}_{2}^{V0C} (#it{q}_{2}^{V0C})", "pl");
           //q2perms->AddEntry(qAepC, "#it{R}_{2}^{V0C} (#it{q}_{2}^{V0A})", "pl");
           q2perms->AddEntry(qMepM, "#it{R}_{2}^{V0M} (#it{q}_{2}^{V0M})", "pl");

  TLegend* Rleg = new TLegend(0.6, 0.2, 0.8, 0.4);
           Rleg->SetBorderSize(0);
           Rleg->SetTextSize(0.045);
           if (R == 2) {
              Rleg->AddEntry((TObject*)0, "Train 8144", "");
              Rleg->AddEntry((TObject*)0, "R = 0.2", "");}
           if (R == 4) {
              Rleg->AddEntry((TObject*)0, "Train 8119", "");
              Rleg->AddEntry((TObject*)0, "R = 0.4", "");}

  TCanvas *c1 = new TCanvas("Q2 vs Cent", "Q2 vs Cent", 500, 500);
  c1->cd();
       gStyle->SetOptStat(0);
       //gStyle->SetPadTickX(1);
       //gStyle->SetPadTickY(1);
       TPad *pad1 = new TPad("pad1", "", 0.0, 0.05, 1.0, 1.0);
          pad1->SetBottomMargin(0.08); // Upper and lower plot are joined
          pad1->SetLeftMargin(0.1); // Upper and lower plot are joined
          pad1->SetRightMargin(0.15);
          pad1->SetTopMargin(0.05);
          pad1->Draw();
          pad1->cd();
          //pad1->SetLogz();
       //q2hist2D->GetYaxis()->SetMaximum(20.0);
       //q2hist2D->GetYaxis()->SetRangeUser(0.0, 20.0);
       q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
       q2hist2D->Draw("same colz");
       q2hist2D->SetTitle("");
       q2hist2D->GetXaxis()->SetTitle("Centrality (%)");
       q2hist2D->GetYaxis()->SetTitle("#it{q}_{2}^{V0M}");
       //rangeup->Draw("same colz");
       perc10->Draw("same");
       perc40->Draw("same");
       perc70->Draw("same");
       perfleg->Draw("same");

  TCanvas *c2 = new TCanvas("EP", "EP", 1200, 500);
  c2->cd();
       gStyle->SetOptStat(0);
       TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 0.33, 1.0);
         pad2->SetBottomMargin(0.15);
         pad2->SetLeftMargin(0.15);
         pad2->Draw();
       TPad *pad3 = new TPad("pad3", "", 0.33, 0.05, 0.67, 1.0);
         pad3->SetBottomMargin(0.15);
         pad3->SetLeftMargin(0.15);
         pad3->Draw();
       TPad *pad4 = new TPad("pad4", "", 0.67, 0.05, 1.0, 1.0);
         pad4->SetBottomMargin(0.15);
         pad4->SetLeftMargin(0.15);
         pad4->Draw();
   pad2->cd();
     EPcalibV0M->SetLineColor(kBlack);
     EPcalibV0M->SetLineWidth(2);
     //EPcalib->SetMarkerColor(kBlack);
     //EPcalib->SetMarkerStyle(20);
     EPcalibV0M->SetMinimum(0);
     EPcalibV0M->SetMaximum(0.006);
     EPcalibV0M->SetYTitle("Normalized Events");
     EPcalibV0M->Draw("same");
     Rleg->Draw("same");
     c2->cd();
  pad3->cd();
     EPcalibV0A->SetLineColor(kBlack);
     EPcalibV0A->SetLineWidth(2);
     //EPcalib->SetMarkerColor(kBlack);
     //EPcalib->SetMarkerStyle(20);
     EPcalibV0A->SetMinimum(0);
     EPcalibV0A->SetMaximum(0.006);
     EPcalibV0A->SetYTitle("");
     EPcalibV0A->Draw("same");
     Rleg->Draw("same");
     c2->cd();
  pad4->cd();
     EPcalibV0C->SetLineColor(kBlack);
     EPcalibV0C->SetLineWidth(2);
     //EPcalib->SetMarkerColor(kBlack);
     //EPcalib->SetMarkerStyle(20);
     EPcalibV0C->SetMinimum(0);
     EPcalibV0C->SetMaximum(0.006);
     EPcalibV0C->SetYTitle("");
     EPcalibV0C->Draw("same");
     Rleg->Draw("same");

 
  TCanvas *c5 = new TCanvas("resolution", "resolution", 600, 600);
  c5->cd();
       gStyle->SetOptStat(0);
       TPad *pad5 = new TPad("pad5", "", 0.0, 0.05, 1.0, 1.0);
       pad5->SetBottomMargin(0.15);
       pad5->SetLeftMargin(0.15);
       pad5->Draw();
       pad5->cd();
       qMepM->SetTitle("");
       qMepM->GetYaxis()->SetTitle("#it{R}_{2}");
       qMepM->GetXaxis()->SetTitle("#it{q}_{2}");
       qMepM->GetXaxis()->SetRangeUser(0.0, 10.0);
       qMepM->SetMinimum(0.01);
       qMepM->SetMaximum(0.99);
       //qAepC->Draw("same");
       //qCepC->Draw("same");
       //qAepA->Draw("same");  
       //qCepA->Draw("same");  
       qMepM->Draw("same");
       q2perms->Draw("same");
       for (int i = 0; i < 9; i++)  lines[i]->Draw("same");
       cout << "Integral: " << qMepM->Integral(qMepM->GetXaxis()->FindBin(cindq2test[0]), qMepM->GetXaxis()->FindBin(cindq2test[3]), "width") <<"\n";
       cout << "Diff: " << cindq2test[3]-cindq2test[0] <<"\n";
       cout << "R2(10-40): " << qMepM->Integral(qMepM->FindBin(cindq2test[0]), qMepM->FindBin(cindq2test[3]), "width")/(cindq2test[3]-cindq2test[0]) <<"\n";

     // refqCepC->lec->Draw("same colz");

/*  TCanvas *c4 = new TCanvas("resolution1", "resolution1", 500, 500);
  c4->cd();
       gStyle->SetOptStat(0);
       TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
       pad4->SetBottomMargin(0.15);
       pad4->SetLeftMargin(0.15);
       pad4->SetLogz();
       pad4->Draw();
       pad4->cd();
       hAClow->SetLineColor(kBlue);
         hAClow->SetMarkerColor(kBlue);
         hAClow->SetMarkerStyle(20);
         hAClow->SetMarkerSize(0.6);
       hAChigh->SetLineColor(kAzure);
         hAChigh->SetMarkerColor(kAzure);
         hAChigh->SetMarkerStyle(20);
         hAChigh->SetMarkerSize(0.6);
       hCTlow->SetLineColor(kRed);
         hCTlow->SetMarkerColor(kRed);
         hCTlow->SetMarkerStyle(20);
         hCTlow->SetMarkerSize(0.6);
       hCThigh->SetLineColor(kPink);
         hCThigh->SetMarkerColor(kPink);
         hCThigh->SetMarkerStyle(20);
         hCThigh->SetMarkerSize(0.6);
       hATlow->SetLineColor(kGreen);
         hATlow->SetMarkerColor(kGreen);
         hATlow->SetMarkerStyle(20);
         hATlow->SetMarkerSize(0.6); 
       hAThigh->SetLineColor(kYellow);
         hAThigh->SetMarkerColor(kYellow);
         hAThigh->SetMarkerStyle(20);
         hAThigh->SetMarkerSize(0.6); 
       //hAClow->Scale(1.0/hAClow->Integral());
       //hAChigh->Scale(1.0/hAChigh->Integral());
       //hCTlow->Scale(1.0/hCTlow->Integral());
       //hCThigh->Scale(1.0/hCThigh->Integral());
       //hATlow->Scale(1.0/hATlow->Integral());
       //hAThigh->Scale(1.0/hAThigh->Integral());

    hahaha->Draw("same");

    cout << "Mean: " << hahaha->GetMean() <<"\n";
   hAChigh->Draw("same");
   hAClow->Draw("same");

   hAClow->ClearUnderflowAndOverflow();

   cout << "hAClow: " << hAClow->GetMean() <<"\n";
   cout << "hAChigh: " << hAChigh->GetMean() <<"\n";

   hCThigh->Draw("same");
   hCTlow->Draw("same");
   cout << "hCTlow: " << hCTlow->GetMean() <<"\n";
   cout << "hCThigh: " << hCThigh->GetMean() <<"\n";

   hAThigh->Draw("same");
   hATlow->Draw("same");
   cout << "hATlow: " << hATlow->GetMean() <<"\n";
   cout << "hAThigh: " << hAThigh->GetMean() <<"\n";
*/
//cout << "Low Mean: " << reso_lq->GetMean() <<"\n";
//cout << "High Mean: " << reso_hq->GetMean() <<"\n";

}
