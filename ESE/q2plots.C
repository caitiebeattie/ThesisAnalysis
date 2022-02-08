#include "TRandom.h"
#include "TRandom3.h"
  
void q2plots ()  {

   //gStyle->SetPalette(kCool);  
   TRandom3* rand = new TRandom3();   

   double centl = 30.0;
   double centr = 50.0;

   TString cFiles2="files15.txt";

  //============================================
  //==========   Get Raw Spectra    ============
  //============================================
   //Raw Spectra from Trains
   //string rootfilename = "../Spectra/AnalysisResults7283.root";   // semi-central
    string rootfilename = "AnalysisResults7812_7857.root";             // R=0.2, 30-50%, pass1
    const char *lerootfile = rootfilename.c_str(); 

    //Load File/Tree into system
    TFile *lefile = new TFile(lerootfile);
    TTree *T = nullptr;
    lefile->GetObject("JetTree_AliAnalysisTaskJetExtractor_Jet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_Jet", T);
    float jetpT, centra, qvec, ep, N;
    T->SetBranchAddress("Jet_Pt", &jetpT);
    T->SetBranchAddress("Event_Centrality", &centra);
    T->SetBranchAddress("Event_Q2Vector", &qvec);    
    T->SetBranchAddress("Jet_EPangle", &ep);    
    T->SetBranchAddress("Jet_Eta", &N);

      //Calculate q2 percentiles
      TDirectoryFile *qnfile = (TDirectoryFile*)lefile->Get("PWGJE_QnVectorTender");
      TList *qnlist = (TList*)qnfile->Get("coutputQnVectorTender_Q2V0C_EPV0C");
      TH2F *q2hist2D = (TH2F*)qnlist->FindObject("fHistqnVsCentrV0C");
      //initialize variables used to calculate percentiles
      long int percentileticker = 0;                                   //number of entries, starting from lowest q2 and counting upward
      int percentilen = 1;                                             //number percentile we're on
      vector<long double> percentiles(0);                              //vector that stores 10th, 20th etc percentile per cent bin
      vector<vector<long double>>  allPercent(0);                      //vector that stores percentiles vectors for all centralities
      double centdiff = (centr-centl);
      int leftc, rightc;
      //get percentiles in bins of 1% centrality
      for (int j = 0; j < (double)centdiff; j++) {
        leftc = (int)centl + j; 
        rightc = (int)centl + j + 1;    
        q2hist2D->GetXaxis()->SetRangeUser(leftc, rightc);       //restrict 2D hist to current 1% centrality bin
        TH1F *q2hist = (TH1F*)q2hist2D->ProjectionY();           
        q2hist->Rebin(100);
        long double percentilemarker = (long double)q2hist->GetEntries()/10.0;      //number of entries that compose 10% of sample
        for (int i = 1; i <= 100; i++)   {
          //double split = rand->Rndm();    
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




  TH1D *perc20 = new TH1D("perc20", "perc20", centdiff, centl, centr);
  TH1D *perc50 = new TH1D("perc50", "perc50", centdiff, centl, centr);
  TH1D *perc80 = new TH1D("perc80", "perc80", centdiff, centl, centr);
  for (int i = 1; i <= centdiff; i++) {
      perc20->SetBinContent(i, allPercent[i-1][2]);
      perc50->SetBinContent(i, allPercent[i-1][4]);
      perc80->SetBinContent(i, allPercent[i-1][6]);
  }
  perc20->SetLineColor(kRed);
    perc20->SetLineWidth(2);
  perc80->SetLineColor(kRed);
    perc80->SetLineWidth(2);
  perc50->SetLineColor(kBlack);
    perc50->SetLineWidth(2);

    std::vector<double> kBinsMeasured = {30.0, 40.0, 50.0, 60.0, 70.0, 85.0, 100.0, 120.0};   //semi-central
    std::vector<double> kq2binsn = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0}; 
    std::vector<double> kBinsEPangle = {0.0, sqrt(1)/2.0, sqrt(2)/2.0, sqrt(3)/2.0, 1.0};
   
    TH2D* jetvq = new TH2D("jetvq", "jetvq", kBinsMeasured.size()-1, kBinsMeasured.data(), kq2binsn.size()-1, kq2binsn.data());
    TH1D* epreal = new TH1D("epreal", "epreal", 10, 0.0, 1.0);  

    int entries = T->GetEntries();
    for (int i = 0; i < entries; i++)   {
       T->GetEntry(i);     
       if (centra < centl || centra > centr)    continue;  
       if (jetpT < 30.0 || jetpT > 120.0)       continue;
       //if (ese == 0)  {if (qvec > allPercent[(int)centra-30][1]) continue;} //20th percentile
       //if (ese == 1)  {if (qvec < allPercent[(int)centra-30][7]) continue;} //80th percentile
       //if (abs(cos(ep)) < sqrt(2)/2.0) continue;  //in-plane
       //if (abs(cos(ep)) > sqrt(2)/2.0) continue;  //out-plane
       jetvq->Fill(jetpT, qvec);
      if (abs(N) > (0.7))                continue;  //TPC fiducial cut
       epreal->Fill(abs(cos(ep)));

    }
   float ptJet,ptPar, ptdet, centrality, matchradius, area, eta, q2, epangle, epPar;
  
   int nEv=0;; 
 
   ifstream infile2;
   infile2.open(cFiles2.Data());
   char filename2[300];


   TH1F* theepangle = new TH1F("theepangle", "the epangle", 10, 0.0, 1.0);
   TH2F* jeteptru = new TH2F("jeteptru", "jeteptru", kBinsMeasured.size()-1, kBinsMeasured.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
   TH2F* jetephyb = new TH2F("jetephyb", "jetephyb", kBinsMeasured.size()-1, kBinsMeasured.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
   TH2F* ptresp = new TH2F("ptresp", "ptresp", kBinsMeasured.size()-1, kBinsMeasured.data(), kBinsMeasured.size()-1, kBinsMeasured.data());
   TH2F* epresp = new TH2F("epresp", "epresp", kBinsEPangle.size()-1, kBinsEPangle.data(), kBinsEPangle.size()-1, kBinsEPangle.data());
 
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
      double xsect = hcross->GetBinContent(htrials->GetMean() + 1);
      double xscale = xsect*hcross->GetEntries();
      double trials = htrials->GetBinContent(htrials->GetMean() + 1);
      double pTHardscalefactor=xscale/(n_event*trials);  

      std::cout << "pT Hard Bin: " << pthardbin-1 << ", Scale Factor:  " << pTHardscalefactor << std::endl;
    TFile *input2=TFile::Open(filename2);
    // get the mc tree
    TTree *mc = nullptr;
    mc=(TTree*)input2->Get("JetTree_AliAnalysisTaskJetExtractor_hybJet_AKTChargedR020_tracks_pT0150_pt_scheme_Rho_hybJet"); 
    int nEv=mc->GetEntries(); 

    int numTracks;
    // set the branches of the mc tree
    mc->SetBranchAddress("Jet_Pt", &ptJet); 
    mc->SetBranchAddress("Jet_Eta", &eta);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_Pt", &ptPar);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Pt", &ptdet);
    float jetTrackPt[400];
    //mc->SetBranchAddress("Jet_Track_Pt", &jetTrackPt);
    mc->SetBranchAddress("Event_Centrality", &centrality);
    mc->SetBranchAddress("Jet_Area", &area);
    mc->SetBranchAddress("Jet_MC_MatchedDetLevelJet_Distance", &matchradius); 
    mc->SetBranchAddress("Event_Q2Vector", &q2);
    mc->SetBranchAddress("Jet_EPangle", &epangle);
    mc->SetBranchAddress("Jet_MC_MatchedPartLevelJet_EPangle", &epPar);
    
  
  int countm=0;
    for(int iEntry=0; iEntry< nEv; iEntry++){
      mc->GetEntry(iEntry);
      //****** Get the scaling factor
      double scalefactor = pTHardscalefactor;
      double EBscale = 1.;
      if(ptJet >= -20. && ptJet < 10.)        EBscale = 10.;
      else if(ptJet >= 10. && ptJet < 20.)    EBscale = 1.0/0.01;
      else if(ptJet >= 20. && ptJet < 30.)    EBscale = 1.0/0.03;
      else if(ptJet >= 30. && ptJet < 40.)    EBscale = 1.0/0.15;
      else if(ptJet >= 40. && ptJet < 60.)    EBscale = 1.0/0.20;
      else if(ptJet >= 60. && ptJet < 80.)    EBscale = 1.0/0.20;
      else if(ptJet >= 80. && ptJet < 100.)   EBscale = 1.0/0.20;
      else if(ptJet >= 100. && ptJet < 140.)  EBscale = 1.0/0.15;
      else if(ptJet >= 140. && ptJet < 200.)  EBscale = 1.0/0.5;
      scalefactor*=EBscale;
      // *********

      if (centrality < centl || centrality > centr)  continue;
      //if (ese == 0)  {if (q2 > allPercent[(int)centrality-30][1]) continue;} //20th percentile
      //if (ese == 1)  {if (q2 < allPercent[(int)centrality-30][7]) continue;} //80th percentile

      if (abs(eta) > (0.7))                continue;  //TPC fiducial cut
      if (ptPar > 200.)                    continue;    
      if (ptPar < 10.0 || ptJet < 10.0 )   continue; 
       if (ptJet < 30.0 || ptJet > 120.0)       continue;

	  theepangle->Fill(abs(cos(epangle)), scalefactor);
          jeteptru->Fill(ptPar, abs(cos(epPar)), scalefactor);
          jetephyb->Fill(ptJet, abs(cos(epangle)), scalefactor);
          ptresp->Fill(ptPar, ptJet, scalefactor);
          epresp->Fill(epPar, epangle, scalefactor);
    }}


   theepangle->SetLineColor(kRed);
   theepangle->SetMarkerColor(kRed);
   theepangle->SetMarkerStyle(20);
   theepangle->SetMarkerSize(0.6);
   theepangle->Scale(1.0/theepangle->Integral());


   epreal->SetLineColor(kBlack);
   epreal->SetMarkerColor(kBlack);
   epreal->SetMarkerStyle(20);
   epreal->SetMarkerSize(0.7);
   epreal->Scale(1.0/epreal->Integral());


   TH1F* lerat = (TH1F*)epreal->Clone("lerat");
   lerat->Divide(theepangle);
   lerat->SetLineColor(kBlack);
   lerat->SetMarkerColor(kBlack);
   lerat->SetMarkerStyle(20);
   lerat->SetMarkerSize(0.7);

   auto specleg = new TLegend(0.15, 0.6, 0.4, 0.85);
        specleg->SetTextSize(0.055);
        specleg->SetBorderSize(0);
        specleg->AddEntry(epreal, "Data", "pl");
        specleg->AddEntry(theepangle, "MC", "pl");  

  TLine *line = new TLine(0.0, 1, 1.0, 1);
  line->SetLineStyle(2);


  //--------------------------------------------------------------------------


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
       q2hist2D->GetXaxis()->SetRangeUser(centl, centr);
       q2hist2D->GetXaxis()->SetTitle("Centrality (%)");
       q2hist2D->GetYaxis()->SetTitle("q_{2}^{V0C}");
       q2hist2D->Draw("same colz");
       perc20->Draw("same");
       perc80->Draw("same");

  TCanvas *c2 = new TCanvas("Q2 vs Jet Spec", "Q2 vs Jet Spec", 500, 500);
  c2->cd();
       TPad *pad2 = new TPad("pad2", "", 0.0, 0.05, 1.0, 1.0);
          pad2->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad2->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad2->SetRightMargin(0.15);
          pad2->Draw();
          pad2->cd();
          pad2->SetLogz();
          jetvq->Draw("same colz");


  TCanvas *c3 = new TCanvas("epan", "epan", 500, 500);
  c3->cd();
       TPad *pad3a = new TPad("pad3a", "", 0.0, 0.4, 1.0, 1.0);
          pad3a->SetBottomMargin(0.0); // Upper and lower plot are joined
          pad3a->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad3a->Draw();
       TPad *pad3b = new TPad("pad3b", "", 0.0, 0.05, 1.0, 0.4);
          pad3b->SetBottomMargin(0.2); // Upper and lower plot are joined
          pad3b->SetTopMargin(0.0); // Upper and lower plot are joined
          pad3b->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad3b->Draw();
          pad3a->cd();
     epreal->SetMinimum(0.01);
     epreal->SetYTitle("Probability");
     epreal->Draw("same");      
     theepangle->Draw("same");
     specleg->Draw("same");
    c3->cd();
      pad3b->cd();
      lerat->SetTitle("");
      lerat->SetMinimum(0.01);
      lerat->SetMaximum(1.99);
      lerat->SetXTitle("|cos(#varphi_{jet} - #varphi_{event-plane})|");
      lerat->SetYTitle("Ratio");
      lerat->GetXaxis()->SetLabelSize(0.08);
      lerat->GetYaxis()->SetLabelSize(0.08);
      lerat->GetXaxis()->SetTitleSize(0.08);
      lerat->GetYaxis()->SetTitleSize(0.08);
      lerat->Draw("same");
      line->Draw("same");

  TCanvas *c4 = new TCanvas("truspec", "truspec", 500, 500);
  c4->cd();
       TPad *pad4 = new TPad("pad4", "", 0.0, 0.05, 1.0, 1.0);
          pad4->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad4->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad4->SetRightMargin(0.15); // Upper and lower plot are joined
          pad4->Draw();
          pad4->cd();
          pad4->SetLogz();
       jeteptru->SetXTitle("#it{p}_{T}^{true}");
       jeteptru->SetYTitle("|cos(#Delta#varphi)^{true}|");   
       jeteptru->Draw("colz same");

  TCanvas *c5 = new TCanvas("hybspec", "hybspec", 500, 500);
  c5->cd();
       TPad *pad5 = new TPad("pad5", "", 0.0, 0.05, 1.0, 1.0);
          pad5->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad5->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad5->SetRightMargin(0.15); // Upper and lower plot are joined
          pad5->Draw();
          pad5->cd();
          pad5->SetLogz();
       jetephyb->SetXTitle("#it{p}_{T}^{hybrid}");
       jetephyb->SetYTitle("|cos(#Delta#varphi)^{hybrid}|");   
       jetephyb->Draw("colz same");
     
  TCanvas *c6 = new TCanvas("pt", "pt", 500, 500);
  c6->cd();
       TPad *pad6 = new TPad("pad6", "", 0.0, 0.05, 1.0, 1.0);
          pad6->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad6->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad6->SetRightMargin(0.15); // Upper and lower plot are joined
          pad6->Draw();
          pad6->cd();
          pad6->SetLogz();
       ptresp->SetXTitle("#it{p}_{T}^{true} (GeV/#it{c})");
       ptresp->SetYTitle("#it{p}_{T}^{hybrid} (GeV/#it{c})");   
       ptresp->SetTitle("");
       ptresp->Draw("colz same");
  

  TCanvas *c7 = new TCanvas("epr", "epr", 500, 500);
  c7->cd();
       TPad *pad7 = new TPad("pad7", "", 0.0, 0.05, 1.0, 1.0);
          pad7->SetBottomMargin(0.15); // Upper and lower plot are joined
          pad7->SetLeftMargin(0.15); // Upper and lower plot are joined
          pad7->SetRightMargin(0.15); // Upper and lower plot are joined
          pad7->Draw();
          pad7->cd();
          pad7->SetLogz();
       epresp->SetXTitle("|cos(#Delta#varphi)^{true}");
       epresp->SetYTitle("|cos(#Delta#varphi)^{hybrid}|");   
       epresp->SetTitle(""); 
       epresp->Draw("colz same");
}
