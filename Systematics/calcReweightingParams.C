void calcReweightingParams(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inFile = TFile::Open("UnfoldingData_2D_non_R02_3050_lo30_Feb4.root");

  //get hist
  TH1D* data = (TH1D*)inFile->Get("specraw1");
  data->SetName("data");
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  data->SetLineColor(kBlack);
  TH1D* mc =(TH1D*)inFile->Get("raw");
  mc->SetName("mc");
  mc->SetMarkerColor(kRed);
  mc->SetMarkerStyle(20);
  mc->SetLineColor(kRed);
  //data->Rebin(5);
  //mc->Rebin(5);
  Int_t bin1 = data->GetXaxis()->FindBin(30.);
  //Int_t bin1 = data->GetXaxis()->FindBin(50.);
  Int_t bin2 = data->GetXaxis()->FindBin(120.);
  //std::cout << data->Integral(bin1, bin2) << std::endl;
  //std::cout << mc->Integral(bin1, bin2) << std::endl;
  data->Scale(1./ data->Integral(bin1, bin2));
  mc->Scale(1./mc->Integral(bin1, bin2));
  data->Divide(mc);

  //fit polynomial
  TF1* linFit = new TF1("linFit", "pol2", 30,120);
  //TF1* linFit = new TF1("linFit", "pol2", 50,200);
  //linFit->SetParameters(0.5, 0.5, 0.9);
  linFit->SetRange(30, 120);
  TLine *line = new TLine(30, 1, 120, 1);
  line->SetLineStyle(2);
  //linFit->SetRange(50, 80);
  data->Fit(linFit, "N R I", "");
  data->GetYaxis()->SetRangeUser(0., 2.);
  data->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  data->SetYTitle("Data/MC");
  data->Draw();
  linFit->Draw("same");
  line->Draw("same");
}
