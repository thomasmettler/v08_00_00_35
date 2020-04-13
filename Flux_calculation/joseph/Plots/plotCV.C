void plotCV(){

  TFile* gsF = new TFile("../output_gsimple.root");
  TFile* dkF = new TFile("../output.root");

  TH1D* gsnumu = (TH1D*)gsF->Get("nue/nue_CV_AV_TPC");
  TH1D* dknumu = (TH1D*)dkF->Get("nue/nue_unweighted_AV_TPC");//"numu/numu_CV_AV_TPC");
  //  TH1D* gsWinnumu = (TH1D*)gsF->Get("numu/numu_CV_Window");
  //TH1D* dkWinnumu = (TH1D*)dkF->Get("numu/numu_CV_Window");

  dknumu->SetLineColor(kRed);
  gsnumu->SetLineColor(kBlue);

  gsnumu->Scale(1./(1.00001e+07*779.));
  dknumu->Scale((1.)/(500000.));
                   
  gsnumu->Draw("hist");
  dknumu->Draw("hist same");

  //  gsnumu->Divide(dknumu);
  // gsnumu->Draw("hist");



  //gsnumu->DrawNormalized("hist",1);
  // dknumu->DrawNormalized("hist same",1);
  //  gsWinnumu->DrawNormalized("hist same",1);
  //dkWinnumu->DrawNormalized("hist same",1);

}
