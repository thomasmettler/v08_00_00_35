void plotSpread(){

  TFile* dkF = new TFile("../output_one_dk2nu.root");
  
  std::vector<string> syst = {"PPFXMIPPKaon","PPFXMIPPPion","PPFXOther","PPFXTargAtten","PPFXThinKaon","PPFXThinMeson","PPFXThinNeutron","PPFXThinNucA","PPFXThinNuc","PPFXThinPion","PPFXTotAbsorp","Total"};
  std::vector<string> flav = {"numu","nue"};

  std::vector< TH1D* > CV;
  std::vector < std::vector< std::vector<TH1D*> > > Syst;

  CV.resize(flav.size());
  Syst.resize(flav.size());

  for(int i = 0; i < flav.size(); i++){

    CV[i] = (TH1D*)dkF->Get(Form("%s/%s_CV_AV_TPC",
				 flav[i].c_str(),
				 flav[i].c_str()));
    CV[i]->SetLineWidth(3);
    CV[i]->SetLineColor(kBlue);
    
    Syst[i].resize(syst.size());
    

    for(int j = 0; j < syst.size(); j++){
      Syst[i][j].resize(1000);
      for(int k = 0; k < 1000; k++){
	Syst[i][j][k] = 
	  (TH1D*)dkF->Get(Form("%s/%s/Active_TPC_Volume/%s_%s_Uni_%d_AV_TPC",
			       flav[i].c_str(),
			       syst[j].c_str(),
			       flav[i].c_str(),
			       syst[j].c_str(),k));
	Syst[i][j][k]->SetLineColor(kGray);
      }//uni      
    }//syst
  }//flav

  std::vector<TCanvas*> c;
  c.resize(flav.size());
for(int i = 0; i < flav.size(); i++){
  c[i] = new TCanvas();
  c[i]->cd();  
  CV[i]->Draw("");
    for(int j = 0; j < syst.size(); j++){
      for(int k = 0; k < 1000; k++){
	Syst[i][j][k]->Draw("same");
      }
    }
    CV[i]->Draw("same");
 }
}
