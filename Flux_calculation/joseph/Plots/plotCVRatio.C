void plotCVRatio(){

  TFile* gsF = new TFile("../output_gsimple.root");
  TFile* dkF = new TFile("../output_one_dk2nu.root");

  std::vector<string> flav = {"numu","nue"};


  std::vector<TH1D*> gs; gs.resize(flav.size());
  std::vector<TH1D*> dk; dk.resize(flav.size());
  std::vector<TCanvas*> c; c.resize(flav.size());
  for(int i = 0; i < flav.size(); i++){
    c[i] = new TCanvas();
    c[i]->cd();
    gs[i] = (TH1D*)gsF->Get(Form("%s/%s_CV_AV_TPC",flav[i].c_str(),flav[i].c_str()));
    dk[i] = (TH1D*)dkF->Get(Form("%s/%s_CV_AV_TPC",flav[i].c_str(),flav[i].c_str()));

    gs[i]->Scale(3e20/(10000068.*779));
    dk[i]->Scale(3e20/500000.);

    dk[i]->Divide(gs[i]);
    dk[i]->Draw();
  }
  
}
