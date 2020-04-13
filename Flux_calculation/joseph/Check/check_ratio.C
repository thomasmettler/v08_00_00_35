void check_ratio(){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0011);
  gStyle->SetStatH(0.05);
  gStyle->SetStatW(0.09);
  gStyle->SetStatX(0.9);
  gStyle->SetOptStat(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetLabelSize(.045, "XY");
  gStyle->SetTitleSize(.045, "XY");
  gStyle->SetTitleFont(62, "XY");
  gStyle->SetHistLineColor(kBlack);
  gStyle->SetHistLineWidth(2);
  gStyle->SetMarkerColor(kBlack);
  gStyle->SetMarkerStyle(33);
  gStyle->SetMarkerSize(1.3);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.05);
    
  TFile* Nf = new TFile("../output_newPos.root");
  TFile* Of = new TFile("../output_oldPost.root");

  std::vector<std::string> flav = {"numu","nue"};
  
  std::vector<TH1D*> N; N.resize(flav.size());
  std::vector<TH1D*> O; O.resize(flav.size());

  std::vector<TCanvas*> c; c.resize(flav.size());

  std::vector<TF1*> f; f.resize(flav.size());

  double fDetLength = 1036.8;
  double fDetHalfWidth = 1.28175;
  double fDetHalfHeight = 1.165;

  for(int i = 0; i < flav.size(); i++){

    N[i] = (TH1D*)Nf->Get(Form("%s/%s_CV_Window",flav[i].c_str(),flav[i].c_str()));
    O[i] = (TH1D*)Of->Get(Form("%s/%s_CV_Window",flav[i].c_str(),flav[i].c_str()));
    
    //N[i]->Scale(1./(fDetHalfWidth*fDetHalfHeight*2.*2.));
    N[i]->Scale(1./197.);
    O[i]->Scale(1./200.);

    N[i]->Divide(O[i]);
    
    c[i] = new TCanvas();
    c[i]->cd();


    f[i] = new TF1(Form("fit_%s",flav[i].c_str()),"pol0", 0.250, 3);
    
    N[i]->Fit(Form("fit_%s",flav[i].c_str()),"R");

    N[i]->GetYaxis()->SetTitle(Form("Ratio New / Old for %s",flav[i].c_str()));    
    N[i]->GetXaxis()->SetTitle("True Neutrino Energy [GeV]");
    N[i]->GetYaxis()->SetRangeUser(0.7,1.3);
    N[i]->Draw();

  }
  

}
