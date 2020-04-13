#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>

// root stuff
#include "TInterpreter.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"

// gallery stuff
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"

// data-products
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "larcoreobj/SummaryData/POTSummary.h"
//#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "MCEventWeight.h"

#include "joseph/GeoVector.h"
#include "joseph/GeoAABox.h"
#include "joseph/GeoHalfLine.h"
#include "joseph/GeoAlgo.h"

// associations 
#include "canvas/Persistency/Common/FindMany.h"

using namespace art;
using namespace std;

int main(int argv, char** argc) {
  vector<string> filename;
  for (int i = 1; i < argv; ++i) { 
    filename.push_back(argc[i]);  
  }
  InputTag mctruths_tag { "flux" };
  InputTag  evtwght_tag { "eventweight" };
  
  //active volume
  double fDetLength = 1036.8; 
  double fDetHalfWidth = 128.175;
  double fDetHalfHeight = 116.5;

  std::pair<float, float>  _xRange;
  std::pair<float, float>  _yRange;
  std::pair<float, float>  _zRange;

  _xRange.first  = 0;
  _xRange.second = 2*fDetHalfWidth;
  _yRange.first  = -1*fDetHalfHeight;
  _yRange.second = fDetHalfHeight;
  _zRange.first =  0;
  _zRange.second = fDetLength;
  
  // fiducial volume
  std::pair<float, float>  _xRange_fid;
  std::pair<float, float>  _yRange_fid;
  std::pair<float, float>  _zRange_fid;

  _xRange_fid.first  = 8.45;
  _xRange_fid.second = 244.8;
  _yRange_fid.first  = -105.53;
  _yRange_fid.second = 107.47;
  _zRange_fid.first =  9.9;
  _zRange_fid.second = 986.9;
  
  
  geoalgo::GeoAlgo const _geo_algo_instance;
  geoalgo::AABox volAVTPC( _xRange.first, _yRange.first, _zRange.first, _xRange.second, _yRange.second, _zRange.second);
  geoalgo::AABox volFVTPC( _xRange_fid.first, _yRange_fid.first, _zRange_fid.first, _xRange_fid.second, _yRange_fid.second, _zRange_fid.second);
  
  
  // TTree stuff here //////////////////////////////
  std::cout << "Initialize variables and histograms for event tree" << std::endl;
  //TFile* output = new TFile("output.root","RECREATE");
  //TTree *my_event_ = new TTree("flux","neutrino flux");
  
  //double nu_energy_;
  //int nu_pdg_;
  //int intercept_;
  //std::vector<double> para[100];
  
  //my_event_->Branch("Energy", &nu_energy_, "Neutrino energy/D");
  //my_event_->Branch("PDG", &nu_pdg_, "Neutrino PDG/I");
  //my_event_->Branch("Intercept", &intercept_, "Hit TPC?/I");
  //my_event_->Branch("Energy", &nu_energy_, "Neutrino energy [GeV]/D");
  
  //h_numu_CV       = new TH1D("h_numu_CV",     "h_numu_CV", 200, 0, 10);
  //h_anti_numu_CV  = new TH1D("h_anti_numu_CV","h_anti_numu_CV", 200, 0, 10);
  //h_nue_CV        = new TH1D("h_nue_CV",      "h_nue_CV", 200, 0, 10);
  //h_anit_nue_CV   = new TH1D("h_anit_nue_CV", "h_anit_nue_CV", 200, 0, 10);
  
  std::vector< TH1D* > Enu_CV_Window;
  std::vector< TH1D* > Enu_CV_AV_TPC;
  std::vector< TH1D* > Enu_CV_FV_TPC;
  std::vector< TH1D* > Enu_UW_Window;
  std::vector< TH1D* > Enu_UW_AV_TPC;
  std::vector< TH1D* > Enu_UW_FV_TPC;

  //flavors - systematic - universe 
  std::vector< std::vector < std::vector < TH1D* > > > Enu_Syst_Window;
  std::vector< std::vector < std::vector < TH1D* > > > Enu_Syst_AV_TPC;
  std::vector< std::vector < std::vector < TH1D* > > > Enu_Syst_FV_TPC;
  
  std::cout << "Initialize TH1D vectors" << std::endl;

  std::vector< double > TotWeight;
  TotWeight.resize(1000);
  
  std::cout << "Resized TotWeights" << std::endl;
  
  std::vector< string > flav;
  flav = {"numu","nue","numubar","nuebar"};
  
  std::vector< string > labels;
  labels = { "expskin_FluxUnisim", "horncurrent_FluxUnisim" "kminus_PrimaryHadronNormalization", "kplus_PrimaryHadronFeynmanScaling",
    "kzero_PrimaryHadronSanfordWang", "nucleoninexsec_FluxUnisim", "nucleonqexsec_FluxUnisim", "nucleontotxsec_FluxUnisim",
    "piminus_PrimaryHadronSWCentralSplineVariation", "pioninexsec_FluxUnisim", "pionqexsec_FluxUnisim", "piontotxsec_FluxUnisim",
    "piplus_PrimaryHadronSWCentralSplineVariation" };
  //labels = {"Total"};
  
  std::cout << "Defined flavours and labels" << std::endl;
  
  //systematic - universe 
  std::vector< std::vector< double > > Weights; 
  
  Weights.resize(labels.size());
  for(int i = 0; i < labels.size(); i++){
    Weights[i].resize(1000);
  }
  
  Enu_CV_Window.resize(4);
  Enu_CV_AV_TPC.resize(4);
  Enu_CV_FV_TPC.resize(4);
  Enu_UW_Window.resize(4);
  Enu_UW_AV_TPC.resize(4);
  Enu_UW_FV_TPC.resize(4);
  Enu_Syst_Window.resize(4);
  Enu_Syst_AV_TPC.resize(4);
  Enu_Syst_FV_TPC.resize(4);
  
  std::cout << "Start inizializing TH1D vectors histo" << std::endl;
  
  for(int i = 0; i < int(flav.size()); i++){//flav
    std::cout << "Loop over flavour: " << flav[i].c_str() << std::endl;
    Enu_CV_Window[i] = new TH1D(Form("%s_CV_Window",flav[i].c_str()),"", 200, 0, 10);
    Enu_CV_AV_TPC[i] = new TH1D(Form("%s_CV_AV_TPC",flav[i].c_str()),"", 200, 0, 10);
    Enu_CV_FV_TPC[i] = new TH1D(Form("%s_CV_FV_TPC",flav[i].c_str()),"", 200, 0, 10);

    Enu_UW_Window[i] = new TH1D(Form("%s_unweighted_Window",flav[i].c_str()),"", 200, 0, 10);
    Enu_UW_AV_TPC[i] = new TH1D(Form("%s_unweighted_AV_TPC",flav[i].c_str()),"", 200, 0, 10);
    Enu_UW_FV_TPC[i] = new TH1D(Form("%s_unweighted_FV_TPC",flav[i].c_str()),"", 200, 0, 10);

    Enu_Syst_Window[i].resize(labels.size());
    Enu_Syst_AV_TPC[i].resize(labels.size());
    Enu_Syst_FV_TPC[i].resize(labels.size());

    for(int j = 0; j < int(labels.size()); j++){//labels
      std::cout << "Loop over Label: " << labels[j].c_str() << std::endl;
      Enu_Syst_Window[i][j].resize(1000);
      Enu_Syst_AV_TPC[i][j].resize(1000);
      Enu_Syst_FV_TPC[i][j].resize(1000);
           
      for(int k = 0; k < 1000; k++){//unis
        Enu_Syst_Window[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k),"", 200, 0, 10);
        Enu_Syst_AV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_AV_TPC",flav[i].c_str(), labels[j].c_str(), k),"", 200, 0, 10);
        Enu_Syst_FV_TPC[i][j][k] =  new TH1D(Form("%s_%s_Uni_%d_FV_TPC",flav[i].c_str(), labels[j].c_str(), k),"", 200, 0, 10);
      }//iterate over universes
    }//iterate over systematics
  }//iterate over flavors 
  std::cout << "End inizializing TH1D vectors histo" << std::endl;
  
  
  ///End TTree stuff  ///////////////////////////////

  //Let's Do Science! 
  int n = 0;
  
  std::cout << "Start loop over events" << std::endl;
  
  for (gallery::Event ev(filename) ; !ev.atEnd(); ev.next()) {
    n++;
    if(n%10000 == 0){
      std::cout << "Event number: " << n << std::endl;
    }
    /*
    if( n == 1){
      auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);
      evwgh::MCEventWeight evtwght;
      evtwght = evtwghts.at(0);
      std::map<std::string, std::vector< double > > Weights;
      for(auto last : evtwght.fWeight){
        std::vector<double> vec(int(last.second.size()));
        for(int i = 0; i < int(last.second.size()); i++){	      
          vec[i] = last.second.at(i);
        }
        Weights[last.first]=vec;
      }
      int counter_para = 0;
      for (auto wgh : Weights) {
        std::cout << "branching parameter: " << counter_para << " name: " << wgh.first << " length: " << wgh.second.size() << std::endl;
        //my_event_->Branch((wgh.first).c_str(),             &para[counter_para]);
        counter_para++;
	      //cout << " Weight "<<wgh.first<<"\t"<<wgh.second.at(0)<<", "<<wgh.second.at(1)<<endl;
      }
    }
    */
    //Next we want to grab from the event the data-product that you want
    auto const& mctruths = *ev.getValidHandle<vector<simb::MCTruth>>(mctruths_tag);   
    auto const& mcfluxs = *ev.getValidHandle<vector<simb::MCFlux>>(mctruths_tag);   

    //std::vector<evwgh::MCEventWeight> evtwghts;
    auto const& evtwghts = *ev.getValidHandle<vector<evwgh::MCEventWeight>>(evtwght_tag);  
    
    if (mctruths.empty() || evtwghts.empty()) {
      continue;    
    }

    //Now we'll iterate through these 
    //std::cout << "mctruths size: " << mctruths.size() << std::endl;
    for (size_t i = 0; i < mctruths.size(); i++) {
      auto const& mctruth = mctruths.at(i);
      auto const& mcflux = mcfluxs.at(i);
      //cout << mctruth << endl;
      //cout << mcflux.Flux(14) << endl;
      evwgh::MCEventWeight evtwght;
      evtwght = evtwghts.at(i);
      
      int pdg = 9; 
      if(mctruth.GetNeutrino().Nu().PdgCode() == 14) {pdg = 0;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() == 12) {pdg = 1;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-14) {pdg = 2;}
      else if(mctruth.GetNeutrino().Nu().PdgCode() ==-12) {pdg = 3;}
      else {std::cout << "shits cray " << mctruth.GetNeutrino().Nu().PdgCode() << std::endl; continue;}

      //nu_pdg_ = mctruth.GetNeutrino().Nu().PdgCode();
      //nu_energy_ = mctruth.GetNeutrino().Nu().E();
      
      /*std::cout << "Neutrino vetex: " << mctruth.GetNeutrino().Nu().Vx() << " - " <<
			    mctruth.GetNeutrino().Nu().Vy() << " - " <<
			    mctruth.GetNeutrino().Nu().Vz() <<  " - " <<
			    mctruth.GetNeutrino().Nu().Px() << " - " <<
			    mctruth.GetNeutrino().Nu().Py() << " - " <<
			    mctruth.GetNeutrino().Nu().Pz() << std::endl;
      */
      geoalgo::HalfLine ray(mctruth.GetNeutrino().Nu().Vx(),
			    mctruth.GetNeutrino().Nu().Vy(),
			    mctruth.GetNeutrino().Nu().Vz(),
			    mctruth.GetNeutrino().Nu().Px(),
			    mctruth.GetNeutrino().Nu().Py(),
			    mctruth.GetNeutrino().Nu().Pz());
      
      auto vec = _geo_algo_instance.Intersection(volAVTPC, ray);
      auto vec_fid = _geo_algo_instance.Intersection(volFVTPC, ray);
     
      //std::cout << "# of intersections : " << vec.size()  << std::endl;

      int intercept = 0;
      int intercept_fid = 0;

      if(vec.size() == 0){ intercept = 0; }
      if(vec.size() == 2){ intercept = 1; }
      if(vec.size() != 2 && vec.size() != 0){ std::cout << "you dum" << std::endl;}
      
      if(vec_fid.size() == 0){ intercept_fid = 0; }
      if(vec_fid.size() == 2){ intercept_fid = 1; }
      if(vec_fid.size() != 2 && vec_fid.size() != 0){ std::cout << "you dum" << std::endl;}
      
      //std::cout << "Neutrino intercepts?  " << intercept << std::endl;
      
      double cv_weight = 1; 
      //intercept_ = intercept;
      
      for(int l = 0; l < int(labels.size()); l++){
        std::fill(Weights[l].begin(), Weights[l].end(), 1.0);
      }
      
      for(auto last : evtwght.fWeight){
        for(int l = 0; l < int(labels.size()); l++){
          if(last.first.find(labels[l].c_str()) != std::string::npos){
            for(int i = 0; i < int(last.second.size()); i++){	      
              Weights[l][i] *= last.second.at(i);
            }
          }
        }
      }
      
      Enu_CV_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
      Enu_UW_Window[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
      
      if(intercept){
        Enu_CV_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
        Enu_UW_AV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
      }
      if(intercept_fid){
        Enu_CV_FV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
        Enu_UW_FV_TPC[pdg]->Fill(mctruth.GetNeutrino().Nu().E());
      }
      
      std::fill(TotWeight.begin(), TotWeight.end(), 1.0);
      if(pdg==0){
        for(int l = 0; l < int(labels.size()); l++){
          for(unsigned int i = 0; i < Weights[l].size(); i++){
            TotWeight[i] *= Weights[l][i];
            //Enu_Syst_Window[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(),TotWeight[i]);

            if(intercept){	  
              Enu_Syst_AV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(),TotWeight[i]);
            }
            if(intercept_fid){	  
              Enu_Syst_FV_TPC[pdg][l][i]->Fill(mctruth.GetNeutrino().Nu().E(),TotWeight[i]);
            }
          }
        }
      }/*
      int full_lable = labels.size()-1;
      for(unsigned int i = 0; i < Weights[0].size(); i++){
        Enu_Syst_Window[pdg][full_lable][i]->Fill(mctruth.GetNeutrino().Nu().E(),TotWeight[i]);

        if(intercept){	  
          Enu_Syst_AV_TPC[pdg][full_lable][i]->Fill(mctruth.GetNeutrino().Nu().E(),TotWeight[i]);
        }	
      }*/
        
    }//Iterate through neutrino interactions
    //std::cout << "End of loop over neutrino thruth interactions" << std::endl;
  }// Iterate through events
  std::cout << "End of loop over events" << std::endl;
      
      //We can also write the output to a output root file like this:
  //TFile* output = new TFile("output.root","RECREATE"); //filename Form("%s_%s_Uni_%d_Window",flav[i].c_str(), labels[j].c_str(), k)
  TFile* output = new TFile(Form("%s_histo.root",filename.at(0).c_str() ),"RECREATE");
  
  TDirectory *savdir = gDirectory;

  std::vector< std::vector< std::vector< TDirectory* > > > subdir; //flav //syst //cont 
  subdir.resize(4);
  for(int i = 0; i < 4; i++){
    subdir[i].resize(labels.size()+1);
    for(int j = 0; j < labels.size()+1; j++){
      subdir[i][j].resize(3);
    }
  }

  std::vector<string> cont; cont = {"Window","Active_TPC_Volume",""};

  for(int f = 0; f < 4; f++){
    std::cout << flav[f] << std::endl;
    subdir[f][0][0] = savdir->mkdir(Form("%s",flav[f].c_str()));
    subdir[f][0][0]->cd();
    
    Enu_CV_Window[f]->Write();      
    Enu_CV_AV_TPC[f]->Write();
    Enu_CV_FV_TPC[f]->Write();
    Enu_UW_Window[f]->Write();      
    Enu_UW_AV_TPC[f]->Write();
    Enu_UW_FV_TPC[f]->Write();

    for(int s = 1; s < labels.size()+1; s++){
      std::cout << labels[s-1] << std::endl;
      subdir[f][s][0] = subdir[f][0][0]->mkdir(Form("%s",labels[s-1].c_str()));
      subdir[f][s][0]->cd();
      for(int c = 1; c < 3; c++){
        std::cout << cont[c-1] << std::endl;
        subdir[f][s][c] = subdir[f][s][0]->mkdir(Form("%s",cont[c-1].c_str()));
        subdir[f][s][c]->cd();
        
        if(c == 1 && f==0){
          for(int i = 0; i < 1000; i++){
            //Enu_Syst_Window[f][s-1][i]->Write();
          }	
        }

        if(c == 2 && f==0){
          for(int i = 0; i < 1000; i++){
            Enu_Syst_AV_TPC[f][s-1][i]->Write();
            Enu_Syst_FV_TPC[f][s-1][i]->Write();
          }	
        }
        
        }//cont
      
    }//systs
    savdir->cd();
  }//flavs

  output->Close();
  
  return 1;

    
}//main 

      
      
      
      
      
      
      
/*      
      

      std::map<std::string, std::vector< double > > Weights;   
      int num_para = 0;
      for(auto last : evtwght.fWeight){
        para[num_para] = last.second;
        num_para++;
        std::vector<double> vec(int(last.second.size()));
        for(int j = 0; j < int(last.second.size()); j++){	      
          vec[j] = last.second.at(j);
        }
        Weights[last.first]=vec;
      }
      
      cout <<"Event "<<n<<" neutrino "<<i<<endl;
      //example dumping first and second universe weight for all weights
      for (auto wgh : Weights) {
	      //cout << " Weight "<<wgh.first<<"\t"<<wgh.second.at(0)<<", "<<wgh.second.at(1)<<endl;
      }
      my_event_->Fill();
    }//end of loop over mcthruths
    //my_event_->Fill();
 }//end of loop over events
  
  output->Close();
  
  return 1;

}
*/