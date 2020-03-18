///////////////////////////////////////////////////////////////////////////////////////////////////////
//    generate one entry per event
//    Written by Thomas Mettler
///////////////////////////////////////////////////////////////////////////////////////////////////////

#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <time.h>
#include <unistd.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>


void groupby_point_V33(const char* filename, int event_size, int is_data_ = 1, int with_genie = 0){
  FILE *data = 0;
  FILE *text = 0; 
  int has_wire = 0;
  int has_space_point = 0;
  int verbose_ = 0;
  //int is_data_ = 1;
  
  TFile *f_in = TFile::Open(filename,"READ");
  if (!f_in) { return; }
  
  //TTree *t_in; 
  //f_in->GetObject("NuCCana/event",t_in);
  TTree * t_in=(TTree*)f_in->Get("numuCCAna/event");
  
  //std::vector<int> *track_key_ = 0;
  //int muon_candidate_key_ = 0;
  //TSetBranchAddress *bvpx = 0;
  //t_in->SetSetBranchAddressAddress("track_key",&track_key_);
  //t_in->SetSetBranchAddressAddress("muon_candidate_key",&muon_candidate_key_);
  
  
  int event_counter = 0;
  uint32_t fEvtNum;                //Number of current event                       
  uint32_t frunNum;                //Run Number taken from event  
  uint32_t fsubRunNum;             //Subrun Number taken from event 
  
  double EventWeight = 1.0;
  double TunedCentralValue_Genie_out = 1.0;
  std::vector<double> *  TunedCentralValue_Genie_ = 0;
  std::vector<double> *  All_Genie = 0;
  
  std::vector<double> *  expskin_FluxUnisim_ = 0;
  std::vector<double> *  horncurrent_FluxUnisim_ = 0;
  std::vector<double> *  kminus_PrimaryHadronNormalization_ = 0;
  std::vector<double> *  kplus_PrimaryHadronFeynmanScaling_ = 0;
  std::vector<double> *  kzero_PrimaryHadronSanfordWang_ = 0;
  std::vector<double> *  nucleoninexsec_FluxUnisim_ = 0;
  std::vector<double> *  nucleonqexsec_FluxUnisim_ = 0;
  std::vector<double> *  nucleontotxsec_FluxUnisim_ = 0;
  std::vector<double> *  piminus_PrimaryHadronSWCentralSplineVariation_ = 0;
  std::vector<double> *  pioninexsec_FluxUnisim_ = 0;
  std::vector<double> *  pionqexsec_FluxUnisim_ = 0;
  std::vector<double> *  piontotxsec_FluxUnisim_ = 0;
  std::vector<double> *  piplus_PrimaryHadronSWCentralSplineVariation_ = 0;
  
  std::vector<double> *  AxFFCCQEshape_Genie_ = 0;
  std::vector<double> *  DecayAngMEC_Genie_ = 0;
  std::vector<double> *  MaNCRES_Genie_ = 0;
  std::vector<double> *  Theta_Delta2Npi_Genie_ = 0;
  std::vector<double> *  VecFFCCQEshape_Genie_ = 0;
  
  int has_neutrino_ = -1;
  double NuScore_ = -1;
  double FlashScore_ = -1;
  double FlashScoreTime_ = -1;
  int NuPDG_ = 0;
  int NumPfp_ = -1;
  double Nu_Vx_=-999,Nu_Vy_=-999,Nu_Vz_=-999;
  
  int NuTracks_ = 0;
  int NuShowers_ = 0;
  
  // vector for all tracks
  
  //std::vector<int> * track_key_;
  std::vector<int> *track_key_ = 0;
  std::vector<int> * TrackPID_ = 0;

  std::vector<double> * a_crthit_ts0 = 0;
  std::vector<double> * a_crthit_ts1 = 0;
  std::vector<int> * a_adc_length = 0;
  std::vector<double> * a_crt_adc = 0;
  std::vector<int> * a_t0_counter = 0;

  std::vector<double> * crtt0_time_ = 0;
  std::vector<int> * crtt0_trig_ = 0;
  std::vector<double> * crtt0_DCA_ = 0;
  std::vector<int> * crtt0_plane_ = 0;

  std::vector<double> * VtxDistance_ = 0;
  std::vector<double> * Vx_ = 0;
  std::vector<double> * Vy_ = 0;
  std::vector<double> * Vz_ = 0;
  std::vector<double> * TrackScore_ = 0;
  std::vector<double> * TrackScoreGlobal_ = 0;
  std::vector<double> * TrackLength_ = 0;
  std::vector<double> * TrackPID_chiproton_ = 0;
  std::vector<double> * TrackPID_chimuon_ = 0;
  std::vector<double> * TrackPID_chipion_ = 0;
  std::vector<double> * TrackPID_chikaon_ = 0;

  std::vector<float> * TrackMomRange_p_ = 0;
  std::vector<float> * TrackMomRange_mu_ = 0;
  std::vector<float> * TrackMomMCS_mom_ = 0;
  std::vector<float> * TrackMomMCS_err_ = 0;
  std::vector<float> * TrackMomMCS_ll_ = 0;

  std::vector<double> * TrackStart_x_ = 0;
  std::vector<double> * TrackStart_y_ = 0;
  std::vector<double> * TrackStart_z_ = 0;
  std::vector<double> * TrackEnd_x_ = 0;
  std::vector<double> * TrackEnd_y_ = 0;
  std::vector<double> * TrackEnd_z_ = 0;
  std::vector<double> * TrackDir_x_ = 0;
  std::vector<double> * TrackDir_y_ = 0;
  std::vector<double> * TrackDir_z_ = 0;
  
  std::vector<double> * AllTrack_point_x_ = 0;
  std::vector<double> * AllTrack_point_y_ = 0;
  std::vector<double> * AllTrack_point_z_ = 0;
  
  std::vector<unsigned int> * Wire_id_ = 0;
  std::vector<unsigned int> * Wire_plane_ = 0;
  
  std::vector<double> * SpacePoint_x_ = 0;
  std::vector<double> * SpacePoint_y_ = 0;
  std::vector<double> * SpacePoint_z_ = 0;

  std::vector<double> * TrackTheta_ = 0;
  std::vector<double> * TrackPhi_ = 0;
  
   std::vector<double> * Vx_sce_ = 0;
   std::vector<double> * Vy_sce_ = 0;
   std::vector<double> * Vz_sce_ = 0;
   std::vector<double> * VtxDistance_sce_ = 0;
   std::vector<int> * TrackPfp_ = 0;
   std::vector<int> * isShowerTrack_ = 0;
   std::vector<double> * TrackStart_x_sce_ = 0;
   std::vector<double> * TrackStart_y_sce_ = 0;
   std::vector<double> * TrackStart_z_sce_ = 0;
   std::vector<double> * TrackEnd_x_sce_ = 0;
   std::vector<double> * TrackEnd_y_sce_ = 0;
   std::vector<double> * TrackEnd_z_sce_ = 0;
   std::vector<int> * TrackNHitsU_ = 0;
   std::vector<int> * TrackNHitsV_ = 0;
   std::vector<int> * TrackNHitsY_ = 0;
   std::vector<float> * TrackCaloU_ = 0;
   std::vector<float> * TrackCaloY_ = 0;
   std::vector<float> * TrackCaloV_ = 0;
   std::vector<float> * TrackCaloAll_ = 0;
   std::vector<double> * crthit_x_ = 0;
   std::vector<double> * crthit_y_ = 0;
   std::vector<double> * crthit_z_ = 0;
   std::vector<double> * crt_plane_ = 0;
   std::vector<double> * a_crthit_x_ = 0;
   std::vector<double> * a_crthit_y_ = 0;
   std::vector<double> * a_crthit_z_ = 0;  
   std::vector<double> * a_crt_plane_ = 0;  
  
  std::vector<int> * ShowerPfp_ = 0;
  std::vector<double> * ShowerScore_ = 0;
  std::vector<double> * ShowerDir_x_ = 0;
  std::vector<double> * ShowerDir_y_ = 0;
  std::vector<double> * ShowerDir_z_ = 0;
  std::vector<int> * ShowerNHitsU_ = 0;
  std::vector<int> * ShowerNHitsV_ = 0;
  std::vector<int> * ShowerNHitsY_ = 0;
  std::vector<double> * ShowerCaloU_ = 0;
  std::vector<double> * ShowerCaloV_ = 0;
  std::vector<double> * ShowerCaloY_ = 0;
  // end track variables
  int muon_candidate_key = 0;
  int muon_candidate_pfp = 0;
  
  double TriTim_sec_ = 0;          //event trigger time sec
  double TriTim_nsec_ = 0;          //event trigger time ns
  
  // CRT in beam variables
  
  int nr_crthit_ = -1; // # crt hits assigned to a tpc track
  std::vector<double> * crthit_ts0_ = 0;
  std::vector<double> * crthit_ts1_ = 0;
  std::vector<int> * adc_length_ = 0;
  std::vector<double> * crt_adc_ = 0;
  std::vector<int> * crtbeam_hit_nr_ = 0;
  
  double TimFla_ = -99;
  double flash_PE_ = -99;
  double flash_y_ = -999;
  double flash_z_ = -999;
  
  
  //mc
  // MC neutrino info
  uint NuMCnu = 0; // number of MC neutrinos in event, only one gets saved!
  int MCNu_Interaction = 0;
  int MCNu_CCNC = 0;
  int MCNu_PDG = 0;
  float MCNu_Energy = 0;
  float MCNu_leptonPx, MCNu_leptonPy, MCNu_leptonPz = 0;
  float MCNu_LeptonEnergy = 0;
  float MCNu_Px, MCNu_Py, MCNu_Pz = 0;
  float MCNu_leptonTheta = 0;
  //float MCNu_leptonPhi = 0;
  float MCNu_time = 0; // time of the true neutrino interaction
  float MCNu_Vx, MCNu_Vy, MCNu_Vz = 0;
  float MCNu_VxSce, MCNu_VySce, MCNu_VzSce = 0;
  float MCNu_vertexDistance = 0;
  
  // Matched MCParticle info
  bool MCNU_matched = 0;
  bool MCCosmic_matched = 0;
  
  std::vector<int> *  MCle_pfp_ = 0;
  std::vector<int> *  MCle_key_ = 0;
  std::vector<int> *  MCle_PDG_ = 0;
  std::vector<float> * MCle_Energy_ = 0;
  std::vector<float> * MCle_Px_ = 0;
  std::vector<float> * MCle_Py_ = 0;
  std::vector<float> * MCle_Pz_ = 0;
  //std::vector<float> * MCle_purity_ = 0;
  std::vector<float> * MCle_Vx_ = 0;
  std::vector<float> * MCle_Vy_ = 0;
  std::vector<float> * MCle_Vz_ = 0;
  std::vector<float> * MCle_Endx_ = 0;
  std::vector<float> * MCle_Endy_ = 0;
  std::vector<float> * MCle_Endz_ = 0;
  std::vector<float> * MCle_length_ = 0;
  std::vector<float> * MCle_VxSce_ = 0;
  std::vector<float> * MCle_VySce_ = 0;
  std::vector<float> * MCle_VzSce_ = 0;
  std::vector<float> * MCle_Theta_ = 0;
  std::vector<float> * MCle_Phi_ = 0;
  
  std::vector<double> * MCTrackPurity_ = 0;
  std::vector<double> * MCTrackPDG_ = 0;
  std::vector<double> * MCTrackEnergy_ = 0;
  std::vector<double> * MCTrackMomentum_ = 0;
  std::vector<double> * MCTrackTheta_ = 0;
  std::vector<double> * MCTrackPhi_ = 0;
  std::vector<double> * MCTrackLength_ = 0;
  
  std::vector<double> * MCTrackStart_x_ = 0;
  std::vector<double> * MCTrackStart_y_ = 0;
  std::vector<double> * MCTrackStart_z_ = 0;
  std::vector<double> * MCTrackEnd_x_ = 0;
  std::vector<double> * MCTrackEnd_y_ = 0;
  std::vector<double> * MCTrackEnd_z_ = 0;
  
  double Genie_Q2 = -1;
  double Genie_q2 = -1;
  double Genie_W = -1;
  double Genie_T = -1;
  double Genie_X = -1;
  double Genie_Y = -1;
  int Genie_nNeutron_preFSI = -1;
  int Genie_nProton_preFSI = -1;
  int Genie_nPi0_preFSI = -1;
  int Genie_nPiPlus_preFSI = -1;
  int Genie_nPiMinus_preFSI = -1;
  
  double Nu_Vx_sce_ = -999, Nu_Vy_sce_ = -999, Nu_Vz_sce_ = -999;
  int Nu_NhitsU_ = -1, Nu_NhitsV_ = -1, Nu_NhitsY_ = -1;
  float Nu_CaloU_ = -1, Nu_CaloV_ = -1, Nu_CaloY_ = -1;
  
  float MCNu_Vx_sce_ = -999, MCNu_Vy_sce_ = -999, MCNu_Vz_sce_ = -999;
  std::vector<float> *  MCle_Vx_sce_;
  std::vector<float> *  MCle_Vy_sce_;
  std::vector<float> *  MCle_Vz_sce_;
  
  double crt_trig_corr_mean = 0;
  double crt_trig_corr_med = 0;
  if(verbose_ !=0) printf("Before t_in branching\n");
  
  t_in->SetBranchAddress("event_counter", &event_counter);
  t_in->SetBranchAddress("frunNum", &frunNum);
  t_in->SetBranchAddress("fsubRunNum", &fsubRunNum);
  t_in->SetBranchAddress("fEvtNum", &fEvtNum);
  t_in->SetBranchAddress("EventWeight", &EventWeight);
  if(with_genie != -1)t_in->SetBranchAddress("TunedCentralValue_Genie", &TunedCentralValue_Genie_);
  if(with_genie == 1) t_in->SetBranchAddress("All_Genie", &All_Genie);
  
  if(with_genie == 10) t_in->SetBranchAddress("expskin_FluxUnisim", &expskin_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("horncurrent_FluxUnisim", &horncurrent_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("kminus_PrimaryHadronNormalization", &kminus_PrimaryHadronNormalization_);
  if(with_genie == 10) t_in->SetBranchAddress("kplus_PrimaryHadronFeynmanScaling", &kplus_PrimaryHadronFeynmanScaling_);
  if(with_genie == 10) t_in->SetBranchAddress("kzero_PrimaryHadronSanfordWang", &kzero_PrimaryHadronSanfordWang_);
  if(with_genie == 10) t_in->SetBranchAddress("nucleoninexsec_FluxUnisim", &nucleoninexsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("nucleonqexsec_FluxUnisim", &nucleonqexsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("nucleontotxsec_FluxUnisim", &nucleontotxsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("piminus_PrimaryHadronSWCentralSplineVariation", &piminus_PrimaryHadronSWCentralSplineVariation_);
  if(with_genie == 10) t_in->SetBranchAddress("pioninexsec_FluxUnisim", &pioninexsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("pionqexsec_FluxUnisim", &pionqexsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("piontotxsec_FluxUnisim", &piontotxsec_FluxUnisim_);
  if(with_genie == 10) t_in->SetBranchAddress("piplus_PrimaryHadronSWCentralSplineVariation", &piplus_PrimaryHadronSWCentralSplineVariation_);
  
  if(with_genie == 10) t_in->SetBranchAddress("AxFFCCQEshape_Genie", &AxFFCCQEshape_Genie_);
  if(with_genie == 10) t_in->SetBranchAddress("DecayAngMEC_Genie", &DecayAngMEC_Genie_);
  if(with_genie == 10) t_in->SetBranchAddress("MaNCRES_Genie", &MaNCRES_Genie_);
  if(with_genie == 10) t_in->SetBranchAddress("Theta_Delta2Npi_Genie", &Theta_Delta2Npi_Genie_);
  if(with_genie == 10) t_in->SetBranchAddress("VecFFCCQEshape_Genie", &VecFFCCQEshape_Genie_);

  t_in->SetBranchAddress("TriTim_sec",        &TriTim_sec_);
  t_in->SetBranchAddress("TriTim_nsec",       &TriTim_nsec_);
  
  t_in->SetBranchAddress("NuScore",           &NuScore_);
  t_in->SetBranchAddress("FlashScore",        &FlashScore_);
  t_in->SetBranchAddress("FlashScoreTime",    &FlashScoreTime_);
  t_in->SetBranchAddress("NuPDG",             &NuPDG_);
  t_in->SetBranchAddress("NumPfp",            &NumPfp_);
  t_in->SetBranchAddress("Nu_Vx",             &Nu_Vx_);
  t_in->SetBranchAddress("Nu_Vy",             &Nu_Vy_);
  t_in->SetBranchAddress("Nu_Vz",             &Nu_Vz_);
  t_in->SetBranchAddress("Nu_Vx_sce",             &Nu_Vx_sce_);
  t_in->SetBranchAddress("Nu_Vy_sce",             &Nu_Vy_sce_);
  t_in->SetBranchAddress("Nu_Vz_sce",             &Nu_Vz_sce_);
  
  t_in->SetBranchAddress("NuTracks",          &NuTracks_);
  t_in->SetBranchAddress("NuShowers",         &NuShowers_);
  
  t_in->SetBranchAddress("Nu_NhitsU",         &Nu_NhitsU_);
  t_in->SetBranchAddress("Nu_NhitsV",         &Nu_NhitsV_);
  t_in->SetBranchAddress("Nu_NhitsY",         &Nu_NhitsY_);
  t_in->SetBranchAddress("Nu_CaloU",         &Nu_CaloU_);
  t_in->SetBranchAddress("Nu_CaloV",         &Nu_CaloV_);
  t_in->SetBranchAddress("Nu_CaloY",         &Nu_CaloY_);
  
  t_in->SetBranchAddress("muon_candidate_key",&muon_candidate_key);
  t_in->SetBranchAddress("muon_candidate_pfp",&muon_candidate_pfp);
  
  t_in->SetBranchAddress("track_key",        &track_key_);
  t_in->SetBranchAddress("TrackPID",         &TrackPID_);
  
  t_in->SetBranchAddress("Vx",                &Vx_);
  t_in->SetBranchAddress("Vy",                &Vy_);
  t_in->SetBranchAddress("Vz",                &Vz_);
  t_in->SetBranchAddress("Vx_sce",                &Vx_sce_);
  t_in->SetBranchAddress("Vy_sce",                &Vy_sce_);
  t_in->SetBranchAddress("Vz_sce",                &Vz_sce_);
  
  t_in->SetBranchAddress("TrackScore",        &TrackScore_);
  t_in->SetBranchAddress("ShowerScore",      &ShowerScore_);
  t_in->SetBranchAddress("TrackScoreGlobal",        &TrackScoreGlobal_);
  t_in->SetBranchAddress("VtxDistance",       &VtxDistance_);
  t_in->SetBranchAddress("VtxDistance_sce",       &VtxDistance_sce_);
  t_in->SetBranchAddress("TrackLength",       &TrackLength_);
  t_in->SetBranchAddress("TrackPfp",       &TrackPfp_);
  t_in->SetBranchAddress("isShowerTrack",       &isShowerTrack_);
  
  t_in->SetBranchAddress("TrackMomRange_p",   &TrackMomRange_p_);
  t_in->SetBranchAddress("TrackMomRange_mu",  &TrackMomRange_mu_);
  t_in->SetBranchAddress("TrackMomMCS_mom",   &TrackMomMCS_mom_);
  t_in->SetBranchAddress("TrackMomMCS_err",   &TrackMomMCS_err_);
  t_in->SetBranchAddress("TrackMomMCS_ll",    &TrackMomMCS_ll_);
  
  t_in->SetBranchAddress("TrackStart_x",      &TrackStart_x_);
  t_in->SetBranchAddress("TrackStart_y",      &TrackStart_y_);
  t_in->SetBranchAddress("TrackStart_z",      &TrackStart_z_);
  t_in->SetBranchAddress("TrackStart_x_sce",      &TrackStart_x_sce_);
  t_in->SetBranchAddress("TrackStart_y_sce",      &TrackStart_y_sce_);
  t_in->SetBranchAddress("TrackStart_z_sce",      &TrackStart_z_sce_);
  t_in->SetBranchAddress("TrackEnd_x",        &TrackEnd_x_);
  t_in->SetBranchAddress("TrackEnd_y",        &TrackEnd_y_);
  t_in->SetBranchAddress("TrackEnd_z",        &TrackEnd_z_);
  t_in->SetBranchAddress("TrackEnd_x_sce",        &TrackEnd_x_sce_);
  t_in->SetBranchAddress("TrackEnd_y_sce",        &TrackEnd_y_sce_);
  t_in->SetBranchAddress("TrackEnd_z_sce",        &TrackEnd_z_sce_);
  t_in->SetBranchAddress("TrackDir_x",        &TrackDir_x_);
  t_in->SetBranchAddress("TrackDir_y",        &TrackDir_y_);
  t_in->SetBranchAddress("TrackDir_z",        &TrackDir_z_);
  t_in->SetBranchAddress("TrackTheta",        &TrackTheta_);
  t_in->SetBranchAddress("TrackPhi",          &TrackPhi_);
  
  if(has_space_point){
  t_in->SetBranchAddress("AllTrack_point_x",          &AllTrack_point_x_);
  t_in->SetBranchAddress("AllTrack_point_y",          &AllTrack_point_y_);
  t_in->SetBranchAddress("AllTrack_point_z",          &AllTrack_point_z_);
  t_in->SetBranchAddress("SpacePoint_x",          &SpacePoint_x_);
  t_in->SetBranchAddress("SpacePoint_y",          &SpacePoint_y_);
  t_in->SetBranchAddress("SpacePoint_z",          &SpacePoint_z_);
  }
  if(has_wire){
  t_in->SetBranchAddress("Wire_id",          &Wire_id_);
  t_in->SetBranchAddress("Wire_plane",          &Wire_plane_);
  }
  t_in->SetBranchAddress("TrackNHitsU",       &TrackNHitsU_);
  t_in->SetBranchAddress("TrackNHitsV",       &TrackNHitsV_);
  t_in->SetBranchAddress("TrackNHitsY",       &TrackNHitsY_);
  t_in->SetBranchAddress("TrackCaloU",        &TrackCaloU_);
  t_in->SetBranchAddress("TrackCaloV",        &TrackCaloV_);
  t_in->SetBranchAddress("TrackCaloY",        &TrackCaloY_);
  //t_in->SetBranchAddress("TrackCaloAll",        &TrackCaloAll_);
  
  t_in->SetBranchAddress("ShowerPfp",      &ShowerPfp_);
  //t_in->SetBranchAddress("ShowerLength",      &ShowerLength_);
  //t_in->SetBranchAddress("ShowerOpenAngle",   &ShowerOpenAngle_);
  //t_in->SetBranchAddress("ShowerEnergy",   &ShowerEnergy_);
  //t_in->SetBranchAddress("ShowerMIPEnergy",   &ShowerMIPEnergy_);
  t_in->SetBranchAddress("ShowerDir_x",       &ShowerDir_x_);
  t_in->SetBranchAddress("ShowerDir_y",       &ShowerDir_y_);
  t_in->SetBranchAddress("ShowerDir_z",       &ShowerDir_z_);
  
  //t_in->SetBranchAddress("Shower_dEdxU",      &Shower_dEdxU_);
  //t_in->SetBranchAddress("Shower_dEdxV",      &Shower_dEdxV_);
  //t_in->SetBranchAddress("Shower_dEdxY",      &Shower_dEdxY_);
  
  //t_in->SetBranchAddress("Shower_dEdxHitsU",  &Shower_dEdxHitsU_);
  //t_in->SetBranchAddress("Shower_dEdxHitsV",  &Shower_dEdxHitsV_);
  //t_in->SetBranchAddress("Shower_dEdxHitsY",  &Shower_dEdxHitsY_);
  
  //t_in->SetBranchAddress("Shower_dEdxPitchU", &Shower_dEdxPitchU_);
  //t_in->SetBranchAddress("Shower_dEdxPitchV", &Shower_dEdxPitchV_);
  //t_in->SetBranchAddress("Shower_dEdxPitchY", &Shower_dEdxPitchY_);
  
  t_in->SetBranchAddress("ShowerNHitsU",      &ShowerNHitsU_);
  t_in->SetBranchAddress("ShowerNHitsV",      &ShowerNHitsV_);
  t_in->SetBranchAddress("ShowerNHitsY",      &ShowerNHitsY_);
  
  t_in->SetBranchAddress("ShowerCaloU",       &ShowerCaloU_);
  t_in->SetBranchAddress("ShowerCaloV",       &ShowerCaloV_);
  t_in->SetBranchAddress("ShowerCaloY",       &ShowerCaloY_);
  
  t_in->SetBranchAddress("crtt0_time",        &crtt0_time_);
  t_in->SetBranchAddress("crtt0_trig",        &crtt0_trig_);
  t_in->SetBranchAddress("crtt0_DCA",         &crtt0_DCA_);
  t_in->SetBranchAddress("crtt0_plane",       &crtt0_plane_);
  
  t_in->SetBranchAddress("TrackPID_chiproton",&TrackPID_chiproton_);
  t_in->SetBranchAddress("TrackPID_chipion",  &TrackPID_chipion_);
  t_in->SetBranchAddress("TrackPID_chikaon",  &TrackPID_chikaon_);
  t_in->SetBranchAddress("TrackPID_chimuon",  &TrackPID_chimuon_);

  t_in->SetBranchAddress("nr_crthit",         &nr_crthit_);
  t_in->SetBranchAddress("crthit_ts0",        &crthit_ts0_);
  t_in->SetBranchAddress("crthit_ts1",        &crthit_ts1_);
  t_in->SetBranchAddress("crthit_x",        &crthit_x_);
  t_in->SetBranchAddress("crthit_y",        &crthit_y_);
  t_in->SetBranchAddress("crthit_z",        &crthit_z_);
  t_in->SetBranchAddress("crt_plane",        &crt_plane_);

  t_in->SetBranchAddress("adc_length",        &adc_length_);
  t_in->SetBranchAddress("crt_adc",           &crt_adc_);
  t_in->SetBranchAddress("crtbeam_hit_nr",   &crtbeam_hit_nr_);

  t_in->SetBranchAddress("a_crthit_ts0",      &a_crthit_ts0);
  t_in->SetBranchAddress("a_crthit_ts1",      &a_crthit_ts1);
  t_in->SetBranchAddress("a_crthit_x",        &a_crthit_x_);
  t_in->SetBranchAddress("a_crthit_y",        &a_crthit_y_);
  t_in->SetBranchAddress("a_crthit_z",        &a_crthit_z_);
  t_in->SetBranchAddress("a_crt_plane",        &a_crt_plane_);
  t_in->SetBranchAddress("a_adc_length",      &a_adc_length);
  t_in->SetBranchAddress("a_crt_adc",         &a_crt_adc);
  t_in->SetBranchAddress("a_t0_counter",      &a_t0_counter);
  
  t_in->SetBranchAddress("TimFla",            &TimFla_);
  t_in->SetBranchAddress("flash_PE",          &flash_PE_);
  t_in->SetBranchAddress("flash_y",           &flash_y_);
  t_in->SetBranchAddress("flash_z",           &flash_z_);
  
  t_in->SetBranchAddress("crt_trig_corr_mean",&crt_trig_corr_mean);
  t_in->SetBranchAddress("crt_trig_corr_med", &crt_trig_corr_med);

  
  if(!is_data_){
    t_in->SetBranchAddress("MCNu_Interaction",&MCNu_Interaction);
    t_in->SetBranchAddress("MCNu_CCNC",           &MCNu_CCNC);
    t_in->SetBranchAddress("MCNu_PDG",            &MCNu_PDG);
    t_in->SetBranchAddress("MCNu_Energy",         &MCNu_Energy);
    t_in->SetBranchAddress("MCNu_leptonPx",       &MCNu_leptonPx);
    t_in->SetBranchAddress("MCNu_leptonPy",       &MCNu_leptonPy);
    t_in->SetBranchAddress("MCNu_leptonPz",       &MCNu_leptonPz);
    t_in->SetBranchAddress("MCNu_LeptonEnergy",   &MCNu_LeptonEnergy);
    t_in->SetBranchAddress("MCNu_Px",             &MCNu_Px);
    t_in->SetBranchAddress("MCNu_Py",             &MCNu_Py);
    t_in->SetBranchAddress("MCNu_Pz",             &MCNu_Pz);
    t_in->SetBranchAddress("MCNu_leptonTheta",    &MCNu_leptonTheta);
    //t_in->SetBranchAddress("MCNu_leptonPhi",      &MCNu_leptonPhi,      "MCNu_leptonPhi/F");
    t_in->SetBranchAddress("MCNu_time",           &MCNu_time);
    t_in->SetBranchAddress("MCNu_Vx",             &MCNu_Vx);
    t_in->SetBranchAddress("MCNu_Vy",             &MCNu_Vy);
    t_in->SetBranchAddress("MCNu_Vz",             &MCNu_Vz);
    t_in->SetBranchAddress("MCNu_VxSce",          &MCNu_VxSce);
    t_in->SetBranchAddress("MCNu_VySce",          &MCNu_VySce);
    t_in->SetBranchAddress("MCNu_VzSce",          &MCNu_VzSce);
    t_in->SetBranchAddress("MCNu_vertexDistance",    &MCNu_vertexDistance);
    
    t_in->SetBranchAddress("MCle_pfp",               &MCle_pfp_);
    t_in->SetBranchAddress("MCle_key",               &MCle_key_);
    t_in->SetBranchAddress("MCle_PDG",               &MCle_PDG_);
    //t_in->SetBranchAddress("MCle_purity",            &MCle_purity_);
    
    t_in->SetBranchAddress("MCle_Energy",            &MCle_Energy_);
    
    t_in->SetBranchAddress("MCle_Px",                &MCle_Px_);
    t_in->SetBranchAddress("MCle_Py",                &MCle_Py_);
    t_in->SetBranchAddress("MCle_Pz",                &MCle_Pz_);
    t_in->SetBranchAddress("MCle_Vx",                &MCle_Vx_);
    t_in->SetBranchAddress("MCle_Vy",                &MCle_Vy_);
    t_in->SetBranchAddress("MCle_Vz",                &MCle_Vz_);
    
    
    t_in->SetBranchAddress("MCle_Endx",                &MCle_Endx_);
    t_in->SetBranchAddress("MCle_Endy",                &MCle_Endy_);
    t_in->SetBranchAddress("MCle_Endz",                &MCle_Endz_);
    t_in->SetBranchAddress("MCle_length",            &MCle_length_);
    t_in->SetBranchAddress("MCle_VxSce",             &MCle_VxSce_);
    t_in->SetBranchAddress("MCle_VySce",             &MCle_VySce_);
    t_in->SetBranchAddress("MCle_VzSce",             &MCle_VzSce_);
    t_in->SetBranchAddress("MCle_Theta",             &MCle_Theta_);
    t_in->SetBranchAddress("MCle_Phi",             &MCle_Phi_);  
    
    
    t_in->SetBranchAddress("MCTrackPurity",             &MCTrackPurity_);
    
    t_in->SetBranchAddress("MCTrackPDG",             &MCTrackPDG_);
    t_in->SetBranchAddress("MCTrackEnergy",             &MCTrackEnergy_);
    t_in->SetBranchAddress("MCTrackMomentum",             &MCTrackMomentum_);
    t_in->SetBranchAddress("MCTrackTheta",             &MCTrackTheta_);
    t_in->SetBranchAddress("MCTrackPhi",             &MCTrackPhi_);
    t_in->SetBranchAddress("MCTrackLength",             &MCTrackLength_);
    
    t_in->SetBranchAddress("MCTrackStart_x",             &MCTrackStart_x_);
    t_in->SetBranchAddress("MCTrackStart_y",             &MCTrackStart_y_);
    t_in->SetBranchAddress("MCTrackStart_z",             &MCTrackStart_z_);
    t_in->SetBranchAddress("MCTrackEnd_x",             &MCTrackEnd_x_);
    t_in->SetBranchAddress("MCTrackEnd_y",             &MCTrackEnd_y_);
    t_in->SetBranchAddress("MCTrackEnd_z",             &MCTrackEnd_z_);
    
    t_in->SetBranchAddress("Genie_Q2",&Genie_Q2);
    t_in->SetBranchAddress("Genie_q2",&Genie_q2);
    t_in->SetBranchAddress("Genie_W",&Genie_W);
    t_in->SetBranchAddress("Genie_T",&Genie_T);
    t_in->SetBranchAddress("Genie_X",&Genie_X);
    t_in->SetBranchAddress("Genie_Y",&Genie_Y);
    t_in->SetBranchAddress("Genie_nNeutron_preFSI",&Genie_nNeutron_preFSI);
    t_in->SetBranchAddress("Genie_nProton_preFSI",&Genie_nProton_preFSI);
    t_in->SetBranchAddress("Genie_nPi0_preFSI",&Genie_nPi0_preFSI);
    t_in->SetBranchAddress("Genie_nPiPlus_preFSI",&Genie_nPiPlus_preFSI);
    t_in->SetBranchAddress("Genie_nPiMinus_preFSI",&Genie_nPiMinus_preFSI);
    
  }
  
  if(verbose_ !=0) printf("After t_in branching\n");
  //data=fopen(filename,"r");
  char name[120] = "";
  //name = filename;
  strcat(name,filename);
  if(with_genie == 1) strcat(name,"out33wG.root");
  else strcat(name,"out33.root");
  //strcat(name,"out33.root");
  TFile f_out(name,"RECREATE");
  //TFile f_out(name,"UPDATE");
  //TFile f("ProdRun20170101_011005-crt04.1.pairs.root","RECREATE");
  
  TTree * t_out = new TTree("t_out","Single");
  int track_key_out_ = 0;
  
  int track_key_out = 0;
  int  TrackPID_out = 0;

  double a_crthit_ts0out = 0;
  double a_crthit_ts1out = 0;
  int  a_adc_lengthout = 0;
  double a_crt_adcout = 0;
  int  a_t0_counterout = 0;

  double crtt0_time_out = 0;
  int  crtt0_trig_out = 0;
  double crtt0_DCA_out = 0;
  int  crtt0_plane_out = 0;

  double VtxDistance_out = 0;
  double Vx_out = 0;
  double Vy_out = 0;
  double Vz_out = 0;
  double TrackScore_out = 0;
  double TrackScoreGlobal_out = 0;
  double TrackLength_out = 0;
  double TrackPID_chiproton_out = 0;
  double TrackPID_chimuon_out = 0;
  double TrackPID_chipion_out = 0;
  double TrackPID_chikaon_out = 0;

  float TrackMomRange_p_out = 0;
  float TrackMomRange_mu_out = 0;
  float TrackMomMCS_mom_out = 0;
  float TrackMomMCS_err_out = 0;
  float TrackMomMCS_ll_out = 0;

  double TrackStart_x_out = 0;
  double TrackStart_y_out = 0;
  double TrackStart_z_out = 0;
  double TrackEnd_x_out = 0;
  double TrackEnd_y_out = 0;
  double TrackEnd_z_out = 0;
  double TrackDir_x_out = 0;
  double TrackDir_y_out = 0;
  double TrackDir_z_out = 0;

  double TrackTheta_out = 0;
  double TrackPhi_out = 0;
  
  int ShowerPfp_out = 0;
  double ShowerDir_x_out = 0;
  double ShowerDir_y_out = 0;
  double ShowerDir_z_out = 0;
  int ShowerNHitsU_out = 0;
  int ShowerNHitsV_out = 0;
  int ShowerNHitsY_out = 0;
  double ShowerCaloU_out = 0;
  double ShowerCaloV_out = 0;
  double ShowerCaloY_out = 0;
  double ShowerCaloAll_out = 0;
  int ShowerNHitsAll_out = 0;
  int NumShowers_corr = 0;
  double ShowerScore_out = 0;
  
  double crthit_ts0_out = 0;
  double crthit_ts1_out = 0;
  int  adc_length_out = 0;
  double crt_adc_out = 0;
  int  crtbeam_hit_nr_out = 0;
  
  double Vx_sce_out = 0, Vy_sce_out = 0, Vz_sce_out = 0;
  double VtxDistance_sce_out = 0;
  int TrackPfp_out = 0;
  int isShowerTrack_out = 0;
  
  double TrackStart_x_sce_out = 0 ,TrackStart_y_sce_out = 0 ,TrackStart_z_sce_out = 0;
  double TrackEnd_x_sce_out = 0 ,TrackEnd_y_sce_out = 0 ,TrackEnd_z_sce_out = 0;
  int TrackNHitsU_out = 0, TrackNHitsV_out = 0,TrackNHitsY_out = 0;
  double TrackCaloU_out = -1,TrackCaloV_out = -1,TrackCaloY_out = -1;
  double TrackCaloAll_out = 0;
  int TrackNHitsAll_out = 0;
  
  double crthit_x_out = 0, crthit_y_out = -9999, crthit_z_out = 9999, crthit_plane_out = -1;
  double a_crthit_x_out = 0, a_crthit_y_out = 0, a_crthit_z_out = 0, a_crthit_plane_out = -1;
  
  int top_y_layer = 0;
  int y_box = 0;
  int z_dead = 0;
  int front_strip = 0;
  int z_dead_y = 0;
  int front_strip_y = 0;
  int front_z_layer = 0;
  
  int lowest_Track_point = 0;
  
  int MCle_pfp_out = 0;
  int MCle_key_out = 0;
  int MCle_PDG_out = 0;
  float MCle_Energy_out = 0;
  float MCle_Px_out = 0;
  float MCle_Py_out = 0;
  float MCle_Pz_out = 0;
  //float MCle_purity_out = 0;
  float MCle_Vx_out, MCle_Vy_out, MCle_Vz_out = 0;
  float MCle_Endx_out = 0, MCle_Endy_out = 0, MCle_Endz_out = 0;
  float MCle_length_out = 0;
  float MCle_Theta_out = 0;
  float MCle_Phi_out = 0;
  float MCle_VxSce_out, MCle_VySce_out, MCle_VzSce_out = 0;
  float MCle_Vx_sce_out, MCle_Vy_sce_out, MCle_Vz_sce_out = 0;
  
  double MCTrackPurity_out = 0, MCTrackPDG_out = 0, MCTrackEnergy_out = 0;
  double MCTrackMomentum_out = 0, MCTrackTheta_out = 0, MCTrackPhi_out = 0;
  double MCTrackLength_out = 0;
  double MCTrackStart_x_out = 0, MCTrackStart_y_out = 0, MCTrackStart_z_out = 0;
  double MCTrackEnd_x_out = 0, MCTrackEnd_y_out = 0, MCTrackEnd_z_out = 0;
  //initalize the variavles
  //t_out->Branch("muon_candidate_key",&muon_candidate_key,"muon_candidate_key /I");
  //t_out->Branch("track_key",&track_key_out_,"track_key /I");
  if(verbose_ !=0) printf("After t_out branching\n");
  
  t_out->Branch("event_counter", &event_counter, "event_counter/I");
  t_out->Branch("frunNum", &frunNum, "Run Number/i");
  t_out->Branch("fsubRunNum", &fsubRunNum, "SubRun Number/i");
  t_out->Branch("fEvtNum", &fEvtNum, "Event Number/i");
  t_out->Branch("EventWeight", &EventWeight, "EventWeight/D");
  t_out->Branch("TunedCentralValue_Genie", &TunedCentralValue_Genie_out,"TunedCentralValue_Genie/D");
  if(with_genie == 1) t_out->Branch("All_Genie", &All_Genie);
  
  if(with_genie == 10) t_out->Branch("expskin_FluxUnisim", &expskin_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("horncurrent_FluxUnisim", &horncurrent_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("kminus_PrimaryHadronNormalization", &kminus_PrimaryHadronNormalization_);
  if(with_genie == 10) t_out->Branch("kplus_PrimaryHadronFeynmanScaling", &kplus_PrimaryHadronFeynmanScaling_);
  if(with_genie == 10) t_out->Branch("kzero_PrimaryHadronSanfordWang", &kzero_PrimaryHadronSanfordWang_);
  if(with_genie == 10) t_out->Branch("nucleoninexsec_FluxUnisim", &nucleoninexsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("nucleonqexsec_FluxUnisim", &nucleonqexsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("nucleontotxsec_FluxUnisim", &nucleontotxsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("piminus_PrimaryHadronSWCentralSplineVariation", &piminus_PrimaryHadronSWCentralSplineVariation_);
  if(with_genie == 10) t_out->Branch("pioninexsec_FluxUnisim", &pioninexsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("pionqexsec_FluxUnisim", &pionqexsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("piontotxsec_FluxUnisim", &piontotxsec_FluxUnisim_);
  if(with_genie == 10) t_out->Branch("piplus_PrimaryHadronSWCentralSplineVariation", &piplus_PrimaryHadronSWCentralSplineVariation_);
  
  if(with_genie == 10) t_out->Branch("AxFFCCQEshape_Genie", &AxFFCCQEshape_Genie_);
  if(with_genie == 10) t_out->Branch("DecayAngMEC_Genie", &DecayAngMEC_Genie_);
  if(with_genie == 10) t_out->Branch("MaNCRES_Genie", &MaNCRES_Genie_);
  if(with_genie == 10) t_out->Branch("Theta_Delta2Npi_Genie", &Theta_Delta2Npi_Genie_);
  if(with_genie == 10) t_out->Branch("VecFFCCQEshape_Genie", &VecFFCCQEshape_Genie_);
    
  t_out->Branch("TriTim_sec",        &TriTim_sec_,            "TriTim_sec/D");
  t_out->Branch("TriTim_nsec",       &TriTim_nsec_,           "TriTim_nsec/D");
  
  t_out->Branch("NuScore",           &NuScore_,               "NuScore/D");
  t_out->Branch("FlashScore",        &FlashScore_,            "FlashScore/D");
  t_out->Branch("FlashScoreTime",    &FlashScoreTime_,        "FlashScoreTime/D");
  t_out->Branch("NuPDG",             &NuPDG_,                 "NuPDG/I");
  t_out->Branch("NumPfp",            &NumPfp_,                "NumPfp/I");
  t_out->Branch("Nu_Vx",             &Nu_Vx_,                 "Nu_Vx/D");
  t_out->Branch("Nu_Vy",             &Nu_Vy_,                 "Nu_Vy/D");
  t_out->Branch("Nu_Vz",             &Nu_Vz_,                 "Nu_Vz/D");
  t_out->Branch("Nu_Vx_sce",             &Nu_Vx_sce_,                 "Nu_Vx_sce/D");
  t_out->Branch("Nu_Vy_sce",             &Nu_Vy_sce_,                 "Nu_Vy_sce/D");
  t_out->Branch("Nu_Vz_sce",             &Nu_Vz_sce_,                 "Nu_Vz_sce/D");
  
  
  t_out->Branch("NuTracks",          &NuTracks_,              "NuTracks/I");
  t_out->Branch("NuShowers",         &NuShowers_,             "NuShowers/I");
  
  t_out->Branch("Nu_NhitsU",         &Nu_NhitsU_,             "Nu_NhitsU/I");
  t_out->Branch("Nu_NhitsV",         &Nu_NhitsV_,             "Nu_NhitsV/I");
  t_out->Branch("Nu_NhitsY",         &Nu_NhitsY_,             "Nu_NhitsY/I");
  t_out->Branch("Nu_CaloU",         &Nu_CaloU_,             "Nu_CaloU/F");
  t_out->Branch("Nu_CaloV",         &Nu_CaloV_,             "Nu_CaloV/F");
  t_out->Branch("Nu_CaloY",         &Nu_CaloY_,             "Nu_CaloY/F");
  
  t_out->Branch("muon_candidate_key",&muon_candidate_key,     "muon_candidate_key/I");
  t_out->Branch("muon_candidate_pfp",&muon_candidate_pfp,     "muon_candidate_pfp/I");
  /*
  t_out->Branch("track_key",        &track_key_out,              "track_key_out/I");
  t_out->Branch("TrackPID",        &TrackPID_out,              "TrackPID_out/I");
  
  t_out->Branch("Vx",                &Vx_out,              "Vx_out/D");
  t_out->Branch("Vy",                &Vy_out,              "Vy_out/D");
  t_out->Branch("Vz",                &Vz_out,              "Vz_out/D");
  t_out->Branch("Vx_sce",                &Vx_sce_out,          "Vx_sce_out/D");
  t_out->Branch("Vy_sce",                &Vy_sce_out,          "Vy_sce_out/D");
  t_out->Branch("Vz_sce",                &Vz_sce_out,          "Vz_sce_out/D");
  
  t_out->Branch("TrackScore",        &TrackScore_out,              "TrackScore_out/D");
  t_out->Branch("TrackScoreGlobal",        &TrackScoreGlobal_out,              "TrackScoreGlobal_out/D");
  t_out->Branch("VtxDistance",       &VtxDistance_out,              "VtxDistance_out/D");
  t_out->Branch("VtxDistance_sce",       &VtxDistance_sce_out,              "VtxDistance_sce_out/D");
  t_out->Branch("TrackLength",       &TrackLength_out,              "TrackLength_out/D");
  t_out->Branch("TrackPfp",       &TrackPfp_out,                "TrackPfp_out/I");
  t_out->Branch("isShowerTrack",       &isShowerTrack_out,      "isShowerTrack_out/I");
  
  t_out->Branch("TrackMomRange_p",   &TrackMomRange_p_out,              "NuTracks/F");
  t_out->Branch("TrackMomRange_mu",  &TrackMomRange_mu_out,              "NuTracks/F");
  t_out->Branch("TrackMomMCS_mom",   &TrackMomMCS_mom_out,              "NuTracks/F");
  t_out->Branch("TrackMomMCS_err",   &TrackMomMCS_err_out,              "NuTracks/F");
  t_out->Branch("TrackMomMCS_ll",    &TrackMomMCS_ll_out,              "NuTracks/F");
  
  t_out->Branch("TrackStart_x",      &TrackStart_x_out,              "NuTracks/D");
  t_out->Branch("TrackStart_y",      &TrackStart_y_out,              "NuTracks/D");
  t_out->Branch("TrackStart_z",      &TrackStart_z_out,              "NuTracks/D");
  t_out->Branch("TrackStart_x_sce",      &TrackStart_x_sce_out,              "TrackStart_x_sce_out/D");
  t_out->Branch("TrackStart_y_sce",      &TrackStart_y_sce_out,              "TrackStart_x_sce_out/D");
  t_out->Branch("TrackStart_z_sce",      &TrackStart_z_sce_out,              "TrackStart_x_sce_out/D");
  t_out->Branch("TrackEnd_x",        &TrackEnd_x_out,              "NuTracks/D");
  t_out->Branch("TrackEnd_y",        &TrackEnd_y_out,              "NuTracks/D");
  t_out->Branch("TrackEnd_z",        &TrackEnd_z_out,              "NuTracks/D");
  t_out->Branch("TrackEnd_x_sce",        &TrackEnd_x_sce_out,              "TrackEnd_x_sce_out/D");
  t_out->Branch("TrackEnd_y_sce",        &TrackEnd_y_sce_out,              "TrackEnd_x_sce_out/D");
  t_out->Branch("TrackEnd_z_sce",        &TrackEnd_z_sce_out,              "TrackEnd_x_sce_out/D");
  t_out->Branch("TrackDir_x",        &TrackDir_x_out,              "NuTracks/D");
  t_out->Branch("TrackDir_y",        &TrackDir_y_out,              "NuTracks/D");
  t_out->Branch("TrackDir_z",        &TrackDir_z_out,              "NuTracks/D");
  t_out->Branch("TrackTheta",        &TrackTheta_out,              "NuTracks/D");
  t_out->Branch("TrackPhi",          &TrackPhi_out,              "NuTracks/D");
  */
  t_out->Branch("track_key",        &track_key_);
  t_out->Branch("TrackPID",        &TrackPID_);
  
  t_out->Branch("Vx",                &Vx_);
  t_out->Branch("Vy",                &Vy_);
  t_out->Branch("Vz",                &Vz_);
  t_out->Branch("Vx_sce",                &Vx_sce_);
  t_out->Branch("Vy_sce",                &Vy_sce_);
  t_out->Branch("Vz_sce",                &Vz_sce_);
  
  t_out->Branch("TrackScore",        &TrackScore_);
  t_out->Branch("TrackScoreGlobal",        &TrackScoreGlobal_);
  t_out->Branch("VtxDistance",       &VtxDistance_);
  t_out->Branch("VtxDistance_sce",       &VtxDistance_sce_);
  t_out->Branch("TrackLength",       &TrackLength_);
  t_out->Branch("TrackPfp",       &TrackPfp_);
  t_out->Branch("isShowerTrack",       &isShowerTrack_);
  
  t_out->Branch("TrackMomRange_p",   &TrackMomRange_p_);
  t_out->Branch("TrackMomRange_mu",  &TrackMomRange_mu_);
  t_out->Branch("TrackMomMCS_mom",   &TrackMomMCS_mom_);
  t_out->Branch("TrackMomMCS_err",   &TrackMomMCS_err_);
  t_out->Branch("TrackMomMCS_ll",    &TrackMomMCS_ll_);
  
  t_out->Branch("TrackStart_x",      &TrackStart_x_);
  t_out->Branch("TrackStart_y",      &TrackStart_y_);
  t_out->Branch("TrackStart_z",      &TrackStart_z_);
  t_out->Branch("TrackStart_x_sce",      &TrackStart_x_sce_);
  t_out->Branch("TrackStart_y_sce",      &TrackStart_y_sce_);
  t_out->Branch("TrackStart_z_sce",      &TrackStart_z_sce_);
  t_out->Branch("TrackEnd_x",        &TrackEnd_x_);
  t_out->Branch("TrackEnd_y",        &TrackEnd_y_);
  t_out->Branch("TrackEnd_z",        &TrackEnd_z_);
  t_out->Branch("TrackEnd_x_sce",        &TrackEnd_x_sce_);
  t_out->Branch("TrackEnd_y_sce",        &TrackEnd_y_sce_);
  t_out->Branch("TrackEnd_z_sce",        &TrackEnd_z_sce_);
  t_out->Branch("TrackDir_x",        &TrackDir_x_);
  t_out->Branch("TrackDir_y",        &TrackDir_y_);
  t_out->Branch("TrackDir_z",        &TrackDir_z_);
  t_out->Branch("TrackTheta",        &TrackTheta_);
  t_out->Branch("TrackPhi",          &TrackPhi_);
  
  
  t_out->Branch("AllTrack_point_x",          &AllTrack_point_x_);
  t_out->Branch("AllTrack_point_y",          &AllTrack_point_y_);
  t_out->Branch("AllTrack_point_z",          &AllTrack_point_z_);
  t_out->Branch("SpacePoint_x",          &SpacePoint_x_);
  t_out->Branch("SpacePoint_y",          &SpacePoint_y_);
  t_out->Branch("SpacePoint_z",          &SpacePoint_z_);
  
  t_out->Branch("TrackNHitsU",       &TrackNHitsU_);
  t_out->Branch("TrackNHitsV",       &TrackNHitsV_);
  t_out->Branch("TrackNHitsY",       &TrackNHitsY_);
  t_out->Branch("TrackCaloU",        &TrackCaloU_);
  t_out->Branch("TrackCaloV",        &TrackCaloV_);
  t_out->Branch("TrackCaloY",        &TrackCaloY_);
  t_out->Branch("TrackCaloAll",        &TrackCaloAll_);
  
  /*
  t_out->Branch("TrackNHitsU",       &TrackNHitsU_out,    "NuTracks/I");
  t_out->Branch("TrackNHitsV",       &TrackNHitsV_out,    "NuTracks/I");
  t_out->Branch("TrackNHitsY",       &TrackNHitsY_out,      "NuTracks/I");
  t_out->Branch("TrackNHitsAll",       &TrackNHitsAll_out,      "TrackNHitsAll/I");
  t_out->Branch("TrackCaloU",        &TrackCaloU_out,       "NuTracks/D");
  t_out->Branch("TrackCaloV",        &TrackCaloV_out,       "NuTracks/D");
  t_out->Branch("TrackCaloY",        &TrackCaloY_out,       "NuTracks/D");
  t_out->Branch("TrackCaloAll",        &TrackCaloAll_out,       "TrackCaloAll/D");
  */
  t_out->Branch("crtt0_time",        &crtt0_time_);
  t_out->Branch("crtt0_trig",        &crtt0_trig_);
  t_out->Branch("crtt0_DCA",         &crtt0_DCA_);
  t_out->Branch("crtt0_plane",       &crtt0_plane_);
  
  t_out->Branch("TrackPID_chiproton",&TrackPID_chiproton_);
  t_out->Branch("TrackPID_chipion",  &TrackPID_chipion_);
  t_out->Branch("TrackPID_chikaon",  &TrackPID_chikaon_);
  t_out->Branch("TrackPID_chimuon",  &TrackPID_chimuon_);
  /*
  t_out->Branch("TrackPID_chiproton",&TrackPID_chiproton_out,              "NuTracks/D");
  t_out->Branch("TrackPID_chipion",  &TrackPID_chipion_out,              "NuTracks/D");
  t_out->Branch("TrackPID_chikaon",  &TrackPID_chikaon_out,              "NuTracks/D");
  t_out->Branch("TrackPID_chimuon",  &TrackPID_chimuon_out,              "NuTracks/D");
  */
  t_out->Branch("ShowerPfp",       &ShowerPfp_out,    "ShowerPfp_out/I");
  t_out->Branch("ShowerDir_x",       &ShowerDir_x_out,    "ShowerDir_x/D");
  t_out->Branch("ShowerDir_y",       &ShowerDir_y_out,    "ShowerDir_y/D");
  t_out->Branch("ShowerDir_z",       &ShowerDir_z_out,    "ShowerDir_z/D");
  t_out->Branch("ShowerNHitsU",       &ShowerNHitsU_out,      "ShowerNHitsU/I");
  t_out->Branch("ShowerNHitsV",       &ShowerNHitsV_out,      "ShowerNHitsV/I");
  t_out->Branch("ShowerNHitsY",       &ShowerNHitsY_out,      "ShowerNHitsY/I");
  t_out->Branch("ShowerCaloU",        &ShowerCaloU_out,       "ShowerCaloU/D");
  t_out->Branch("ShowerCaloV",        &ShowerCaloV_out,       "ShowerCaloV/D");
  t_out->Branch("ShowerCaloY",        &ShowerCaloY_out,       "ShowerCaloY/D");
  
  t_out->Branch("ShowerCaloAll",        &ShowerCaloAll_out,       "ShowerCaloAll/D");
  t_out->Branch("ShowerNHitsAll",        &ShowerNHitsAll_out,       "ShowerNHitsAll/I");
  t_out->Branch("NumShowers_corr",        &NumShowers_corr,       "NumShowers_corr/I");
  t_out->Branch("ShowerScore",        &ShowerScore_out,       "ShowerScore/D");

  t_out->Branch("nr_crthit",         &nr_crthit_,              "nr_crthit_/I");
  t_out->Branch("crthit_ts0",        &crthit_ts0_);
  t_out->Branch("crthit_ts1",        &crthit_ts1_);
  t_out->Branch("crthit_x",        &crthit_x_);
  t_out->Branch("crthit_y",        &crthit_y_);
  t_out->Branch("crthit_z",        &crthit_z_);
  t_out->Branch("crt_plane",        &crt_plane_);
  
  t_out->Branch("adc_length",        &adc_length_);
  t_out->Branch("crt_adc",           &crt_adc_);
  t_out->Branch("crtbeam_hit_nr",   &crtbeam_hit_nr_);
  
  t_out->Branch("a_crthit_ts0",      &a_crthit_ts0);
  t_out->Branch("a_crthit_ts1",      &a_crthit_ts1);
  t_out->Branch("a_crthit_x",        &a_crthit_x_);
  t_out->Branch("a_crthit_y",        &a_crthit_y_);
  t_out->Branch("a_crthit_z",        &a_crthit_z_);
  t_out->Branch("a_crt_plane",        &a_crt_plane_);
  t_out->Branch("a_adc_length",      &a_adc_length);
  t_out->Branch("a_crt_adc",         &a_crt_adc);
  t_out->Branch("a_t0_counter",      &a_t0_counter);
  
  /*
  t_out->Branch("nr_crthit",         &nr_crthit_,              "nr_crthit_/I");
  t_out->Branch("crthit_ts0",        &crthit_ts0_out,              "NuTracks/D");
  t_out->Branch("crthit_ts1",        &crthit_ts1_out,              "NuTracks/D");
  t_out->Branch("crthit_x",        &crthit_x_out,              "NuTracks/D");
  t_out->Branch("crthit_y",        &crthit_y_out,              "NuTracks/D");
  t_out->Branch("crthit_z",        &crthit_z_out,              "NuTracks/D");
  t_out->Branch("crthit_plane",        &crthit_plane_out,              "NuTracks/D");
  t_out->Branch("adc_length",        &adc_length_out,              "NuTracks/I");
  t_out->Branch("crt_adc",           &crt_adc_out,              "NuTracks/D");
  t_out->Branch("crtbeam_hit_nr",   &crtbeam_hit_nr_out,              "NuTracks/I");

  t_out->Branch("a_crthit_ts0",      &a_crthit_ts0out,              "NuTracks/D");
  t_out->Branch("a_crthit_ts1",      &a_crthit_ts1out,              "NuTracks/D");
  t_out->Branch("a_crthit_x",      &a_crthit_x_out,              "NuTracks/D");
  t_out->Branch("a_crthit_y",      &a_crthit_y_out,              "NuTracks/D");
  t_out->Branch("a_crthit_z",      &a_crthit_z_out,              "NuTracks/D");
  t_out->Branch("a_crthit_plane",      &a_crthit_plane_out,              "NuTracks/D");

  t_out->Branch("a_adc_length",      &a_adc_lengthout,              "NuTracks/I");
  t_out->Branch("a_crt_adc",         &a_crt_adcout,              "NuTracks/D");
  t_out->Branch("a_t0_counter",      &a_t0_counterout,              "NuTracks/I");
  */
  
  t_out->Branch("TimFla",            &TimFla_,                "TimFla/D");
  t_out->Branch("flash_PE",          &flash_PE_,              "flash_PE/D");
  t_out->Branch("flash_y",           &flash_y_,               "flash_y/D");
  t_out->Branch("flash_z",           &flash_z_,               "flash_z/D");
  
  t_out->Branch("crt_trig_corr_mean",&crt_trig_corr_mean,     "crt_trig_corr_mean/D");
  t_out->Branch("crt_trig_corr_med", &crt_trig_corr_med,      "crt_trig_corr_med/D");
  
  int nr_crthit_beam_ = 0;
  int nr_crthit_beam_tres_ = 0;
  int nr_track_asso_out_ = 0;
  int nr_track_asso_ = 0;
  int nr_track_asso_uncon_ = 0;
  int nr_track_asso_out_uncon_ = 0;
  int nr_tracks_uncontained_ = 0;
  int track_uncontained_ = 0;
  
  int nr_crthit_top_ = 0;
  int crthit_vertex_zcut_ = 0;
  double crthit_vertex_z_ = 0;
  double crthit_vertex_z_adc = 0;
  
  int has_matched_muon_ = 0;
  int has_matched_muonHigh_ = 0;
  
  int small_wire = 0;
  int wire_plane = 100;
  
  double PID_muprotdiff = -999;
  int key_muprotdiff = -999;
  double PID_muprotratio = -999;
  int key_muprotratio = -999;
  double PID_mupiondiff = -999;
  int key_mupiondiff = -999;
  double PID_mupionratio = -999;
  int key_mupionratio = -999;
  double longest_track = -999;
  int key_longest_track = -999;
  int num_mc_muon = 0;
  
  t_out->Branch("nr_crthit_beam",         &nr_crthit_beam_,              "nr_crthit_beam_/I");
  t_out->Branch("nr_crthit_beam_tres",         &nr_crthit_beam_tres_,              "nr_crthit_beam_tres_/I");
  t_out->Branch("nr_track_asso_out",         &nr_track_asso_out_,              "nr_track_asso_out/I");
  t_out->Branch("nr_track_asso",         &nr_track_asso_,              "nr_track_asso/I");
  t_out->Branch("nr_track_asso_uncon",         &nr_track_asso_uncon_,              "nr_track_asso_uncon/I");
  t_out->Branch("nr_track_asso_out_uncon",         &nr_track_asso_out_uncon_,              "nr_track_asso_out_uncon/I");
  t_out->Branch("nr_tracks_uncontained",         &nr_tracks_uncontained_,              "nr_tracks_uncontained/I");
  t_out->Branch("track_uncontained",         &track_uncontained_,              "track_uncontained/I");
  
  t_out->Branch("nr_crthit_top",         &nr_crthit_top_,              "nr_crthit_top/I");
  t_out->Branch("crthit_vertex_zcut",         &crthit_vertex_zcut_,              "crthit_vertex_zcut/I");
  t_out->Branch("crthit_vertex_z",         &crthit_vertex_z_,              "crthit_vertex_z/D");
  t_out->Branch("crthit_vertex_z_adc",         &crthit_vertex_z_adc,              "crthit_vertex_z_adc/D");
  
  t_out->Branch("has_matched_muon",         &has_matched_muon_,              "has_matched_muon/I");
  t_out->Branch("has_matched_muonHigh",         &has_matched_muonHigh_,              "has_matched_muonHigh/I");
  
  t_out->Branch("top_y_layer",         &top_y_layer,              "top_y_layer/I");
  t_out->Branch("y_box",         &y_box,              "y_box/I");
  t_out->Branch("z_dead",         &z_dead,              "z_dead/I");
  t_out->Branch("front_strip",         &front_strip,              "front_strip/I");
  t_out->Branch("z_dead_y",         &z_dead_y,              "z_dead_y/I");
  t_out->Branch("front_strip_y",         &front_strip_y,              "front_strip_y/I");
  t_out->Branch("front_z_layer",         &front_z_layer,              "front_z_layer/I");
  
  t_out->Branch("lowest_Track_point",         &lowest_Track_point,              "lowest_Track_point/I");
  
  t_out->Branch("small_wire",         &small_wire,              "small_wire/I");
  t_out->Branch("wire_plane",         &wire_plane,              "wire_plane/I");
  
  t_out->Branch("PID_muprotdiff",         &PID_muprotdiff,              "PID_muprotdiff/D");
  t_out->Branch("key_muprotdiff",         &key_muprotdiff,              "key_muprotdiff/I");
  t_out->Branch("PID_muprotratio",         &PID_muprotratio,              "PID_muprotratio/D");
  t_out->Branch("key_muprotratio",         &key_muprotratio,              "key_muprotratio/I");
  t_out->Branch("PID_mupiondiff",         &PID_mupiondiff,              "PID_mupiondiff/D");
  t_out->Branch("key_mupiondiff",         &key_mupiondiff,              "key_mupiondiff/I");
  t_out->Branch("PID_mupionratio",         &PID_mupionratio,              "PID_mupionratio/D");
  t_out->Branch("key_mupionratio",         &key_mupionratio,              "key_mupionratio/I");
  t_out->Branch("longest_track",         &longest_track,              "longest_track/D");
  t_out->Branch("key_longest_track",         &key_longest_track,              "key_longest_track/I");
  t_out->Branch("num_mc_muon",         &num_mc_muon,              "num_mc_muon/I");
  int key_muon = 0;
  t_out->Branch("key_muon",         &key_muon,              "key_muon/I");
  
  
  if(!is_data_){
    t_out->Branch("MCNu_Interaction",&MCNu_Interaction,"MCNu_Interaction/I");
    t_out->Branch("MCNu_CCNC",           &MCNu_CCNC,           "MCNu_CCNC/I");
    t_out->Branch("MCNu_PDG",            &MCNu_PDG,            "MCNu_PDG/I");
    t_out->Branch("MCNu_Energy",         &MCNu_Energy,         "MCNu_Energy/F");
    t_out->Branch("MCNu_leptonPx",       &MCNu_leptonPx,       "MCNu_leptonPx/F");
    t_out->Branch("MCNu_leptonPy",       &MCNu_leptonPy,       "MCNu_leptonPy/F");
    t_out->Branch("MCNu_leptonPz",       &MCNu_leptonPz,       "MCNu_leptonPz/F");
    t_out->Branch("MCNu_LeptonEnergy",   &MCNu_LeptonEnergy,   "MCNu_LeptonEnergy/F");
    t_out->Branch("MCNu_Px",             &MCNu_Px,             "MCNu_Px/F");
    t_out->Branch("MCNu_Py",             &MCNu_Py,             "MCNu_Py/F");
    t_out->Branch("MCNu_Pz",             &MCNu_Pz,             "MCNu_Pz/F");
    t_out->Branch("MCNu_leptonTheta",    &MCNu_leptonTheta,    "MCNu_leptonTheta/F");
    //t_out->Branch("MCNu_leptonPhi",      &MCNu_leptonPhi,      "MCNu_leptonPhi/F");
    t_out->Branch("MCNu_time",           &MCNu_time,           "MCNu_time/F");
    t_out->Branch("MCNu_Vx",             &MCNu_Vx,             "MCNu_Vx/F");
    t_out->Branch("MCNu_Vy",             &MCNu_Vy,             "MCNu_Vy/F");
    t_out->Branch("MCNu_Vz",             &MCNu_Vz,             "MCNu_Vz/F");
    t_out->Branch("MCNu_VxSce",          &MCNu_VxSce,          "MCNu_VxSce/F");
    t_out->Branch("MCNu_VySce",          &MCNu_VySce,          "MCNu_VySce/F");
    t_out->Branch("MCNu_VzSce",          &MCNu_VzSce,          "MCNu_VzSce/F");
    t_out->Branch("MCNu_Vx_sce",          &MCNu_Vx_sce_,          "MCNu_Vx_sce/F");
    t_out->Branch("MCNu_Vy_sce",          &MCNu_Vy_sce_,          "MCNu_Vy_sce/F");
    t_out->Branch("MCNu_Vz_sce",          &MCNu_Vz_sce_,          "MCNu_Vz_sce/F");
    t_out->Branch("MCNu_vertexDistance",    &MCNu_vertexDistance,    "MCNu_vertexDistance/F");
    
    t_out->Branch("MCle_key",               &MCle_key_out,               "MCle_key/I");
    t_out->Branch("MCle_pfp",               &MCle_pfp_out,               "MCle_pfp/I");
    t_out->Branch("MCle_PDG",               &MCle_PDG_out,               "MCle_PDG/I");
    //t_out->Branch("MCle_purity",            &MCle_purity_out,            "MCle_purity/F");
    t_out->Branch("MCle_Energy",            &MCle_Energy_out,            "MCle_Energy/F");
    t_out->Branch("MCle_Px",                &MCle_Px_out,                "MCle_Px/F");
    t_out->Branch("MCle_Py",                &MCle_Py_out,                "MCle_Py/F");
    t_out->Branch("MCle_Pz",                &MCle_Pz_out,                "MCle_Pz/F");
    t_out->Branch("MCle_Vx",                &MCle_Vx_out,                "MCle_Vx/F");
    t_out->Branch("MCle_Vy",                &MCle_Vy_out,                "MCle_Vy/F");
    t_out->Branch("MCle_Vz",                &MCle_Vz_out,                "MCle_Vz/F");
    t_out->Branch("MCle_Endx",                &MCle_Endx_out,                "MCle_Endx_/F");
    t_out->Branch("MCle_Endx",                &MCle_Endy_out,                "MCle_Endy/F");
    t_out->Branch("MCle_Endx",                &MCle_Endz_out,                "MCle_Endz/F");
    t_out->Branch("MCle_length",            &MCle_length_out,            "MCle_length/F");
    t_out->Branch("MCle_VxSce",             &MCle_VxSce_out,             "MCle_VxSce/F");
    t_out->Branch("MCle_VySce",             &MCle_VySce_out,             "MCle_VySce/F");
    t_out->Branch("MCle_VzSce",             &MCle_VzSce_out,             "MCle_VzSce/F");
    t_out->Branch("MCle_Theta",             &MCle_Theta_out,             "MCle_Theta/F");
    t_out->Branch("MCle_Phi",             &MCle_Phi_out,             "MCle_Phi_/F");
    /*
    t_out->Branch("MCTrackPurity",             &MCTrackPurity_out,             "MCTrackPurity/D");
    t_out->Branch("MCTrackPDG",             &MCTrackPDG_out,             "MCTrackPDG/D");
    t_out->Branch("MCTrackEnergy",             &MCTrackEnergy_out,             "MCTrackEnergy/D");
    t_out->Branch("MCTrackMomentum",             &MCTrackMomentum_out,             "MCTrackMomentum/D");
    t_out->Branch("MCTrackTheta",             &MCTrackTheta_out,             "MCTrackTheta/D");
    t_out->Branch("MCTrackPhi",             &MCTrackPhi_out,             "MCTrackPhi/D");
    t_out->Branch("MCTrackLength",             &MCTrackLength_out,             "MCTrackLength/D");
    
    t_out->Branch("MCTrackStart_x",             &MCTrackStart_x_out,             "MCTrackStart_x/D");
    t_out->Branch("MCTrackStart_y",             &MCTrackStart_y_out,             "MCTrackStart_y/D");
    t_out->Branch("MCTrackStart_z",             &MCTrackStart_z_out,             "MCTrackStart_z/D");
    t_out->Branch("MCTrackEnd_x",             &MCTrackEnd_x_out,             "MCTrackEnd_x/D");
    t_out->Branch("MCTrackEnd_y",             &MCTrackEnd_y_out,             "MCTrackEnd_y/D");
    t_out->Branch("MCTrackEnd_z",             &MCTrackEnd_z_out,             "MCTrackEnd_z/D");
    */
    t_out->Branch("MCTrackPurity",             &MCTrackPurity_);
    t_out->Branch("MCTrackPDG",             &MCTrackPDG_);
    t_out->Branch("MCTrackEnergy",             &MCTrackEnergy_);
    t_out->Branch("MCTrackMomentum",             &MCTrackMomentum_);
    t_out->Branch("MCTrackTheta",             &MCTrackTheta_);
    t_out->Branch("MCTrackPhi",             &MCTrackPhi_);
    t_out->Branch("MCTrackLength",             &MCTrackLength_);
    
    t_out->Branch("MCTrackStart_x",             &MCTrackStart_x_);
    t_out->Branch("MCTrackStart_y",             &MCTrackStart_y_);
    t_out->Branch("MCTrackStart_z",             &MCTrackStart_z_);
    t_out->Branch("MCTrackEnd_x",             &MCTrackEnd_x_);
    t_out->Branch("MCTrackEnd_y",             &MCTrackEnd_y_);
    t_out->Branch("MCTrackEnd_z",             &MCTrackEnd_z_);
    
    t_out->Branch("Genie_Q2",&Genie_Q2,"Genie_Q2/D");
    t_out->Branch("Genie_q2",&Genie_q2,"Genie_q2/D");
    t_out->Branch("Genie_W",&Genie_W,"Genie_W/D");
    t_out->Branch("Genie_T",&Genie_T,"Genie_T/D");
    t_out->Branch("Genie_X",&Genie_X,"Genie_X/D");
    t_out->Branch("Genie_Y",&Genie_Y,"Genie_Y/D");
    t_out->Branch("Genie_nNeutron_preFSI",&Genie_nNeutron_preFSI,"Genie_nNeutron_preFSI/I");
    t_out->Branch("Genie_nProton_preFSI",&Genie_nProton_preFSI,"Genie_nProton_preFSI/I");
    t_out->Branch("Genie_nPi0_preFSI",&Genie_nPi0_preFSI,"Genie_nPi0_preFSI/I");
    t_out->Branch("Genie_nPiPlus_preFSI",&Genie_nPiPlus_preFSI,"Genie_nPiPlus_preFSI/I");
    t_out->Branch("Genie_nPiMinus_preFSI",&Genie_nPiMinus_preFSI,"Genie_nPiMinus_preFSI/I");

  }
  
  
  long int tot_ev=t_in->GetEntries()-1;
  
  long counter_old=0;
  int wire_plane1 = 0, wire_plane2 = 0;
  //double min_muprotdiff = 999;
  //double min_muprotration = 999;
  //double max_tracklength = -999;
  
  
  for(long i = 0; i< tot_ev; i++){
    t_in->GetEntry(i);
    if(i%1000 == 0){
      printf("At entry: %ld from %ld  (%ld\%) \n",i, tot_ev, i*100/tot_ev);
    }
    PID_muprotdiff = 999;
    PID_muprotratio = 999;
    PID_mupiondiff = 999;
    PID_mupionratio = 999;
    longest_track = -999;
    num_mc_muon = 0;
    TunedCentralValue_Genie_out = 1.0;
    if(with_genie != -1) TunedCentralValue_Genie_out = TunedCentralValue_Genie_->at(0);
    
    for(unsigned int j = 0; j<track_key_->size(); j++){
      if(TrackPID_chimuon_->at(j) != -999 && TrackPID_chiproton_->at(j) != 0 && TrackLength_->at(j)>8.0){
        if(PID_muprotdiff > (TrackPID_chimuon_->at(j) - TrackPID_chiproton_->at(j))){
          PID_muprotdiff = (TrackPID_chimuon_->at(j) - TrackPID_chiproton_->at(j));
          key_muprotdiff = track_key_->at(j);
        }
        if(PID_muprotratio > (TrackPID_chimuon_->at(j)/TrackPID_chiproton_->at(j))){
          PID_muprotratio = (TrackPID_chimuon_->at(j)/TrackPID_chiproton_->at(j));
          key_muprotratio = track_key_->at(j);
        }
        if(PID_mupiondiff > (TrackPID_chimuon_->at(j) - TrackPID_chipion_->at(j))){
          PID_mupiondiff = (TrackPID_chimuon_->at(j) - TrackPID_chipion_->at(j));
          key_mupiondiff = track_key_->at(j);
        }
        if(PID_mupionratio > (TrackPID_chimuon_->at(j)/TrackPID_chipion_->at(j))){
          PID_mupionratio = (TrackPID_chimuon_->at(j)/TrackPID_chipion_->at(j));
          key_mupionratio = track_key_->at(j);
        }
        if(longest_track < (TrackLength_->at(j))){
          longest_track = TrackLength_->at(j);
          key_longest_track = track_key_->at(j);
        }
        if(!is_data_){
          if(MCTrackPDG_->at(j)==13){
            num_mc_muon++;
          }
        }
      }
    }
    std::vector<int> pass_proton_cut;
    std::vector<int> pass_pion_cut;
    //get muon candidate with proton, pion procedure:
    for(unsigned int j = 0; j<track_key_->size(); j++){
      if(track_key_->size()==1){
        key_muon = track_key_->at(j);
        //break;
      }
      if(TrackPID_chimuon_->at(j)!=-1 && (TrackPID_chimuon_->at(j)/TrackPID_chiproton_->at(j))<0.168 ){
        pass_proton_cut.push_back(j);
        if(TrackPID_chimuon_->at(j)!=-1 && (TrackPID_chimuon_->at(j)/TrackPID_chipion_->at(j))<1.06 ){
          pass_pion_cut.push_back(j);
        }
      }
      //if(TrackPID_chimuon_->at(j)!=-1 && (TrackPID_chimuon_->at(j)/TrackPID_chiproton_->at(j))<0.16 ){
      //  pass_proton_cut.push_back(j);
      //  if(TrackPID_chimuon_->at(j)!=-1 && (TrackPID_chimuon_->at(j)/TrackPID_chipion_->at(j))<0.75 ){
      //    pass_pion_cut.push_back(j);
      //  }
      //}
    }
    if(key_muon==-1){
      if(pass_proton_cut.size()<1) key_muon = key_longest_track;
      if(pass_proton_cut.size()==1) key_muon = track_key_->at(pass_proton_cut.at(0));
      else{
        if(pass_pion_cut.size()==1) key_muon = track_key_->at(pass_pion_cut.at(0));
        if(pass_pion_cut.size()<1){ //longest of pass proton
          double longest = -1;
          for(unsigned int k = 0; k<pass_proton_cut.size(); k++){
            if(longest<TrackLength_->at(pass_proton_cut.at(k))){
              longest=TrackLength_->at(pass_proton_cut.at(k));
              key_muon = track_key_->at(pass_proton_cut.at(k));
            }
          }
        }
        if(pass_pion_cut.size()>1){//longest of pion pass
          double longest = -1;
          for(unsigned int k = 0; k<pass_pion_cut.size(); k++){
            if(longest<TrackLength_->at(pass_pion_cut.at(k))){
              longest=TrackLength_->at(pass_pion_cut.at(k));
              key_muon = track_key_->at(pass_pion_cut.at(k));
            }
          }
          
        }
        //for(unsigned int j = 0; j<pass_proton_cut->size(); j++){
          
        
      //}
     //for(unsigned int j = 0; j<track_key_->size(); j++){
       
     }
    }
    
    
    //printf("At entry: %ld\n",i);
    for(unsigned int j = 0; j<track_key_->size(); j++){
      //printf("Entry: %ld - %d, muoncandidate: %d, track_key: %d \n", i, j, muon_candidate_key, track_key_->at(j));
      if(track_key_->at(j) == muon_candidate_key){
        //track_key_out_ = track_key_->at(j);
        printf("Fill moun candidate, %ld, %d\n",i,j);
  
        track_key_out = track_key_->at(j);
        TrackPID_out = TrackPID_->at(j);
        
        VtxDistance_out = VtxDistance_->at(j);
        VtxDistance_sce_out = VtxDistance_sce_->at(j);
        Vx_out = Vx_->at(j);
        Vy_out = Vy_->at(j);
        Vz_out = Vz_->at(j);
        Vx_sce_out = Vx_sce_->at(j);
        Vy_sce_out = Vy_sce_->at(j);
        Vz_sce_out = Vz_sce_->at(j);
        TrackScore_out = TrackScore_->at(j);
        TrackLength_out = TrackLength_->at(j);
        TrackPfp_out = TrackPfp_->at(j);
        
        if(TrackPID_chikaon_->size() > j){
          TrackPID_chiproton_out = TrackPID_chiproton_->at(j);
          TrackPID_chimuon_out = TrackPID_chimuon_->at(j);
          TrackPID_chipion_out = TrackPID_chipion_->at(j);
          TrackPID_chikaon_out = TrackPID_chikaon_->at(j);
        }
        else{
          //printf("Missed CHi2!!!!!!!!!!!!!!!!!!!!!\n");
          TrackPID_chiproton_out = -1;
          TrackPID_chimuon_out = -1;
          TrackPID_chipion_out = -1;
          TrackPID_chikaon_out = -1;
        }
        TrackMomRange_p_out = TrackMomRange_p_->at(j);
        TrackMomRange_mu_out = TrackMomRange_mu_->at(j);
        TrackMomMCS_mom_out = TrackMomMCS_mom_->at(j);
        TrackMomMCS_err_out = TrackMomMCS_err_->at(j);
        TrackMomMCS_ll_out = TrackMomMCS_ll_->at(j);
        
        TrackStart_x_out = TrackStart_x_->at(j);
        TrackStart_y_out = TrackStart_y_->at(j);
        TrackStart_z_out = TrackStart_z_->at(j);
        TrackStart_x_sce_out = TrackStart_x_sce_->at(j);
        TrackStart_y_sce_out = TrackStart_y_sce_->at(j);
        TrackStart_z_sce_out = TrackStart_z_sce_->at(j);
        TrackEnd_x_out = TrackEnd_x_->at(j);
        TrackEnd_y_out = TrackEnd_y_->at(j);
        TrackEnd_z_out = TrackEnd_z_->at(j);
        TrackEnd_x_sce_out = TrackEnd_x_sce_->at(j);
        TrackEnd_y_sce_out = TrackEnd_y_sce_->at(j);
        TrackEnd_z_sce_out = TrackEnd_z_sce_->at(j);
        TrackDir_x_out = TrackDir_x_->at(j);
        TrackDir_y_out = TrackDir_y_->at(j);
        TrackDir_z_out = TrackDir_z_->at(j);

        TrackTheta_out = TrackTheta_->at(j);
        TrackPhi_out = TrackPhi_->at(j);
        TrackNHitsU_out = TrackNHitsU_->at(j);
        TrackNHitsV_out = TrackNHitsV_->at(j);
        TrackNHitsY_out = TrackNHitsY_->at(j);
        TrackCaloU_out = TrackCaloU_->at(j);
        
        TrackNHitsAll_out = TrackNHitsU_out+TrackNHitsV_out+TrackNHitsY_out;
        TrackCaloAll_out = 0;
        
        a_crthit_ts0out = a_crthit_ts0->at(j);
        a_crthit_ts1out = a_crthit_ts1->at(j);
        if(a_crt_plane_->size()!=0){
          a_crthit_x_out = a_crthit_x_->at(0);
          a_crthit_y_out = a_crthit_y_->at(0);
          a_crthit_z_out = a_crthit_z_->at(0);
          a_crthit_plane_out = a_crt_plane_->at(0);
        }
         a_adc_lengthout = a_adc_length->at(j);
         a_crt_adcout = a_crt_adc->at(j);
         a_t0_counterout = a_t0_counter->at(j);
        printf("CRT assigned: %lf - %lf\n",a_adc_lengthout,a_adc_length->at(j));
        crtt0_time_out = crtt0_time_->at(j);
         crtt0_trig_out = crtt0_trig_->at(j);
        crtt0_DCA_out = crtt0_DCA_->at(j);
         crtt0_plane_out = crtt0_plane_->at(j);
        printf("CRT T0: %lf - %lf\n",crtt0_time_out,crtt0_time_->at(j));
        if(TrackEnd_x_->at(j)>(254.8-5) || TrackEnd_x_->at(j)< (-1.55+5) || TrackEnd_y_->at(j)>(117.47-5) || TrackEnd_y_->at(j)<(-115.53+5) || TrackEnd_z_->at(j)>(1039.9-5) || TrackEnd_z_->at(j)<(0.1+5)){
          track_uncontained_++;
        }
        /*crthit_ts0_out = crthit_ts0_->at(0);
        crthit_ts1_out = crthit_ts1_->at(0);
         adc_length_out = adc_length_->at(0);
        crt_adc_out = crt_adc_->at(0);
         crtbeam_hit_nr_out = crtbeam_hit_nr_->at(0);
        */
        //for(int j = 0; j<a_crthit_ts0->size(); j++){
          //printf("Entry: %ld - %d, muoncandidate: %d, crthit_ts0_: %d \n", i, j, muon_candidate_key, crthit_ts0_->at(j));
          //if(track_key_->at(j) == muon_candidate_key){}

        //}
        /*
        crthit_ts0_out = crthit_ts0_->at(j);
        crthit_ts1_out = crthit_ts1_->at(j);
         adc_length_out = adc_length_->at(j);
        crt_adc_out = crt_adc_->at(j);
         crtbeam_hit_nr_out = crtbeam_hit_nr_->at(j);
        
        a_crthit_ts0out = a_crthit_ts0->at(j);
        a_crthit_ts1out = a_crthit_ts1->at(j);
         a_adc_lengthout = a_adc_length->at(j);
        a_crt_adcout = a_crt_adc->at(j);
         a_t0_counterout = a_t0_counter->at(j);

        crtt0_time_out = crtt0_time_->at(j);
         crtt0_trig_out = crtt0_trig_->at(j);
        crtt0_DCA_out = crtt0_DCA_->at(j);
         crtt0_plane_out = crtt0_plane_->at(j);*/
      }
      if(TrackEnd_x_->at(j)>(254.8-5) || TrackEnd_x_->at(j)< (-1.55+5) || TrackEnd_y_->at(j)>(117.47-5) || TrackEnd_y_->at(j)<(-115.53+5) || TrackEnd_z_->at(j)>(1039.9-5) || TrackEnd_z_->at(j)<(0.1+5)){
        nr_tracks_uncontained_++;
      }
      //printf("Check asso...\n");
      if(crtt0_time_->at(j)>3.15 && crtt0_time_->at(j)<4.75){
        nr_track_asso_++;
        //printf("End: %lf - %lf - %lf\n",TrackEnd_x_->at(j),TrackEnd_y_->at(j),TrackEnd_z_->at(j));
        if(TrackEnd_x_->at(j)>(254.8-5) || TrackEnd_x_->at(j)< (-1.55+5) || TrackEnd_y_->at(j)>(117.47-5) || TrackEnd_y_->at(j)<(-115.53+5) || TrackEnd_z_->at(j)>(1039.9-5) || TrackEnd_z_->at(j)<(0.1+5)){
          nr_track_asso_uncon_++;
        }
      }
      else if(crtt0_time_->at(j)!=-1){
        nr_track_asso_out_++;
        //printf("End: %lf - %lf - %lf\n",TrackEnd_x_->at(j),TrackEnd_y_->at(j),TrackEnd_z_->at(j));
        if(TrackEnd_x_->at(j)>(254.8-5) || TrackEnd_x_->at(j)< (-1.55+5) || TrackEnd_y_->at(j)>(117.47-5) || TrackEnd_y_->at(j)<(-115.53+5) || TrackEnd_z_->at(j)>(1039.9-5) || TrackEnd_z_->at(j)<(0.1+5)){
          nr_track_asso_out_uncon_++;
        }
      }
      
    }
    double maxShowerCaloU = 0;
    int shower_index = -1;
    NumShowers_corr = 0;
    for(unsigned int j = 0; j<ShowerPfp_->size(); j++){
      //printf("In Loop: i = %d - SHowerScore = %lf\n",j,ShowerScore_->at(j));
      if(ShowerScore_->at(j)<0.5){
        NumShowers_corr++;
        if(ShowerCaloU_->at(j)>maxShowerCaloU){
          maxShowerCaloU = ShowerCaloU_->at(j);
          shower_index = j; 
        }
      }
    }
    if(shower_index != -1){
      ShowerPfp_out = ShowerPfp_->at(shower_index);
      //ShowerDir_x_out = ShowerDir_x_->at(shower_index);
      //ShowerDir_y_out = ShowerDir_y_->at(shower_index);
      //ShowerDir_z_out = ShowerDir_z_->at(shower_index);
      ShowerNHitsU_out = ShowerNHitsU_->at(shower_index);
      ShowerNHitsV_out = ShowerNHitsV_->at(shower_index);
      ShowerNHitsY_out = ShowerNHitsY_->at(shower_index);
      ShowerCaloU_out = ShowerCaloU_->at(shower_index);
      ShowerCaloV_out = ShowerCaloV_->at(shower_index);
      ShowerCaloY_out = ShowerCaloY_->at(shower_index);
      ShowerCaloAll_out = ShowerCaloU_out+ShowerCaloV_out+ShowerCaloY_out;
      ShowerNHitsAll_out = ShowerNHitsU_out+ShowerNHitsV_out+ShowerNHitsY_out;
      ShowerScore_out = ShowerScore_->at(shower_index);
      //printf("i = %d - SHowerScore = %lf\n",shower_index,ShowerScore_out);
    }
    if(!is_data_){
      for(unsigned int j = 0; j<MCle_pfp_->size(); j++){
        if(MCle_pfp_->at(j) == muon_candidate_pfp){
          //fill all MCle_out variables
          MCle_key_out = MCle_key_->at(j);
          MCle_pfp_out = MCle_pfp_->at(j);
          MCle_PDG_out = MCle_PDG_->at(j);
          //MCle_purity_out = MCle_purity_->at(j);
          MCle_Energy_out = MCle_Energy_->at(j);
          MCle_Px_out = MCle_Px_->at(j);
          MCle_Py_out = MCle_Py_->at(j);
          MCle_Pz_out = MCle_Pz_->at(j);
          //MCle_Vx_out = MCle_Vx_->at(j);
          //MCle_Vy_out = MCle_Vy_->at(j);
          //MCle_Vz_out = MCle_Vz_->at(j);
          //MCle_length_out = MCle_length_->at(j);
          //MCle_VxSce_out = MCle_VxSce_->at(j);
          //MCle_VySce_out = MCle_VySce_->at(j);
          //MCle_VzSce_out = MCle_VzSce_->at(j);
          //MCle_Vx_sce_out = MCle_Vx_sce_->at(j);
          //MCle_Vy_sce_out = MCle_Vy_sce_->at(j);
          //MCle_Vz_sce_out = MCle_Vz_sce_->at(j);
        }
        if(MCle_PDG_->at(j) == 13){
          has_matched_muon_++;
          //if(MCle_purity_->at(j) > 0.5) has_matched_muonHigh_++;
        }
      }
      //printf("Check0\n");
      if(MCNu_PDG==14 && MCNu_CCNC==0){
        //printf("Check1 size: %d\n",MCle_PDG_->size());
        for(unsigned int j = 0; j<MCle_PDG_->size(); j++){
          //printf("Check2 PDG: %d\n",MCle_PDG_->at(j));
          if(MCle_PDG_->at(j) == 13){
            //printf("Check2\n");
            //fill all MCle_out variables
            //MCle_key_out = MCle_key_->at(j);
            //MCle_pfp_out = MCle_pfp_->at(j);
            MCle_PDG_out = MCle_PDG_->at(j);
            //MCle_purity_out = MCle_purity_->at(j);
            MCle_Energy_out = MCle_Energy_->at(j);
            MCle_Px_out = MCle_Px_->at(j);
            MCle_Py_out = MCle_Py_->at(j);
            MCle_Pz_out = MCle_Pz_->at(j);
            MCle_Vx_out = MCle_Vx_->at(j);
            MCle_Vy_out = MCle_Vy_->at(j);
            MCle_Vz_out = MCle_Vz_->at(j);
            MCle_length_out = MCle_length_->at(j);
            //MCle_VxSce_out = MCle_VxSce_->at(j);
            //MCle_VySce_out = MCle_VySce_->at(j);
            //MCle_VzSce_out = MCle_VzSce_->at(j);
            //MCle_Vx_sce_out = MCle_Vx_sce_->at(j);
            //MCle_Vy_sce_out = MCle_Vy_sce_->at(j);
            //MCle_Vz_sce_out = MCle_Vz_sce_->at(j);
          }
        }
      }
    }
    double max_crtz_diff = 9999;
    for(unsigned int j = 0; j<crthit_ts0_->size(); j++){
      //printf("Entry: %ld - %d, muoncandidate: %d, crthit_ts0_: %lf \n", i, j, muon_candidate_key, crthit_ts0_->at(j));
      if(crthit_ts0_->at(j)<4.75 && crthit_ts0_->at(j)>3.15){
        //printf("Found CRT hit in beam!\n");
        nr_crthit_beam_++;
        if(crt_adc_->at(j)>70){
          nr_crthit_beam_tres_++;
          if(crt_plane_->at(j) == 3) nr_crthit_top_++;
          if( crthit_z_->at(j)< crthit_z_out) crthit_z_out = crthit_z_->at(j);
          if( crthit_y_->at(j)> crthit_y_out) crthit_y_out = crthit_y_->at(j);
          if( (crthit_z_->at(j) - Nu_Vz_sce_)<0 ){
            crthit_vertex_zcut_++;
          }
          if(max_crtz_diff>(crthit_z_->at(j) - Nu_Vz_sce_)){
            crthit_vertex_z_ = (crthit_z_->at(j) - Nu_Vz_sce_);
            max_crtz_diff = crthit_vertex_z_;
            crthit_vertex_z_adc = adc_length_->at(j);
          }
        }
      }
      
    }
    front_z_layer = 2000;
    for(unsigned int j = 0; j<AllTrack_point_z_->size(); j++){
      //printf("In Loop: i = %d - SHowerScore = %lf\n",j,ShowerScore_->at(j));
      if(AllTrack_point_y_->at(j)>110){
        top_y_layer = 1;
      }
      if(AllTrack_point_z_->at(j)<front_z_layer && AllTrack_point_z_->at(j)>-100){
        front_z_layer = AllTrack_point_z_->at(j);
      }
      if(AllTrack_point_y_->at(j)>100 && (AllTrack_point_z_->at(j)>500 && AllTrack_point_z_->at(j)<570)){
        y_box = 1;
      }
      if(AllTrack_point_z_->at(j)>700 && AllTrack_point_z_->at(j)<740){
        z_dead = 1;
      }
      if( abs(170.0/300*AllTrack_point_z_->at(j)-111-AllTrack_point_y_->at(j))<5 ){
        front_strip = 1;
      }
      if((AllTrack_point_z_->at(j)>700 && AllTrack_point_z_->at(j)<740) && (AllTrack_point_y_->at(j)>50)){
        z_dead_y = 1;
      }
      if(( abs(170.0/300*AllTrack_point_z_->at(j)-111-AllTrack_point_y_->at(j))<5) && (AllTrack_point_y_->at(j)>50)){
        front_strip_y = 1;
      }
    }
    if(has_wire){
      for(unsigned int j = 0; j<Wire_id_->size(); j++){
        if( Wire_id_->at(j)<10){
          small_wire++;
          if(wire_plane>Wire_id_->at(j)) wire_plane=Wire_id_->at(j);
          //small_wire+=10*Wire_plane_->at(j);
          //if(Wire_plane_->at(j)==0){
          //  wire_plane = 1;
          //}
          //if(Wire_plane_->at(j)==1){
          //  wire_plane1 = 1;
          //}
          //if(Wire_plane_->at(j)==2){
          //  wire_plane2 = 1;
          //}
        }
      }
    }
    //wire_plane += wire_plane1*10 + wire_plane2*100;
    if( abs(TunedCentralValue_Genie_out)<10000 ) t_out->Fill();
    nr_crthit_beam_ = 0;
    nr_crthit_beam_tres_ = 0;
    nr_track_asso_ = 0;
    nr_track_asso_out_=0;
    nr_track_asso_out_uncon_ = 0;
    nr_track_asso_uncon_ = 0;
    nr_tracks_uncontained_ = 0;
    track_uncontained_ = 0;
    
    track_key_out = -1;
    TrackPID_out = -1;

    VtxDistance_out = -1;
    VtxDistance_sce_out = -1;
    Vx_out = -9999;
    Vy_out = -9999;
    Vz_out = -9999;
    Vx_sce_out = -9999;
    Vy_sce_out = -9999;
    Vz_sce_out = -9999;
    TrackScore_out = -1;
    TrackLength_out = -1;
    TrackPfp_out = -1;
    isShowerTrack_out = -1;
    TrackMomRange_p_out = -1;
    TrackMomRange_mu_out = -1;
    TrackMomMCS_mom_out = -1;
    TrackMomMCS_err_out = -1;
    TrackMomMCS_ll_out = -1;

    TrackStart_x_out = -9999;
    TrackStart_y_out = -9999;
    TrackStart_z_out = -9999;
    TrackStart_x_sce_out = -9999;
    TrackStart_y_sce_out = -9999;
    TrackStart_z_sce_out = -9999;
    TrackEnd_x_out = -9999;
    TrackEnd_y_out = -9999;
    TrackEnd_z_out = -9999;
    TrackEnd_x_sce_out = -9999;
    TrackEnd_y_sce_out = -9999;
    TrackEnd_z_sce_out = -9999;
    TrackDir_x_out = -9999;
    TrackDir_y_out = -9999;
    TrackDir_z_out = -9999;

    TrackTheta_out = -9999;
    TrackPhi_out = -9999;
    
    TrackNHitsU_out = 0;
    TrackNHitsV_out = 0;
    TrackNHitsY_out = 0;
    TrackCaloU_out = 0;
    TrackCaloAll_out = 0;
    TrackNHitsAll_out = 0;
    
    ShowerPfp_out = 0;
    ShowerDir_x_out = 0;
    ShowerDir_y_out = 0;
    ShowerDir_z_out = 0;
    ShowerNHitsU_out = 0;
    ShowerNHitsV_out = 0;
    ShowerNHitsY_out = 0;
    ShowerCaloU_out = 0;
    ShowerCaloV_out = 0;
    ShowerCaloY_out = 0;
    ShowerCaloAll_out = 0;
    ShowerNHitsAll_out = 0;
    NumShowers_corr = 0;
    ShowerScore_out = 0;
    //TrackNHitsU_out = TrackNHitsU_->at(j);
    //TrackNHitsV_out = TrackNHitsV_->at(j);
    //TrackNHitsY_out = TrackNHitsY_->at(j);
    //TrackCaloU_out = TrackCaloU_->at(j);

    a_crthit_ts0out = -9999;
    a_crthit_ts1out = -9999;
    a_adc_lengthout = 0;
    a_crt_adcout = 0;
    a_t0_counterout = 0;
             
    crthit_z_out = 9999;
    crthit_y_out = -9999;

    crtt0_time_out = -1;
    crtt0_trig_out = -1;
    crtt0_DCA_out = -1;
    crtt0_plane_out = -1;
    
    nr_crthit_top_ = 0;
    crthit_vertex_zcut_ = 0;
    crthit_vertex_z_ = 9999;
    crthit_vertex_z_adc = -1;
    
    has_matched_muon_ = 0;
    has_matched_muonHigh_ = 0;
    
    key_muon = -1;
    
    top_y_layer = 0;
    y_box = 0;
    z_dead = 0;
    front_strip = 0;
    z_dead_y = 0;
    front_strip_y = 0;
    front_z_layer = 0;
    
    small_wire = 0;
    wire_plane = 100;
    wire_plane1 = 0;
    wire_plane2 = 0;
    
    MCle_key_out = -1;
    MCle_pfp_out =-1;
    MCle_PDG_out = -1;
    //MCle_purity_out = -1;
    MCle_Energy_out =-1;
    MCle_Px_out = -1;
    MCle_Py_out = -1;
    MCle_Pz_out = -1;
    //MCle_Vx_out = -9999;
    //MCle_Vy_out = -9999;
    //MCle_Vz_out = -9999;
    //MCle_length_out =  -9999;
    //MCle_VxSce_out =  -9999;
    //MCle_VySce_out = -9999;
    //MCle_VzSce_out =  -9999;
    //MCle_Vx_sce_out = -9999;
    //MCle_Vy_sce_out = -9999;
    //MCle_Vz_sce_out =  -9999;
    
    if(i> event_size) break;
  }
  
  printf("\n");
  t_out->Write("",TObject::kOverwrite);
  
  f_in->Close();
  f_out.Close();
  printf("finished\n");
  //fclose(text);
}