////////////////////////////////////////////////////////////////////////
// Class:       NuCC
// Plugin Type: analyzer (art v2_11_03)
// File:        NuCC_module.cc
//
// Generated by Wouter Van de pontseele using cetskelgen
////////////////////////////////////////////////////////////////////////

#ifndef NUCC_H
#define NUCC_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcore/Geometry/Geometry.h"


#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "helpers/PandoraInterfaceHelper.h"
#include "helpers/EnergyHelper.h"
#include "helpers/TrackHelper.h"

#include "TTree.h"

class NuCC;

class NuCC : public art::EDAnalyzer
{
public:
    explicit NuCC(fhicl::ParameterSet const &p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    NuCC(NuCC const &) = delete;
    NuCC(NuCC &&) = delete;
    NuCC &operator=(NuCC const &) = delete;
    NuCC &operator=(NuCC &&) = delete;

    // Required functions.
    void analyze(art::Event const &e) override;
    void reconfigure(fhicl::ParameterSet const &p);
    void clearEvent();
    void clearDaughter();

    /**
     *  @brief  Collect and fill the reco-truth matching information.
     *
     *  @param  e Art event
     */
    void FillReconTruthMatching(art::Event const &e);

    /**
     *  @brief  Collect and fill the reconstructed information.
     *
     *  @param  e Art event
     */
    void FillReconstructed(art::Event const &e);

    /**
     *  @brief  Collect and fill the MC based neutrino information.
     *
     *  @param  e Art event
     */
    void FillTrueNu(art::Event const &e);

    /**
     *  @brief  Collect and fill the MC neutrino daughter information.
     *
     *  @param  e Art event
     */
    void FillTrueNuDaughters(art::Event const &e);

    /**
     *  @brief  Fill the tree for every daughter of the neutrino candidate.
     *
     *  @param  pfparticle ptr The Pfparticle corresponding to the daughter.
     *  @return 1 if succesful, 0 if failure.
     */
    bool FillDaughters(const art::Ptr<recob::PFParticle> &pfp,
                       const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                       const art::FindManyP<anab::ParticleID> &trackPIDAssn);

    /**
     *  @brief  Fill the information about the matching of the reconstructed daughter
     *
     *  @param  pfparticle ptr The Pfparticle corresponding to the daughter.
     *  @return 1 if succesfully matched, 0 if not matched.
     */
    bool MatchDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp);

    /**
     *  @brief  Tag a daughter as a muon candidate
     *
     *  @return 1 if succesfully, 0 if not.
     */
    bool IsMuonCandidate();

    /**
     *  @brief  Tag event as NuMuCC, if so, create association
     *
     *  @return 1 if succesfully, 0 if not.
     */
    bool IsNuMuCC();

    /**
     *  @brief  Returns if point is inside a fiducial volume
     *
     *  @param fiducial volume tolerance: -x,+x,-y,+y,-z,+z
     *  @return 1 if succesfully, 0 if not.
     */
    bool IsContained(float x, float y, float z, const std::vector<float> &borders) const;


    void endSubRun(const art::SubRun &subrun);

private:
    // Fields needed for the analyser
    std::string m_pfp_producer;
    std::string m_hitfinder_producer;
    std::string m_geant_producer;
    std::string m_hit_mcp_producer;

    float m_vtx_fid_x_start;
    float m_vtx_fid_y_start;
    float m_vtx_fid_z_start;
    float m_vtx_fid_x_end;
    float m_vtx_fid_y_end;
    float m_vtx_fid_z_end;
    float m_pfp_start_border;

    bool m_hasMCNeutrino;
    bool m_isData;

    float m_muon_cut_trackscore;
    float m_muon_cut_vtxdistance;
    float m_muon_cut_protonchi2;
    float m_muon_cut_muonchi2;
    float m_muon_cut_chiratio;
    float m_muon_cut_length;

    float m_event_cut_flashchi2;
    float m_event_cut_nuscore_soft;
    float m_event_cut_nuscore_hard;
    float m_event_cut_flashchi2_ratio;
    float m_event_cut_length;

    PandoraInterfaceHelper pandoraInterfaceHelper;
    EnergyHelper energyHelper;
    TrackHelper trackHelper;

    // Store the pfps that qualify as muon candidates
    std::vector<art::Ptr<recob::PFParticle>> m_muon_candidates;

    // LAr Pandora Helper fields
    lar_pandora::LArPandoraHelper larpandora;
    lar_pandora::PFParticleVector pfparticles;
    lar_pandora::PFParticleVector pfneutrinos;
    lar_pandora::PFParticleVector pfdaughters;
    lar_pandora::ShowerVector pfshowers;
    lar_pandora::TrackVector pftracks;
    lar_pandora::PFParticleMap particleMap;
    lar_pandora::PFParticlesToMetadata particlesToMetadata;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::PFParticlesToClusters particlesToClusters;
    lar_pandora::PFParticlesToShowers particlesToShowers;
    lar_pandora::PFParticlesToTracks particlesToTracks;
    lar_pandora::PFParticlesToSpacePoints particlesToSpacePoints;
    lar_pandora::ClustersToHits clustersToHits;
    lar_pandora::HitsToSpacePoints hitsToSpacePoints;
    lar_pandora::SpacePointsToHits spacePointsToHits;
    // Used for reco truth matching
    lar_pandora::PFParticlesToMCParticles matchedParticles;
    std::set<art::Ptr<simb::MCParticle>> matchedMCParticles;
    std::map<art::Ptr<recob::PFParticle>, float> matchedHitFractions;

    //// Tree for every event
    TTree *fEventTree;
    uint fRun, fSubrun, fEvent;
    UInt_t fTimeHigh, fTimeLow;
    uint fNumPfp;
    // MC neutrino info
    uint fNumNu; // number of MC neutrinos in event, only one gets saved!
    int fTrueNu_InteractionType;
    int fTrueNu_CCNC;
    int fTrueNu_PDG;
    float fTrueNu_Energy;
    float fTrueNu_LeptonPx, fTrueNu_LeptonPy, fTrueNu_LeptonPz;
    float fTrueNu_LeptonEnergy;
    float fTrueNu_Px, fTrueNu_Py, fTrueNu_Pz;
    float fTrueNu_LeptonTheta;
    float fTrueNu_Time; // time of the true neutrino interaction
    float fTrueNu_Vx, fTrueNu_Vy, fTrueNu_Vz;
    float fTrueNu_VxSce, fTrueNu_VySce, fTrueNu_VzSce;
    float fTrueNu_VtxDistance;
    // MC neutrino daughter
    std::vector<int> fTrueNu_DaughterPDG;
    std::vector<float> fTrueNu_DaughterE;
    std::vector<bool> fTrueNu_DaughterMatched;
    // Reco candidate info
    int fNu_PDG;
    float fNu_Score;
    uint fNu_SliceIndex;
    float fNu_Vx, fNu_Vy, fNu_Vz;
    bool fNu_Contained;
    bool fDaughtersStartContained;
    uint fNu_NhitsU, fNu_NhitsV, fNu_NhitsY;
    float fNu_CaloU, fNu_CaloV, fNu_CaloY;
    uint fNu_NSpacepoints;
    float fNu_FlashChi2;
    float fBestObviousCosmic_FlashChi2;
    uint fNumPrimaryDaughters;
    uint fNumDaughters;
    uint fNumMatchedDaughters;
    uint fNumShowers;
    uint fNumTracks;
    bool fDaughtersStored; // if all the neutrino daughters were stored correctly
    bool fCosmicMatched;   // One of the daughters is not matched to neutrino origin
    bool fIsNuMuCC;

    //// Tree for every daughter
    TTree *fNueDaughtersTree;
    // Reco candidate info
    uint fGeneration;
    float fTrackScore;
    bool fIsTrack;  // if the pfparticle was associated with a track
    bool fIsShower; // if the pfparticle was associated with a shower
    bool fHasShowerDaughter;
    bool fIsTrackDaughter;
    float fVx, fVy, fVz;
    float fStartContained;
    float fVtxDistance;
    uint fNhitsU, fNhitsV, fNhitsY;
    float fCaloU, fCaloV, fCaloY;
    uint fNSpacepoints;
    // Track info
    float fTrackLength;
    float fTrackDirX;
    float fTrackDirY;
    float fTrackDirZ;
    float fTrackEndX;
    float fTrackEndY;
    float fTrackEndZ;
    float fTrackRange_mom_p;
    float fTrackRange_mom_mu;
    float fTrackMCS_mom;
    float fTrackMCS_err;
    float fTrackMCS_ll;
    float fTrackPID_chiproton;
    float fTrackPID_chimuon;
    bool fIsMuonCandidate;
    // Shower info
    float fShowerLength;
    float fShowerOpenAngle;
    float fShowerDirX;
    float fShowerDirY;
    float fShowerDirZ;
    float fDedxU;
    float fDedxV;
    float fDedxY;
    uint fDedxHitsU;
    uint fDedxHitsV;
    uint fDedxHitsY;
    float fDedxPitchU;
    float fDedxPitchV;
    float fDedxPitchY;
    // Matched MCParticle info
    bool fMatchedNeutrino;
    float fTrueHitFraction;
    int fTruePDG;
    float fTrueEnergy;
    float fTrueVx, fTrueVy, fTrueVz;
    float fTrueLength;
    float fTrueVxSce, fTrueVySce, fTrueVzSce;
    float fTruePx, fTruePy, fTruePz;

    //// Tree for the POT subrun info
    TTree *fSubrunTree;
    uint m_run, m_subrun;
    float m_pot;
};

void NuCC::reconfigure(fhicl::ParameterSet const &p)
{
    m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraConsolidated");
    m_hitfinder_producer = p.get<std::string>("hitfinder_producer", "gaushit");
    m_geant_producer = p.get<std::string>("geant_producer", "largeant");
    m_hit_mcp_producer = p.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");

    m_vtx_fid_x_start = p.get<float>("vtx_fid_x_start", 10);
    m_vtx_fid_y_start = p.get<float>("vtx_fid_y_start", 10);
    m_vtx_fid_z_start = p.get<float>("vtx_fid_z_start", 10);
    m_vtx_fid_x_end = p.get<float>("vtx_fid_x_end", 10);
    m_vtx_fid_y_end = p.get<float>("vtx_fid_y_end", 10);
    m_vtx_fid_z_end = p.get<float>("vtx_fid_z_end", 50);
    m_pfp_start_border = p.get<float>("pfp_start_border", 10);

    m_isData = p.get<bool>("is_data", false);
    m_hasMCNeutrino = p.get<bool>("has_MC_neutrino", false);

    m_muon_cut_trackscore  = p.get<float>("muon_cut_trackscore", 0.8);
    m_muon_cut_vtxdistance  = p.get<float>("muon_cut_vtxdistance", 4.0);
    m_muon_cut_protonchi2 = p.get<float>("muon_cut_protonchi2", 60);
    m_muon_cut_muonchi2 = p.get<float>("muon_cut_muonchi2", 30);
    m_muon_cut_chiratio = p.get<float>("muon_cut_chiratio", 7);
    m_muon_cut_length = p.get<float>("muon_cut_length", 5);

    m_event_cut_flashchi2 = p.get<float>("event_cut_flashchi2", 10);
    m_event_cut_nuscore_soft = p.get<float>("event_cut_nuscore_soft", 0.25);
    m_event_cut_nuscore_hard = p.get<float>("event_cut_nuscore_hard", 0.06);
    m_event_cut_flashchi2_ratio = p.get<float>("event_cut_flashchi2_ratio", 5);
    m_event_cut_length = p.get<float>("event_cut_length", 20);

    energyHelper.reconfigure(p);
}

NuCC::NuCC(fhicl::ParameterSet const &p)
    : EDAnalyzer(p)
{
    art::ServiceHandle<art::TFileService> tfs;
    this->reconfigure(p);

    //// Check if things are set up properly:
    std::cout << std::endl;
    std::cout << "[NuCC constructor] Checking set-up" << std::endl;
    std::cout << "[NuCC constructor] pfp_producer: " << m_pfp_producer << std::endl;
    std::cout << "[NuCC constructor] is_data: " << m_isData << std::endl;
    std::cout << "[NuCC constructor] has_MC_neutrino: " << m_hasMCNeutrino << std::endl;

    //// Tree for every event
    fEventTree = tfs->make<TTree>("Event", "Event Tree");
    fEventTree->Branch("event", &fEvent, "event/i");
    fEventTree->Branch("run", &fRun, "run/i");
    fEventTree->Branch("subrun", &fSubrun, "subrun/i");
    fEventTree->Branch("evt_time_sec", &fTimeHigh, "evt_time_sec/i");
    fEventTree->Branch("evt_time_nsec", &fTimeLow, "evt_time_nsec/i");
    fEventTree->Branch("numpfp", &fNumPfp, "numpfp/i");
    fEventTree->Branch("hitsU", &fNu_NhitsU, "hitsU/i");
    fEventTree->Branch("hitsV", &fNu_NhitsV, "hitsV/i");
    fEventTree->Branch("hitsY", &fNu_NhitsY, "hitsY/i");
    fEventTree->Branch("caloU", &fNu_CaloU, "caloU/F");
    fEventTree->Branch("caloV", &fNu_CaloV, "caloV/F");
    fEventTree->Branch("caloY", &fNu_CaloY, "caloY/F");
    fEventTree->Branch("hitsSps", &fNu_NSpacepoints, "hitsSps/i");
    fEventTree->Branch("num_primary_daughters", &fNumPrimaryDaughters, "num_primary_daughters/i");
    fEventTree->Branch("num_daughters", &fNumDaughters, "num_daughters/i");
    fEventTree->Branch("num_showers", &fNumShowers, "num_showers/i");
    fEventTree->Branch("num_tracks", &fNumTracks, "num_tracks/i");
    fEventTree->Branch("daughters_stored", &fDaughtersStored, "daughters_stored/O");
    fEventTree->Branch("nu_vx", &fNu_Vx, "nu_vx/F");
    fEventTree->Branch("nu_vy", &fNu_Vy, "nu_vy/F");
    fEventTree->Branch("nu_vz", &fNu_Vz, "nu_vz/F");
    fEventTree->Branch("nu_contained", &fNu_Contained, "nu_contained/O");
    fEventTree->Branch("nu_pdg", &fNu_PDG, "nu_pdg/I");
    fEventTree->Branch("nu_flash_chi2", &fNu_FlashChi2, "nu_flash_chi2/F");
    fEventTree->Branch("obvious_cosmic_chi2", &fBestObviousCosmic_FlashChi2, "obvious_cosmic_chi2/F");
    fEventTree->Branch("nu_mu_cc_selected", &fIsNuMuCC, "nu_mu_cc_selected/0");


    if (m_hasMCNeutrino && !m_isData)
    {
        fEventTree->Branch("num_neutrinos", &fNumNu, "num_neutrinos/i");
        fEventTree->Branch("mc_nu_vx", &fTrueNu_Vx, "mc_nu_vx/F");
        fEventTree->Branch("mc_nu_vy", &fTrueNu_Vy, "mc_nu_vy/F");
        fEventTree->Branch("mc_nu_vz", &fTrueNu_Vz, "mc_nu_vz/F");
        fEventTree->Branch("mc_nu_vx_sce", &fTrueNu_VxSce, "mc_nu_vx_sce/F");
        fEventTree->Branch("mc_nu_vy_sce", &fTrueNu_VySce, "mc_nu_vy_sce/F");
        fEventTree->Branch("mc_nu_vz_sce", &fTrueNu_VzSce, "mc_nu_vz_sce/F");
        fEventTree->Branch("mc_nu_energy", &fTrueNu_Energy, "mc_nu_energy/F");
        fEventTree->Branch("mc_nu_px", &fTrueNu_Px, "mc_nu_px/F");
        fEventTree->Branch("mc_nu_py", &fTrueNu_Py, "mc_nu_py/F");
        fEventTree->Branch("mc_nu_pz", &fTrueNu_Pz, "mc_nu_pz/F");
        fEventTree->Branch("mc_nu_lepton_energy", &fTrueNu_LeptonEnergy, "mc_nu_lepton_energy/F");
        fEventTree->Branch("mc_nu_lepton_px", &fTrueNu_LeptonPx, "mc_nu_lepton_px/F");
        fEventTree->Branch("mc_nu_lepton_py", &fTrueNu_LeptonPy, "mc_nu_lepton_py/F");
        fEventTree->Branch("mc_nu_lepton_pz", &fTrueNu_LeptonPz, "mc_nu_lepton_pz/F");
        fEventTree->Branch("mc_nu_lepton_theta", &fTrueNu_LeptonTheta, "mc_nu_lepton_theta/F");
        fEventTree->Branch("mc_nu_time", &fTrueNu_Time, "mc_nu_time/F");
        fEventTree->Branch("mc_nu_pdg", &fTrueNu_PDG, "mc_nu_pdg/I");
        fEventTree->Branch("mc_nu_interaction_type", &fTrueNu_InteractionType, "mc_nu_interaction_type/I");
        fEventTree->Branch("mc_nu_ccnc", &fTrueNu_CCNC, "mc_nu_ccnc/O");
        fEventTree->Branch("mc_nu_vtx_distance", &fTrueNu_VtxDistance, "mc_nu_vtx_distance/F");
        fEventTree->Branch("num_matched_daughters", &fNumMatchedDaughters, "num_matched_daughters/i");
        fEventTree->Branch("cosmic_matched", &fCosmicMatched, "cosmic_matched/O");

        fEventTree->Branch("mc_nu_daughter_matched", "std::vector< bool >", &fTrueNu_DaughterMatched);
        fEventTree->Branch("mc_nu_daughter_pdg", "std::vector< int >", &fTrueNu_DaughterPDG);
        fEventTree->Branch("mc_nu_daughter_energy", "std::vector< float >", &fTrueNu_DaughterE);
    }

    //// Tree for every daughter
    fNueDaughtersTree = tfs->make<TTree>("Daughters", "Daughters Tree");
    fNueDaughtersTree->Branch("event", &fEvent, "event/i");
    fNueDaughtersTree->Branch("run", &fRun, "run/i");
    fNueDaughtersTree->Branch("subrun", &fSubrun, "subrun/i");
    fNueDaughtersTree->Branch("evt_time_sec", &fTimeHigh, "evt_time_sec/i");
    fNueDaughtersTree->Branch("evt_time_nsec", &fTimeLow, "evt_time_nsec/i");
    fNueDaughtersTree->Branch("hitsU", &fNhitsU, "hitsU/i");
    fNueDaughtersTree->Branch("hitsV", &fNhitsV, "hitsV/i");
    fNueDaughtersTree->Branch("hitsY", &fNhitsY, "hitsY/i");
    fNueDaughtersTree->Branch("caloU", &fCaloU, "caloU/F");
    fNueDaughtersTree->Branch("caloV", &fCaloV, "caloV/F");
    fNueDaughtersTree->Branch("caloY", &fCaloY, "caloY/F");
    fNueDaughtersTree->Branch("hitsSps", &fNSpacepoints, "hitsSps/i");
    fNueDaughtersTree->Branch("generation", &fGeneration, "generation/i");
    fNueDaughtersTree->Branch("track_score", &fTrackScore, "track_score/F");
    fNueDaughtersTree->Branch("is_shower", &fIsShower, "is_shower/O");
    fNueDaughtersTree->Branch("is_track", &fIsTrack, "is_track/O");
    fNueDaughtersTree->Branch("has_shower_daughter", &fHasShowerDaughter, "has_shower_daughter/O");
    fNueDaughtersTree->Branch("is_track_daughter", &fIsTrackDaughter, "is_track_daughter/O");
    fNueDaughtersTree->Branch("vx", &fVx, "vx/F");
    fNueDaughtersTree->Branch("vy", &fVy, "vy/F");
    fNueDaughtersTree->Branch("vz", &fVz, "vz/F");
    fNueDaughtersTree->Branch("start_contained", &fStartContained, "start_contained/O");
    fNueDaughtersTree->Branch("vtx_distance", &fVtxDistance, "vtx_distance/F");

    fNueDaughtersTree->Branch("track_length", &fTrackLength, "track_length/F");
    fNueDaughtersTree->Branch("track_endx", &fTrackEndX, "track_endx/F");
    fNueDaughtersTree->Branch("track_endy", &fTrackEndY, "track_endy/F");
    fNueDaughtersTree->Branch("track_endz", &fTrackEndZ, "track_endz/F");
    fNueDaughtersTree->Branch("track_dirx", &fTrackDirX, "track_dirx/F");
    fNueDaughtersTree->Branch("track_diry", &fTrackDirY, "track_diry/F");
    fNueDaughtersTree->Branch("track_dirz", &fTrackDirZ, "track_dirz/F");
    fNueDaughtersTree->Branch("track_is_muon_candidate", &fIsMuonCandidate, "track_is_muon_candidate/O");
    fNueDaughtersTree->Branch("track_range_mom_p", &fTrackRange_mom_p, "track_range_mom_p/F");
    fNueDaughtersTree->Branch("track_range_mom_mu", &fTrackRange_mom_mu, "track_range_mom_mu/F");
    fNueDaughtersTree->Branch("track_mcs_mom", &fTrackMCS_mom, "track_mcs_mom/F");
    fNueDaughtersTree->Branch("track_mcs_err", &fTrackMCS_err, "track_mcs_err/F");
    fNueDaughtersTree->Branch("track_mcs_ll", &fTrackMCS_ll, "track_mcs_ll/F");
    fNueDaughtersTree->Branch("track_chi2_proton", &fTrackPID_chiproton, "track_chi2_proton/F");
    fNueDaughtersTree->Branch("track_chi2_muon", &fTrackPID_chimuon, "track_chi2_muon/F");
    fNueDaughtersTree->Branch("shower_length", &fShowerLength, "shower_length/F");
    fNueDaughtersTree->Branch("shower_openangle", &fShowerOpenAngle, "shower_openangle/F");
    fNueDaughtersTree->Branch("shower_dirx", &fShowerDirX, "shower_dirx/F");
    fNueDaughtersTree->Branch("shower_diry", &fShowerDirY, "shower_diry/F");
    fNueDaughtersTree->Branch("shower_dirz", &fShowerDirZ, "shower_dirz/F");
    fNueDaughtersTree->Branch("start_dedxU", &fDedxU, "start_dedxU/F");
    fNueDaughtersTree->Branch("start_dedxV", &fDedxV, "start_dedxV/F");
    fNueDaughtersTree->Branch("start_dedxY", &fDedxY, "start_dedxY/F");
    fNueDaughtersTree->Branch("start_hitsU", &fDedxHitsU, "start_hitsU/I");
    fNueDaughtersTree->Branch("start_hitsV", &fDedxHitsV, "start_hitsV/I");
    fNueDaughtersTree->Branch("start_hitsY", &fDedxHitsY, "start_hitsY/I");
    fNueDaughtersTree->Branch("start_pitchU", &fDedxPitchU, "start_pitchU/F");
    fNueDaughtersTree->Branch("start_pitchV", &fDedxPitchV, "start_pitchV/F");
    fNueDaughtersTree->Branch("start_pitchY", &fDedxPitchY, "start_pitchY/F");

    if (!m_isData)
    {
        fNueDaughtersTree->Branch("mc_neutrino", &fMatchedNeutrino, "mc_neutrino/O");
        fNueDaughtersTree->Branch("mc_pdg", &fTruePDG, "mc_pdg/I");
        fNueDaughtersTree->Branch("mc_hit_fraction", &fTrueHitFraction, "mc_hit_fraction/F");
        fNueDaughtersTree->Branch("mc_vx", &fTrueVx, "mc_vx/F");
        fNueDaughtersTree->Branch("mc_vy", &fTrueVy, "mc_vy/F");
        fNueDaughtersTree->Branch("mc_vz", &fTrueVz, "mc_vz/F");
        fNueDaughtersTree->Branch("mc_vx_sce", &fTrueVxSce, "mc_vx_sce/F");
        fNueDaughtersTree->Branch("mc_vy_sce", &fTrueVySce, "mc_vy_sce/F");
        fNueDaughtersTree->Branch("mc_vz_sce", &fTrueVzSce, "mc_vz_sce/F");
        fNueDaughtersTree->Branch("mc_length", &fTrueLength, "mc_length/F");
        fNueDaughtersTree->Branch("mc_energy", &fTrueEnergy, "mc_energy/F");
        fNueDaughtersTree->Branch("mc_px", &fTruePx, "mc_px/F");
        fNueDaughtersTree->Branch("mc_py", &fTruePy, "mc_py/F");
        fNueDaughtersTree->Branch("mc_pz", &fTruePz, "mc_pz/F");
    }

    fSubrunTree = tfs->make<TTree>("subruns", "SubRun Tree");
    fSubrunTree->Branch("run", &m_run, "run/i");
    fSubrunTree->Branch("subRun", &m_subrun, "subRun/i");

    if (!m_isData)
        fSubrunTree->Branch("pot", &m_pot, "pot/F");
}

void NuCC::clearEvent()
{
    fNu_PDG = 0; // if 0, no neutrinocandidate was found, only look at truth information.
    fDaughtersStored = true;
    fCosmicMatched = false;
    fIsNuMuCC = false;
    fDaughtersStartContained = true;
    fNumMatchedDaughters = 0;
    fNumShowers = 0;
    fNumTracks = 0;
    fNu_NhitsU = 0;
    fNu_NhitsV = 0;
    fNu_NhitsY = 0;
    fNu_CaloU = 0;
    fNu_CaloV = 0;
    fNu_CaloY = 0;
    fNu_NSpacepoints = 0;
    fNumNu = 0;
    fNu_FlashChi2 = 0;
    fBestObviousCosmic_FlashChi2 = 0;
    fTrueNu_DaughterPDG.clear();
    fTrueNu_DaughterE.clear();
    fTrueNu_DaughterMatched.clear();

    pfparticles.clear();
    pfneutrinos.clear();
    pfdaughters.clear();
    pfshowers.clear();
    pftracks.clear();
    particleMap.clear();
    particlesToMetadata.clear();
    particlesToVertices.clear();
    particlesToClusters.clear();
    particlesToSpacePoints.clear();
    particlesToShowers.clear();
    particlesToTracks.clear();
    clustersToHits.clear();
    hitsToSpacePoints.clear();
    spacePointsToHits.clear();

    matchedParticles.clear();
    matchedHitFractions.clear();
    matchedMCParticles.clear();

    m_muon_candidates.clear();
}

void NuCC::clearDaughter()
{
    fIsShower = false;
    fIsTrack = false;
    fHasShowerDaughter = false;
    fIsTrackDaughter = false;

    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fCaloU = 0;
    fCaloV = 0;
    fCaloY = 0;
    fNSpacepoints = 0;

    // Track info
    fTrackLength = 0;
    fTrackDirX = 0;
    fTrackDirY = 0;
    fTrackDirZ = 0;
    fTrackEndX = 0;
    fTrackEndY = 0;
    fTrackEndZ = 0;
    fTrackRange_mom_p = 0;
    fTrackRange_mom_mu = 0;
    fTrackMCS_mom = 0;
    fTrackMCS_err = 0;
    fTrackMCS_ll = 0;
    fTrackPID_chiproton = 0;
    fTrackPID_chimuon = 0;
    fIsMuonCandidate = false;
    // Shower info
    fShowerLength = 0;
    fShowerOpenAngle = 0;
    fShowerDirX = 0;
    fShowerDirY = 0;
    fShowerDirZ = 0;
    fDedxU = 0;
    fDedxV = 0;
    fDedxY = 0;
    fDedxHitsU = 0;
    fDedxHitsV = 0;
    fDedxHitsY = 0;
    fDedxPitchU = 0;
    fDedxPitchV = 0;
    fDedxPitchY = 0;
}

DEFINE_ART_MODULE(NuCC)
#endif // NUCC_H