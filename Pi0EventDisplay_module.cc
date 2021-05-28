////////////////////////////////////////////////////////////////////////
// Class:       pi0EventDisplay
// Plugin Type: analyzer (art v2_08_03)
// File:        pi0EventDisplay_module.cc
//
// Generated at Fri Oct 20 15:42:43 2017 by Lorena Escudero Sanchez using cetskelgen
// from cetlib version v3_01_01.
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nug4/ParticleNavigation/ParticleList.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/Simulation/sim.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "fhiclcpp/ParameterSet.h"
#include "larreco/RecoAlg/PMAlg/Utilities.h"
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardata/Utilities/DatabaseUtil.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "canvas/Utilities/InputTag.h"

#include <vector>
#include <map>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "larevt/Filters/ChannelFilter.h"

#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

namespace test {
  class pi0EventDisplay;
  struct RecoTrueMatch;
  std::pair<std::pair<double,double>,std::vector<const recob::Hit*>> getMCParticleHits(const simb::MCParticle &mcpart, const art::Event &evt, std::string hitModule);
  std::vector<double> estimateEnergyFromHitCharge(const art::Event &evt, const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg);
  std::vector<RecoTrueMatch> GetMCParticleListFromShowerHits (const std::vector<const recob::Hit*>& hitVec, const art::Event &evt, std::string hitModule);
  int GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule);
  const std::vector<const recob::Hit*> GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule);
  std::vector<std::vector<double>> SimIDEXYZ_Finder (const art::Event &evt, const std::vector<const recob::Hit*> &hitVec);
  void ShwrSpcPntXYZ_Finder (const art::Event &evt, std::string pfpTag, std::string showerTag, std::vector<std::vector<double>> &recoX, std::vector<std::vector<double>> &recoY, std::vector<std::vector<double>> &recoZ);
}

struct test::RecoTrueMatch {
  const simb::MCParticle* part;
  double weight;
  unsigned int matched_hits;
};

//------------------------------------------------------------------------------------------------------------------------------------------

class test::pi0EventDisplay : public art::EDAnalyzer 
{
public:
  explicit pi0EventDisplay(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  pi0EventDisplay(pi0EventDisplay const &) = delete;
  pi0EventDisplay(pi0EventDisplay &&) = delete;
  pi0EventDisplay & operator = (pi0EventDisplay const &) = delete;
  pi0EventDisplay & operator = (pi0EventDisplay &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;
  virtual void endSubRun(const art::SubRun& sr) override;
  
  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // fcl parameters
  std::string fCalorimetryTag;
  std::string fTrackerTag;
  std::string fShowerTag;
  std::string fPFParticleTag;
  std::string fGeneratorTag;
  std::string fHitTag;
  std::string fSpcPntTag;
  std::string fSimTag;
  bool fVerbose;

  art::ServiceHandle<cheat::ParticleInventoryService> particleinventory;
  calo::CalorimetryAlg fCalorimetryAlg;
  
  //////////////////
  ///// Create trees
  //////////////////
  TTree *fPiZeroShowerTree;
  TTree *fPOTTree;
  
  ///////////////////////
  ///// declare variables
  ///////////////////////
  // event info
  unsigned int fEvent;
  unsigned long long fTime;
  unsigned int nu_in_fv;  
  int fIsCC, fNuPDG, fNuPDGunosc, fMode, fLepPDG; 
  double fEv, fQ2, fW, fX, fY, fNuMomX, fNuMomY, fNuMomZ, fLepMomX, fLepMomY, fLepMomZ, fLepE, fLepNuAngle;
  double vtx_x, vtx_y, vtx_z, vtx_X;
  double fPOT;
  
  // pi0 info
  std::vector<unsigned int> fPi0CandidateID;
  std::vector<int>          fPi0CandidatePDG;
  std::vector<unsigned int> fPi0CandidateNumDaughters;
  std::vector<unsigned int> fPi0CandidateNumPhotonDaughters;
  std::vector<double>       fPi0CandidateStartX;
  std::vector<double>       fPi0CandidateStartY;
  std::vector<double>       fPi0CandidateStartZ;
  std::vector<double>       fPi0CandidateEndX;
  std::vector<double>       fPi0CandidateEndY;
  std::vector<double>       fPi0CandidateEndZ;

  // photon info
  std::vector<unsigned int> fPhotonID;
  std::vector<double>       fPhotonStartX;
  std::vector<double>       fPhotonStartY;
  std::vector<double>       fPhotonStartZ;
  std::vector<double>       fPhotonEndX;
  std::vector<double>       fPhotonEndY;
  std::vector<double>       fPhotonEndZ;
  std::vector<double>       fPhotonE;
  std::vector<double>       fPhotonIDEEnergy;
  std::vector<double>       fPhotonIDEEFrac;
  std::vector<double>       fPhotonEfromHits;
  std::vector<unsigned int> fPhotonNumHits;

  // recob::Shower info
  std::vector<int>          fPFShowerBestMCParentPDG;
  std::vector<unsigned int> fPFShowerBestMCParentID;
  std::vector<double>       fCheatedShowerEnergy;
  std::vector<double>       fPFShowerEnergy0;
  std::vector<double>       fPFShowerEnergy1;
  std::vector<double>       fPFShowerEnergyCollection;
  std::vector<double>       fPFShowerBestMCParentWeight;
  std::vector<double>       fPFShowerDirectionX;
  std::vector<double>       fPFShowerDirectionY;
  std::vector<double>       fPFShowerDirectionZ;
  std::vector<double>       fPFShowerStartX;
  std::vector<double>       fPFShowerStartY;
  std::vector<double>       fPFShowerStartZ;
  std::vector<unsigned int> fPFShowerNumHits;
  std::vector<double>       fPFShowerBestMCParentCompleteness;
  std::vector<double>       fPFShowerBestMCParentPurity;
  std::vector<double>       fPFShowerBestMCParentMatchedHits;

  // spacial info
  // sim::IDEs from true photons
  std::vector<std::vector<double>> fTruePhotonSimIDEX;
  std::vector<std::vector<double>> fTruePhotonSimIDEY;
  std::vector<std::vector<double>> fTruePhotonSimIDEZ;

  std::vector<std::vector<double>> fTruePhotonSimIDEE;

  // sim::IDEs backtracked from recob::Hits in recob::Showers
  std::vector<std::vector<double>> fShowerSimIDEX;
  std::vector<std::vector<double>> fShowerSimIDEY;
  std::vector<std::vector<double>> fShowerSimIDEZ;
  // sim::IDEs backtracked from all recob::Hits in event
  std::vector<std::vector<double>> fTrueIDEX;
  std::vector<std::vector<double>> fTrueIDEY;
  std::vector<std::vector<double>> fTrueIDEZ;
  // recob::SpacePoints from recob::Showers
  std::vector<std::vector<double>> fShowerSpacePointX;
  std::vector<std::vector<double>> fShowerSpacePointY;
  std::vector<std::vector<double>> fShowerSpacePointZ;
  // all recob::SpacePoints from event
  std::vector<std::vector<double>> fEventSpacePointX;
  std::vector<std::vector<double>> fEventSpacePointY;
  std::vector<std::vector<double>> fEventSpacePointZ;
  // hit info (all hits)
  std::vector<std::vector<double>> fHitPeakTime0;
  std::vector<std::vector<double>> fHitWire0;
  std::vector<std::vector<double>> fHitPlane0;
  std::vector<std::vector<double>> fHitPeakTime1;
  std::vector<std::vector<double>> fHitWire1;
  std::vector<std::vector<double>> fHitPlane1;
  std::vector<std::vector<double>> fHitPeakTime2;
  std::vector<std::vector<double>> fHitWire2;
  std::vector<std::vector<double>> fHitPlane2;
  // hit info (showers)
  std::vector<std::vector<double>> fShowerHitPeakTime0;
  std::vector<std::vector<double>> fShowerHitWire0;
  std::vector<std::vector<double>> fShowerHitPlane0;
  std::vector<std::vector<double>> fShowerHitPeakTime1;
  std::vector<std::vector<double>> fShowerHitWire1;
  std::vector<std::vector<double>> fShowerHitPlane1;
  std::vector<std::vector<double>> fShowerHitPeakTime2;
  std::vector<std::vector<double>> fShowerHitWire2;
  std::vector<std::vector<double>> fShowerHitPlane2;
  
};

//------------------------------------------------------------------------------------------------------------------------------------------

////////////////////
///// fcl parameters
////////////////////
test::pi0EventDisplay::pi0EventDisplay(fhicl::ParameterSet const & p) :
  EDAnalyzer(p),
  fCalorimetryTag(p.get<std::string>("CalorimetryTag")),
  fTrackerTag(p.get<std::string>("TrackerTag")), // This was set to ShowerTag initially? Not sure why
  fShowerTag(p.get<std::string>("ShowerTag")),
  fPFParticleTag(p.get<std::string>("PFParticleTag")),
  fGeneratorTag(p.get<std::string>("GeneratorTag")),
  fHitTag(p.get<std::string>("HitTag")),
  fSpcPntTag(p.get<std::string>("SpcPntTag")),
  fSimTag(p.get<std::string>("SimTag")),
  fVerbose(p.get<bool>("Verbose")),
  fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
                                                                                                          
}

//------------------------------------------------------------------------------------------------------------------------------------------

void test::pi0EventDisplay::beginJob()
{                                                                                                                                 
  art::ServiceHandle<art::TFileService> tfs;
  fPiZeroShowerTree = tfs->make<TTree>("PiZeroShower","PiZeroShower");
  fPOTTree = tfs->make<TTree>("POT", "POT");

  fPiZeroShowerTree->Branch("event",                                &fEvent, "event/i");
  fPiZeroShowerTree->Branch("nu_in_fd",                             &nu_in_fv, "nu_in_fv/i");  
  fPiZeroShowerTree->Branch("isCC",                                 &fIsCC,       "isCC/I");
  fPiZeroShowerTree->Branch("nuPDG",                                &fNuPDG,        "nuPDG/I");
  fPiZeroShowerTree->Branch("nuPDGunosc",                           &fNuPDGunosc,   "nuPDGunosc/I");
  fPiZeroShowerTree->Branch("NuMomX",                               &fNuMomX,       "NuMomX/D");
  fPiZeroShowerTree->Branch("NuMomY",                               &fNuMomY,       "NuMomY/D");
  fPiZeroShowerTree->Branch("NuMomZ",                               &fNuMomZ,       "NuMomZ/D");
  fPiZeroShowerTree->Branch("Ev",                                   &fEv,           "Ev/D");
  fPiZeroShowerTree->Branch("mode",                                 &fMode,         "mode/I");
  fPiZeroShowerTree->Branch("LepPDG",                               &fLepPDG,       "LepPDG/I");
  fPiZeroShowerTree->Branch("LepMomX",                              &fLepMomX,      "LepMomX/D");
  fPiZeroShowerTree->Branch("LepMomY",                              &fLepMomY,      "LepMomY/D");
  fPiZeroShowerTree->Branch("LepMomZ",                              &fLepMomZ,      "LepMomZ/D");
  fPiZeroShowerTree->Branch("LepE",                                 &fLepE,         "LepE/D");
  fPiZeroShowerTree->Branch("LepNuAngle",                           &fLepNuAngle,   "LepNuAngle/D");
  fPiZeroShowerTree->Branch("Q2",                                   &fQ2,           "Q2/D");
  fPiZeroShowerTree->Branch("W",                                    &fW,            "W/D");
  fPiZeroShowerTree->Branch("X",                                    &fX,            "X/D");
  fPiZeroShowerTree->Branch("Y",                                    &fY,            "Y/D");
  fPiZeroShowerTree->Branch("vtx_x",                                &vtx_x,    "vtx_x/D");
  fPiZeroShowerTree->Branch("vtx_y",                                &vtx_y,    "vtx_y/D");
  fPiZeroShowerTree->Branch("vtx_z",                                &vtx_z,    "vtx_z/D");
  fPiZeroShowerTree->Branch("vtx_X",                                &vtx_X,    "vtx_X/D");
  fPOTTree->Branch("POT",                                           &fPOT, "POT/D");
  fPOT=0;
                              
  fPiZeroShowerTree->Branch("Pi0CandidateID",                       &fPi0CandidateID);
  fPiZeroShowerTree->Branch("Pi0CandidatePDG",                      &fPi0CandidatePDG);
  fPiZeroShowerTree->Branch("Pi0CadidateNumDaughters",              &fPi0CandidateNumDaughters);
  fPiZeroShowerTree->Branch("Pi0CandidateNumPhotonDaughters",       &fPi0CandidateNumPhotonDaughters);
  fPiZeroShowerTree->Branch("Pi0CandidateStartX",                   &fPi0CandidateStartX);
  fPiZeroShowerTree->Branch("Pi0CandidateStartY",                   &fPi0CandidateStartY);
  fPiZeroShowerTree->Branch("Pi0CandidateStartZ",                   &fPi0CandidateStartZ);
  fPiZeroShowerTree->Branch("Pi0CandidateEndX",                     &fPi0CandidateEndX);
  fPiZeroShowerTree->Branch("Pi0CandidateEndY",                     &fPi0CandidateEndY);
  fPiZeroShowerTree->Branch("Pi0CandidateEndZ",                     &fPi0CandidateEndZ);

  fPiZeroShowerTree->Branch("PhotonID",                             &fPhotonID);
  fPiZeroShowerTree->Branch("PhotonStartX",                         &fPhotonStartX);
  fPiZeroShowerTree->Branch("PhotonStartY",                         &fPhotonStartY);
  fPiZeroShowerTree->Branch("PhotonStartZ",                         &fPhotonStartZ );
  fPiZeroShowerTree->Branch("PhotonEndX",                           &fPhotonEndX);
  fPiZeroShowerTree->Branch("PhotonEndY",                           &fPhotonEndY);
  fPiZeroShowerTree->Branch("PhotonEndZ",                           &fPhotonEndZ);
  fPiZeroShowerTree->Branch("PhotonE",                              &fPhotonE);
  fPiZeroShowerTree->Branch("PhotonIDEEnergy",                      &fPhotonIDEEnergy);
  fPiZeroShowerTree->Branch("PhotonIDEEFrac",                       &fPhotonIDEEFrac);
  fPiZeroShowerTree->Branch("PhotonEfromHits",                      &fPhotonEfromHits);
  fPiZeroShowerTree->Branch("PhotonNumHits",                        &fPhotonNumHits);

  fPiZeroShowerTree->Branch("PFShowerBestMCParentPDG",              &fPFShowerBestMCParentPDG);
  fPiZeroShowerTree->Branch("PFShowerBestMCParentID",               &fPFShowerBestMCParentID);
  fPiZeroShowerTree->Branch("CheatedShowerEnergy",                  &fCheatedShowerEnergy);
  fPiZeroShowerTree->Branch("PFShowerEnergy0",                      &fPFShowerEnergy0);
  fPiZeroShowerTree->Branch("PFShowerEnergy1",                      &fPFShowerEnergy1);
  fPiZeroShowerTree->Branch("PFShowerEnergyCollection",             &fPFShowerEnergyCollection);
  fPiZeroShowerTree->Branch("PFShowerBestMCParentWeight",           &fPFShowerBestMCParentWeight);
  fPiZeroShowerTree->Branch("PFShowerDirectionX",                   &fPFShowerDirectionX);
  fPiZeroShowerTree->Branch("PFShowerDirectionY",                   &fPFShowerDirectionY);
  fPiZeroShowerTree->Branch("PFShowerDirectionZ",                   &fPFShowerDirectionZ);
  fPiZeroShowerTree->Branch("PFShowerStartX",                       &fPFShowerStartX);
  fPiZeroShowerTree->Branch("PFShowerStartY",                       &fPFShowerStartY);
  fPiZeroShowerTree->Branch("PFShowerStartZ",                       &fPFShowerStartZ);
  fPiZeroShowerTree->Branch("PFShowerNumHits",                      &fPFShowerNumHits);
  fPiZeroShowerTree->Branch("PFShowerBestMCParentCompleteness",     &fPFShowerBestMCParentCompleteness);
  fPiZeroShowerTree->Branch("PFShowerBestMCParentPurity",           &fPFShowerBestMCParentPurity);
  fPiZeroShowerTree->Branch("PFShowerBestMCParentMatchedHits",      &fPFShowerBestMCParentMatchedHits);
  
  fPiZeroShowerTree->Branch("ShowerSpacePointX",                    &fShowerSpacePointX);
  fPiZeroShowerTree->Branch("ShowerSpacePointY",                    &fShowerSpacePointY);
  fPiZeroShowerTree->Branch("ShowerSpacePointZ",                    &fShowerSpacePointZ);
  fPiZeroShowerTree->Branch("ShowerSimIDEX",                        &fShowerSimIDEX);
  fPiZeroShowerTree->Branch("ShowerSimIDEY",                        &fShowerSimIDEY);
  fPiZeroShowerTree->Branch("ShowerSimIDEZ",                        &fShowerSimIDEZ);
  fPiZeroShowerTree->Branch("EventSpacePointX",                     &fEventSpacePointX);
  fPiZeroShowerTree->Branch("EventSpacePointY",                     &fEventSpacePointY);
  fPiZeroShowerTree->Branch("EventSpacePointZ",                     &fEventSpacePointZ);
  fPiZeroShowerTree->Branch("TrueIDEX",                             &fTrueIDEX);
  fPiZeroShowerTree->Branch("TrueIDEY",                             &fTrueIDEY);
  fPiZeroShowerTree->Branch("TrueIDEZ",                             &fTrueIDEZ);
  fPiZeroShowerTree->Branch("TruePhotonSimIDEX",                    &fTruePhotonSimIDEX);
  fPiZeroShowerTree->Branch("TruePhotonSimIDEY",                    &fTruePhotonSimIDEY);
  fPiZeroShowerTree->Branch("TruePhotonSimIDEZ",                    &fTruePhotonSimIDEZ);

  fPiZeroShowerTree->Branch("TruePhotonSimIDEE",                    &fTruePhotonSimIDEE);

  fPiZeroShowerTree->Branch("HitPeakTime0",                         &fHitPeakTime0);
  fPiZeroShowerTree->Branch("HitWire0",                             &fHitWire0);
  fPiZeroShowerTree->Branch("HitPlane0",                            &fHitPlane0);
  fPiZeroShowerTree->Branch("HitPeakTime1",                         &fHitPeakTime1);
  fPiZeroShowerTree->Branch("HitWire1",                             &fHitWire1);
  fPiZeroShowerTree->Branch("HitPlane1",                            &fHitPlane1);
  fPiZeroShowerTree->Branch("HitPeakTime2",                         &fHitPeakTime2);
  fPiZeroShowerTree->Branch("HitWire2",                             &fHitWire2);
  fPiZeroShowerTree->Branch("HitPlane2",                            &fHitPlane2);
  fPiZeroShowerTree->Branch("ShowerHitPeakTime0",                   &fShowerHitPeakTime0);
  fPiZeroShowerTree->Branch("ShowerHitWire0",                       &fShowerHitWire0);
  fPiZeroShowerTree->Branch("ShowerHitPlane0",                      &fShowerHitPlane0);
  fPiZeroShowerTree->Branch("ShowerHitPeakTime1",                   &fShowerHitPeakTime1);
  fPiZeroShowerTree->Branch("ShowerHitWire1",                       &fShowerHitWire1);
  fPiZeroShowerTree->Branch("ShowerHitPlane1",                      &fShowerHitPlane1);
  fPiZeroShowerTree->Branch("ShowerHitPeakTime2",                   &fShowerHitPeakTime2);
  fPiZeroShowerTree->Branch("ShowerHitWire2",                       &fShowerHitWire2);
  fPiZeroShowerTree->Branch("ShowerHitPlane2",                      &fShowerHitPlane2);

}

//-----------------------------------------------------------------------------------------------------------------------------------------

void test::pi0EventDisplay::endSubRun(const art::SubRun& sr)
{
  art::Handle< sumdata::POTSummary > pots;
  if( sr.getByLabel("generator",pots) ) fPOT += pots->totpot;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Functions used for analysis

///////////////////////////////////////////////////////////////
///// For each event, find recob::Hits belonging to MCParticles
///////////////////////////////////////////////////////////////
std::pair<std::pair<double,double>, std::vector<const recob::Hit*>> test::getMCParticleHits(
  const simb::MCParticle &mcpart, const art::Event &evt, std::string hitModule) {

  std::vector<const recob::Hit*> outVec;
  
  art::Handle<std::vector<recob::Hit>> hitHandle;
  if(!evt.getByLabel(hitModule, hitHandle)) {
    std::cout << "GetMCParticleHits: could not find hits in event.\n";
    return std::make_pair(std::make_pair(0,0),outVec);
  }
  double energy = 0;
  double energyFrac = 0;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // Backtrack all hits to verify whether they belong to the current MCParticle.
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  for(const recob::Hit& hit : *hitHandle) {
    for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(clockData,hit)) {
      if(pi_serv->TrackIdToParticle_P(ide.trackID) == 
         pi_serv->TrackIdToParticle_P(mcpart.TrackId())) {
        outVec.push_back(&hit);
        if (hit.WireID().Plane != 2) break;
        energy+=ide.energy;
        energyFrac += (ide.energy*ide.energyFrac);
        break;
      }
    }
  }
  energyFrac /= energy;

  return std::make_pair(std::make_pair(energy,energyFrac),outVec);
}

////////////////////////////////
///// Get index of recob::Shower
////////////////////////////////
int test::GetShowerIndex(const recob::Shower &shower, art::Event const &evt, const std::string showerModule){

  if(shower.ID() != -999) return shower.ID();

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);

  // Iterate through all showers to find the matching one to our shower
  int actualIndex = shower.ID();
  if(shower.ID() < 0){
    for(unsigned int s = 0; s < recoShowers->size(); ++s){
      const recob::Shower thisShower = (*recoShowers)[s];
      // Can't compare actual objects so look at a property
      if(fabs(thisShower.Length() - shower.Length()) < 1e-5){
        actualIndex = s;
        continue;
      }
    }
  }

  return actualIndex;

}

/////////////////////////////////////////////////
///// Get recob::Hits belonging to recob::Showers
/////////////////////////////////////////////////
const std::vector<const recob::Hit*> test::GetRecoShowerHits(const recob::Shower &shower, art::Event const &evt, const std::string showerModule) {

  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower> >(showerModule);
  art::FindManyP<recob::Hit> findHits(recoShowers,evt,showerModule);

  // Shower.ID is sometimes at a default value - make sure we get the correct one
  int actualIndex = GetShowerIndex(shower,evt,showerModule);

  std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(actualIndex);
  std::vector<const recob::Hit*> showerHits;

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  for(const art::Ptr<recob::Hit> hit : inputHits)
    {
      showerHits.push_back(hit.get());	  	  
    }
  return showerHits;  

}

//////////////////////////////////////
///// Truth/Reco matching happens here
//////////////////////////////////////
std::vector<test::RecoTrueMatch>
test::GetMCParticleListFromShowerHits
(const std::vector<const recob::Hit*>& hitVec, const art::Event &evt, std::string hitModule) {

  using weightedMCPair = std::pair<const simb::MCParticle*, double>;
  std::vector<weightedMCPair> outVec;
  std::vector<RecoTrueMatch> matches;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

  // Loop over all hits in the input vector and record the contributing MCParticles.                                                                         
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::unordered_map<const simb::MCParticle*, double> mcEMap;
  std::unordered_map<const simb::MCParticle*, unsigned int> mcHMap;
  double hitTotalE = 0;
  for(const recob::Hit* hitp : hitVec) {
    for(const sim::TrackIDE& ide : bt_serv->HitToEveTrackIDEs(clockData,*hitp)) {
      const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
      mcEMap[curr_part] += ide.energy;
      mcHMap[curr_part] += 1;
      hitTotalE += ide.energy;
    }
  }

  // Fill and sort the output vector                                                                               
  for (weightedMCPair const& p : mcEMap) {
    outVec.push_back(p);
  }
  std::sort(outVec.begin(), outVec.end(),
	    [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});

  // Normalise the weights by the total track energy.         
  if (hitTotalE < 1e-5) { hitTotalE = 1; } // Protect against zero division                                                                                                                                                                 
  for (weightedMCPair& p : outVec) {
    double wei = p.second / hitTotalE;
    RecoTrueMatch m {p.first, wei, mcHMap[p.first]};
    matches.push_back(m);
  }

  return matches;
}

//////////////////////////////////////
///// Approximate recob::Shower energy
//////////////////////////////////////
std::vector<double> test::estimateEnergyFromHitCharge(const art::Event &evt, const std::vector<const recob::Hit*> &hits, calo::CalorimetryAlg caloAlg)
{
  double kGeVtoElectrons { 4.237e7 }; // obtained from utils class.. Copied for now, should use class (although this is a physical constant, so hopefully doesn't change).
  double recombination   { 1/0.63 };

  std::vector<double> showerEnergy = {0,0,0};

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob();

  // Find the total charge on each plane
  for ( size_t h{0} ; h < hits.size() ; h++ ) {
    const recob::Hit* hit = hits[h];
    const int plane = hit->WireID().Plane;
    showerEnergy[ plane ] += ( caloAlg.ElectronsFromADCArea( hit->Integral(), plane) * caloAlg.LifetimeCorrection(clockData,detProp,hit->PeakTime(),clockData.TriggerTime()) ) / kGeVtoElectrons;
  }

  showerEnergy[0] *= recombination;
  showerEnergy[1] *= recombination;
  showerEnergy[2] *= recombination;

  // caloAlg.ElectronsFromADCArea( hit->Integral(), plane) -> Does hit->Integral()/AreaConstants(plane)
  // AreaConstants(plane) is defined in calorimetry_pdune.fcl. Although that fcl file has a typo.
  // These probably need tuning for protodune data.

  return showerEnergy;
}

/////////////////////////////////////////////////////
///// recob::Hits --> sim::IDEs using the backtracker
/////////////////////////////////////////////////////
std::vector<std::vector<double>> test::SimIDEXYZ_Finder (const art::Event &evt, const std::vector<const recob::Hit*> &hitVec) {

  std::vector<std::vector<double>> outVec;
  std::vector<double>              xVec;
  std::vector<double>              yVec;
  std::vector<double>              zVec;

  auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  for (const auto &hit : hitVec)
    {
      std::vector< const sim::IDE * > ides = bt_serv->HitToSimIDEs_Ps(clockData,*hit);
      for (auto& ide: ides)
	{
	  xVec.emplace_back(ide->x);
	  yVec.emplace_back(ide->y);
	  zVec.emplace_back(ide->z);
	}
    }
  
  outVec.emplace_back(xVec);
  outVec.emplace_back(yVec);
  outVec.emplace_back(zVec);
  return outVec;
}

////////////////////////////////////////////////////////////
///// Find recob::SpacePoints associated with recob::Showers
////////////////////////////////////////////////////////////
void test::ShwrSpcPntXYZ_Finder (const art::Event &evt, std::string pfpTag, std::string showerTag, std::vector<std::vector<double>> &recoX, std::vector<std::vector<double>> &recoY, std::vector<std::vector<double>> &recoZ) {

  std::vector<double> xVec;
  std::vector<double> yVec;
  std::vector<double> zVec;

  // Get all the showers in the event                                                                                                                  
  art::Handle<std::vector<recob::Shower>> showerListHandle;
  std::vector<art::Ptr<recob::Shower>>    recoShowers;
  if (evt.getByLabel(showerTag, showerListHandle)) {
    art::fill_ptr_vector(recoShowers, showerListHandle);
  }

  // Get all the pfparticles in the event                                                                                                              
  art::Handle<std::vector<recob::PFParticle>> pfpListHandle;
  std::vector<art::Ptr<recob::PFParticle>>    pfps;
  if (evt.getByLabel(pfpTag, pfpListHandle)) {
    art::fill_ptr_vector(pfps, pfpListHandle);
  }

  // Association between shower and pfparticle                                                                                                         
  art::FindManyP<recob::PFParticle> fmpf(showerListHandle, evt, showerTag);
  // Association between pfparticle and space point                                                                                                    
  art::FindManyP<recob::SpacePoint> fmsp(pfpListHandle, evt, pfpTag);


  // Loop through all showers in the event and find the associated pfparticle and then the associated space points through that                        
  art::Ptr<recob::PFParticle>              pfpart;
  std::vector<art::Ptr<recob::SpacePoint>> spcPnts;

  for (const auto &s : recoShowers) {

    pfpart  = fmpf.at(s.key()).at(0);  // Note: there is only one(?) associated pfparticle with each shower
    std::cout <<fmpf.at(s.key()).size();                                     
    spcPnts = fmsp.at(pfpart.key());

    // Store all of the spacepoints x y z information into temporary vectors                                                                           
    for (const auto &sp : spcPnts) {
      const double *xyz = sp.get()->XYZ();
      xVec.emplace_back(xyz[0]);
      yVec.emplace_back(xyz[1]);
      zVec.emplace_back(xyz[2]);
    }

    // emplace the vectors into place and clear the vectors to store the next shower                                                                   
    recoX.emplace_back(xVec);
    recoY.emplace_back(yVec);
    recoZ.emplace_back(zVec);
    xVec.clear();
    yVec.clear();
    zVec.clear();

  }

}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void test::pi0EventDisplay::analyze(art::Event const & evt)
{
  // clear all variables that are being saved
  nu_in_fv = 0;
  fPi0CandidateID.clear();
  fPi0CandidatePDG.clear();
  fPi0CandidateNumDaughters.clear();
  fPi0CandidateNumPhotonDaughters.clear();
  fPi0CandidateStartX.clear();
  fPi0CandidateStartY.clear();
  fPi0CandidateStartZ.clear();
  fPi0CandidateEndX.clear();
  fPi0CandidateEndY.clear();
  fPi0CandidateEndZ.clear();
  
  fPhotonID.clear();
  fPhotonStartX.clear();
  fPhotonStartY.clear();
  fPhotonStartZ.clear();
  fPhotonEndX.clear();
  fPhotonEndY.clear();
  fPhotonEndZ.clear();
  fPhotonE.clear();
  fPhotonIDEEnergy.clear();
  fPhotonIDEEFrac.clear();
  fPhotonEfromHits.clear();
  fPhotonNumHits.clear();

  fPFShowerBestMCParentPDG.clear();
  fPFShowerBestMCParentID.clear();
  fCheatedShowerEnergy.clear();
  fPFShowerEnergy0.clear();
  fPFShowerEnergy1.clear();
  fPFShowerEnergyCollection.clear();
  fPFShowerBestMCParentWeight.clear();
  fPFShowerDirectionX.clear();
  fPFShowerDirectionY.clear();
  fPFShowerDirectionZ.clear();
  fPFShowerStartX.clear();
  fPFShowerStartY.clear();
  fPFShowerStartZ.clear();
  fPFShowerNumHits.clear();
  fPFShowerBestMCParentCompleteness.clear();
  fPFShowerBestMCParentPurity.clear();
  fPFShowerBestMCParentMatchedHits.clear();
  
  fShowerSpacePointX.clear();
  fShowerSpacePointY.clear();
  fShowerSpacePointZ.clear();
  fShowerSimIDEX.clear();
  fShowerSimIDEY.clear();
  fShowerSimIDEZ.clear();
  fEventSpacePointX.clear();
  fEventSpacePointY.clear();
  fEventSpacePointZ.clear();
  fTrueIDEX.clear();
  fTrueIDEY.clear();
  fTrueIDEZ.clear();
  fTruePhotonSimIDEX.clear();
  fTruePhotonSimIDEY.clear();
  fTruePhotonSimIDEZ.clear();

  fTruePhotonSimIDEE.clear();

  fHitPeakTime0.clear();
  fHitWire0.clear();
  fHitPlane0.clear();
  fHitPeakTime1.clear();
  fHitWire1.clear();
  fHitPlane1.clear();
  fHitPeakTime2.clear();
  fHitWire2.clear();
  fHitPlane2.clear();

  fShowerHitPeakTime0.clear();
  fShowerHitWire0.clear();
  fShowerHitPlane0.clear();
  fShowerHitPeakTime1.clear();
  fShowerHitWire1.clear();
  fShowerHitPlane1.clear();
  fShowerHitPeakTime2.clear();
  fShowerHitWire2.clear();
  fShowerHitPlane2.clear();

  ////////////////////
  ///// Begin analysis
  ////////////////////
  
  // get the event number that is selected
  fEvent = evt.id().event();
  fTime = evt.time().value();
  bool beamTriggerEvent = false;

  std::map<int, std::pair<unsigned int, int>> MCRecoMap;
  MCRecoMap[0] = std::make_pair(0,0);

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  // find beam particle if event is MC
  if(!evt.isRealData())
    {
    art::Handle< std::vector<simb::MCTruth> > mct;
    std::vector< art::Ptr<simb::MCTruth> > truth;
    if( evt.getByLabel("generator", mct) )
      art::fill_ptr_vector(truth, mct);
    else
      mf::LogWarning("CAFMaker") << "No MCTruth.";

    art::Handle< std::vector<simb::MCFlux> > mcf;
    std::vector< art::Ptr<simb::MCFlux> > flux;
    if( evt.getByLabel("generator", mcf) )
      art::fill_ptr_vector(flux, mcf);
    else
      mf::LogWarning("CAFMaker") << "No MCFlux.";

    for(size_t i=0; i<truth.size(); i++){
      if(i>1){
        mf::LogWarning("CAFMaker") << "Skipping MC truth index " << i;
        continue;
      }

      fIsCC     = !(truth[i]->GetNeutrino().CCNC());  // ccnc is 0=CC 1=NC
      fNuPDG    = truth[i]->GetNeutrino().Nu().PdgCode();
      fNuPDGunosc = flux[i]->fntype;
      fMode     = truth[i]->GetNeutrino().Mode(); //0=QE/El, 1=RES, 2=DIS, 3=Coherent production; this is different than mode in ND
      fEv       = truth[i]->GetNeutrino().Nu().E();
      fQ2       = truth[i]->GetNeutrino().QSqr();
      fW        = truth[i]->GetNeutrino().W();
      fX        = truth[i]->GetNeutrino().X();
      fY        = truth[i]->GetNeutrino().Y();
      fNuMomX   = truth[i]->GetNeutrino().Nu().Momentum().X();
      fNuMomY   = truth[i]->GetNeutrino().Nu().Momentum().Y();
      fNuMomZ   = truth[i]->GetNeutrino().Nu().Momentum().Z();
      vtx_X     = truth[i]->GetNeutrino().Nu().Vx();
      vtx_x     = truth[i]->GetNeutrino().Lepton().Vx();
      vtx_y     = truth[i]->GetNeutrino().Lepton().Vy();
      vtx_z     = truth[i]->GetNeutrino().Lepton().Vz();
      //Lepton stuff
      fLepPDG     = truth[i]->GetNeutrino().Lepton().PdgCode();
      fLepMomX    = truth[i]->GetNeutrino().Lepton().Momentum().X();
      fLepMomY    = truth[i]->GetNeutrino().Lepton().Momentum().Y();
      fLepMomZ    = truth[i]->GetNeutrino().Lepton().Momentum().Z();
      fLepE       = truth[i]->GetNeutrino().Lepton().Momentum().T();
      fLepNuAngle = truth[i]->GetNeutrino().Nu().Momentum().Vect().Angle(truth[i]->GetNeutrino().Lepton().Momentum().Vect());
    }

    const sim::ParticleList& trueParticles = particleinventory->ParticleList();
    
    for (auto particle = trueParticles.begin(); particle != trueParticles.end(); ++particle){
      //get primary particle information
      auto geantGoodParticle = particle->second;
						
      unsigned int PrimaryNumDaughters = geantGoodParticle->NumberDaughters();
      std::vector<simb::MCParticle> PhotonDaughters;

      //get info on every daughter photon
      for(unsigned int i=0; i<PrimaryNumDaughters; ++i){
        int PrimaryDaughterId = geantGoodParticle->Daughter(i);
	bool includedTrack = trueParticles.HasParticle(PrimaryDaughterId);
	if(includedTrack == 1){
	const simb::MCParticle* DaughterParticle = particleinventory->TrackIdToParticle_P(PrimaryDaughterId);//get the daughter particle
        if(DaughterParticle->PdgCode() == 22) PhotonDaughters.push_back(*DaughterParticle);
	}
      }  
      if(geantGoodParticle->PdgCode() == 111){//particle qualifies as a pi0 candidate and daughter info is recorded       
        std::cout<<"Pi0 Candidate Found!"<<std::endl;
        fPi0CandidateID.push_back(geantGoodParticle->TrackId());
        fPi0CandidatePDG.push_back(geantGoodParticle->PdgCode());
        fPi0CandidateNumDaughters.push_back(PrimaryNumDaughters);
        fPi0CandidateNumPhotonDaughters.push_back(PhotonDaughters.size());
        fPi0CandidateStartX.push_back(geantGoodParticle->Position().X());
        fPi0CandidateStartY.push_back(geantGoodParticle->Position().Y());
        fPi0CandidateStartZ.push_back(geantGoodParticle->Position().Z());
        fPi0CandidateEndX.push_back(geantGoodParticle->EndPosition().X());
        fPi0CandidateEndY.push_back(geantGoodParticle->EndPosition().Y());
        fPi0CandidateEndZ.push_back(geantGoodParticle->EndPosition().Z());
        
        for (auto part : PhotonDaughters) {
          MCRecoMap[part.TrackId()] = std::make_pair(fPi0CandidateID.back(), fPi0CandidatePDG.back());
          fPhotonID        .push_back(part.TrackId());

	  // Get true photon sim::IDEs
	  std::vector<std::vector<double>> outVec;
	  //std::vector<double>              outVecE;
	  std::vector<double>              eVec;
	  std::vector<double>              xVec;
	  std::vector<double>              yVec;
	  std::vector<double>              zVec;
	  std::vector< const sim::IDE * > ides = bt_serv->TrackIdToSimIDEs_Ps(part.TrackId());
       
	  for (auto& ide: ides)
	    {
	      //eVec.emplace_back(ide->energy);
	      //if (ide->energy > 0.02)
		
	      xVec.emplace_back(ide->x);
	      yVec.emplace_back(ide->y);
	      zVec.emplace_back(ide->z);
		
	    }

	  outVec.emplace_back(xVec);
	  outVec.emplace_back(yVec);
	  outVec.emplace_back(zVec);
	  
	  fTruePhotonSimIDEE.emplace_back(eVec);

	  fTruePhotonSimIDEX.push_back(outVec.at(0));
	  fTruePhotonSimIDEY.push_back(outVec.at(1));
	  fTruePhotonSimIDEZ.push_back(outVec.at(2));


	  // other photon quantities
          fPhotonStartX    .push_back(part.Position().X());
          fPhotonStartY    .push_back(part.Position().Y());
          fPhotonStartZ    .push_back(part.Position().Z());
          fPhotonEndX      .push_back(part.EndPosition().X());
          fPhotonEndY      .push_back(part.EndPosition().Y());
          fPhotonEndZ      .push_back(part.EndPosition().Z());
          fPhotonE         .push_back(part.Momentum().E());
          auto phits = getMCParticleHits(part, evt, fHitTag);
          fPhotonNumHits.push_back(phits.second.size());
          fPhotonIDEEnergy.push_back(phits.first.first);
          fPhotonIDEEFrac.push_back(phits.first.second);
          fPhotonEfromHits .push_back(estimateEnergyFromHitCharge(evt, phits.second, fCalorimetryAlg)[2]);
        }
      }
      PhotonDaughters.clear();
    }
    }
  else
    {
      if(beamTriggerEvent)
	{
	  std::cout << "This data event has a beam trigger" << std::endl;
	}
    }

  //////////////////////////////////////////////////////
  ///// Get info about event hits, SpacePoints, and IDEs
  //////////////////////////////////////////////////////
  std::vector<double> xVec;
  std::vector<double> yVec;
  std::vector<double> zVec;

  std::vector<double> hit_peaktime_temp0;
  std::vector<double> hit_wire_temp0;
  std::vector<double> hit_plane_temp0;
  std::vector<double> hit_peaktime_temp1;
  std::vector<double> hit_wire_temp1;
  std::vector<double> hit_plane_temp1;
  std::vector<double> hit_peaktime_temp2;
  std::vector<double> hit_wire_temp2;
  std::vector<double> hit_plane_temp2;

  //std::vector<const recob::Hit*> outVec;

  // hit handle                                                                                                                                              
  auto hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitTag);

  // space point handle                                                                                                                                      
  auto spcPntHandle = evt.getValidHandle<std::vector<recob::SpacePoint>>(fSpcPntTag);

  // all hits                                                                                                                                                
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hitHandle);

  std::vector<const recob::Hit*> eventHits;

  // Get info from all hits in event
  for(art::Ptr<recob::Hit> hit : hits)
    {
      eventHits.push_back(hit.get());
      auto hit_peaktime = hit->PeakTime();
      auto hit_wire = hit->WireID().Wire;
      auto hit_plane = hit->WireID().Plane;
      if(hit_plane == 0)
	{
          hit_peaktime_temp0.emplace_back(hit_peaktime);
          hit_wire_temp0.emplace_back(hit_wire);
          hit_plane_temp0.emplace_back(hit_plane);
	}
      if(hit_plane == 1)
	{
          hit_peaktime_temp1.emplace_back(hit_peaktime);
          hit_wire_temp1.emplace_back(hit_wire);
          hit_plane_temp1.emplace_back(hit_plane);
	}
      if(hit_plane == 2)
	{
	  hit_peaktime_temp2.emplace_back(hit_peaktime);
	  hit_wire_temp2.emplace_back(hit_wire);
	  hit_plane_temp2.emplace_back(hit_plane);
	}
    } 

  fHitPeakTime0.emplace_back(hit_peaktime_temp0);
  fHitWire0.emplace_back(hit_wire_temp0);
  fHitPlane0.emplace_back(hit_plane_temp0);
  fHitPeakTime1.emplace_back(hit_peaktime_temp1);
  fHitWire1.emplace_back(hit_wire_temp1);
  fHitPlane1.emplace_back(hit_plane_temp1);
  fHitPeakTime2.emplace_back(hit_peaktime_temp2);
  fHitWire2.emplace_back(hit_wire_temp2);
  fHitPlane2.emplace_back(hit_plane_temp2);
  hit_peaktime_temp0.clear();
  hit_wire_temp0.clear();
  hit_plane_temp0.clear();
  hit_peaktime_temp1.clear();
  hit_wire_temp1.clear();
  hit_plane_temp1.clear();
  hit_peaktime_temp2.clear();
  hit_wire_temp2.clear();
  hit_plane_temp2.clear();

  // Find locations of all ides in event
  auto xyzVec = SimIDEXYZ_Finder(evt,eventHits);
  fTrueIDEX.emplace_back(xyzVec.at(0));                                                                                          
  fTrueIDEY.emplace_back(xyzVec.at(1));
  fTrueIDEZ.emplace_back(xyzVec.at(2));


  // Get all SpacePoints in event                                                                   
  std::vector<art::Ptr<recob::SpacePoint>> spcPnts;
  art::fill_ptr_vector(spcPnts, spcPntHandle);

  for(const auto& sp : spcPnts)
    {
      const double *xyz = sp.get()->XYZ();
      xVec.emplace_back(xyz[0]);
      yVec.emplace_back(xyz[1]);
      zVec.emplace_back(xyz[2]);

    }

  fEventSpacePointX.emplace_back(xVec);
  fEventSpacePointY.emplace_back(yVec);
  fEventSpacePointZ.emplace_back(zVec);
  xVec.clear();
  yVec.clear();
  zVec.clear();

  // End of hit/SpacePoint/IDE section 
  
  //pandora shower info
  ShwrSpcPntXYZ_Finder(evt, fPFParticleTag, fShowerTag, fShowerSpacePointX, fShowerSpacePointY, fShowerSpacePointZ); // recob::SpacePoints
  auto recoShowers = evt.getValidHandle<std::vector<recob::Shower>>(fShowerTag);
  for ( auto s : *recoShowers) {
    auto showerHits = GetRecoShowerHits(s, evt, fShowerTag);
    fPFShowerNumHits.push_back(showerHits.size());   

    // sim::IDEs (backtracked)
    auto xyzVec = SimIDEXYZ_Finder(evt,showerHits);
    fShowerSimIDEX.emplace_back(xyzVec.at(0));
    fShowerSimIDEY.emplace_back(xyzVec.at(1));
    fShowerSimIDEZ.emplace_back(xyzVec.at(2));

    // hit info
    std::vector<double> showerHit_peaktime_temp0;
    std::vector<double> showerHit_wire_temp0;
    std::vector<double> showerHit_plane_temp0;
    std::vector<double> showerHit_peaktime_temp1;
    std::vector<double> showerHit_wire_temp1;
    std::vector<double> showerHit_plane_temp1;
    std::vector<double> showerHit_peaktime_temp2;
    std::vector<double> showerHit_wire_temp2;
    std::vector<double> showerHit_plane_temp2;

    for(const recob::Hit* hit: showerHits)
      {
	auto hit_peaktime = hit->PeakTime();
	auto hit_wire = hit->WireID().Wire;
	auto hit_plane = hit->WireID().Plane;

	if(hit_plane == 0)
	  {
	    showerHit_peaktime_temp0.emplace_back(hit_peaktime);
	    showerHit_wire_temp0.emplace_back(hit_wire);
	    showerHit_plane_temp0.emplace_back(hit_plane);
	  }
	if(hit_plane == 1)
	  {
	    showerHit_peaktime_temp1.emplace_back(hit_peaktime);
	    showerHit_wire_temp1.emplace_back(hit_wire);
	    showerHit_plane_temp1.emplace_back(hit_plane);
	  }
	if(hit_plane == 2)
	  {
	    showerHit_peaktime_temp2.emplace_back(hit_peaktime);
	    showerHit_wire_temp2.emplace_back(hit_wire);
	    showerHit_plane_temp2.emplace_back(hit_plane);
	  }	
      }

    fShowerHitPeakTime0.emplace_back(showerHit_peaktime_temp0);
    fShowerHitWire0.emplace_back(showerHit_wire_temp0);
    fShowerHitPlane0.emplace_back(showerHit_plane_temp0);
    fShowerHitPeakTime1.emplace_back(showerHit_peaktime_temp1);
    fShowerHitWire1.emplace_back(showerHit_wire_temp1);
    fShowerHitPlane1.emplace_back(showerHit_plane_temp1);
    fShowerHitPeakTime2.emplace_back(showerHit_peaktime_temp2);
    fShowerHitWire2.emplace_back(showerHit_wire_temp2);
    fShowerHitPlane2.emplace_back(showerHit_plane_temp2);
    showerHit_peaktime_temp0.clear();
    showerHit_wire_temp0.clear();
    showerHit_plane_temp0.clear();
    showerHit_peaktime_temp1.clear();
    showerHit_wire_temp1.clear();
    showerHit_plane_temp1.clear();
    showerHit_peaktime_temp2.clear();
    showerHit_wire_temp2.clear();
    showerHit_plane_temp2.clear();

    // other shower quantities
    auto direction = s.Direction();
    fPFShowerDirectionX.push_back(direction.X());
    fPFShowerDirectionY.push_back(direction.Y());
    fPFShowerDirectionZ.push_back(direction.Z());
    auto start = s.ShowerStart();
    fPFShowerStartX.push_back(start.X());
    fPFShowerStartY.push_back(start.Y());
    fPFShowerStartZ.push_back(start.Z());
    
    auto match = GetMCParticleListFromShowerHits(showerHits, evt, fHitTag)[0];
    fPFShowerBestMCParentMatchedHits.push_back(match.matched_hits);
    const simb::MCParticle* matchedShowerParticle = match.part;
    auto true_hits = getMCParticleHits(*matchedShowerParticle, evt, fHitTag).second;
    auto completeness = (true_hits.size() ? double(match.matched_hits) / double(true_hits.size()) : 0);
    auto purity = (showerHits.size() ? double(match.matched_hits) / double(showerHits.size()) : 0);
    fPFShowerBestMCParentPDG.push_back(matchedShowerParticle->PdgCode());
    fPFShowerBestMCParentID.push_back(matchedShowerParticle->TrackId());
    fPFShowerBestMCParentWeight.push_back(match.weight);
    fPFShowerBestMCParentCompleteness.push_back(completeness);
    fPFShowerBestMCParentPurity.push_back(purity);

    //the EstimateEnergyFromHitCharge function is only an approximation and should be replaced in the future with a function that gets both energy and momentum of the shower
    std::vector<double> energyvec = estimateEnergyFromHitCharge(evt, showerHits ,fCalorimetryAlg);    
    // I think i=2 is the collection plane?
    fPFShowerEnergy0.push_back(energyvec[0]);
    fPFShowerEnergy1.push_back(energyvec[1]);
    fPFShowerEnergyCollection.push_back(energyvec[2]);       
  }


  fPiZeroShowerTree->Fill();
  std::cout<<"=============================================================="<<std::endl;
}

void test::pi0EventDisplay::endJob()
{
  fPOTTree->Fill();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DEFINE_ART_MODULE(test::pi0EventDisplay)
