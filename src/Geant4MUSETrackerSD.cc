// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrackerSD.cc
/// \brief Implementation of the Geant4MUSETrackerSD class

#include "Geant4MUSETrackerSD.hh"
#include "Geant4MUSEAnalysis.hh"
#include "Geant4MUSEPrimaryGeneratorAction.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include <iomanip>
#include <G4VProcess.hh>
#include <G4VUserTrackInformation.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSETrackerSD::Geant4MUSETrackerSD(const G4String& name, const G4String& hitsCollectionName) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSETrackerSD::~Geant4MUSETrackerSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSETrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection 
    = new Geant4MUSETrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool Geant4MUSETrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    G4String DetectorName           = SensitiveDetectorName;
    Geant4MUSETrackerHit* newHit    = new Geant4MUSETrackerHit();
    G4int           Track_ID        = aStep -> GetTrack() -> GetTrackID();
    G4double        HitTime         = aStep -> GetPreStepPoint() -> GetGlobalTime();
    G4int           particle_ID     = aStep -> GetTrack() -> GetDefinition() -> GetPDGEncoding();
    G4String        particle_Name   = aStep -> GetTrack() -> GetDefinition() -> GetParticleName();
    G4double        E_kin           = aStep -> GetPostStepPoint() -> GetKineticEnergy();
    G4double        E_tot           = aStep -> GetPostStepPoint() -> GetTotalEnergy();
    G4double        p               = aStep -> GetPostStepPoint() -> GetMomentum() . mag();
    G4double        Edep            = aStep -> GetTotalEnergyDeposit();
    G4ThreeVector   Position        = aStep -> GetPreStepPoint() -> GetPosition();
    G4double        Theta           = (360./CLHEP::twopi)*(aStep -> GetPostStepPoint() -> GetMomentum()).getTheta();
    G4double        Velocity        = aStep -> GetPreStepPoint() -> GetVelocity();
    G4int           ParentTrack_ID  = aStep -> GetTrack() -> GetParentID();
    G4ThreeVector   VertexPosition  = aStep -> GetTrack() -> GetVertexPosition();
    G4String        VertexLVName    = aStep -> GetTrack() -> GetLogicalVolumeAtVertex() -> GetName();
    G4String        CreatorProcessName = (ParentTrack_ID != 0) ? aStep -> GetTrack() -> GetCreatorProcess() -> GetProcessName() : "Primary";
    G4double        PreTheta        = (360./CLHEP::twopi)*(aStep -> GetPreStepPoint() -> GetMomentum()).getTheta();
    G4double        PostTheta       = (360./CLHEP::twopi)*(aStep -> GetPostStepPoint() -> GetMomentum()).getTheta();
    G4double        PrePhi          = (360./CLHEP::twopi)*(aStep -> GetPreStepPoint() -> GetMomentum()).getPhi();
    G4double        PostPhi         = (360./CLHEP::twopi)*(aStep -> GetPostStepPoint() -> GetMomentum()).getPhi();
    
    
    newHit -> SetDetectorName(DetectorName);
    newHit -> SetParticleCode(particle_ID);
    newHit -> SetParentID(ParentTrack_ID);
    newHit -> SetParticleName(particle_Name);
    newHit -> SetKineticEnergy(E_kin);
    newHit -> SetMomentum(p);
    newHit -> SetTotalEnergy(E_tot);
    newHit -> SetTrackID(Track_ID);
    newHit -> SetEdep(Edep);
    newHit -> SetPos(Position);
    newHit -> SetTheta(Theta);
    newHit -> SetTime(HitTime);
    newHit -> SetVelocity(Velocity);
    newHit -> SetVertexPosition (VertexPosition);
    newHit -> SetVertexLVName (VertexLVName);
    newHit -> SetCreatorProcessName (CreatorProcessName);
    newHit -> SetPreTheta(PreTheta);
    newHit -> SetPostTheta(PostTheta);
    newHit -> SetPrePhi(PrePhi);
    newHit -> SetPostPhi(PostPhi);
    
    fHitsCollection->insert( newHit );
    
    
    if (Track_ID != 1)
        return false;
    else
        return true;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSETrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
//     G4int nofHits = fHitsCollection->entries();
//     G4cout << "\n-------->Hits Collection: in this event they are " << nofHits
//            << " hits in the tracker chambers: " << G4endl;
//     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
