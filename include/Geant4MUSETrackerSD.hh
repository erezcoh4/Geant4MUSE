// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrackerSD.hh
/// \brief Definition of the Geant4MUSETrackerSD class

#ifndef Geant4MUSETrackerSD_h
#define Geant4MUSETrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "Geant4MUSETrackerHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Geant4MUSETracker sensitive detector class
///
/// The hits are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step. A hit is created with each step with non zero 
/// energy deposit.

class Geant4MUSETrackerSD : public G4VSensitiveDetector
{
  public:
    Geant4MUSETrackerSD(const G4String& name, 
                const G4String& hitsCollectionName);
    virtual ~Geant4MUSETrackerSD();
  
    // methods from base class
    virtual void   Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

  private:
    Geant4MUSETrackerHitsCollection* fHitsCollection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
