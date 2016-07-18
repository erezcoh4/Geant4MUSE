// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEDetectorMessenger.hh
/// \brief Definition of the Geant4MUSEDetectorMessenger class

#ifndef Geant4MUSEDetectorMessenger_h
#define Geant4MUSEDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Geant4MUSEDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Messenger class that defines commands for Geant4MUSEDetectorConstruction.
///
/// It implements commands:
/// - /Geant4MUSE/det/setTargetMaterial name
/// - /Geant4MUSE/det/setChamberMaterial name
/// - /Geant4MUSE/det/stepMax value unit

class Geant4MUSEDetectorMessenger: public G4UImessenger
{
  public:
    Geant4MUSEDetectorMessenger(Geant4MUSEDetectorConstruction* );
    virtual ~Geant4MUSEDetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Geant4MUSEDetectorConstruction*  fDetectorConstruction;

    G4UIdirectory*           fGeant4MUSEDirectory;
    G4UIdirectory*           fDetDirectory;

    G4UIcmdWithAString*      fTargMatCmd;
    G4UIcmdWithAString*      fChamMatCmd;

    G4UIcmdWithADoubleAndUnit* fStepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
