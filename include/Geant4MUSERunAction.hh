// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSERunAction.hh
/// \brief Definition of the Geant4MUSERunAction class

#ifndef Geant4MUSERunAction_h
#define Geant4MUSERunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "Geant4MUSEParameters.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// Run action class

class Geant4MUSERunAction : public G4UserRunAction
{
private:

public:
    Geant4MUSEParameters muse_par;
   
    Geant4MUSERunAction();
    virtual ~Geant4MUSERunAction();

    virtual void BeginOfRunAction(const G4Run* run);
    virtual void EndOfRunAction(const G4Run* run);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
