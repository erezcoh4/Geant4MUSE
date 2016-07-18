// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEActionInitialization.cc
/// \brief Implementation of the Geant4MUSEActionInitialization class

#include "Geant4MUSEActionInitialization.hh"
#include "Geant4MUSEPrimaryGeneratorAction.hh"
#include "Geant4MUSERunAction.hh"
#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSEParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEActionInitialization::Geant4MUSEActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEActionInitialization::~Geant4MUSEActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEActionInitialization::BuildForMaster() const
{
  SetUserAction(new Geant4MUSERunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEActionInitialization::Build() const
{
    SetUserAction(new Geant4MUSEPrimaryGeneratorAction);
    SetUserAction(new Geant4MUSERunAction);
    SetUserAction(new Geant4MUSEEventAction);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
