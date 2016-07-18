// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEPrimaryGeneratorAction.hh
/// \brief Definition of the Geant4MUSEPrimaryGeneratorAction class

#ifndef Geant4MUSEPrimaryGeneratorAction_h
#define Geant4MUSEPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "Geant4MUSEParameters.hh"
#include <G4strstreambuf.hh>

class G4ParticleGun;
class G4Event;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the Tracker 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class Geant4MUSEPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:

    Geant4MUSEParameters muse_par;

    Geant4MUSEPrimaryGeneratorAction();
    virtual ~Geant4MUSEPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* );

    G4ParticleGun* GetParticleGun() {return fParticleGun;}
  
    // Set methods
    void SetRandomFlag(G4bool );
    
    G4double ParticleRate( G4int PID , G4double momentum );
    G4double pin;
    
    G4int pGenerated[3];
    G4double x0 , y0 , z0;
    G4double xp , yp;
    
    
    void SetParticleGunMomentumAndPosition();
    void  GeneratePrimaryVertexes(G4Event* anEvent);

  private:

    G4ParticleGun*          fParticleGun; // G4 particle gun
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
