// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEActionInitialization.hh
/// \brief Definition of the Geant4MUSEActionInitialization class

#ifndef Geant4MUSEActionInitialization_h
#define Geant4MUSEActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class B4DetectorConstruction;

/// Action initialization class.
///

class Geant4MUSEActionInitialization : public G4VUserActionInitialization
{
  public:
    Geant4MUSEActionInitialization();
    virtual ~Geant4MUSEActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;
};

#endif

    
