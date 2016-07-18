// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrack.hh
/// \brief Definition of the Geant4MUSETrack class

#ifndef Geant4MUSETrack_h
#define Geant4MUSETrack_h 1


#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSERunAction.hh"
#include "Geant4MUSEAnalysis.hh"
#include "Geant4MUSETrackerHit.hh"
#include "Geant4MUSEEventAction.hh"

#include "Geant4MUSEPrimaryGeneratorAction.hh"

#include "G4UserEventAction.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "TRandom2.h"

#include <vector>
#include "G4UserRunAction.hh"
#include "G4Event.hh"
#include "globals.hh"
#include "G4Track.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// My own clas to control a track running along MUSE beam line
class Geant4MUSETrack : public G4Track
{
    
public:
    Geant4MUSETrack();
    virtual ~Geant4MUSETrack();
    
    
    // Set methods
    void SetTrackNumber     (G4int track)       { fTrackNumber = track; };
    void SetTrackID         (G4int track_id)    { fTrackID = track_id; };
    void SetGenTrackID      (G4int track_id_g)  { fGenTrackID = track_id_g; };
    void SetParentTrackID   (G4int par_track_id){ fParentTrackID = par_track_id; };
    void SetVertex          (G4ThreeVector ver) { fVertex = ver; };
    void SetVertexLVName    (G4String ver_lv)   { fVertexLVName = ver_lv; };
    void SetProcessName     (G4String process)  { fCreatorProcess = process; };
    void SetParticleID      (G4int par_id)      { fParticleID = par_id; };
    void SetPrimaryGenerator(G4int prim_gen)    { fPrimaryGenerator = prim_gen; };
    void SetGEMTimeWindow   (G4int time_win)    { fGEMTimeWindow = time_win; };
    
    void InitializeTrack    ();
//    void SetSciFiData       ( int, G4int, G4int, G4ThreeVector, G4double , G4double , G4double , G4double , G4double );
    void SetSciHodData       ( int, G4int, G4int, G4ThreeVector, G4double , G4double , G4double , G4double , G4double , G4int);
    void SetGEMData         ( int, G4int, G4ThreeVector, G4double , G4int);
    void SetTargetData      ( G4int, G4ThreeVector, G4double );
    void AddBeamFlightTime  ( G4double );
    
    
    // Get methods
    G4int           GetTrackNumber()        const { return fTrackNumber; }
    G4int           GetTrackID()            const { return fTrackID; }
    G4int           GetGenTrackID()         const { return fGenTrackID; }
    G4int           GetParentTrackID()      const { return fParentTrackID; }
    G4ThreeVector   GetVertex()             const { return fVertex; }
    G4String        GetVertexLVName()       const { return fVertexLVName; }
    G4String        GetProcessName()        const { return fCreatorProcess; }
    G4int           GetParticleID()         const { return fParticleID; }
    G4int           GetPrimaryGenerator()   const { return fPrimaryGenerator; }
    G4int           GetGEMTimeWindow()      const { return fGEMTimeWindow; }

    
    
    
    G4int           fGEMTimeWindow;
    G4int           fTrackNumber;
    G4int           fTrackID;
    G4int           fGenTrackID;
    G4int           fParentTrackID;
    G4ThreeVector   fVertex;
    G4String        fVertexLVName;
    G4String        fCreatorProcess;
    G4int           fParticleID;
    G4int           fPrimaryGenerator;
    
    // SciHod
    static const int Nplanes = 3;
    G4int           fSciHodHitFired      [Nplanes];
    G4int           fSciHodBarFired      [Nplanes];
    G4ThreeVector   fSciHodHitPosition   [Nplanes];
    G4double        fSciHodHitTime       [Nplanes];
    G4double        fSciHodPreTheta      [Nplanes];
    G4double        fSciHodPostTheta     [Nplanes];
    G4double        fSciHodPrePhi        [Nplanes];
    G4double        fSciHodPostPhi       [Nplanes];
    G4int           fSciHodPID           [Nplanes];
    
    // GEM
    static const int Ngem = 3;
    G4int           fGEMHitFired      [Ngem];
    G4ThreeVector   fGEMHitPosition   [Ngem];
    G4double        fGEMHitTime       [Ngem];
    G4int           fGEMPID           [Ngem];
    

    // Target
    G4int           fTargetHitFired;
    G4ThreeVector   fTargetHitPosition;
    G4double        fTargetHitTime;
    
    
    bool            TrackDetected       ();
    
    // prints
    void            PrintTrackData      ();
    void            PrintSciFiData      ( int );
    void            PrintSciHodData      ( int );
    void            PrintGEMData        ( int );
    void            PrintTargetData     ( );

    
    // Fill root tree
    void            FillTrackDataToROOT ();
    
    
private:
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


