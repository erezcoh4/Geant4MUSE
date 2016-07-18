// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEEventAction.hh
/// \brief Definition of the Geant4MUSEEventAction class

#ifndef Geant4MUSEEventAction_h
#define Geant4MUSEEventAction_h 1

#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSERunAction.hh"
#include "Geant4MUSEAnalysis.hh"
#include "Geant4MUSETrackerHit.hh"
#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSETrack.hh"

#include "Geant4MUSEPrimaryGeneratorAction.hh"

#include "G4UserEventAction.hh"

#include "G4UnitsTable.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "TRandom2.h"

#include <vector>

#include "G4UserLimits.hh"

#include "G4SystemOfUnits.hh"
#include "Geant4MUSEParameters.hh"

/// Event action class

class Geant4MUSEEventAction : public G4UserEventAction
{

private:
    

public:


    Geant4MUSEParameters muse_par;
    
    Geant4MUSEEventAction();
    virtual ~Geant4MUSEEventAction();

    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );

    static const G4double GEMTimeWindowLength = 100.*ns;

    
    
    //global variables
    bool        PrintOutData;
    G4double    Threshold_Edep;
    
    
    // Generated Particle
    G4PrimaryParticle * primary_particle[100] ;
    G4int       primary_TrackID         [100] ;
    G4int       particle_ID_g           [100] ;
    G4double    primary_velocity        [100] ;
    G4double    p_g                     [100] ;
    G4double    E_kin_g                 [100] ;
    G4double    x_g                     [100] ;
    G4double    y_g                     [100] ;
    G4String    PrtclName_g             [100] ;
    G4double    pGenerated[3];
    
    

    
    
    // Detected Particle
    bool        ParticleHitDetector;
    bool        TrackIDWasDetected;
    G4String    DetectorName;
    G4String    PrtclName;
    G4double    GenerationTime;        // time that particle was deflected from sinchotron
    int         NPrimaryParticles;
    
    G4int    gem_time_window;
    G4int    ParentID;
    G4int    TrackID;
    G4double particle_ID;
    G4double Edep;
    G4double E_kin;
    G4double p;
    G4double PreTheta , PostTheta;
    G4double PrePhi , PostPhi;
    G4double x;
    G4double y;
    G4double z;
    G4double HitTime;
    G4ThreeVector VertexPosition;
    G4String VertexLVName;
    G4String CreatorProcessName;
    G4double PrimaryVelocity;
    

    // SciHod
    static const int Nplanes = 3;//muse_par.Nplanes;
    static const int Nbars   = 20;//muse_par.Nbars;
    G4double    SciHod_Edep          [Nplanes][100];
    G4int       SciHod_TrackID       [Nplanes][100];
    G4int       SciHod_HitFired      [Nplanes][100];
    G4int       SciHod_BarFired      [Nplanes][100];
    G4double    SciHod_particle_ID   [Nplanes][100];
    G4double    SciHod_p             [Nplanes][100];
    G4double    SciHod_E_kin         [Nplanes][100];
    G4ThreeVector SciHod_HitPos      [Nplanes][100];
    G4double    SciHod_x             [Nplanes][100];
    G4double    SciHod_y             [Nplanes][100];
    G4double    SciHod_z             [Nplanes][100];
    G4double    SciHod_HitTime       [Nplanes][100];
    G4double    SciHod_PreTheta      [Nplanes][100];
    G4double    SciHod_PostTheta     [Nplanes][100];
    G4double    SciHod_PrePhi        [Nplanes][100];
    G4double    SciHod_PostPhi       [Nplanes][100];

    
    

    // GEM
    static const int Ngem = 3;//muse_par.Ngem;
    G4double    GEM_Edep        [Ngem][100];
    G4int       GEM_TrackID     [Ngem][100];
    G4int       GEM_HitFired    [Ngem][100];
    G4double    GEM_particle_ID [Ngem][100];
    G4double    GEM_p           [Ngem][100];
    G4double    GEM_E_kin       [Ngem][100];
    G4ThreeVector GEM_HitPos    [Ngem][100];
    G4double    GEM_x           [Ngem][100];
    G4double    GEM_y           [Ngem][100];
    G4double    GEM_z           [Ngem][100];
    G4double    GEM_HitTime     [Ngem][100];

    
    
    // Target
    G4double    TARGET_Edep         [100];
    G4int       TARGET_TrackID      [100];
    G4int       TARGET_HitFired     [100];
    G4double    TARGET_particle_ID  [100];
    G4double    TARGET_p            [100];
    G4double    TARGET_E_kin        [100];
    G4ThreeVector TARGET_HitPos     [100];
    G4double    TARGET_x            [100];
    G4double    TARGET_y            [100];
    G4double    TARGET_z            [100];
    G4double    TARGET_HitTime      [100];

    
    
    
    // track counter
    size_t      MaxNumebrOfTracksDetected;
    int         track_ctr;                              // counter running over tracks detected in detector X
    size_t      track;                                  // a counter running over all tracks detected in either of the detectors
    TRandom2 *  rand2;

    std::vector<int>    DetectedTracksIDs;              // a vector that holds the trackID - s detected by SciFi
    std::vector<int>    DetectedTracksIDsSciHod[Nplanes]; // a vector that holds the trackID - s detected by SciHod
    std::vector<int>    DetectedTracksIDsGEM[Ngem];     // a vector that holds the trackID - s detected by GEM
    std::vector<int>    DetectedTracksIDsTarget;        // a vector that holds the trackID - s detected by Target

    
    
    
    // Methods
    void            InitializeEvent( const G4Event* );
    void               InitializeHC( );
    void              InitializeHit( );
    void           GetPrimariesData( const G4Event* );
    void     GetPrimaryParticleData( const G4Event* , int );
    void              GetTracksData( const G4Event* );
    void                 GetHitData( const Geant4MUSETrackerHit* );
    bool       PlugDataIntoDetector( G4String , int );
    
    bool           TrackWasDetected( G4String , int );
    void         SawTrackInDetector( G4String , int );
    bool     WasTrackSeenInDetector( G4String , int );
    void                   SawTrack( int );
    bool               WasTrackSeen( int );
    G4int       GetGeneratedTrackID( G4int );
    
    // prints
    void   PrintPrimaryParticleData( int );
    void        PrintTracksDetected( );
    void               PrintHitData();
    G4String              FiberName(int ,int );
    G4String              SciHodName(int ,int );
    G4String                GEMName(int);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
