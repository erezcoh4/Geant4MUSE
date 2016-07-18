// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrackerHit.hh
/// \brief Definition of the Geant4MUSETrackerHit class

#ifndef Geant4MUSETrackerHit_h
#define Geant4MUSETrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class Geant4MUSETrackerHit : public G4VHit
{
public:
    Geant4MUSETrackerHit();
    Geant4MUSETrackerHit(const Geant4MUSETrackerHit&);
    virtual ~Geant4MUSETrackerHit();
    
    // operators
    const Geant4MUSETrackerHit& operator=(const Geant4MUSETrackerHit&);
    G4int operator==(const Geant4MUSETrackerHit&) const;
    
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
    // methods from base class
    virtual void Draw();
    virtual void Print();
    
    // Set methods
    void SetDetectorName    (G4String det_name) { fDetectorName = det_name; };
    void SetTrackID         (G4int track)       { fTrackID = track; };
    void SetParentID        (G4int parent_id)   { fParentID = parent_id; };
    void SetChamberNb       (G4int chamb)       { fChamberNb = chamb; };
    void SetEdep            (G4double de)       { fEdep = de; };
    void SetPos             (G4ThreeVector xyz) { fPos = xyz; };
    void SetTotalEnergy     (G4double e_tot)    { fTotalEnergy = e_tot; };
    void SetParticleCode    (G4int p_code)      { fParticleCode = p_code; };
    void SetParticleName    (G4String p_name)   { fParticleName = p_name; };
    void SetKineticEnergy   (G4int e_kin)       { fKineticEnergy = e_kin; };
    void SetMomentum        (G4double p_momen)  { fMomentum = p_momen; };
    void SetTheta           (G4double p_theta)  { fTheta = p_theta; };
    void SetTime            (G4double hit_time) { fTime = hit_time; };
    void SetVelocity        (G4double vel)      { fVelocity = vel; };
    void SetVertexPosition  (G4ThreeVector ver) { fVertexPosition = ver; };
    void SetVertexLVName    (G4String lvname)   { fVertexLVName = lvname; };
    void SetCreatorProcessName(G4String cpname) { fCreatorProcessName = cpname; };
    void SetPreTheta        (G4double p_theta)  { fPreTheta = p_theta; };
    void SetPostTheta       (G4double p_theta)  { fPostTheta = p_theta; };
    void SetPrePhi          (G4double p_phi)    { fPrePhi = p_phi; };
    void SetPostPhi         (G4double p_phi)    { fPostPhi = p_phi; };
    
    // Get methods
    G4String GetDetectorName()  const { return fDetectorName; };
    G4int GetTrackID()          const { return fTrackID; };
    G4int GetParentID()         const { return fParentID; };
    G4int GetChamberNb()        const { return fChamberNb; };
    G4double GetEdep()          const { return fEdep; };
    G4ThreeVector GetPos()      const { return fPos; };
    G4double GetTotalEnergy()   const { return fTotalEnergy; };
    G4double GetKineticEnergy() const { return fKineticEnergy; };
    G4double GetMomentum()      const { return fMomentum; };
    G4double GetTheta()         const { return fTheta; };
    G4double GetX()             const { return fPos.x(); };
    G4double GetY()             const { return fPos.y(); };
    G4double GetZ()             const { return fPos.z(); };
    G4String GetParticleName()  const { return fParticleName; };
    G4double GetTime()          const { return fTime; };
    G4int GetParticleCode()     const { return fParticleCode; };
    G4double GetVelocity()      const { return fVelocity;};
    G4ThreeVector GetVertexPosition() const { return fVertexPosition;};
    G4String GetVertexLVName()  const { return fVertexLVName; };
    G4String GetCreatorProcessName()  const { return fCreatorProcessName; };
    G4double GetPreTheta()      const { return fPreTheta; };
    G4double GetPostTheta()     const { return fPostTheta; };
    G4double GetPrePhi()        const { return fPrePhi; };
    G4double GetPostPhi()       const { return fPostPhi; };
    
    private:

    G4String      fDetectorName;
    G4int         fTrackID;
    G4int         fParentID;
    G4double      fEdep;
    G4ThreeVector fPos;
    G4double      fTotalEnergy;
    G4double      fKineticEnergy;
    G4double      fMomentum;
    G4int         fParticleCode;
    G4String      fParticleName;
    G4int         fChamberNb;
    G4double      fTheta;
    G4double      fX;
    G4double      fY;
    G4double      fZ;
    G4double      fTime;
    G4double      fVelocity;
    G4ThreeVector fVertexPosition;
    G4String      fVertexLVName;
    G4String      fCreatorProcessName;
    G4double      fPreTheta;
    G4double      fPostTheta;
    G4double      fPrePhi;
    G4double      fPostPhi;
    };

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    typedef G4THitsCollection<Geant4MUSETrackerHit> Geant4MUSETrackerHitsCollection;
    
    extern G4ThreadLocal G4Allocator<Geant4MUSETrackerHit>* Geant4MUSETrackerHitAllocator;
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    inline void* Geant4MUSETrackerHit::operator new(size_t)
    {
        if(!Geant4MUSETrackerHitAllocator)
            Geant4MUSETrackerHitAllocator = new G4Allocator<Geant4MUSETrackerHit>;
        return (void *) Geant4MUSETrackerHitAllocator->MallocSingle();
    }
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    inline void Geant4MUSETrackerHit::operator delete(void *hit)
    {
        Geant4MUSETrackerHitAllocator->FreeSingle((Geant4MUSETrackerHit*) hit);
    }
    
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
#endif
