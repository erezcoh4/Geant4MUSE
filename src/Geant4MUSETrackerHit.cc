// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrackerHit.cc
/// \brief Implementation of the Geant4MUSETrackerHit class

#include "Geant4MUSETrackerHit.hh"
#include "Geant4MUSEAnalysis.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4ThreadLocal G4Allocator<Geant4MUSETrackerHit>* Geant4MUSETrackerHitAllocator=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSETrackerHit::Geant4MUSETrackerHit()
: G4VHit(),
fTrackID(-10000),
fParentID(-10000),
fEdep(-10000.),
fPos(G4ThreeVector(-10000.,-10000.,-10000.)),
fTotalEnergy(-10000.),
fKineticEnergy(-10000.),
fMomentum(-10000.),
fParticleCode(-10000),
fTheta(-10000.),
fX(-10000.),
fY(-10000.),
fZ(-10000.),
fVertexPosition(G4ThreeVector(-10000.,-10000.,-10000.))
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSETrackerHit::~Geant4MUSETrackerHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSETrackerHit::Geant4MUSETrackerHit(const Geant4MUSETrackerHit& right)
: G4VHit()
{
    fTrackID        = right.fTrackID;
    fParentID       = right.fParentID;
    fEdep           = right.fEdep;
    fPos            = right.fPos;
    fTotalEnergy    = right.fTotalEnergy;
    fKineticEnergy  = right.fKineticEnergy;
    fMomentum       = right.fMomentum;
    fParticleCode   = right.fParticleCode;
    fParticleName   = right.fParticleName;
    fTheta          = right.fTheta;
    fX              = right.fPos.x();
    fY              = right.fPos.y();
    fZ              = right.fPos.z();
    fTime           = right.fTime;
    fVelocity       = right.fVelocity;
    fVertexPosition = right.fVertexPosition;
    fVertexLVName   = right.fVertexLVName ;
    //fCreatorProcessName = right.fCreatorProcessName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const Geant4MUSETrackerHit& Geant4MUSETrackerHit::operator=(const Geant4MUSETrackerHit& right)
{
    fTrackID        = right.fTrackID;
    fParentID       = right.fParentID;
    fEdep           = right.fEdep;
    fPos            = right.fPos;
    fTotalEnergy    = right.fTotalEnergy;
    fKineticEnergy  = right.fKineticEnergy;
    fMomentum       = right.fMomentum;
    fParticleCode   = right.fParticleCode;
    fParticleName   = right.fParticleName;
    fTheta          = right.fTheta;
    fX              = right.fPos.x();
    fY              = right.fPos.y();
    fZ              = right.fPos.z();
    fTime           = right.fTime;
    fVelocity       = right.fVelocity;
    fVertexPosition = right.fVertexPosition;
    fVertexLVName   = right.fVertexLVName ;
    //fCreatorProcessName = right.fCreatorProcessName;
    
    return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int Geant4MUSETrackerHit::operator==(const Geant4MUSETrackerHit& right) const
{
    return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSETrackerHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
        G4Circle circle(fPos);
        circle.SetScreenSize(4.);
        circle.SetFillStyle(G4Circle::filled);
        G4Colour colour(1.,0.,0.);
        G4VisAttributes attribs(colour);
        circle.SetVisAttributes(attribs);
        pVVisManager->Draw(circle);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSETrackerHit::Print()
{
    G4cout << "\n---------------- Hit - start --------------------"
    << "\n Particle Code: " << fParticleCode
    << "\n Particle Name: " << fParticleName
    << "\n Momentum: "  << G4BestUnit(fMomentum , "Energy" )
    << "\n Kinetic Energy: " << std::setw(7) << G4BestUnit( fKineticEnergy , "Energy" )
    << "\n Theta =  " << fTheta
    << "\n X =  " << fX
    << "\n Y =  " << fY
    << "\n Z =  " << fZ
    << "\n-------------------- Hit - End ---------------------\n"
    << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
