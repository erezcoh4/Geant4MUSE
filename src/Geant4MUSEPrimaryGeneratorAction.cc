
// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEPrimaryGeneratorAction.cc
/// \brief Implementation of the Geant4MUSEPrimaryGeneratorAction class

#include "Geant4MUSEPrimaryGeneratorAction.hh"
#include "Geant4MUSETrackerSD.hh"
#include "Geant4MUSEAnalysis.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Poisson.hh"

#include "Randomize.hh"
#include "G4UnitsTable.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEPrimaryGeneratorAction::Geant4MUSEPrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction()
{
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);
    for (int p = 0 ; p < 3 ; p++)
        pGenerated[p] = 0;
    pin = muse_par.momentum * (1 + CLHEP::RandFlat::shoot(-0.5 , 0.5)*muse_par.momentumSpread);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEPrimaryGeneratorAction::~Geant4MUSEPrimaryGeneratorAction()
{
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{   // This function is called at the begining of event
    G4double worldZHalfLength = 0;
    G4LogicalVolume* worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    G4Box* worldBox = NULL;
    if ( worldLV ) worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
    if ( worldBox ) worldZHalfLength = worldBox->GetZHalfLength();
    else
        G4cerr << "World volume of box not found." << "Perhaps you have changed geometry." << "The gun will be place in the center." << G4endl;
    
    long seeds[2] = {(long)(1e6*G4UniformRand()),(long)(1e6*G4UniformRand())};
    CLHEP::HepRandom::setTheSeeds(seeds);
    SetParticleGunMomentumAndPosition();
    GeneratePrimaryVertexes(anEvent);
}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEPrimaryGeneratorAction::SetParticleGunMomentumAndPosition(){
    G4double z0_beamline = -150.*cm;
    
    
    if (muse_par.BeamType == "Point beam") {
        x0 = 0; xp=0;
        y0 = 0; yp=0;
        z0 = -10.*cm;
    }

    
    else if (muse_par.BeamType == "beam 2013") {
        z0 = z0_beamline;
        G4double sig_xp = 0;    G4double sig_yp = 0;
        G4double sig_x0 = 0;    G4double sig_y0 = 0;
        G4double fpx =  0.; G4double fpy =  0.;
        if (pin < 130.0 * MeV) {
            sig_xp = 16./1000.; sig_yp = 10./1000.;
            sig_x0 = 0.2 * cm;  sig_y0 = 0.7 * cm;
            fpx =  50. * cm;    fpy =  30. * cm;
        } else if (pin < 180.0 * MeV) {
            sig_xp = 14./1000.; sig_yp = 10./1000.;
            sig_x0 = 0.2 * cm;  sig_y0 = 0.8 * cm;
            fpx =  50. * cm;    fpy =  0. * cm;
        } else {
            sig_xp = 12./1000.; sig_yp = 10./1000.;
            sig_x0 = 0.2 * cm;  sig_y0 = 0.7 * cm;
            fpx =  30. * cm;    fpy =  0. * cm;
        }
        xp = CLHEP::RandGauss::shoot(0, sig_xp);
        yp = CLHEP::RandGauss::shoot(0, sig_yp);
        x0 = CLHEP::RandGauss::shoot(0, sig_x0) + (z0 - fpx) * xp;
        y0 = CLHEP::RandGauss::shoot(0, sig_y0) + (z0 - fpy) * yp;
    }
    
    
    
    else if (muse_par.BeamType == "beam 2015") {
        z0 = z0_beamline;
        G4double sig_xp     = 0;
        G4double sig_yp     = 0;
        G4double sig_x0     = 0;
        G4double sig_y0     = 0;
        G4double mean_xp    = 0;
        G4double mean_yp    = 0;
        G4double mean_x0    = 0;
        G4double mean_y0    = 0;
        G4double fpx        = 0.;
        G4double fpy        = 0.;
        G4double charge     = -1 ; //particle->GetPDGCharge();
        if (charge <= 0) {// negatively charged particles, or neutrals, geantino
            if (pin < 130.0 * MeV) {
                mean_xp = -8./1000.;    mean_yp = +8./1000.;
                mean_x0 = 3. * mm;      mean_y0 = 2. * mm;
                sig_xp = 18./1000.;     sig_yp =  6./1000.;
                sig_x0 = 0. * mm;       sig_y0 = 6. * mm;
                fpx =  0. * mm;         fpy =  -300. * mm;
            } else if (pin < 180.0 * MeV) {
                mean_xp = -8./1000.;    mean_yp = +4./1000.;
                mean_x0 = 2. * mm;      mean_y0 = 1. * mm;
                sig_xp = 15./1000.;     sig_yp = 10./1000.;
                sig_x0 = 2. * mm;       sig_y0 = 5. * mm;
                fpx =  -20. * mm;       fpy =  -300. * mm;
            } else {
                mean_xp = -2./1000.;    mean_yp = +4./1000.;
                mean_x0 = 0. * mm;      mean_y0 = -1. * mm;
                sig_xp = 10./1000.;     sig_yp = 7./1000.;
                sig_x0 = 2. * mm;       sig_y0 = 4. * mm;
                fpx =  -20. * mm;       fpy =  -400. * mm;
            }
        } else { // positively charged particles
            if (pin < 130.0 * MeV) {
                mean_xp = -8./1000.;    mean_yp = +0./1000.;
                mean_x0 = 1.5 * mm;     mean_y0 = 0. * mm;
                sig_xp = 24./1000.;     sig_yp =  6./1000.;
                sig_x0 = 2. * mm;       sig_y0 = 6. * mm;
                fpx =  0. * mm;         fpy =  -300. * mm;
            } else if (pin < 180.0 * MeV) {
                mean_xp = -8./1000.;    mean_yp = +0./1000.;
                mean_x0 = 1. * mm;      mean_y0 = 0. * mm;
                sig_xp = 16./1000.;     sig_yp = 10./1000.;
                sig_x0 = 3. * mm;       sig_y0 = 5. * mm;
                fpx =  -20. * mm;       fpy =  -300. * mm;
            } else {
                mean_xp = +0./1000.;    mean_yp = +2./1000.;
                mean_x0 = -1. * mm;     mean_y0 = -1. * mm;
                sig_xp = 14./1000.;     sig_yp = 7./1000.;
                sig_x0 = 2. * mm;       sig_y0 = 4. * mm;
                fpx =  -20. * mm;       fpy =  -400. * mm;
            }
        }
        xp = CLHEP::RandGauss::shoot(mean_xp, sig_xp);
        yp = CLHEP::RandGauss::shoot(mean_yp, sig_yp);
        x0 = CLHEP::RandGauss::shoot(mean_x0, sig_x0) + (z0 - fpx) * xp;
        y0 = CLHEP::RandGauss::shoot(mean_y0, sig_y0) + (z0 - fpy) * yp;
    }
    
    
    
    // set momentum...
    // ---------------------------------
    fParticleGun -> SetParticlePosition(G4ThreeVector(x0, y0, z0));
    fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(xp, yp, 1.));
    G4strstreambuf* oldBuffer = dynamic_cast<G4strstreambuf*>(G4cout.rdbuf(0));
    fParticleGun -> SetParticleMomentum( pin );
    G4cout.rdbuf(oldBuffer);
}







//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEPrimaryGeneratorAction::GeneratePrimaryVertexes(G4Event* anEvent){
    
    
    if (muse_par.BeamParticleType == "One e- at a time") {
        //   One simple particle (e-) at a time
        // ---------------------------------
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle;
        particle = particleTable -> FindParticle ("e-") ;
        fParticleGun -> SetParticleDefinition(particle);
        fParticleGun -> GeneratePrimaryVertex(anEvent);
        pGenerated[0] ++;
    }
    
    
    
    
    else if (muse_par.BeamParticleType == "TDR 3 particles"){
        //   Mixed beam pion/muon/electron
        // ---------------------------------
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        G4ParticleDefinition* particle;
        if (anEvent->GetEventID()%muse_par.PrintEachEvent==0){
            G4cout << "************************************************************************************"<< G4endl;
            G4cout << " Starting event "<< anEvent->GetEventID() <<"  with " << G4endl;
        }
        
        G4double    BeamRate        = muse_par.fBeamRate;                   // The 50 MHz is fixed by the proton accelerator
        G4int       fNBunch[3];
        for (int ParticleType = 0 ; ParticleType < 3 ; ParticleType++ ){
            switch( ParticleType ) { // generate an e / µ / π bunch
                case 0:     particle = particleTable -> FindParticle ("e-") ;   break;
                case 1:     particle = particleTable -> FindParticle ("mu-");   break;
                default:    particle = particleTable -> FindParticle ("pi-");   break;
            }
            fParticleGun->SetParticleDefinition(particle);
            
            // determine how many particles were in the bunch (fNBunch events in each bunch)
            G4double    mean        = ParticleRate( fParticleGun -> GetParticleDefinition() -> GetPDGEncoding() , pin ) / BeamRate;
            fNBunch[ParticleType]   = G4Poisson( mean );
            pGenerated[ParticleType] += fNBunch[ParticleType];
            
            if (anEvent->GetEventID()%muse_par.PrintEachEvent==0)
                G4cout << " " << fNBunch[ParticleType] << " "  << fParticleGun -> GetParticleDefinition() -> GetParticleName();
            
            for (int particle_ctr = 0 ; particle_ctr < fNBunch[ParticleType] ; particle_ctr++ )
                fParticleGun -> GeneratePrimaryVertex(anEvent);
            
        }
        if (anEvent->GetEventID()%muse_par.PrintEachEvent==0){
            G4cout << "\n in Total: " << pGenerated[0] << " e-,\t" << pGenerated[1] << " mu-,\t" << pGenerated[2] << " pi-" << G4endl;
            G4cout << "*************************************************"<< G4endl;
        }
    }
}








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double Geant4MUSEPrimaryGeneratorAction::ParticleRate( G4int PID , G4double momentum ){
    // particle rate for 5MHz limitation, taken from TDR[Gilman et al. PSI R-12-01.1 (Nov 2014)]
    
    G4double fParticleRate = 0.0*CLHEP::hertz;
    switch ((int)(momentum/CLHEP::MeV)){
        case (115): {// 115 MeV/c
            switch (PID) {
                case 11:    {fParticleRate = 4.93e6*CLHEP::hertz;    break;}      // e
                case 13:    {fParticleRate = 0.05e6*CLHEP::hertz;    break;}      // mu
                default:    {fParticleRate = 0.03e6*CLHEP::hertz;    break;}      // pi
            }
            break;
        }
        case (153): {//153 MeV/c
            switch (PID) {
                case 11:    {fParticleRate = 4.50e6*CLHEP::hertz;    break;}      // e
                case 13:    {fParticleRate = 0.16e6*CLHEP::hertz;    break;}      // mu
                default:    {fParticleRate = 0.34e6*CLHEP::hertz;    break;}      // pi
            }
            break;
        }
        case (160): {// 160 MeV/c
            switch (PID) {
                case 11:    {fParticleRate = 2.2e6*CLHEP::hertz;    break;}      // e
                case 13:    {fParticleRate = 0.7e6*CLHEP::hertz;    break;}      // mu
                default:    {fParticleRate = 2.1e6*CLHEP::hertz;    break;}      // pi
            }
            break;
        }
        default:    {// 210 MeV/c
            switch (PID) {
                case 11:    {fParticleRate = 2.35e6*CLHEP::hertz;    break;}      // e
                case 13:    {fParticleRate = 0.20e6*CLHEP::hertz;    break;}      // mu
                default:    {fParticleRate = 2.45e6*CLHEP::hertz;    break;}      // pi
            }
            break;
        }
    }
    return fParticleRate;
}


