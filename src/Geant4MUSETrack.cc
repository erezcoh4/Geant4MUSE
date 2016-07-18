// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSETrack.cc
/// \brief Implementation of the Geant4MUSETrack class

#include "Geant4MUSETrackerHit.hh"
#include "Geant4MUSEAnalysis.hh"
#include "Geant4MUSETrack.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Geant4MUSETrack::Geant4MUSETrack(){}
Geant4MUSETrack::~Geant4MUSETrack(){}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::InitializeTrack(){
    for ( int p = 0 ; p < Nplanes ; p++ ){
        fSciHodHitFired      [p] = 0;
        fSciHodBarFired      [p] = -1000;
        fSciHodHitPosition   [p] = G4ThreeVector();
        fSciHodHitTime       [p] = 0;
        fSciHodPreTheta      [p] = 0;
        fSciHodPostTheta     [p] = 0;
        fSciHodPrePhi        [p] = 0;
        fSciHodPostPhi       [p] = 0;
    }
    for ( int gem = 0 ; gem < Ngem ; gem++ ){
        fGEMHitFired      [gem] = 0;
        fGEMHitPosition   [gem] = G4ThreeVector();
        fGEMHitTime       [gem] = 0;
    }
    fTrackNumber        = 0;
    fTargetHitFired     = 0;
    fTargetHitPosition  = G4ThreeVector();
    fTargetHitTime      = 0;
    fGEMTimeWindow      = 0;
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::SetSciHodData(int p, G4int hit_fired, G4int bar_fired, G4ThreeVector hit_pos, G4double hit_time, G4double pre_theta, G4double post_theta, G4double pre_phi, G4double post_phi, G4int pID){
    
    fSciHodHitFired      [p] = hit_fired;
    fSciHodBarFired      [p] = bar_fired;
    fSciHodHitPosition   [p] = hit_pos;
    fSciHodHitTime       [p] = hit_time;
    fSciHodPreTheta      [p] = pre_theta;
    fSciHodPostTheta     [p] = post_theta;
    fSciHodPrePhi        [p] = pre_phi;
    fSciHodPostPhi       [p] = post_phi;
    fSciHodPID           [p] = pID;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::SetGEMData ( int gem , G4int hit_fired , G4ThreeVector hit_pos , G4double hit_time, G4int pID ){
    fGEMHitFired[gem]       = hit_fired;
    fGEMHitPosition[gem]    = hit_pos;
    fGEMHitTime[gem]        = hit_time;
    fGEMPID[gem]            = pID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::SetTargetData ( G4int hit_fired , G4ThreeVector hit_pos , G4double hit_time ){
    fTargetHitFired       = hit_fired;
    fTargetHitPosition    = hit_pos;
    fTargetHitTime        = hit_time;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::AddBeamFlightTime  ( G4double PrimaryVelocity ){// Primary flew in ~ 24 m beamline
//    for ( int plane = 0 ; plane < Nplane ; plane++ )
//        fSciFiHitTime[plane]    += 24.*CLHEP::m / PrimaryVelocity;
    for ( int p = 0 ; p < Nplanes ; p++ )
        fSciHodHitTime[p]    += 24.*CLHEP::m / PrimaryVelocity;
    for ( int gem = 0 ; gem < Ngem ; gem++ )
        fGEMHitTime[gem]    += 24.*CLHEP::m / PrimaryVelocity;
    fTargetHitTime          += 24.*CLHEP::m / PrimaryVelocity;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Geant4MUSETrack::TrackDetected(){ // used to delete a track if it didn't deposit enough energy in any detector
//    if (fTargetHitFired) return true; // This is for beam profile monitoring! only to keep all tracks from all events...
    int N_SciHod_Planes_Fired = 0 , N_GEMs_Fired = 0;
    for ( int p = 0 ; p < Nplanes ; p++ )
        N_SciHod_Planes_Fired += fSciHodHitFired[p];
    for ( int gem = 0 ; gem < Ngem ; gem++ )
        N_GEMs_Fired += fGEMHitFired[gem];
    if( (N_SciHod_Planes_Fired == Nplanes) ) // Trigger = all SciHod planes fired! // can also add if desired: && (N_GEMs_Fired==Ngem))
        return true;
    return false;
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::PrintTrackData(){
    G4cout << "-------------------------------------------------------- "
    << "\nTrack "                 << fTrackID << " (# " << fTrackNumber << " this event)"
    << "\n-------------------------------------------------------- \n"
    << "Generated Track "       << fGenTrackID
    << "\t Particle Code "      << fParticleID
    << "\t Parent Track "       << fParentTrackID
    << "\t Primary Generator "  << fPrimaryGenerator
    << "\n Vertex ("            << std::setw(7)     << G4BestUnit( fVertex , "Length"  )    << ") "
    << " at "                   << fVertexLVName
    << " ( "                    << fCreatorProcess  << " )"
    << " GEM Time window "      << fGEMTimeWindow  
    << G4endl;
    for ( int p = 0 ; p < Nplanes ; p++ )
        PrintSciHodData  ( p );
    for ( int gem = 0 ; gem < Ngem ; gem++ )
        PrintGEMData  ( gem );
    PrintTargetData  ( );
    G4cout << "-------------------------------------------------------- " << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::PrintSciHodData  ( int p ){
    if ( fSciHodHitFired[p] == 1)
        G4cout
        << "SciHod plane "   << p+1                  << " hit Bar "                                        << fSciHodBarFired[p]
        << " ("             << std::setw(2)         << G4BestUnit( fSciHodHitPosition[p] , "Length"  )       << ") "
        << " time "         << std::setw(2)         << G4BestUnit( fSciHodHitTime[p]     , "Time"  )
        << " theta: ("      << fSciHodPreTheta[p]    << " -> "                                               << fSciHodPostTheta[p]   << ")"
        << " phi: ("        << fSciHodPrePhi[p]      << " -> "                                               << fSciHodPostPhi[p]     << ")"
        << ", ID = "        << fSciHodPID[p]
        << G4endl;
    else
        G4cout << "SciHod plane " << p+1  << " did not fire "     << G4endl;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::PrintGEMData  ( int gem ){
    if ( fGEMHitFired[gem] == 1)
        G4cout
        << "GEM "       << gem+1            << " hit "
        << "("          << std::setw(2)     << G4BestUnit( fGEMHitPosition[gem] , "Length"  )   << ") "
        << " time "     << std::setw(2)     << G4BestUnit( fGEMHitTime[gem]     , "Time"  )
        << ", ID = "    << fGEMPID[gem]
        << G4endl;
    else
        G4cout << "GEM plane " << gem+1  << " did not fire "         << G4endl;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::PrintTargetData  ( ){
    if ( fTargetHitFired == 1)
        G4cout
        << "Target hit at ("    << std::setw(2)     << G4BestUnit( fTargetHitPosition , "Length"  )   << ") "
        << " time "             << std::setw(2)     << G4BestUnit( fTargetHitTime     , "Time"  )
        << G4endl;
    else
        G4cout << "Target did not fire "         << G4endl;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSETrack::FillTrackDataToROOT(){
    
    
    //----------STREAM DATA INTO DAQ-------------------//
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    // SciHod
    for ( int p = 0 ; p < Nplanes ; p++ ){
        man -> FillNtupleIColumn( 11*p+0    , fSciHodHitFired     [p] );
        man -> FillNtupleIColumn( 11*p+1    , fSciHodBarFired     [p] );
        man -> FillNtupleDColumn( 11*p+2    , fSciHodHitPosition  [p].x()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*p+3    , fSciHodHitPosition  [p].y()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*p+4    , fSciHodHitPosition  [p].z()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*p+5    , fSciHodHitTime      [p]       /CLHEP::ns  );
        man -> FillNtupleDColumn( 11*p+6    , fSciHodPreTheta     [p] );
        man -> FillNtupleDColumn( 11*p+7    , fSciHodPostTheta    [p] );
        man -> FillNtupleDColumn( 11*p+8    , fSciHodPrePhi       [p] );
        man -> FillNtupleDColumn( 11*p+9    , fSciHodPostPhi      [p] );
        man -> FillNtupleIColumn( 11*p+10   , fSciHodPID          [p] );
    }
    // GEM
    for ( int gem = 0 ; gem < Ngem ; gem++ ){
        man -> FillNtupleIColumn( 11*Nplanes+6*gem+0    , fGEMHitFired[gem] );
        man -> FillNtupleDColumn( 11*Nplanes+6*gem+1    , fGEMHitPosition[gem].x()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*Nplanes+6*gem+2    , fGEMHitPosition[gem].y()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*Nplanes+6*gem+3    , fGEMHitPosition[gem].z()   /CLHEP::mm  );
        man -> FillNtupleDColumn( 11*Nplanes+6*gem+4    , fGEMHitTime[gem]           /CLHEP::ns  );
        man -> FillNtupleIColumn( 11*Nplanes+6*gem+5    , fGEMPID[gem]                  );
    }
    // Target
    man -> FillNtupleIColumn( 11*Nplanes+6*Ngem+0    , fTargetHitFired );
    man -> FillNtupleDColumn( 11*Nplanes+6*Ngem+1    , fTargetHitPosition.x()   /CLHEP::mm  );
    man -> FillNtupleDColumn( 11*Nplanes+6*Ngem+2    , fTargetHitPosition.y()   /CLHEP::mm  );
    man -> FillNtupleDColumn( 11*Nplanes+6*Ngem+3    , fTargetHitPosition.z()   /CLHEP::mm  );
    man -> FillNtupleDColumn( 11*Nplanes+6*Ngem+4    , fTargetHitTime           /CLHEP::ns  );
    
    // Global track variables
    man -> FillNtupleIColumn( 11*Nplanes+6*Ngem+5    , fTrackID );
    man -> FillNtupleIColumn( 11*Nplanes+6*Ngem+6    , fGenTrackID );
    man -> FillNtupleIColumn( 11*Nplanes+6*Ngem+7    , fParticleID );
    man -> FillNtupleIColumn( 11*Nplanes+6*Ngem+8    , fGEMTimeWindow );

    
    man -> AddNtupleRow();
}







