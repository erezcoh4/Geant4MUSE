// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEEventAction.cc
/// \brief Implementation of the Geant4MUSEEventAction class

#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSERunAction.hh"
#include "Geant4MUSEAnalysis.hh"
#include "Geant4MUSETrackerHit.hh"
#include "Geant4MUSEEventAction.hh"
#include "Geant4MUSETrack.hh"

#include "Geant4MUSEPrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include <G4VHit.hh>
#include <G4VHitsCollection.hh>

#include "G4UnitsTable.hh"
#include <iomanip>
#include "Randomize.hh"
#include <vector>



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEEventAction::Geant4MUSEEventAction()
: G4UserEventAction()
{   //units: http://www-zeuthen.desy.de/geant4/geant4.9.3.b01/G4SystemOfUnits_8hh-source.html
    Threshold_Edep      = 1.*keV;
    rand2               = new TRandom2();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Geant4MUSEEventAction::~Geant4MUSEEventAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::BeginOfEventAction(const G4Event* event){
    InitializeEvent(event);
    if (PrintOutData)
        G4cout << "\n********** Event: " << event->GetEventID()
        << " Generation time " << GenerationTime << " ns, "
        << " GEM time window " << gem_time_window
        << " ****************" << G4endl;
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::InitializeEvent( const G4Event* event){
    GenerationTime  = (G4double)event->GetEventID()/muse_par.fBeamRate;        // time that particle was deflected from sinchotron
    PrintOutData    = ((event->GetEventID()%muse_par.PrintEachEvent==0) && (event->GetNumberOfPrimaryVertex()>0)) ? true : false; // print every x events
    gem_time_window = ((G4int)GenerationTime/(G4int)GEMTimeWindowLength);
    for (int p = 0 ; p < Nplanes ; p++ )
        if (!DetectedTracksIDsSciHod[p].empty())
            DetectedTracksIDsSciHod[p].erase(DetectedTracksIDsSciHod[p].begin(),DetectedTracksIDsSciHod[p].end());
    for (int gem = 0 ; gem < Ngem ; gem++)
        if (!DetectedTracksIDsGEM[gem].empty())
            DetectedTracksIDsGEM[gem].erase(DetectedTracksIDsGEM[gem].begin(),DetectedTracksIDsGEM[gem].end());
    if (!DetectedTracksIDsTarget.empty())
        DetectedTracksIDsTarget.erase(DetectedTracksIDsTarget.begin(),DetectedTracksIDsTarget.end());
    if (!DetectedTracksIDs.empty())
        DetectedTracksIDs.erase(DetectedTracksIDs.begin(),DetectedTracksIDs.end());
    
    for (track = 0 ; track < 100 ; track++){
        for (int p = 0 ; p < Nplanes ; p++){
            SciHod_Edep          [p][track]=-1000;
            SciHod_BarFired      [p][track]=-1000;
            SciHod_TrackID       [p][track]=-1000;
            SciHod_HitFired      [p][track]=-1000;
            SciHod_particle_ID   [p][track]=-1000;
            SciHod_p             [p][track]=-1000;
            SciHod_E_kin         [p][track]=-1000;
            SciHod_HitPos        [p][track]=G4ThreeVector(-1000,-1000,-1000);
            SciHod_x             [p][track]=-1000;
            SciHod_y             [p][track]=-1000;
            SciHod_z             [p][track]=-1000;
            SciHod_HitTime       [p][track]=-1000;
            SciHod_PreTheta      [p][track]=-1000;
            SciHod_PostTheta     [p][track]=-1000;
            SciHod_PrePhi        [p][track]=-1000;
            SciHod_PostPhi       [p][track]=-1000;
        }
        
        for (int gem = 0 ; gem < Ngem ; gem++){
            GEM_Edep        [gem][track]=-1000;
            GEM_TrackID     [gem][track]=-1000;
            GEM_HitFired    [gem][track]=-1000;
            GEM_particle_ID [gem][track]=-1000;
            GEM_p           [gem][track]=-1000;
            GEM_E_kin       [gem][track]=-1000;
            GEM_HitPos      [gem][track]=G4ThreeVector(-1000,-1000,-1000);
            GEM_x           [gem][track]=-1000;
            GEM_y           [gem][track]=-1000;
            GEM_z           [gem][track]=-1000;
            GEM_HitTime     [gem][track]=-1000;
        }
        TARGET_Edep         [track]=-1000;
        TARGET_TrackID      [track]=-1000;
        TARGET_HitFired     [track]=-1000;
        TARGET_particle_ID  [track]=-1000;
        TARGET_p            [track]=-1000;
        TARGET_E_kin        [track]=-1000;
        TARGET_HitPos       [track]=G4ThreeVector(-1000,-1000,-1000);
        TARGET_x            [track]=-1000;
        TARGET_y            [track]=-1000;
        TARGET_z            [track]=-1000;
        TARGET_HitTime      [track]=-1000;
     }
    track = 0;  // a counter running over all tracks detected in either of the detectors
 }




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::InitializeHC(){
    track_ctr   = 0;  //counter running over detected tracks
    TrackID     = -1;
}






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::InitializeHit(){
    PrtclName = DetectorName   = "";
    ParticleHitDetector = false;
    // Initialize Track Hit Data
    Edep            = 0;
    particle_ID     = p  = E_kin  = -10000;
    x               = y  = z      = HitTime = -10000;
    TrackID         = -1;
}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEEventAction::EndOfEventAction(const G4Event* event)
{
    NPrimaryParticles = event -> GetNumberOfPrimaryVertex();
    if (NPrimaryParticles > 0){
        GetPrimariesData( event );
        GetTracksData( event );
    }
    if ( PrintOutData ){
        G4cout << "************************************************************************************"<< G4endl;
        G4cout << " Finished event "<< event->GetEventID() <<"  with " << G4endl;
        for (int i = 0 ; i < NPrimaryParticles ; i++ )
            G4cout << " " << PrtclName_g[i] << "\t" ;
        G4cout << G4endl;
        G4cout << "in total \t" << pGenerated[0] << " e-,\t" << pGenerated[1] << " mu-,\t" << pGenerated[2] << " pi-" << G4endl;
        G4cout << "\n************************************************************************************\n"<< G4endl;
    }
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::GetPrimariesData( const G4Event* event ){    //get generated event data
    for ( int primary_ctr = 0 ; primary_ctr < NPrimaryParticles ; primary_ctr ++ ){
        GetPrimaryParticleData( event , primary_ctr );
        PrintPrimaryParticleData( primary_ctr );
        if (event->GetEventID()==0)
            pGenerated[0] = pGenerated[1] = pGenerated[2] = 0;
        switch (particle_ID_g[primary_ctr]) {
            case 11: // e
                pGenerated[0]++;
                break;
            case 13: // µ
                pGenerated[1]++;
                break;
            default: // π
                pGenerated[2]++;
                break;
        }
    }
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::GetPrimaryParticleData( const G4Event* event , int primary_ctr ){
    primary_particle    [primary_ctr]       = event -> GetPrimaryVertex(primary_ctr) -> GetPrimary();
    x_g                 [primary_ctr]       = event -> GetPrimaryVertex() -> GetX0();
    y_g                 [primary_ctr]       = event -> GetPrimaryVertex() -> GetY0();
    particle_ID_g       [primary_ctr]       = primary_particle[primary_ctr] -> GetPDGcode();
    primary_TrackID     [primary_ctr]       = primary_particle[primary_ctr] -> GetTrackID();
    E_kin_g             [primary_ctr]       = primary_particle[primary_ctr] -> GetKineticEnergy();
    PrtclName_g         [primary_ctr]       = primary_particle[primary_ctr] -> GetG4code() -> GetParticleName();
    p_g                 [primary_ctr]       = (primary_particle[primary_ctr]-> GetMomentum()).mag();
    G4double mass                           = primary_particle[primary_ctr] -> GetMass();
    G4double T                              = E_kin_g[primary_ctr] / mass;
    primary_velocity    [primary_ctr]       = CLHEP::c_light*std::sqrt(T*(T+2.))/(T+1.0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::PrintPrimaryParticleData( int primary_ctr ){
    if (PrintOutData) {
        G4cout
        << "\n Primary:"
        << "\n Track ID: "                  << primary_TrackID[primary_ctr]
        << "\n Particle Code: "             << particle_ID_g[primary_ctr]
        << "\n Particle Name: "             << PrtclName_g[primary_ctr]
        << "\n Primary velocity: "          << primary_velocity[primary_ctr] / (CLHEP::m/CLHEP::s) << " m/s "
        << "\n Generated Momentum: "        << G4BestUnit( p_g[primary_ctr] , "Energy" )
        << "\n Generated Kinetic Energy: "  << G4BestUnit( E_kin_g[primary_ctr] , "Energy" )
        << "\n Generated x: "               << G4BestUnit( x_g[primary_ctr] , "Length" )
        << "\n Generated y: "               << G4BestUnit( y_g[primary_ctr] , "Length" )
        << G4endl;
    }
}








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::GetTracksData(const G4Event* event){
    
    // Detected Tracks
    Geant4MUSETrack MUSE_track[100];
    for ( int j = 0 ; j < 100 ; j++ )
        MUSE_track[j].InitializeTrack();

    for (int hitsCollection = 0 ; hitsCollection < event->GetHCofThisEvent()->GetNumberOfCollections() ; hitsCollection++ ){
        G4VHitsCollection* hc = event -> GetHCofThisEvent() -> GetHC( hitsCollection );

        if (hitsCollection <= Ngem                  // gem planes
            || hitsCollection==Ngem+1
            || hitsCollection==Ngem+Nbars
            || hitsCollection==Ngem+2*Nbars
            || hitsCollection==Ngem+3*Nbars         //SS planes
            || hitsCollection==Ngem+Nplanes*Nbars ) // target detector
            InitializeHC();
        
        if (PrintOutData && hc->GetSize()>0)
            G4cout << "\n------\n" << hc -> GetName() << "\t hits collection" << ", Size = " << hc -> GetSize() << G4endl;
        
        for (size_t hit = 0 ; hit < hc -> GetSize() ; hit++){
            InitializeHit();
            Geant4MUSETrackerHit* MUSE_Hit = (Geant4MUSETrackerHit*)(event->GetHCofThisEvent()->GetHC(hitsCollection)->GetHit(hit));
            GetHitData( MUSE_Hit );
            
            // Plug data into detectors - each hit collection is streamed into the corresponding detector
            if (!WasTrackSeenInDetector( DetectorName , TrackID )){//{!TrackWasDetected( DetectorName , TrackID )
                if (PlugDataIntoDetector( DetectorName , track_ctr )){
                    SawTrackInDetector( DetectorName , TrackID );
                    track_ctr++ ;
                }
                if (PrintOutData)
                    PrintHitData();
            }
            
            // determine all the tracks that were detected in this event
            if (!WasTrackSeen( TrackID )){
                SawTrack( TrackID );
                MUSE_track[track].SetTrackNumber    ( track );
                MUSE_track[track].SetTrackID        ( TrackID );
                MUSE_track[track].SetGenTrackID     ( GetGeneratedTrackID(TrackID) );
                MUSE_track[track].SetParentTrackID  ( ParentID );
                MUSE_track[track].SetVertex         ( VertexPosition );
                MUSE_track[track].SetVertexLVName   ( VertexLVName );
                MUSE_track[track].SetProcessName    ( CreatorProcessName );
                MUSE_track[track].SetParticleID     ( particle_ID );
                MUSE_track[track].SetGEMTimeWindow  ( gem_time_window );
                track++;
            }
            
        }
    }
    
    if (PrintOutData) PrintTracksDetected();
    
    // Plug data into tracks - each track gets the track detection data in all detectors
    for ( track = 0 ; track < DetectedTracksIDs.size() ; track++){
        for (int p = 0 ; p < Nplanes ; p++){
            for (size_t SS = 0 ; SS < DetectedTracksIDsSciHod[p].size() ; SS ++){
                if ( SciHod_TrackID[p][SS] == MUSE_track[track].GetTrackID() ){
                    MUSE_track[track].SetSciHodData( p
                                                   , SciHod_HitFired     [p][SS]
                                                   , SciHod_BarFired     [p][SS]
                                                   , SciHod_HitPos       [p][SS]
                                                   , SciHod_HitTime      [p][SS]
                                                   , SciHod_PreTheta     [p][SS]
                                                   , SciHod_PostTheta    [p][SS]
                                                   , SciHod_PrePhi       [p][SS]
                                                   , SciHod_PostPhi      [p][SS]
                                                   , SciHod_particle_ID  [p][SS]);
                }
            }
        }
        for (int gem = 0 ; gem < Ngem ; gem++)
            for (size_t GM = 0 ; GM < DetectedTracksIDsGEM[gem].size() ; GM ++)
                if ( GEM_TrackID[gem][GM] == MUSE_track[track].GetTrackID() )
                    MUSE_track[track].SetGEMData( gem
                                                 , GEM_HitFired     [gem][GM]
                                                 , GEM_HitPos       [gem][GM]
                                                 , GEM_HitTime      [gem][GM]
                                                 , GEM_particle_ID  [gem][GM]);
        
        for (size_t TGT = 0 ; TGT < DetectedTracksIDsTarget.size() ; TGT ++){
            if ( TARGET_TrackID[TGT] == MUSE_track[track].GetTrackID() )
                MUSE_track[track].SetTargetData( TARGET_HitFired    [TGT]
                                                , TARGET_HitPos     [TGT]
                                                , TARGET_HitTime    [TGT]);
        }
     }
    
    // find out who is the primary generator of each track
    for ( track = 0 ; track < DetectedTracksIDs.size() ; track++){
        G4int tmp           = MUSE_track[track].GetTrackNumber();
        G4int p_gen         = MUSE_track[track].GetTrackID();
        while ( (MUSE_track[tmp].GetTrackID()==p_gen) && (MUSE_track[tmp].GetParentTrackID()!=0) ){
            p_gen       = MUSE_track[tmp].GetParentTrackID();
            for (size_t trk = 0 ; trk < DetectedTracksIDs.size() ; trk++ )
                if ( MUSE_track[trk].GetTrackID() == p_gen )
                    tmp = MUSE_track[trk].GetTrackNumber();
        }
        MUSE_track[track].SetPrimaryGenerator( p_gen );
        
        // add each track the corresponding beam flight time in the beam line of the responsible primary generator
        for (int primary = 0 ; primary < NPrimaryParticles ; primary++){
            if  ( primary_TrackID[primary] == p_gen )
                MUSE_track[track].AddBeamFlightTime( primary_velocity[primary] );

        }
    }
    
    // Print the track data and stream it into ROOT ntuple
    for ( track = 0 ; track < DetectedTracksIDs.size() ; track++){
        if ( ! MUSE_track[track].TrackDetected() ){
            if (PrintOutData)
                G4cout << "Track "<< MUSE_track[track].GetTrackID() << " (particle " << MUSE_track[track].GetParticleID() <<") did not trigger so it is not filled into ROOT..." << G4endl;
//            continue;
        } else {
            if (PrintOutData){
                MUSE_track[track].PrintTrackData();
                G4cout << "Filling track number " << track << ", which is track " << MUSE_track[track].GetTrackID() << ", of particle "<< MUSE_track[track].GetParticleID() << " into ROOT..." << G4endl;
            }
            MUSE_track[track].FillTrackDataToROOT();
        }
    }
}








//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::GetHitData(const Geant4MUSETrackerHit* MUSE_Hit){
    if ( !ParticleHitDetector ){
        // take the hit iformation from the first primary particle hit
        ParentID        = MUSE_Hit -> GetParentID();
        TrackID         = MUSE_Hit -> GetTrackID();
        PrtclName       = MUSE_Hit -> GetParticleName();
        DetectorName    = MUSE_Hit -> GetDetectorName();
        x               = MUSE_Hit -> GetX();
        y               = MUSE_Hit -> GetY();
        z               = MUSE_Hit -> GetZ();
        PrimaryVelocity = (ParentID==0) ? MUSE_Hit->GetVelocity() : -1 ;
        HitTime         = MUSE_Hit -> GetTime() + GenerationTime;
        particle_ID     = MUSE_Hit -> GetParticleCode();
        VertexPosition  = MUSE_Hit -> GetVertexPosition();
        VertexLVName    = MUSE_Hit -> GetVertexLVName();
        CreatorProcessName = MUSE_Hit -> GetCreatorProcessName();
        PreTheta        = MUSE_Hit -> GetPreTheta();                 // angle when entering the detector
        PrePhi          = MUSE_Hit -> GetPrePhi();                   // angle when entering the detector

        ParticleHitDetector = true;
    }
    // take the energy, momentum, and energy deposition information from the traveling within the detector
    Edep                += MUSE_Hit -> GetEdep();
    E_kin               = MUSE_Hit -> GetKineticEnergy();
    p                   = MUSE_Hit -> GetMomentum();
    PostTheta           = MUSE_Hit -> GetPostTheta();                   // angle when leaving the detector
    PostPhi             = MUSE_Hit -> GetPostPhi();                     // angle when leaving the detector
}










//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Geant4MUSEEventAction::PlugDataIntoDetector( G4String detector_name , int trk_ctr ){
    // return a boolean that determines if this hit crossed the threshold or not...
    // SciHod
    for (int p = 0 ; p < Nplanes ; p++ ){
        for (int bar = 0 ; bar < Nbars ; bar++ ){
            if (detector_name == SciHodName(p,bar) ){
                if (Edep > Threshold_Edep){
                    SciHod_HitFired      [p][trk_ctr] = 1;
                    SciHod_BarFired      [p][trk_ctr] = bar+1;
                    SciHod_Edep          [p][trk_ctr] = Edep;
                    SciHod_TrackID       [p][trk_ctr] = TrackID;
                    SciHod_particle_ID   [p][trk_ctr] = particle_ID;
                    SciHod_E_kin         [p][trk_ctr] = E_kin;
                    SciHod_p             [p][trk_ctr] = p;
                    SciHod_HitPos        [p][trk_ctr] = G4ThreeVector( x , y , z );
                    SciHod_HitTime       [p][trk_ctr] = rand2 -> Gaus(HitTime, muse_par.sigma_SciHod);             //bar time resolution
                    SciHod_PreTheta      [p][trk_ctr] = PreTheta;
                    SciHod_PostTheta     [p][trk_ctr] = PostTheta;
                    SciHod_PrePhi        [p][trk_ctr] = PrePhi;
                    SciHod_PostPhi       [p][trk_ctr] = PostPhi;
                    return true;
                } else {
                    if (PrintOutData)
                        G4cout << "TrackID " << TrackID << " hit SciHod " << p+1 <<" but did not deposit enough energy: "<<G4endl  ;
                    return false;
                }
            }
        }
    }
    
    // GEM
    for (int gem = 0 ; gem < Ngem ; gem++ ){
        if (detector_name == GEMName(gem)){
            if (Edep > 100*CLHEP::eV){  // minimal threshold just to supress photons...
                GEM_Edep          [gem][trk_ctr] = Edep;
                GEM_TrackID       [gem][trk_ctr] = TrackID;
                GEM_HitFired      [gem][trk_ctr] = 1;
                GEM_particle_ID   [gem][trk_ctr] = particle_ID;
                GEM_E_kin         [gem][trk_ctr] = E_kin;
                GEM_p             [gem][trk_ctr] = p;
                GEM_HitPos        [gem][trk_ctr] = G4ThreeVector( x , y , z);
                GEM_HitTime       [gem][trk_ctr] = HitTime;
                return true;
            } else {
                if (PrintOutData)
                    G4cout << "TrackID " << TrackID << " hit GEM " << gem+1 <<" but did not deposit enough energy...."<<G4endl  ;
                return false;
            }
        }
    }
    
    // Target
    if (detector_name == (G4String)("TargetDetector")){
        if (Edep > Threshold_Edep){
            TARGET_Edep         [trk_ctr] = Edep;
            TARGET_TrackID      [trk_ctr] = TrackID;
            TARGET_HitFired     [trk_ctr] = 1;
            TARGET_particle_ID  [trk_ctr] = particle_ID;
            TARGET_E_kin        [trk_ctr] = E_kin;
            TARGET_p            [trk_ctr] = p;
            TARGET_HitPos       [trk_ctr] = G4ThreeVector( x , y , z);
            TARGET_HitTime      [trk_ctr] = HitTime;
            return true;
        }
    }
    return false;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String Geant4MUSEEventAction::SciHodName(int p,int bar){
    return (G4String)("SciHod_p_"+G4UIcommand::ConvertToString(p+1)+"_bar_"+G4UIcommand::ConvertToString(bar+1));    }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String Geant4MUSEEventAction::GEMName(int gem){    return (G4String)("GEM_"+G4UIcommand::ConvertToString(gem+1)); }




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::SawTrackInDetector( G4String detector_name , int fTrackID ){

    for (int p = 0 ; p < Nplanes ; p++ ){
        for ( int bar = 0 ; bar < Nbars ; bar++ )
            if (detector_name == SciHodName(p,bar)){
                DetectedTracksIDsSciHod[p].push_back( fTrackID );
                return;
            }
    }
    for (int gem = 0 ; gem < Ngem ; gem++ ){
        if (detector_name == GEMName(gem)){
            DetectedTracksIDsGEM[gem].push_back( fTrackID );
            return;
        }
    }
    if (detector_name == (G4String)("TargetDetector")){
        DetectedTracksIDsTarget.push_back( fTrackID );
        return;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Geant4MUSEEventAction::WasTrackSeenInDetector( G4String detector_name , int fTrackID ){
    for (int p = 0 ; p < Nplanes ; p++ )
        for ( int bar = 0 ; bar < Nbars ; bar++ )
            if (detector_name == SciHodName(p,bar))
                for (size_t i = 0 ; i < DetectedTracksIDsSciHod[p].size() ; i++ )
                    if ( fTrackID == DetectedTracksIDsSciHod[p].at(i) )
                        return true;
    for (int gem = 0 ; gem < Ngem ; gem++ )
        if (detector_name == GEMName(gem))
            for (size_t i = 0 ; i < DetectedTracksIDsGEM[gem].size() ; i++ )
                if ( fTrackID == DetectedTracksIDsGEM[gem].at(i) )
                    return true;
    if (detector_name == (G4String)("TargetDetector"))
        for (size_t i = 0 ; i < DetectedTracksIDsTarget.size() ; i++ )
            if ( fTrackID == DetectedTracksIDsTarget.at(i) )
                return true;
    return false;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::SawTrack( int fTrackID ){ DetectedTracksIDs.push_back( fTrackID ); }



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Geant4MUSEEventAction::WasTrackSeen( int fTrackID ){
    for (size_t i = 0 ; i < DetectedTracksIDs.size() ; i++ ){
       if ( fTrackID == DetectedTracksIDs.at(i) ){
            //            printf("Track %d was detected...\n",fTrackID);
            return true;
        }
    }
    return false;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool Geant4MUSEEventAction::TrackWasDetected(  G4String detector_name , int fTrackID ){
    for (int p = 0 ; p < Nplanes ; p++ )
        if (detector_name == (G4String)("SciHod_p_"+G4UIcommand::ConvertToString(p+1)))
            for (size_t i = 0 ; i < DetectedTracksIDsSciHod[p].size() ; i++ )
                if ( fTrackID == DetectedTracksIDsSciHod[p].at(i) ){
                    //                    printf("Track %d allready detected by SciHod %d...\n",fTrackID,p+1);
                    return true;
                }
    for (int gem = 0 ; gem < Ngem ; gem++ )
        if (detector_name == (G4String)("GEM_"+G4UIcommand::ConvertToString(gem+1)))
            for (size_t i = 0 ; i < DetectedTracksIDsGEM[gem].size() ; i++ )
                if ( fTrackID == DetectedTracksIDsGEM[gem].at(i) ) {
                    //                    printf("Track %d allready detected by GEM %d...\n",fTrackID,gem+1);
                    return true;
                }
    if (detector_name == (G4String)("TargetDetector"))
        for (size_t i = 0 ; i < DetectedTracksIDsTarget.size() ; i++ )
            if ( fTrackID == DetectedTracksIDsTarget.at(i) ) {
                //                printf("Track %d allready detected by Target...\n",fTrackID);
                return true;
            }
    return false;
}









//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Geant4MUSEEventAction::GetGeneratedTrackID( G4int track_id ){
    for ( int i = 0 ; i < NPrimaryParticles ; i ++ ){
        if (primary_TrackID[i] == track_id){
            return track_id;
        }
    }
    return (-1);
}







//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::PrintTracksDetected(){
    G4cout << "-------------------------------------------------------\n";
    G4cout << "The Detected Tracks IDs In this event were:";
    for (std::vector<int>::iterator it=DetectedTracksIDs.begin(); it!=DetectedTracksIDs.end(); ++it)
        G4cout << ' ' << *it;
    G4cout << '\n';
    for (int p = 0 ; p < Nplanes ; p++){
        G4cout << "by SciHod p " << p+1 << '\t' ;
        for (std::vector<int>::iterator it=DetectedTracksIDsSciHod[p].begin(); it!=DetectedTracksIDsSciHod[p].end(); ++it)
            G4cout << ' ' << *it;
        G4cout << '\n';
    }
    for (int gem = 0 ; gem < Ngem ; gem++ ){
        G4cout << "by GEM " << gem+1 << '\t';
        for (std::vector<int>::iterator it=DetectedTracksIDsGEM[gem].begin(); it!=DetectedTracksIDsGEM[gem].end(); ++it)
            G4cout << ' ' << *it;
        G4cout << '\n';
    }
    G4cout << "by Target " << '\t';
    for (std::vector<int>::iterator it=DetectedTracksIDsTarget.begin(); it!=DetectedTracksIDsTarget.end(); ++it)
        G4cout << ' ' << *it;
    G4cout << "\n-------------------------------------------------------\n";
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Geant4MUSEEventAction::PrintHitData(){
    G4cout
    << PrtclName        << " Hit "                  << DetectorName
    << ", Track ID = "  << TrackID
    << " ("             << std::setw(2)             << G4BestUnit( G4ThreeVector( x , y , z) , "Length"  )   << ") "
    << " time "         << std::setw(2)             << G4BestUnit( HitTime  , "Time"  )
    << " e. deposition "<< std::setw(2)             << G4BestUnit( Edep     , "Energy"  )
    << G4endl;
}




