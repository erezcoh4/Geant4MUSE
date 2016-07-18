// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
// $Id: Geant4MUSERunAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file Geant4MUSERunAction.cc
/// \brief Implementation of the Geant4MUSERunAction class

#include "Geant4MUSERunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Geant4MUSEAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSERunAction::Geant4MUSERunAction()
: G4UserRunAction()
{
//    G4cout << "\n----------Geant4MUSERunAction: begining of run----------" << G4endl;
    // Create an Ntuple. The number of entries would be the number of detected tracks in
    // N(entries) = Sum(over all events) [ Number of tracks in largest hit collection per event ]
    // set printing event number per each 100 events
    G4RunManager::GetRunManager()->SetPrintProgress(muse_par.PrintEachEvent);
    
    // Create Ntuple
    G4AnalysisManager* man = G4AnalysisManager::Instance(); // Open an output file
    man->CreateNtuple("MUSE", "beam detectors data");
    
    // SciHod
    const int Nplanes = 3;
    for(int p = 0 ; p < Nplanes ; p++ ){
        man->CreateNtupleIColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_HitFired");           //11*p+0
        man->CreateNtupleIColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_BarFired");           //11*p+1
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_x");                  //11*p+2
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_y");                  //11*p+3
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_z");                  //11*p+4
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_Time");               //11*p+5
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_PreTheta");           //11*p+6
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_PostTheta");          //11*p+7
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_PrePhi");             //11*p+8
        man->CreateNtupleDColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_PostPhi");            //11*p+9
        man->CreateNtupleIColumn("SciHod_"+G4UIcommand::ConvertToString(p+1)+"_PID");                //11*p+10
    }
    // GEM
    const int Ngem = 3;
    for(int gem = 0 ; gem < Ngem ; gem++ ){
        man->CreateNtupleIColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_HitFired");         //11*Nplaness+6*gem+0
        man->CreateNtupleDColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_x");                //11*Nplanes+6*gem+1
        man->CreateNtupleDColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_y");                //11*Nplanes+6*gem+2
        man->CreateNtupleDColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_z");                //11*Nplanes+6*gem+3
        man->CreateNtupleDColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_Time");             //11*Nplanes+6*gem+4
        man->CreateNtupleIColumn("GEM_"+G4UIcommand::ConvertToString(gem+1)+"_PID");              //11*Nplanes+6*gem+5
    }
    // Target
    man->CreateNtupleIColumn("Target_HitFired");         //11*Nplanes+6*Ngem+0
    man->CreateNtupleDColumn("Target_x");                //11*Nplanes+6*Ngem+1
    man->CreateNtupleDColumn("Target_y");                //11*Nplanes+6*Ngem+2
    man->CreateNtupleDColumn("Target_z");                //11*Nplanes+6*Ngem+3
    man->CreateNtupleDColumn("Target_Time");             //11*Nplanes+6*Ngem+4
    // Global track
    man->CreateNtupleIColumn("DetectedTrackID");              // detected track ID (events in Ntuple are sorted by track ids in each event)
    man->CreateNtupleIColumn("GeneratedTrackID");             // generated track ID (set to -10000 if particle is secondary)
    man->CreateNtupleIColumn("TrackParticleID");              // track particle type code (for wherever it was detected....)
    man->CreateNtupleIColumn("GEMTimeWindow");                // divide the tracks into GEM time windows of 100 ns

    man->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Geant4MUSERunAction::~Geant4MUSERunAction()
{    //    G4cout << "ending Geant4MUSERunAction" << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSERunAction::BeginOfRunAction(const G4Run*)
{
    //inform the runManager to save random number seed
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    // Get analysis manager
    G4AnalysisManager* man = G4AnalysisManager::Instance(); // Open an output file
    man -> OpenFile("Output");
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Geant4MUSERunAction::EndOfRunAction(const G4Run* ){
    
//    G4cout << "\n----------Geant4MUSERunAction: end of run----------" << G4endl;
    G4AnalysisManager * man = G4AnalysisManager::Instance();
    man -> Write();
    man -> CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
