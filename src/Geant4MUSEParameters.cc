// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
// $Id: Geant4MUSEParameters.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file Geant4MUSEParameters.cc
/// \brief Implementation of the Geant4MUSEParameters class

#include "Geant4MUSEParameters.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Geant4MUSEParameters::Geant4MUSEParameters() {
    
    
    // run parameters....
    // ------------------------------
    momentum                    = 210*MeV;
    momentumSpread              = 0.008;
    BeamType                    = "beam 2015";
    BeamParticleType            = "TDR 3 particles";
    
    
    
    
    
    // event action
    // ------------------------------
    sigma_SciHod                = 0.100*(CLHEP::ns);      //time resolution per SciHod bar is 90 ps
    fBeamRate                   = 50.e6*CLHEP::hertz;   // each bunch comes after 20 ns
    PrintEachEvent              = 50000;

    
    

    // beam
    // ------------------------------
    BuildEntranceWindow         = true;
    EntrancePosition            = -140.*cm;
    EntranceVisAttributes       = new G4VisAttributes( G4Colour( 0.4 , 0.4 , 0.4 ) );
    
    
    
    
    
    // target
    // ------------------------------
    BuildTarget                 = true;
    TargetPosition              = 0.*cm;
    z0                          = -7.*cm;
    EntWindDetecVisAttributes   = new G4VisAttributes( G4Colour( 0.4 , 0.4 , 0.4 ) );
    TargetDetecVisAttributes    = new G4VisAttributes( G4Colour( 0.75 , 0.6 , 0.75 ) );
   
    
    
    
    // GEM
    // ------------------------------
    BuildGEM                = true;
    GEMPosition             = -300.*mm + z0; // Upsream GEM plane
    GEM_side                = 10.*cm;
    // Mylar radiation length = 39.95 cm [http://pdg.lbl.gov/2002/atomicrpp.ps], MK: "each GEM radiation length = 0.32% => GEM_length = 1.284*mm
    GEM_thick               = 1.284*mm;
    GEMPlanesSeperation     = 8.4*cm;
    GEMAttributes           = new G4VisAttributes(G4Colour( 0.99 , 0.88 , 0.66 ));
    
    
    
    // SciHod
    // ------------------------------
    BuildSciHod             = true;
    SciLength               = 100.*mm;
    SciWidth                = 5.*mm;
    SciThick                = 2.*mm;
    SprtrWidth              = 14*um;            // 2 x 7 µm aluminized mylar cross-talk seperator between paddles + 0 µm gap
    SciHodPosition          = -410.*mm + z0;    // Scintillator + SiPMs - Upstream plane U
    SciHodPlanesSeperation  = 15*mm;
    SciAttributes           = new G4VisAttributes(G4Colour(0.25,0.25,0.75));
    SprtrAttributes         = new G4VisAttributes(G4Colour(0.90 , 0.88 , 0.88 ));
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String Geant4MUSEParameters::GetSciHodName(G4String name, int p ,int bar){
    return name + "_p_" + G4UIcommand::ConvertToString(p+1) + "_bar_" + G4UIcommand::ConvertToString(bar+1) + "_LV";
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4String Geant4MUSEParameters::GetGEMName(int gem){
    return "GEM_" + G4UIcommand::ConvertToString(gem+1) + "LV";
}










