// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEParameters.hh
/// \brief Definition of the Geant4MUSEParameters class


#ifndef __Geant4MUSEParameters_H__
#define __Geant4MUSEParameters_H__



#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "G4VisAttributes.hh"


#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"


#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class Geant4MUSEParameters {
public:
    
    
    
    //global variables
    G4double    sigma_SciHod;        //time resolution per SciHod bar
    G4double    fBeamRate;
    int         PrintEachEvent;
    G4double    momentum;
    G4double    momentumSpread;
    G4String    BeamType;
    G4String    BeamParticleType;
    
    
    bool BuildEntranceWindow , BuildSciHod , BuildGEM , BuildTarget;

 
    // beam components positions
    G4double z0;
    G4double TargetPosition;
    G4double EntrancePosition;
    G4VisAttributes* EntranceVisAttributes;
    G4VisAttributes* EntWindDetecVisAttributes;
    G4VisAttributes* TargetDetecVisAttributes;
   

    
    // GEMs
    static const int Ngem       = 3;
    G4double GEMPosition;
    G4double GEM_thick;
    G4double GEM_side;
    G4double GEMPlanesSeperation;
    G4VisAttributes* GEMAttributes;
    
    
    
    
    
    // SciHod
    static const int Nplanes    = 3; // Beam ---> U (-510 mm) ---> Y (-495 mm) ---> X (-480 mm)
    static const int Nbars      = 20;
    G4double SciHodPosition;
    G4double SciWidth;
    G4double SciLength;
    G4double SciThick;
    G4double SprtrWidth;
    G4double SciHodPlanesSeperation;
    G4VisAttributes* SciAttributes;
    G4VisAttributes* SprtrAttributes;
    
    
    

    
    
    // constructors
    Geant4MUSEParameters();
    ~Geant4MUSEParameters(){};
    
    
    
    
    G4String GetSciHodName  (G4String , int , int );
    G4String GetGEMName     (int);
    
    
    // get methods
    int GetNSciBar      ()      {return Nbars;};
    int GetNSciPlanes   ()      {return Nplanes;};
    int GetNgem         ()      {return Ngem;};
    
    
private:
    
    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


