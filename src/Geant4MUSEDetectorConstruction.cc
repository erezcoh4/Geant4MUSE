// ********************************************************************
// Geant4MUSE simulation
// ********************************************************************
/// \file Geant4MUSEDetectorConstruction.cc
/// \brief Implementation of the Geant4MUSEDetectorConstruction class

#include "Geant4MUSEDetectorConstruction.hh"
#include "Geant4MUSEDetectorMessenger.hh"
#include "Geant4MUSETrackerSD.hh"
#include "Geant4MUSEPrimaryGeneratorAction.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* Geant4MUSEDetectorConstruction::fMagFieldMessenger = 0;

Geant4MUSEDetectorConstruction::Geant4MUSEDetectorConstruction()
:G4VUserDetectorConstruction(),
fNbOfChambers(0),
fLogicTarget(NULL), fLogicChamber(NULL),
fTargetMaterial(NULL), fChamberMaterial(NULL),
fStepLimit(NULL),
fCheckOverlaps(true)
{
    fMessenger = new Geant4MUSEDetectorMessenger(this);
    fNbOfChambers = 5;
    fLogicChamber = new G4LogicalVolume*[fNbOfChambers];

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Geant4MUSEDetectorConstruction::~Geant4MUSEDetectorConstruction()
{
    delete [] fLogicChamber;
    delete fStepLimit;
    delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Geant4MUSEDetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();
    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::DefineMaterials()
{
    // Material definition
    G4NistManager* man = G4NistManager::Instance();
    
    G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
    G4Material* vacuum = new G4Material("Vacuum", 1.e-5*g/cm3, 1, kStateGas, 2.e-2*bar, CLHEP::STP_Temperature);
    vacuum -> AddMaterial(Air, 1.);
    
    G4Element* elO = new G4Element("Oxygen","O2",8.,16.*g/mole);
    G4Element* elH = new G4Element("Hydrogen"  ,"H2" , 1., 1.01*g/mole );
    G4Element* elC = new G4Element("Carbon", "C", 6., 12.01*g/mole);
    G4Element* elN = new G4Element("Nitrogen", "N2", 7., 14.01*g/mole);
    G4Material* Kapton = new G4Material("Kapton", 1.42*g/cm3, 4);
    Kapton->AddElement(elH, 0.0273);
    Kapton->AddElement(elC, 0.7213);
    Kapton->AddElement(elN, 0.0765);
    Kapton->AddElement(elO, 0.1749);

    G4double density;
    std::vector<G4int> natoms;
    std::vector<G4String> elements;
    G4NistManager* fNistMan = G4NistManager::Instance();
    elements.push_back("C");     natoms.push_back(8);
    elements.push_back("H");     natoms.push_back(8);
    density = 1.050*g/cm3;
    fPolystyrene = fNistMan-> ConstructNewMaterial("Polystyrene", elements, natoms, density);
    elements.clear();
    natoms.clear();
    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Geant4MUSEDetectorConstruction::DefineVolumes()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4Material* vacuum          = nist->FindOrBuildMaterial("Vacuum");
    G4Material* default_mat     = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* Kapton_mat      = nist->FindOrBuildMaterial("Kapton");
    G4Material* Quartz          = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    G4Material* Scintillator    = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    G4Material* MYLAR           = nist->FindOrBuildMaterial("G4_MYLAR");
    G4RotationMatrix    rot_mat = G4RotationMatrix();
    G4Transform3D       transform;

    
    
    //-----------------------------//
    //------- World  --------------//
    //-----------------------------//
    G4double worldXYZ =  300.*cm;
    G4Box* solidWorld = new G4Box("World", 0.5*worldXYZ, 0.5*worldXYZ, 0.5*worldXYZ );
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, default_mat, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld,  "World", 0,  false,  0, fCheckOverlaps);
    
    
    
    
    
    //-----------------------------//
    //---- Entrance Window --------//
    //-----------------------------//
    if (muse_par.BuildEntranceWindow){
        G4Box * solidEntrance = new G4Box("EntranceWindow" , 10.*cm/2. , 10.*cm/2. , 180*um/2. );
        G4LogicalVolume* logicEntrance = new G4LogicalVolume(solidEntrance,MYLAR,"EntranceWindow");
        logicEntrance -> SetVisAttributes(muse_par.EntranceVisAttributes);
        new G4PVPlacement( 0 , G4ThreeVector( 0 , 0 , muse_par.EntrancePosition ) , logicEntrance , "EntranceWindow" , logicWorld, false, 0);
    }
    
    
    
    
    
    
    
    
    
    
    
    //------- SciHod ---------------//
    if (muse_par.BuildSciHod){
        G4Box*              SciHodBar    [muse_par.Nplanes][muse_par.Nbars];
        G4LogicalVolume*    SciHodLog    [muse_par.Nplanes][muse_par.Nbars];
        G4VPhysicalVolume*  SciHodPhys   [muse_par.Nplanes][muse_par.Nbars];
        G4Box*              SprtrBar    [muse_par.Nplanes][muse_par.Nbars];
        G4LogicalVolume*    SprtrLog    [muse_par.Nplanes][muse_par.Nbars];
        G4VPhysicalVolume*  SprtrPhys   [muse_par.Nplanes][muse_par.Nbars];
        G4ThreeVector       BarPosition;
        G4ThreeVector       SprtrPosition;
        G4double            Bar_x , Bar_y , Sprtr_x , Sprtr_y , Plane_z ;
        for (int p = 0 ; p < muse_par.Nplanes ; p++ ) {
            rot_mat = G4RotationMatrix();
            switch (p) {
                case 0://**********  U - plane ********//
                    rot_mat.rotateZ(135*deg);
                    break;
                case 1://**********  Y - plane ********//
                    rot_mat.rotateZ(0*deg);
                    break;
                default: //**********  X - plane ********//
                    rot_mat.rotateZ(90*deg);
                    break;
            }
            Plane_z = muse_par.SciHodPosition + p * muse_par.SciHodPlanesSeperation;
            for (int bar = 0 ; bar < muse_par.Nbars ; bar++ ) {
                SciHodBar   [p][bar] = new G4Box("SciHodBar", muse_par.SciLength/2. , muse_par.SciWidth/2. , muse_par.SciThick/2.  );
                SciHodLog   [p][bar] = new G4LogicalVolume(SciHodBar[p][bar] , fPolystyrene  , muse_par.GetSciHodName("SciHod",p,bar) );
                SciHodLog   [p][bar] -> SetVisAttributes(muse_par.SciAttributes);
                switch (p) {
                    case 0://**********  U - plane ********//
                        Bar_x   = std::cos(45) * (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2);
                        Bar_y   = std::sin(45) * (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2);
                        Sprtr_x = std::cos(45) * (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2 + 1/2.);
                        Sprtr_y = std::sin(45) * (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2 + 1/2.);
                        break;
                    case 1://**********  Y - plane ********//
                        Bar_x   = 0;
                        Bar_y   = (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2);
                        Sprtr_x = 0;
                        Sprtr_y = (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2 + 1/2.);
                        break;
                    default: //**********  X - plane ********//
                        Bar_x   = (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2);
                        Bar_y   = 0;
                        Sprtr_x = (muse_par.SciWidth + muse_par.SprtrWidth) * (bar - muse_par.Nbars/2 + 1/2.);
                        Sprtr_y = 0;
                        break;
                }
                BarPosition         = G4ThreeVector( Bar_x      , Bar_y     , Plane_z );
                transform           = G4Transform3D( rot_mat , BarPosition );
                SciHodPhys  [p][bar] = new G4PVPlacement(transform ,SciHodLog[p][bar], "SciHodBar" , logicWorld , false, 0);
                if (bar < muse_par.Nbars - 1 ) {
                    SprtrBar   [p][bar] = new G4Box("SprtrBar", muse_par.SciLength/2. , muse_par.SprtrWidth/2. , muse_par.SciThick/2.  );
                    SprtrLog   [p][bar] = new G4LogicalVolume(SprtrBar[p][bar] , MYLAR , muse_par.GetSciHodName("Sprtr",p,bar) );
                    SprtrLog   [p][bar] -> SetVisAttributes(muse_par.SprtrAttributes);
                    SprtrPosition       = G4ThreeVector( Sprtr_x    , Sprtr_y   , Plane_z );
                    transform           = G4Transform3D( rot_mat , SprtrPosition );
                    SprtrPhys  [p][bar] = new G4PVPlacement(transform ,SprtrLog[p][bar], "SprtrBar" , logicWorld , false, 0);
                }
            }
        }
    }
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    //------- GEM -----------------//
    if(muse_par.BuildGEM){
        G4Box * solidGEM[3];
        G4LogicalVolume * logicGEM[3];
        G4VPhysicalVolume *  physGEM[3];
        for (int gem = 0 ; gem < 3 ; gem++ ){
            solidGEM[gem] = new G4Box( "GEM" , muse_par.GEM_side/2. , muse_par.GEM_side/2. , muse_par.GEM_thick/2. );
            logicGEM[gem] = new G4LogicalVolume(solidGEM[gem],MYLAR, muse_par.GetGEMName(gem) );
            logicGEM[gem] -> SetVisAttributes(muse_par.GEMAttributes);
            physGEM[gem] = new G4PVPlacement( 0 , G4ThreeVector( 0 , 0 , muse_par.GEMPosition + gem * muse_par.GEMPlanesSeperation )
                                             , logicGEM[gem], "GEM" , logicWorld, false, 0);
        }
    }
    
    
    
    
    
    
    
    
    
    
    //-----------------------------//
    //----  Target              ---//
    //-----------------------------//
    if (muse_par.BuildTarget){
        // simplified entracne window
        G4double EntWind_Length      =   200.*um; // entrance window thickness is 200 um
        G4double EntWind_Radius      =   10.*cm;  // diameter should find out! this is just a guess
        
        G4Tubs* solidEntWindDetector = new G4Tubs("EntWindDetector", 0 , EntWind_Radius , EntWind_Length/2. , 0., CLHEP::twopi);
        G4LogicalVolume* logicEntWindDetector = new G4LogicalVolume(solidEntWindDetector,Kapton_mat,"EntWindDetectorLV");
        new G4PVPlacement(0 , G4ThreeVector(0,0,muse_par.z0),logicEntWindDetector, "EntWindDetector",logicWorld, false,  0);
        logicEntWindDetector -> SetVisAttributes(muse_par.EntWindDetecVisAttributes);
        // simplified target chamber - just as a sensitive detector
        G4double Target_Length      =   1.*mm;
        G4double Target_Radius      =   30.*cm;  // target diameter 6 cm (Ron Gilman Aug-16)
        G4Tubs* solidTargetDetector = new G4Tubs("TargetDetector" , 0 , Target_Radius, Target_Length/2., 0., CLHEP::twopi);
        G4LogicalVolume* logicTargetDetector = new G4LogicalVolume(solidTargetDetector,Scintillator,"TargetDetectorLV");
        new G4PVPlacement(0 , G4ThreeVector(0,0,muse_par.TargetPosition),logicTargetDetector, "TargetDetector",logicWorld, false,  0);
        logicTargetDetector -> SetVisAttributes( muse_par.TargetDetecVisAttributes);
    }
    
    
    
    
    // Always return the physical world
    return physWorld;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::ConstructSDandField()
{   // Sensitive detectors
    

    //GEM active detector
    if (muse_par.BuildGEM){
        Geant4MUSETrackerSD* GEMTrackerSD[muse_par.Ngem];
        for (int gem = 0 ; gem < muse_par.Ngem ; gem++ ){
            G4String GEMSDname = "Detector/GEM_" + G4UIcommand::ConvertToString(gem+1) ;
            GEMTrackerSD[gem] = new Geant4MUSETrackerSD(GEMSDname,"GEM_"+G4UIcommand::ConvertToString(gem+1));
            SetSensitiveDetector(muse_par.GetGEMName(gem), GEMTrackerSD[gem] , true);
        }
    }
    
    
    
    //SciHod active detector
    if (muse_par.BuildSciHod){
        Geant4MUSETrackerSD* SciHodTrackerSD[muse_par.Nplanes][muse_par.Nbars];
        for (int p = 0 ; p < muse_par.Nplanes ; p++ ) {
            for (int bar = 0 ; bar < muse_par.Nbars ; bar++ ) {
                G4String SciHodSDname = "Detector/SciHod_p_" + G4UIcommand::ConvertToString(p+1) + "_bar_" + G4UIcommand::ConvertToString(bar+1)  ;
                SciHodTrackerSD[p][bar] = new Geant4MUSETrackerSD(SciHodSDname,"SciHod_p_"+G4UIcommand::ConvertToString(p+1)+"_bar_"+G4UIcommand::ConvertToString(bar+1) );
                SetSensitiveDetector(muse_par.GetSciHodName("SciHod" , p , bar ), SciHodTrackerSD[p][bar] , true);
            }
        }
    }
    
  
    
    
    
    //Target Detector
    if(muse_par.BuildTarget){
        G4String trackerChamberSDname = "Detector/TargetDetector";
        Geant4MUSETrackerSD* aTrackerSD = new Geant4MUSETrackerSD(trackerChamberSDname,"TargetDetector");
        SetSensitiveDetector("TargetDetectorLV" , aTrackerSD , true);
    }
    
    
    
    
    
    G4ThreeVector fieldValue = G4ThreeVector();
    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::SetTargetMaterial(G4String materialName)
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4Material* pttoMaterial =
    nistManager->FindOrBuildMaterial(materialName);
    
    if (fTargetMaterial != pttoMaterial) {
        if ( pttoMaterial ) {
            fTargetMaterial = pttoMaterial;
            if (fLogicTarget) fLogicTarget->SetMaterial(fTargetMaterial);
            //        G4cout << "\n----> The target is made of " << materialName << G4endl;
        } else {
            //        G4cout << "\n-->  WARNING from SetTargetMaterial : "
            //               << materialName << " not found" << G4endl;
        }
    }
}












//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::SetChamberMaterial(G4String materialName)
{
    G4NistManager* nistManager = G4NistManager::Instance();
    
    G4Material* pttoMaterial =
    nistManager->FindOrBuildMaterial(materialName);
    
    if (fChamberMaterial != pttoMaterial) {
        if ( pttoMaterial ) {
            fChamberMaterial = pttoMaterial;
            for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
                if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->
                    SetMaterial(fChamberMaterial);
            }
            //        G4cout << "\n----> The chambers are made of " << materialName << G4endl;
        } else {
            //        G4cout << "\n-->  WARNING from SetChamberMaterial : "
            //               << materialName << " not found" << G4endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::SetMaxStep(G4double maxStep)
{
    if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Geant4MUSEDetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps)
{
    fCheckOverlaps = checkOverlaps;
}  
