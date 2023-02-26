#include "GanESS.hh"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
//#include <G4SubtractionSolid.hh>
//#include <G4UnionSolid.hh>
//#include <G4MultiUnion.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include <G4UserLimits.hh>
#include <G4SDManager.hh>
#include <G4PhysicalConstants.hh>

#include "nexus/FactoryBase.h"
#include "nexus/UniformElectricDriftField.h"
#include "nexus/OpticalMaterialProperties.h"
#include "nexus/XenonProperties.h"
#include "nexus/IonizationSD.h"

#include <iostream>
#include <cmath>
#include <vector>


using namespace nexus;

REGISTER_CLASS(GanESS, GeometryBase)

GanESS::GanESS():
    GeometryBase(),
    msg_ (nullptr),
    gas_ (nullptr),

    world_rad_         (70*cm),
    gas_rad_ 	       (35*cm),
    gas_length_         (70*cm),
    photoe_prob_       (0.),

    pressure_          (10.* bar),
    temperature_       (293. * kelvin), //esta ok??
    
   //dudas con todo esto
    sc_yield_          (22222 * 1/MeV), // Wsc = 45 eV, fr
    elifetime_         (1e6* ms),
    drift_vel_         (1. * mm/microsecond),
    drift_transv_diff_ (1. * mm/sqrt(cm)),
    drift_long_diff_   (.3 * mm/sqrt(cm)),
    el_field_          (16.0 * kilovolt/cm),
    el_vel_            (3. * mm/microsecond),
    el_transv_diff_    (1. * mm/sqrt(cm)),
    el_long_diff_      (.3 * mm/sqrt(cm)),
    
    specific_vertex_{}
{
  // Messenger
  msg_ = new G4GenericMessenger(this, "/Geometry/GanESS/",
                                "Control commands of the GanESS geometry.");
  // Parametrized dimensions
  DefineConfigurationParameters();
}

GanESS::~GanESS()
{
  delete msg_;

  delete drift_gen_;
  delete sphere_gen_;
  
}

void GanESS::Construct()
{
    //Materials
    steel_ = materials::Steel();
    steel_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
    //Gases: Xe, Kr, Ar
    
    //gas_   = materials::GXe(pressure_, temperature_);
    //gas_->SetMaterialPropertiesTable(opticalprops::GXe(pressure_,
                                                    //  temperature_,
                                                    //  sc_yield_,
                                                    //  elifetime_));

    gas_     = gasses::fXegas(pressure_, temperature_);
    //gas_     = gasses::fArgas(pressure_, temperature_);
    //gas_     = gasses::fKrgas(pressure_, temperature_);
    
    // Print gas information
    G4cout<< "The drift gas is: " << gas_->GetName()<<G4endl;
    G4cout<< "The drift gas pressure (bar) is: "<<gas_->GetPressure()/bar<<G4endl;
    G4cout<< "The drift gas temperature (K) is: "<<gas_->GetTemperature()/kelvin<<G4endl;
    G4cout<< "The drift gas density (kg/m3) is: "<<gas_->GetDensity()/(kg/m3)<<G4endl;
    
    vacuum_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
 
    //Sphere, acting as the world volume
    G4Sphere	    *solid_world_vac = new G4Sphere("World", 0, world_rad_, 0*deg, 360.*deg, 0*deg, 180.*deg);
    G4LogicalVolume *logic_world_vac = new G4LogicalVolume(solid_world_vac, vacuum_, "World");
    this->SetLogicalVolume(logic_world_vac);
    
    //Build inside detector
    BuildTPC(gas_, logic_world_vac);
    
    //Sphere generator                                    
    G4RotationMatrix* rotation_gen_ = new G4RotationMatrix();
    rotation_gen_->rotateX(0*deg);
    rotation_gen_-> rotateY(0*deg);
    rotation_gen_->rotateZ(0*deg);
    
    sphere_gen_  = new SpherePointSampler(world_rad_,0*cm,G4ThreeVector(0.,0.,0.),rotation_gen_,0,twopi,0,pi);
}

G4ThreeVector GanESS::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex;

  if (region == "AD_HOC") {
    vertex = specific_vertex_;
  }

  //// Gas regions
  else if (
    (region == "GasEL") ||
    (region == "GasDrift")) {
    vertex = GenerateVertexGas(region);
  }

  //// Sphere
  if(region == "SPHERE") {vertex = GenerateVertexSphere(region);}

  else {
    G4Exception("[GanESS]", "GenerateVertex()", FatalException,
      "Unknown vertex generation region!");
  }

  return vertex;
}

G4ThreeVector GanESS::GenerateVertexGas(const G4String& region) const
{
    G4ThreeVector vertex;

    if     (region == "GasEL")     {vertex = el_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasDrift")  {vertex = drift_gen_->GenerateVertex("VOLUME");}
    else{G4Exception("[GaP]", "GenerateVertex()", FatalException,
                "Unknown vertex generation region!");}
    return vertex;
}


G4ThreeVector GanESS::GenerateVertexSphere(const G4String& region) const
{
     G4ThreeVector vertex;
     if     (region == "SPHERE")    {vertex = sphere_gen_->GenerateVertex("SURFACE");}
     
     return vertex;
}

void GanESS::DefineConfigurationParameters()
{
  // Gas pressure
  G4GenericMessenger::Command& pressure_cmd =
    msg_->DeclareProperty("pressure", pressure_,
                          "Pressure of the gas.");
  pressure_cmd.SetUnitCategory("Pressure");
  pressure_cmd.SetParameterName("pressure", false);
  pressure_cmd.SetRange("pressure>0.");

  // Gas temperature
  G4GenericMessenger::Command& temperature_cmd =
    msg_->DeclareProperty("temperature", temperature_,
                          "Temperature of the gas.");
  temperature_cmd.SetUnitCategory("Temperature");
  temperature_cmd.SetParameterName("temperature", false);
  temperature_cmd.SetRange("temperature>0.");

  // e- lifetime
  G4GenericMessenger::Command& e_lifetime_cmd =
    msg_->DeclareProperty("elifetime", elifetime_,
                          "Electron lifetime in gas.");
  e_lifetime_cmd.SetParameterName("elifetime", false);
  e_lifetime_cmd.SetUnitCategory("Time");
  e_lifetime_cmd.SetRange("elifetime>0.");

  // Drift velocity in drift region
  G4GenericMessenger::Command& drift_vel_cmd =
    msg_->DeclareProperty("drift_vel", drift_vel_,
                          "Electron drift velocity in the drift region.");
  drift_vel_cmd.SetParameterName("drift_vel", false);
  drift_vel_cmd.SetUnitCategory("Velocity");
  drift_vel_cmd.SetRange("drift_vel>0.");

  // Transverse diffusion in drift region
  new G4UnitDefinition("mm/sqrt(cm)", "mm/sqrt(cm)", "Diffusion", mm/sqrt(cm));
  G4GenericMessenger::Command& drift_transv_diff_cmd =
    msg_->DeclareProperty("drift_transv_diff", drift_transv_diff_,
                          "Tranvsersal diffusion in the drift region");
  drift_transv_diff_cmd.SetParameterName("drift_transv_diff", false);
  drift_transv_diff_cmd.SetUnitCategory("Diffusion");
  drift_transv_diff_cmd.SetRange("drift_transv_diff>0.");

  // Longitudinal diffusion in drift region
  G4GenericMessenger::Command& drift_long_diff_cmd =
    msg_->DeclareProperty("drift_long_diff", drift_long_diff_,
                          "Longitudinal diffusion in the drift region");
  drift_long_diff_cmd.SetParameterName("drift_long_diff", false);
  drift_long_diff_cmd.SetUnitCategory("Diffusion");
  drift_long_diff_cmd.SetRange("drift_long_diff>0.");

  // Scintillation yield (for S1)
  new G4UnitDefinition("1/MeV","1/MeV", "1/Energy", 1/MeV);
  G4GenericMessenger::Command& sc_yield_cmd =
    msg_->DeclareProperty("sc_yield", sc_yield_,
        "Set scintillation yield for gas. It is in photons/MeV");
  sc_yield_cmd.SetParameterName("sc_yield", true);
  sc_yield_cmd.SetUnitCategory("1/Energy");

  // Drift velocity in EL region
  G4GenericMessenger::Command& el_vel_cmd =
    msg_->DeclareProperty("el_vel", el_vel_,
                          "Electron drift velocity in the EL region.");
  el_vel_cmd.SetParameterName("el_vel", false);
  el_vel_cmd.SetUnitCategory("Velocity");
  el_vel_cmd.SetRange("el_vel>0.");

  // Transverse diffusion in EL region
  G4GenericMessenger::Command& el_transv_diff_cmd =
    msg_->DeclareProperty("el_transv_diff", el_transv_diff_,
                          "Tranvsersal diffusion in the EL region");
  el_transv_diff_cmd.SetParameterName("el_transv_diff", false);
  el_transv_diff_cmd.SetUnitCategory("Diffusion");
  el_transv_diff_cmd.SetRange("el_transv_diff>0.");

  // Longitudinal diffusion in EL region
  G4GenericMessenger::Command& el_long_diff_cmd =
    msg_->DeclareProperty("el_long_diff", el_long_diff_,
                          "Longitudinal diffusion in the EL region");
  el_long_diff_cmd.SetParameterName("el_long_diff", false);
  el_long_diff_cmd.SetUnitCategory("Diffusion");
  el_long_diff_cmd.SetRange("el_long_diff>0.");

  // EL field
  new G4UnitDefinition("kilovolt/cm", "kV/cm", "Electric field", kilovolt/cm);
  G4GenericMessenger::Command& el_field_cmd =
    msg_->DeclareProperty("el_field", el_field_,
                          "EL electric field intensity");
  el_field_cmd.SetParameterName("el_field", false);
  el_field_cmd.SetUnitCategory("Electric field");

  // Photoelectric probability
  msg_->DeclareProperty("photoe_prob", photoe_prob_,
                        "Probability of optical photon to ie- conversion");
                        
  // Specific vertex in case region to shoot from is AD_HOC
  msg_->DeclarePropertyWithUnit("specific_vertex", "mm",  specific_vertex_, "Set generation vertex.");
}

void GanESS::BuildTPC(G4Material* gas, G4LogicalVolume* logic_world_vac)
{
    //// Drift
    G4Tubs          *solid_gas_drift = new G4Tubs("GasDrift", 0., gas_rad_, (gas_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_drift = new G4LogicalVolume(solid_gas_drift, gas, "GasDrift");
    
    G4VPhysicalVolume* drift_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logic_gas_drift, "GasDrift", logic_world_vac, false, 0, true);
    drift_gen_  = new CylinderPointSampler2020(drift_phys_);
    
    // Define the drift field
    /*UniformElectricDriftField* drift_field = new UniformElectricDriftField();
    drift_field->SetCathodePosition(drift_z + drift_length_/2);
    drift_field->SetAnodePosition  (drift_z - drift_length_/2);
    drift_field->SetDriftVelocity(drift_vel_);
    drift_field->SetTransverseDiffusion(drift_transv_diff_);
    drift_field->SetLongitudinalDiffusion(drift_long_diff_);*/
    
    G4Region* drift_region = new G4Region("DRIFT");
    //drift_region->SetUserInformation(drift_field);
    drift_region->AddRootLogicalVolume(logic_gas_drift);

    logic_gas_drift->SetUserLimits(new G4UserLimits(1.*mm));
    // Set the DRIFT volume as an ionization sensitive detector
    IonizationSD* active_sd = new IonizationSD("/GanESS/DRIFT");
    logic_gas_drift->SetSensitiveDetector(active_sd);
    G4SDManager::GetSDMpointer()->AddNewDetector(active_sd);
}
