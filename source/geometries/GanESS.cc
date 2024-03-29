#include "GanESS.hh"

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Orb.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UnitsTable.hh"
#include <G4UserLimits.hh>
#include <G4SDManager.hh>
#include <G4PhysicalConstants.hh>
#include <G4GenericMessenger.hh>

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


namespace nexus {
GanESS::GanESS():
    GeometryBase(),
    msg_ (nullptr),
    gas_ (nullptr),
    
    shielding_	       (false),
    world_rad_         (2000*cm),
    sphere_rad_        (70*cm),
    vessel_thickn_     (0.6*cm),
    
    gas_rad_ 	       (35*cm),
    gas_length_         (70*cm),
    
    photoe_prob_       (0.),

    temperature_       (293. * kelvin), 
    
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
  msg_ = new G4GenericMessenger(this, "/Geometry/GanESS/", "Control commands of the GanESS geometry.");
                                
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
    //MATERIALS
    
    //Select gas 
     DefineGas(gas_name_);
  
    //Print gas information
     G4cout<<"\n"<<G4endl;
     G4cout<<"//-----------------------GAS INFO--------------------------//"<<G4endl;
     G4cout<< "The drift gas is: " << gas_->GetName()<<G4endl;
     G4cout<< "The drift gas pressure (bar) is: "<<gas_->GetPressure()/bar<<G4endl;
     G4cout<< "The drift gas temperature (K) is: "<<gas_->GetTemperature()/kelvin<<G4endl;
     G4cout<< "The drift gas density (kg/m3) is: "<<gas_->GetDensity()/(kg/m3)<<G4endl;
     G4cout<<"//---------------------------------------------------------//"<<G4endl;
     G4cout<<"\n"<<G4endl;
    
    steel_ = materials::Steel();
    steel_->SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
                                                   
    vacuum_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    
    steel316ti_ = materials::Steel316Ti();
    steel316ti_ -> SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
    
    water_ = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    
    hdpe_ = materials::HDPE();
    hdpe_ -> SetMaterialPropertiesTable(new G4MaterialPropertiesTable());
    
    //SPHERE, ACTING AS WORLD VOLUME
    G4Sphere	    *solid_world_vac = new G4Sphere("World",0,world_rad_,0*deg,360.*deg,0*deg,180.*deg);
    G4LogicalVolume *logic_world_vac = new G4LogicalVolume(solid_world_vac, vacuum_, "World");
    this->SetLogicalVolume(logic_world_vac);
    
    
    //SHIELDING
    if(shielding_){
    
    	G4double sh_rad_ = gas_rad_ + vessel_thickn_ + sh_thickn_;
    	G4double sh_length_ = gas_length_ + vessel_thickn_*2 + sh_thickn_*2;
    	G4Tubs	    *solid_shield = new G4Tubs("Shielding",0,sh_rad_,sh_length_/2,0.,360.*deg);
    		
    	if(sh_material_ == "Water"){
    		G4LogicalVolume *logic_shield = new G4LogicalVolume(solid_shield,water_,"Shielding");
    		G4VPhysicalVolume *shield_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logic_shield,"Shielding",logic_world_vac,false,0,true);
    		//// Build cylinder inside shielding volume
   		BuildTPC(gas_,steel316ti_,logic_shield);
    	}
    		
    	else if(sh_material_ == "HDPE"){
    		G4LogicalVolume *logic_shield = new G4LogicalVolume(solid_shield,hdpe_,"Shielding");
    		G4VPhysicalVolume *shield_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logic_shield,"Shielding",logic_world_vac,false,0,true);
    		//// Build cylinder inside shielding volume
   		BuildTPC(gas_,steel316ti_,logic_shield);
   	}
   		
   	else{G4Exception("[GanESS]", "Construct()", FatalException,
                "Shielding material has not been specified! Valid options are: Water or HDPE");}
    }
    	
    else{
    
    	////Build cylinder inside world sphere
    	BuildTPC(gas_,steel316ti_,logic_world_vac);
    
    }
    
    
    //Print Shielding information
     G4cout<<"\n"<<G4endl;
     G4cout<<"//-----------------------SHIELDING INFO--------------------------//"<<G4endl;
     G4cout<< "Shielding " << shielding_ <<G4endl;
     if(shielding_)
     {
     	G4cout<<"Shielding material: "<< sh_material_ <<G4endl;
     	G4cout<< "Shielding Thickness (cm):  "<< sh_thickn_/cm<<G4endl;
     	G4cout<<"//----------------------------------------------------------------//"<<G4endl;
     	G4cout<<"\n"<<G4endl;
     }
    
   
    //SPHERE GENERATOR                                   
    G4RotationMatrix* rotation_gen_ = new G4RotationMatrix();
    rotation_gen_->rotateX(0*deg);
    rotation_gen_-> rotateY(0*deg);
    rotation_gen_->rotateZ(0*deg);
    sphere_gen_  = new SpherePointSampler(sphere_rad_,0*cm,G4ThreeVector(0.,0.,0.),rotation_gen_,0,twopi,0,pi);
}

G4ThreeVector GanESS::GenerateVertex(const G4String& region) const
{
  G4ThreeVector vertex;

  if (region == "AD_HOC") {
    vertex = specific_vertex_;
  }
  else {
  	//// Gas regions
  	if ((region == "GasEL") ||(region == "GasDrift")) {vertex = GenerateVertexGas(region);}

  	//// Sphere
  	else{
  		if(region == "SPHERE") {vertex = GenerateVertexSphere(region);}
  		else{
  			   	G4Exception("[GanESS]", "GenerateVertex()", FatalException,
     	 			"Unknown vertex generation region!");
  	            }
  	    }	
	}
  return vertex;
}

G4ThreeVector GanESS::GenerateVertexGas(const G4String& region) const
{
    G4ThreeVector vertex;

    if     (region == "GasEL")     {vertex = el_gen_->GenerateVertex("VOLUME");}
    else if(region == "GasDrift")  {vertex = drift_gen_->GenerateVertex("VOLUME");}
    else{G4Exception("[GanESS]", "GenerateVertex()", FatalException,
                "Unknown vertex generation region!");}
    return vertex;
}


G4ThreeVector GanESS::GenerateVertexSphere(const G4String& region) const
{
     G4ThreeVector vertex;
     if     (region == "SPHERE")    {vertex = sphere_gen_->GenerateVertex("SURFACE");}
     
     return vertex;
}


void GanESS::DefineGas(G4String gasname)
{
  if(gas_name_ == "KRIPTON"){
  		G4cout<< " gas_name_ is " << gasname << G4endl;
  		gas_     = gasses::fKrgas(pressure_, temperature_);
  }
  else{
  	if(gas_name_ == "ARGON"){
  		G4cout<< " gas_name_ is " << gasname << G4endl;
  		gas_     = gasses::fArgas(pressure_, temperature_);
  	}
  	else{
  		if(gas_name_ == "XENON"){
  			G4cout<< " gas_name_ is " << gasname << G4endl;
  			gas_     = gasses::fXegas(pressure_, temperature_);
  		}
  		else{
  			G4Exception("[GanESS]","DefineGas()",FatalException,
  				"Unknown kind of gas, valid options are: XENON, ARGON, KRIPTON, VACUUM");
  		}
  	}
  }
}


void GanESS::DefineConfigurationParameters()
{  
	
   //Shielding configuration
   msg_->DeclareProperty("shielding",shielding_,"Set shielding");
   msg_->DeclareProperty("sh_material",sh_material_,"Shielding material: in case of 1 layer");

   //Shielding thickness
   G4GenericMessenger::Command& sh_thickn_cmd =
   msg_->DeclareProperty("sh_thickn",sh_thickn_,"Outter shielding layer thickness");
   sh_thickn_cmd.SetUnitCategory("Length");
   sh_thickn_cmd.SetParameterName("sh_thickn",false);
   sh_thickn_cmd.SetRange("sh_thickn>0.");
      
  //Gas material
    msg_->DeclareProperty("gas_name",gas_name_,
  			"Gas");	
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
  
  //Sphere point sampler radius 
  G4GenericMessenger::Command& sphere_rad_cmd =
   msg_->DeclareProperty("sphere_rad",sphere_rad_,
   				"Radius of sphere point sampler.");
   sphere_rad_cmd.SetUnitCategory("Length");
   sphere_rad_cmd.SetParameterName("sphere_rad",false);
   sphere_rad_cmd.SetRange("sphere_rad>0.");
  
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


void GanESS::BuildTPC(G4Material* gas,G4Material* vessel_mat, G4LogicalVolume* logic_mother_volume)
{    
    
    //// Vessel
    G4double vessel_rad_ = gas_rad_ + vessel_thickn_;
    G4double vessel_length_ = gas_length_ + vessel_thickn_*2;
    
    G4Tubs	    *solid_vessel = new G4Tubs("Vessel",0,vessel_rad_,vessel_length_/2,0.,360.*deg);
    G4LogicalVolume *logic_vessel = new G4LogicalVolume(solid_vessel,vessel_mat,"Vessel");
    G4VPhysicalVolume *vessel_phys = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logic_vessel,"Vessel",logic_mother_volume,false,0,true);
    
    
    //// Drift
    G4Tubs          *solid_gas_drift = new G4Tubs("GasDrift", 0., gas_rad_, (gas_length_)/2, 0., 360.*deg);
    G4LogicalVolume *logic_gas_drift = new G4LogicalVolume(solid_gas_drift, gas, "GasDrift");
    G4VPhysicalVolume* drift_phys_ = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logic_gas_drift, "GasDrift", logic_vessel, false, 0, true);
    drift_gen_  = new CylinderPointSampler2020(drift_phys_);
    
    
    G4Region* drift_region = new G4Region("DRIFT");
    drift_region->AddRootLogicalVolume(logic_gas_drift);

    logic_gas_drift->SetUserLimits(new G4UserLimits(1.*mm));
    
    // Set the DRIFT volume as an ionization sensitive detector
    //position, time and energy deposition will be stored for
    //each step of any charged particle crossssing the volume
    
    IonizationSD* active_sd = new IonizationSD("/GanESS/DRIFT");
    logic_gas_drift->SetSensitiveDetector(active_sd);
    G4SDManager::GetSDMpointer()->AddNewDetector(active_sd);
    
    
}

}


