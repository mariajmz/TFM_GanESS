#ifndef GanESS_H
#define GanESS_H

#include "nexus/GeometryBase.h"
#include "nexus/CylinderPointSampler2020.h"
#include "nexus/SpherePointSampler.h"
#include "nexus/MaterialsList.h"

#include "GanESSMaterials.hh"
#include "UserGenerator.hh"


#include <G4ThreeVector.hh>
//#include <G4GenericMessenger.hh>


class G4GenericMessenger;

//using namespace nexus;

namespace nexus {
class GanESS : public GeometryBase
{

    public:
    // Constructor
        GanESS();
    //Destructor
        ~GanESS();
        
          
          G4ThreeVector GenerateVertex   (const G4String& region) const;
          G4ThreeVector GenerateVertexGas(const G4String& region) const;
          G4ThreeVector GenerateVertexSphere(const G4String& region) const;
          void DefineGas(G4String gasname);
          void DefineConfigurationParameters();
          void Construct();
	  void BuildTPC(G4Material* gas, G4LogicalVolume* logic_world_vac);

    private:
        G4GenericMessenger* msg_;

        // Materials
        G4String gas_name_; //to define via macro material
        G4Material* gas_;
        G4Material* vacuum_;
        G4Material* steel_;

        // Gas parameters
        G4double pressure_;
        G4double temperature_;
        
        G4double sc_yield_;
        G4double elifetime_;
        G4double drift_vel_;
        G4double drift_transv_diff_;
        G4double drift_long_diff_;

	
	//Dimensions
	G4double world_rad_;
	G4double gas_rad_;
	G4double gas_length_;
	
	/*
	// Vessel parameters
        G4double vessel_out_rad_   ;
        G4double vessel_out_length_;
        G4double vessel_rad_   ;
        G4double vessel_length_;

        // Mesh
        G4double mesh_rad_   ;
        G4double mesh_thickn_;
        G4double mesh_transparency_;

        // Mesh Bracket
        G4double meshBracket_rad_;
        G4double meshBracket_thickn_ ;
        G4double anodeBracket_rad_ ;
        G4double anodeBracket_thickn_ ;*/

        G4double pmt_rad_;
        
        
        // EL field (in kV/cm)
        G4double el_field_;
        G4double el_vel_;
        G4double el_transv_diff_;
        G4double el_long_diff_;

	// photoelectric prob
	G4double photoe_prob_;

        // Vertex generation
        G4ThreeVector specific_vertex_;

        CylinderPointSampler2020* drift_gen_;
        CylinderPointSampler2020* el_gen_;
        SpherePointSampler* sphere_gen_;

};
}
#endif
