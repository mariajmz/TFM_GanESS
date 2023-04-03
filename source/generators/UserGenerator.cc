//-----------------------------------------------------------------------------
// UserGenerator.cc
//   
// This class is the primary generator for particles with a user defined energy 
// histogram      
// 
// Modification based on MuonGenerator from nexus/generator
//-----------------------------------------------------------------------------

#include "UserGenerator.hh"

#include "nexus/MuonGenerator.h"
#include "nexus/DetectorConstruction.h"
#include "nexus/GeometryBase.h"
#include "nexus/AddUserInfoToPV.h"
#include "nexus/FactoryBase.h"
#include "nexus/RandomUtils.h"
#include "nexus/IOUtils.h"

#include <G4Event.hh>
#include <G4GenericMessenger.hh>
#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryVertex.hh>
#include <G4Event.hh>
#include <G4RandomDirection.hh>
#include <Randomize.hh>

#include "CLHEP/Units/SystemOfUnits.h"

using namespace nexus;
REGISTER_CLASS(UserGenerator, G4VPrimaryGenerator)

UserGenerator::UserGenerator():
G4VPrimaryGenerator(), msg_(0), particle_definition_(0),
user_ene_dist_(true),geom_(0),IsDistribution_(false),mono_ene_(false), geom_solid_(0), radius_(70*cm)
{

msg_ = new G4GenericMessenger(this, "/Generator/UserGenerator/",
				"Control commands for user generator.");
				 
 msg_->DeclareProperty("region", region_,
			"Set the region of the geometry where the vertex will be generated.");

msg_->DeclareProperty("user_ene_dist", user_ene_dist_,
			"Distribute particles energies according to file?");
			
msg_->DeclareProperty("mono_ene",mono_ene_,
			"Generate mono energetic particles");
						
msg_->DeclareProperty("ene_file", ene_file_,
			"Name of the file containing angular distribution.");
			
G4GenericMessenger::Command& radius_cmd =
msg_->DeclareProperty("radius",radius_,
   			"Radius of sphere point sampler.");
radius_cmd.SetUnitCategory("Length");
radius_cmd.SetParameterName("radius",false);
radius_cmd.SetRange("radius>0.");


//Print information
G4cout<<"\n"<<G4endl;
G4cout<<"//-----------------------USER GENERATOR INFO--------------------------//"<<G4endl;
G4cout<< "Sphere radius (cm): " << radius_/cm <<G4endl;
if(user_ene_dist_)
{
     	G4cout<<"Using energy file distribution for neutrons."<<G4endl;
}
G4cout<<"//--------------------------------------------------------------------//"<<G4endl;
G4cout<<"\n"<<G4endl;

DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();
				
}

UserGenerator::~UserGenerator()
{
  delete msg_;
}


void UserGenerator::LoadUserDistribution()
{
	//energy_ and energy_bins_ only used
	LoadHistData1D(ene_file_,flux_, energy_, energy_bins_);
	
	//convert vector to array
	auto arr_flux = flux_.data();
	
	//Initialice random number generator based on flux distribution
	fRandomGeneral_ = new G4RandGeneral(arr_flux,flux_.size());
}


void UserGenerator::GeneratePrimaryVertex(G4Event* event)
{

	G4double kinetic_energy, energy, mass;
	if(!IsDistribution_){
		if (user_ene_dist_){
      		//std::cout << "Generating particles using user distribution loaded from file" << std::endl;
      		LoadUserDistribution();
      		
     	 }
		//Set initialisation
		IsDistribution_ = true;
	
	}
	
	if (user_ene_dist_){
      		kinetic_energy = GetEnergy();
     	 }
	
	else{
		//if user does not especified file Generate Uniform Random Energy in [Emin,Emax]
		if(!mono_ene_) {kinetic_energy = UniformRandomInRange(0.002*GeV,0.001*GeV);}
		else {kinetic_energy = 2.*MeV;}	
	}
	
	
	particle_definition_ = 
		G4ParticleTable::GetParticleTable()->FindParticle("neutron");
	
	// Particle properties
  	mass          = particle_definition_->GetPDGMass();
  	energy        = kinetic_energy + mass;
  	
	//Generate uniform random direction in 4pi: RandomDirectionInRange
	G4ThreeVector p_dir, pos;
	p_dir = GetRandomDirection4pi();
	
	//Generate uniform random position in sphere surface 
	G4RotationMatrix* rotation_gen_ = new G4RotationMatrix();
    	rotation_gen_->rotateX(0*deg);
    	rotation_gen_-> rotateY(0*deg);
    	rotation_gen_->rotateZ(0*deg);
    
    	sphere_gen_  = new SpherePointSampler(radius_,0*cm,G4ThreeVector(0.,0.,0.),rotation_gen_,0,twopi,0,pi);
	pos = sphere_gen_->GenerateVertex("SURFACE");
	
	
	//momentum
	G4double pmod = std::sqrt(energy*energy - mass*mass);
	
	
	G4double px = pmod*p_dir[0];
	G4double py = pmod*p_dir[1];
	G4double pz = pmod*p_dir[2];
	

	G4double time = 0.;
	//Create a new vertex
	G4PrimaryVertex *vertex = new G4PrimaryVertex(pos,time);

	//Create the new primary particle
	G4PrimaryParticle* particle =
			new G4PrimaryParticle(particle_definition_,px,py,pz);
	vertex->SetPrimary(particle);
	event->AddPrimaryVertex(vertex);
	
	
	/*G4cout<<"----------------------------------------------------------------"<<G4endl;
	G4cout << "PRIMARY GENERATED. Event ID: " << event->GetEventID() << G4endl;
	G4cout<<"Energy is: "<<energy/MeV<<G4endl;
	G4cout<<"Kinetic energy is: "<<kinetic_energy/MeV<<G4endl;
	G4cout<<"Mass: "<<mass<<G4endl;
	
	G4cout<<"Energy is: "<<vertex->GetPrimary()->GetTotalEnergy()<<G4endl;
	G4cout<<"position " <<pos<<G4endl;
	G4cout<<"position " <<vertex->GetPosition()<<G4endl;
	G4cout<<"momentum is: "<<pmod<<G4endl;
	G4cout<<"momentum components: px = "<<px << " py = "<<py<<" pz = "<<pz<<G4endl;
	G4cout<<"momentum components: " <<vertex->GetPrimary()->GetTotalMomentum()<<G4endl;
	G4cout<<"----------------------------------------------------------------"<<G4endl;*/
	

}


G4double UserGenerator::GetEnergy()
{
	//Generate random index 
	//Scale by the flux and get an index
	
	G4int rndm_index = GetRandBinIndex(fRandomGeneral_, flux_);
	
	//bool: Sample between bins with gaussian interpolation
	G4double energy = Sample(energy_[rndm_index]*MeV,true,energy_bins_[rndm_index]*MeV);
	return energy;
}


G4ThreeVector UserGenerator::GetRandomDirection4pi()
{
	G4double cosTheta = G4UniformRand()*2 - 1.;
	G4double sinTheta2 = 1. - cosTheta*cosTheta;
	if(sinTheta2 < 0.) sinTheta2 = 0.;
	G4double sinTheta = std::sqrt(sinTheta2);
	
	G4double phi = G4UniformRand()*twopi;
	return G4ThreeVector(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta).unit();
}

