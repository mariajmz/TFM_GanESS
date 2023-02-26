//-----------------------------------------------------------------------------
// UserGenerator.cc
//   
// This class is the primary generator for particles with a user defined energy 
// histogram      
//
//-----------------------------------------------------------------------------

#include "DetectorConstruction.h"
#include "GeometryBase.h"
#include "MuonsPointSampler.h"
#include "AddUserInfoToPV.h"
#include "FactoryBase.h"
#include "RandomUtils.h"
#include "IOUtils.h"

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
user_ene_dist_(true), user_dir_{}, energy_min_(0.),
energy_max_(0.),geom_(0), geom_solid_(0)
{

 msg_ = new G4GenericMessenger(this, "/Generator/UserGenerator/",
				"Control commands for user generator.");
				
G4GenericMessenger::Command& min_energy =
    msg_->DeclareProperty("min_energy", energy_min_, "Set minimum kinetic energy of the particle.");
  min_energy.SetUnitCategory("Energy");
  min_energy.SetParameterName("min_energy", false);
  min_energy.SetRange("min_energy>0.");

G4GenericMessenger::Command& max_energy =
    msg_->DeclareProperty("max_energy", energy_max_, "Set maximum kinetic energy of the particle");
  max_energy.SetUnitCategory("Energy");
  max_energy.SetParameterName("max_energy", false);
  max_energy.SetRange("max_energy>0.");
  
msg_->DeclareProperty("region", region_,
			"Set the region of the geometry where the vertex will be generated.");

msg_->DeclareProperty("user_ene_dist", user_ene_dist_,
			"Distribute particles energies according to file?");
			
msg_->DeclareProperty("angle_dist", dist_name_,
			"Name of the angular distribution histogram.");
			
msg_->DeclareProperty("energy_file", ene_file_,
			"Name of the file containing angular distribution.");
			
msg_->DeclareProperty("energy_dist", ene_name_,
			"Name of the angular distribution histogram.");	

DetectorConstruction* detconst = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();
				
}

UserGenerator::~UserGenerator()
{
  delete msg_;
}


Generator::Generator()
 : G4VUserPrimaryGeneratorAction(),fParticleGun(0)
{

  //to run with particle gun    
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4double Rmax = 5*cm;//must contain gas volume
  fRmax3 = Rmax*Rmax*Rmax;
  
  G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particle);
  GetHistogram();
  
}

void UserGenerator::LoadUserDistribution()
{

	LoadDataFromFile(ene_file_,flux_,energy_,energy_bin_);
	
	//convert vector to array
	auto arr_flux = flux.data();
	
	//Initialice random number generator based on flux distribution
	fRandomGeneral_ = new G4RandGeneral(arr_flux,flux_.size());

}


//file format: value=intensity in bin, x=bin centre, x_bin=bin width
void UserGenerator::LoadDatafromFile(std::string filename, std::vector<G4double> &intensity, std::vector<G4double> &x, std::vector<G4double> &x_bin)
{

	//open file
	std::ifstream File_(filename);
	
	//Check if the file is properly opened
	if (!FileIn_.is_open()){
	      G4cout<<"Could not read file"<<G4endl;
    	}
    	
      	// Read the Data from the file as strings
    	std::string s_header, s_intensity, s_x;
    	std::string s_x_bin;

    	// Loop over the lines in the file and add the values to a vector
    	while (FileIn_.peek()!=EOF) {

      		std::getline(FileIn_, s_header, ',');

      		if (s_header == "value"){

       			std::getline(File_, s_instensity, ',');
        		std::getline(File_, s_x, ',');
        		std::getline(File_, s_x_bin, '\n');

        		value.push_back(stod(s_intensity));
        		x.push_back(stod(s_x));
       		 	x_smear.push_back(stod(s_x_bin));
      			}

    	} 

    File_.close();
}


void UserGenerator::GeneratePrimaryVertex(G4Event* event)
{

if (user_ener_dist_){
      std::cout << "Generating particles using user distribution loaded from file" << std::endl;
      LoadUserDistribution();
      }
//Generate uniform random direction in 4pi: RandomDirectionInRange


//Generate uniform random position in sphere surface 


//Create a new vertex
G4VPrimaryVertex *vertex = new G4PrimaryVertex(position, time);

//Create the new primary particle

}


G4double UserGenerator::GetEnergy()
{
	//Generate random index 
	//Scale by the flux and get and index
	G4int rndm_index = GetRandBinIndex(fRandomGeneral_,flux_);
	//Sample between bins with gaussian interpolation
	G4double energy = Sample(energy_[rndm_index]*MeV,true,energy_bin_[rndm_index]*MeV);
	
	return energy;
}

