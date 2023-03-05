/// \file UserGenerator.hh

#ifndef USER_GENERATOR_H
#define USER_GENERATOR_H


#include <G4VPrimaryGenerator.hh>
#include <Randomize.hh>

#include "nexus/SpherePointSampler.h"

#include <vector>
#include <G4SystemOfUnits.hh>

class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;
class G4VSolid;


namespace nexus {

  class GeometryBase;

  class UserGenerator: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    UserGenerator();
    /// Destructor
    ~UserGenerator();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a particle in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:
    // Sample the particle Distribution loaded from file
    G4double GetEnergy();
    
    //G4bool CheckOverlap(const G4threeVector& vtx, const G4ThreeVector& dir);

    /// Load in the Muon Angular/Energy Distribution from CSV file
    /// and initialise the discrete flux distribution
    void LoadUserDistribution();
    G4ThreeVector GetRandomDirection4pi();
    
    //G4int GetRandomBin(G4RandGeneral *fRandomGeneral, std::vector<G4double> intensity);
    //G4double SampleBins (G4double sample, G4bool option,G4double bin_width);

  private:
    G4GenericMessenger* msg_;

    G4ParticleDefinition* particle_definition_;

    G4bool user_ene_dist_; ///< Use energy distribution according to input file
    G4bool mono_ene_; ///< generate primaries with monoenergetic distribution
 
    G4String region_; ///< Name of generator region
    
    G4String ene_file_; ///< Name of file with distribution
 
    G4bool IsDistribution_;  ///< Check if the energy distribution is already loaded (don't want memory leaks!! :) )

    const GeometryBase* geom_; ///< Pointer to the detector geometry

    G4VSolid * geom_solid_;

    std::vector<G4double> flux_, energy_,azimuths_,zeniths_; ///< Values of flux, azimuth and zenith from file
    std::vector<G4double> energy_bins_,azimuth_smear_,zenith_smear_;  ///< List of Energy bin smear values
    G4RandGeneral *fRandomGeneral_; ///< Pointer to the RNG flux distribution
    
    G4double world_rad_;
    SpherePointSampler* sphere_gen_;
  };
}
#endif
