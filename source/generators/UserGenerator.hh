/// \file UserGenerator.hh

#ifndef USER_GENERATOR_H
#define USER_GENERATOR_H


#include <G4VPrimaryGenerator.hh>
#include <Randomize.hh>

#include "nexus/SpherePointSampler.h"
#include "nexus/RandomUtils.h"

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

    /// Load in the Muon Angular/Energy Distribution from CSV file
    /// and initialise the discrete flux distribution
    void LoadDatafromFile(std::string filename, std::vector<G4double> &intensity, std::vector<G4double> &x, std::vector<G4double> &x_bin);
    void LoadUserDistribution();
    G4ThreeVector GetRandomDirection4pi();
    G4int GetRandomBin(G4RandGeneral *fRandomGeneral, std::vector<G4double> intensity);
    G4double SampleBins (G4double sample, G4bool option,G4double bin_width);

  private:
    G4GenericMessenger* msg_;

    G4ParticleDefinition* particle_definition_;

    G4bool user_ene_dist_; ///< Use muon distribution according to input file
 
    G4String region_; ///< Name of generator region
    G4String ene_file_; ///< Name of file with distribution


    const GeometryBase* geom_; ///< Pointer to the detector geometry

    G4VSolid * geom_solid_;

    std::vector<G4double> flux_, energy_; ///< Values of flux, azimuth and zenith from file
    std::vector<G4double> energy_bins_;  ///< List of Energy bin smear values
    G4RandGeneral *fRandomGeneral_; ///< Pointer to the RNG flux distribution
    
    G4double world_rad_;
    SpherePointSampler* sphere_gen_;
  };
}
#endif
