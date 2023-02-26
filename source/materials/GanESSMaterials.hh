// ----------------------------------------------------------------------------
// GanESSMaterials.hh
//
// Definition of materials of common use.
//
// ----------------------------------------------------------------------------

#ifndef GanESSMaterials_h
#define GanESSMaterials_h 1

#include "globals.hh"

#include <vector>

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

class G4Material;

namespace gasses
{

	G4double fArDensity_T293(G4double pressure);
	G4double fKrDensity_T293(G4double pressure);
	G4double fXeDensity_T293(G4double pressure);

	G4Material* fXegas(G4double pressure, G4double temperature);
	G4Material* fArgas(G4double pressure, G4double temperature);
	G4Material* fKrgas(G4double pressure, G4double temperature);
	
	
	
}
#endif
