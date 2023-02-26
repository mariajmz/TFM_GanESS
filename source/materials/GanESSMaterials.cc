// ----------------------------------------------------------------------------
// GanESSMaterials.cc
//
// Definition of materials of common use.
//
// Build material functions
// ----------------------------------------------------------------------------

#include "globals.hh"

#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

namespace gasses
{
  
////------------------------------------------------------------------
///	GAS DENSITY (for T=293K)
////------------------------------------------------------------------
  
  G4double fArDensity_T293(G4double pressure)
  {
     
     //Ar (gas) density for T=293K
     //linear interpolation between a pair of values
     //nist
     //p>50 bar supercritical phase
     
     G4double density;

     const G4int n_points = 11;
     G4double data[n_points][2] = {{  0.0 * bar,  0.0 * kg/m3},
                                   {  5.0 * bar,  8.2268 * kg/m3},
                                   { 10.0 * bar,  16.508 * kg/m3},
                                   { 15.0 * bar,  24.843 * kg/m3},
                                   { 20.0 * bar,  33.231 * kg/m3},
                                   { 25.0 * bar,  41.668 * kg/m3},
                                   { 30.0 * bar,  50.155 * kg/m3},
                                   { 35.0 * bar,  58.689 * kg/m3},
                                   { 40.0 * bar,  67.267 * kg/m3},
                                   { 45.0 * bar,  75.889 * kg/m3},
                                   { 50.0 * bar,  84.551 * kg/m3}};
     G4bool list = false;

     for (G4int i=0; i<n_points-1; ++i) {
    	if  (pressure >= data[i][0] && pressure < data[i+1][0]) {
      		G4double x1 = data[i][0];
      		G4double x2 = data[i+1][0];
      		G4double y1 = data[i][1];
      		G4double y2 = data[i+1][1];
      		density = y1 + (y2-y1)*(pressure-x1)/(x2-x1);
      		list = true;
      		break;
    		}
  	}

     if (!list) {
    	if (pressure == data[n_points-1][0]) {
      	density = data[n_points-1][1];
    	}
     else {throw "Unknown argon density for this pressure!";}
  	}

  return density;
  }


  G4double fXeDensity_T293(G4double pressure)
  {
     
     //Xe (gas) density for T=293K
     //linear interpolation between a pair of values
     //nist
     //p>60 bar supercritical phase
     
     G4double density;

     const G4int n_points = 11;
     G4double data[n_points][2] = {{  0.0 * bar,  0.0 * kg/m3},
                                   {  5.0 * bar,  27.722 * kg/m3},
                                   { 10.0 * bar,  57.160 * kg/m3},
                                   { 15.0 * bar,  88.621 * kg/m3},
                                   { 20.0 * bar,  122.51 * kg/m3},
                                   { 25.0 * bar,  159.35 * kg/m3},
                                   { 30.0 * bar,  199.92 * kg/m3},
                                   { 35.0 * bar,  245.31 * kg/m3},
                                   { 40.0 * bar,  297.25 * kg/m3},
                                   { 45.0 * bar,  358.72 * kg/m3},
                                   { 50.0 * bar,  435.52 * kg/m3}};
     G4bool list = false;

     for (G4int i=0; i<n_points-1; ++i) {
    	if  (pressure >= data[i][0] && pressure < data[i+1][0]) {
      		G4double x1 = data[i][0];
      		G4double x2 = data[i+1][0];
      		G4double y1 = data[i][1];
      		G4double y2 = data[i+1][1];
      		density = y1 + (y2-y1)*(pressure-x1)/(x2-x1);
      		list = true;
      		break;
    		}
  	}

     if (!list) {
    	if (pressure == data[n_points-1][0]) {
      	density = data[n_points-1][1];
    	}
     else {throw "Unknown xenon density for this pressure!";}
  	}

     return density;
  }
  
  G4double fKrDensity_T293(G4double pressure)
  {
     
     //Kr (gas) density for T=293K
     //linear interpolation between a pair of values
     //nist
     //p>60 bar supercritical phase
     
     G4double density;

     const G4int n_points = 11;
     G4double data[n_points][2] = {{  0.0 * bar,  0.0 * kg/m3},
                                   {  5.0 * bar,  17.388 * kg/m3},
                                   { 10.0 * bar,  35.162 * kg/m3},
                                   { 15.0 * bar,  53.336 * kg/m3},
                                   { 20.0 * bar,  71.927 * kg/m3},
                                   { 25.0 * bar,  90.947 * kg/m3},
                                   { 30.0 * bar,  110.41 * kg/m3},
                                   { 35.0 * bar,  130.34 * kg/m3},
                                   { 40.0 * bar,  150.74 * kg/m3},
                                   { 45.0 * bar,  171.63 * kg/m3},
                                   { 50.0 * bar,  193.03 * kg/m3}};
     G4bool list = false;

     for (G4int i=0; i<n_points-1; ++i) {
    	if  (pressure >= data[i][0] && pressure < data[i+1][0]) {
      		G4double x1 = data[i][0];
      		G4double x2 = data[i+1][0];
      		G4double y1 = data[i][1];
      		G4double y2 = data[i+1][1];
      		density = y1 + (y2-y1)*(pressure-x1)/(x2-x1);
      		list = true;
      		break;
    		}
  	}

     if (!list) {
    	if (pressure == data[n_points-1][0]) {
      	density = data[n_points-1][1];
    	}
     else {throw "Unknown argon density for this pressure!";}
  	}

  return density;
  }

////------------------------------------------------------------------
///	GAS DEFINITONS
////------------------------------------------------------------------

  //Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
	

  //pressure, temperature and density??
  
  G4Material* fArgas(G4double pressure, G4double temperature)
  {
   G4String name = "fArgas";
   G4int ncomp = 1;
   G4double density = fArDensity_T293(pressure); 
   G4Material* mat = new G4Material(name, density, ncomp, kStateGas,temperature,pressure);
   //Ar
   G4Isotope* gasAr_36 = new G4Isotope("Ar36",18,36, 35.967*g/mole);
   G4Isotope* gasAr_38 = new G4Isotope("Ar38",18,38, 37.962*g/mole);
   G4Isotope* gasAr_40 = new G4Isotope("Ar40",18,40, 39.962*g/mole);
   
   G4Element* fArgas_Ar = new G4Element("Argon", "Ar", 3);
   fArgas_Ar->AddIsotope(gasAr_36,0.00336);
   fArgas_Ar->AddIsotope(gasAr_38,0.00064);
   fArgas_Ar->AddIsotope(gasAr_40,0.996);
   
   mat->AddElement(fArgas_Ar, 1.0);
   return mat;
  
  }
  
  //pressure, temperature and density??
  G4Material* fKrgas(G4double pressure, G4double temperature)
  {
  G4String name = "fKrgas";
  G4int ncomp = 1;
  G4double density = fKrDensity_T293(pressure); 
  G4Material* mat = new G4Material(name, density, ncomp, kStateGas,temperature,pressure);
  //Kr
  G4Isotope* gasKr_78 = new G4Isotope("Kr78",36,78,77.920*g/mole);
  G4Isotope* gasKr_80 = new G4Isotope("Kr80",36,80,79.916*g/mole);
  G4Isotope* gasKr_82 = new G4Isotope("Kr82",36,82,81.913*g/mole);
  G4Isotope* gasKr_83 = new G4Isotope("Kr83",36,83,82.914*g/mole);
  G4Isotope* gasKr_84 = new G4Isotope("Kr84",36,84,83.911*g/mole);
  G4Isotope* gasKr_86 = new G4Isotope("Kr86",36,86,85.910*g/mole);
  
  G4Element* fKrgas_Kr = new G4Element("Kripton", "Kr", 6);
  fKrgas_Kr->AddIsotope(gasKr_78,0.355*perCent);
  fKrgas_Kr->AddIsotope(gasKr_80,2.286*perCent);
  fKrgas_Kr->AddIsotope(gasKr_82,11.593*perCent);
  fKrgas_Kr->AddIsotope(gasKr_83,11.500*perCent);
  fKrgas_Kr->AddIsotope(gasKr_84,56.987*perCent);
  fKrgas_Kr->AddIsotope(gasKr_86,17.29*perCent);
   
  mat->AddElement(fKrgas_Kr, 1.0);
  return mat;
  
  }
  
  
  //pressure, temperature and density??
  G4Material* fXegas(G4double pressure, G4double temperature)
  {
  
  G4String name = "fXegas";
  G4int ncomp = 1;
  G4double density = fXeDensity_T293(pressure); 
  G4Material* mat = new G4Material(name, density, ncomp, kStateGas,temperature,pressure);
  //G4Material* mat = new G4Material(name, density, ncomp, kStateGas);
  
  //Xe
  G4Isotope* gasXe_124 = new G4Isotope("Xe124",54,124,123.905*g/mole);
  G4Isotope* gasXe_126 = new G4Isotope("Xe126",54,126,125.904*g/mole);
  G4Isotope* gasXe_128 = new G4Isotope("Xe128",54,128,127.903*g/mole);
  G4Isotope* gasXe_129 = new G4Isotope("Xe129",54,129,128.904*g/mole);
  G4Isotope* gasXe_130 = new G4Isotope("Xe130",54,130,129.903*g/mole);
  G4Isotope* gasXe_131 = new G4Isotope("Xe131",54,131,130.905*g/mole);
  G4Isotope* gasXe_132 = new G4Isotope("Xe132",54,132,131.904*g/mole);
  G4Isotope* gasXe_134 = new G4Isotope("Xe134",54,134,133.905*g/mole);
  G4Isotope* gasXe_136 = new G4Isotope("Xe136",54,136,135.907*g/mole);
  
  G4Element* fXegas_Xe = new G4Element("Xenon", "Xe", 9);
  fXegas_Xe->AddIsotope(gasXe_124,0.000952);
  fXegas_Xe->AddIsotope(gasXe_126,0.000890);
  fXegas_Xe->AddIsotope(gasXe_128,0.019102);
  fXegas_Xe->AddIsotope(gasXe_129,0.264006);
  fXegas_Xe->AddIsotope(gasXe_130,0.040710);
  fXegas_Xe->AddIsotope(gasXe_131,0.212324);
  fXegas_Xe->AddIsotope(gasXe_132,0.269086);
  fXegas_Xe->AddIsotope(gasXe_134,0.104357);
  fXegas_Xe->AddIsotope(gasXe_136,0.088573);
   
  mat->AddElement(fXegas_Xe, 1.0);
  
  
  return mat;
  }
  
}
