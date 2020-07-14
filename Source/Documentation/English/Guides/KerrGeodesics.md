{
 "Application" -> "KerrGeodesics",
 "Package" -> "KerrGeodesics",
 "Title" -> "Geodesics in Kerr spacetime",
 "Summary" -> 
   "Functions for working with bound timelike geodesics in Kerr spacetime.",
 "Description" -> 
   {"The KerrGeodesics package provides functions for computing bound timelike geodesics and their properties in Kerr spacetime. It also includes functions for computing any timelike or null future-directed near-horizon geodesic."},
 "Keywords" -> {"Kerr", "Geodesic"},
 "Label" -> "KerrGeodesics Application",
 "Synonyms" -> {"KerrGeodesics", "KerrGeodesics"},
 "URL" -> "" ,
 "Packages" -> {
   {"Title" -> "Kerr geodesics",
    "DetailedFunctions" -> {
	  {"KerrGeoEnergy", "Computes the orbital energy"},
	  {"KerrGeoAngularMomentum", "Computes the orbit angular momentum about the symmetry axes of the spacetime"},
      {"KerrGeoCarterConstant", "Computes the Carter Constant"},
      {"KerrGeoConstantsOfMotion", "Computes the constants of motion (energy, angular momentum and Carter constant)"},
      {"KerrGeoFrequencies", "Compute the orbital frequencies"},
	  {"KerrGeoBoundOrbitQ", "Checks if the given (numerical) orbital parameters corresponds to a bound orbit"},
	  {"KerrGeoPhotonSphereRadius", "Computes the radius of the photon sphere for given black hole spin and inclination"},
	  {"KerrGeoISCO", "Computes the location of the inner-most stable circular orbit (ISCO)"},
	  {"KerrGeoIBSO", "Computes the location of the inner-most bound spherical orbit (IBSO)"},
	  {"KerrGeoISSO", "Computes the location of the inner-most stable spherical orbit (ISSO)"},
	  {"KerrGeoSeparatrix", "Computes the value p at the separatrix between stable and plunging/scattered orbits"},
	  {"KerrGeoOrbit", "Computes the orbit trajectory in Boyer-Lindquist coordinates"}
    }
   },
   {"Title" -> "Near-horizon geodesics",
    "DetailedFunctions" -> {
	  {"NearHorizonGeoOrbit", "Computes the trajectory and properties of a near-horizon geodesic from numerical orbital parameters."},
	  {"NearHorizonGeoOrbitClass", "Computes the symbolic orbit and properties of a chosen near-horizon geodesic."},
	  {"NearHorizonGeoOrbitFunction", "An object for storing the trajectory and orbital parameters of a near-horizon geodesic."}
    }
   }
 },
 "Tutorials" -> {
   "KerrGeodesics", "NearHorizonGeodesics"
 } 
}
