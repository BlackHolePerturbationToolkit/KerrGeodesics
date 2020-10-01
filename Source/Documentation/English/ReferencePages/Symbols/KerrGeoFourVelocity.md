{
  "More Information" -> {
    "Computes the components of the body's four-velocity. Components are normalized such that Subscript[u, \[Alpha]] u^\[Alpha]=-1."
    },
    "Option Descriptions" ->{
  	  {"Covariant" -> "True for the covariant four-velocity, False for the contravariant (default)",
	   "Parametrization"-> "Choose between parametrizing the components with Mino Time or the Darwin Parameter."}
    },
  "Basic Examples" -> {
    "contravariantVelocity = KerrGeoVelocity[0.9, 10, 0.2, 0.5];
contravariantFourVelocity[[2]][10]",
    "covariantVelocity = 
  KerrGeoFourVelocity[0.9, 10, 0.2, 0.5, "Covariant" -> True];
covariantVelocity[[2]][10]",
	"DarwinVelocity = 
  KerrGeoFourVelocity[0.9, 10, 0.2, 1, "Parametrization" -> "Darwin"];
DarwinVelocity[[2]][\[Pi]/2]"
  },
  "See Also" -> {"KerrGeoOrbit","KerrGeoFrequencies", "KerrGeoConstantsOfMotion"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"KerrGeodesics"}
}