{
  "More Information" -> {
    "Computes the components of the body's four-velocity."
    },
    "Option Descriptions" ->{
  	  {"Index"->"Choose between contravariant or covariant components.",
	   "Parametrization"-> "Choose between parametrizing the components with Mino Time or the Darwin Parameter."}
    },
  "Basic Examples" -> {
    "contravariantVelocity = KerrGeoVelocity[0.9, 10, 0.2, 0.5];
contravariantVelocity["ur"][10]",
    "covariantVelocity = 
  KerrGeoVelocity[0.9, 10, 0.2, 0.5, "Index" -> "Covariant"];
covariantVelocity["ur"][10]",
	"DarwinVelocity = 
  KerrGeoVelocity[0.9, 10, 0.2, 1, "Parametrization" -> "Darwin"];
DarwinVelocity["ur"][\[Pi]/2]"
  },
  "See Also" -> {"KerrGeoOrbit","KerrGeoFrequencies", "KerrGeoConstantsOfMotion"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"KerrGeodesics"}
}