{
  "More Information" -> {
    "Computes the components of the body's four-velocity. Components are normalized such that \!\(\*SuperscriptBox[\(u\), \(\[alpha]\)]\)\!\(\*SubscriptBox[\(u\), \(\[alpha]\)]\) = -1."
    },
    "Option Descriptions" ->{
  	  {"Index"->"Choose between contravariant or covariant components.",
	   "Parametrization"-> "Choose between parametrizing the components with Mino Time or the Darwin Parameter."}
    },
  "Basic Examples" -> {
    "contravariantVelocity = KerrGeoVelocity[0.9, 10, 0.2, 0.5];
contravariantFourVelocity[[2]][10]",
    "covariantVelocity = 
  KerrGeoFourVelocity[0.9, 10, 0.2, 0.5, "Index" -> "Covariant"];
covariantVelocity[[2]][10]",
	"DarwinVelocity = 
  KerrGeoFourVelocity[0.9, 10, 0.2, 1, "Parametrization" -> "Darwin"];
DarwinVelocity[[2]][\[Pi]/2]"
  },
  "See Also" -> {"KerrGeoOrbit","KerrGeoFrequencies", "KerrGeoConstantsOfMotion"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"KerrGeodesics"}
}