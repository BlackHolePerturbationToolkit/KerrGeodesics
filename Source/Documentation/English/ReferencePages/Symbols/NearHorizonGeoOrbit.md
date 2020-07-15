{
  "More Information" -> {
    "This function should be used to get numerical trajectory or any other orbital information about a near-horizon geodesic of given orbital parameters. The function determines automatically the radial and polar classes to which the geodesic belongs, and returns a NearHorizonGeoOrbitFunction storing the trajectory and orbital information. See the NearHorizonGeodesics tutorial for more informations and examples."
  },
  "Option Descriptions" ->{
	"Parametrization"->"Parametrization as a function of the Mino type or of the radial coordinate R",
	"RadialMotion"->"Define the initial sign of the radial velocity. Enable the distinction between plunging and outward geodesics",
	"CosTheta"->"Gives the trajectory in terms of z=cos \[Theta] instead of \[Theta]",
	"Numerical"->"Gives the parametrization in terms of numerical expressions",
	"RadialMotion"->"Enable the distinction betwwen plunging and outward geodesics",
	"SimplificationRule"->"Choose the simplification rule used to process the output expressions",
	"ExplicitMass"->"Set the spacetime mass to any desired value"
  },
  "Numerical Evaluation" -> {
    "NearHorizonGeoOrbit["NHEK", 3.3, 0.01, 0.1, 1]"
    },
  "See Also" -> {"NearHorizonGeoOrbitClass","NearHorizonGeoOrbitFunction","KerrGeoOrbit"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"NearHorizonGeodesics"}
}