{
  "More Information" -> {
    "This function should be used to get symbolic trajectory or any other orbital information about a near-horizon geodesic of given radial and polar classes. The function returns a NearHorizonGeoOrbitFunction object, from which any information about the geodesic can be accessed, as described below in the corresponding section. See the NearHorzionGeodesics tutorial for more informations and examples."
  },
  "Option Descriptions" ->{
	"Parametrization"->"Parametrization as a function of the Mino type or of the radial coordinate R",
	"Type"->"Type of the geodesic",
	"Style"->"Typographic style of the output",
	"SimplificationRule"->"Simplification rules applied to the expressions produced",
	"Retrograde"->"Must be set to True for a retrograde geodesic",
	"ReplaceLStar"->"Replace Subscript[\[ScriptCapitalL], *] by its value",
	"ReplaceC"->"Replace \[ScriptCapitalC] by its value",
	"ReplaceLNought"->"Replace Subscript[\[ScriptCapitalL], \[SmallCircle]]by its value",
	"ReplaceCNought"->"Replace Subscript[\[ScriptCapitalC], \[SmallCircle]]by its value",
	"ReplaceRoots"->"Replace the polar roots by their values",
	"ReplaceTurningPoints"->"Replace the number of polar turning points m by its value ",
	"CosTheta"->"Gives the trajectory in terms of z=cos \[Theta] instead of \[Theta]",
	"ReplacePhiTheta"->"Replace Subscript[\[CapitalPhi], \[Theta] ]by its value"
  },
  "Numerical Evaluation" -> {
    "NearHorizonGeoOrbit["NHEK", 3.3, 0.01, 0.1, 1]"
    },
  "See Also" -> {"NearHorizonGeoOrbit","NearHorizonGeoOrbitFunction","KerrGeoOrbit"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"NearHorizonGeodesics"}
}