{
  "More Information" -> {
    "Computes the orbital frequencies"
    },
    "Option Descriptions" ->{
  	  "Time"->"Choose between calculating the frequencies w.r.t BoyerLindquist or Mino time"
    },
  "Basic Examples" -> {
    "KerrGeoFrequencies[0.9`20, 5, 0.7`20, Cos[\[Pi]/4]]",
    "KerrGeoFrequencies[0.9`20, 5, 0.7`20, Cos[\[Pi]/4], Time->\"Mino\"]",
	"KerrGeoFrequencies[0, r0, 0, 1]"
  },
  "See Also" -> {"KerrGeoAngularMomentum","KerrGeoEnergy", "KerrGeoCarterConstant", "KerrGeoConstantsOfMotion"},
  "More About" -> {"KerrGeodesics"},
  "Tutorials" -> {"KerrGeodesics"}
}