(* ::Package:: *)

(* ::Title:: *)
(*NearHorizonGeoOrbit subpackage <for> KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["NearHorizonGeoOrbit`"];

NearHorizonGeoOrbit::usage = "Some stuffs";

Begin["`Private`"];


(* ::Section:: *)
(*NHEK radial parametrizations*)


(* ::Section:: *)
(*Near-NHEK radial parametrizations*)


(* ::Section:: *)
(*Polar parametrizations*)


(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


Options[NearHorizonGeoOrbit] = {"Parametrization" -> "Mino", "Spacetime" -> "NHEK"}
SyntaxInformation[NearHorizonGeoOrbit] = {"ArgumentsPattern"->{_,_,OptionsPattern[]}};


NearHorizonGeoOrbit[e_, l_, Q_, \[Mu]_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, spacetime},

param = OptionValue["Parametrization"];
spacetime = OptionValue["Spacetime"];

If[spacetime == "NHEK",

];

If[spacetime == "Near-NHEK",

];


]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
