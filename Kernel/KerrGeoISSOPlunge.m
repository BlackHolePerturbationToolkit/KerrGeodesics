(* ::Package:: *)

(* ::Title:: *)
(*KerrGeoOrbit subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`KerrGeoISSOPlunge`"];

KerrGeoISSOPlunge::usage = "KerrGeoISSOPlunge[a,rI] returns a KerrGeoISSOPlungeFunction[..] which stores the orbital trajectory and parameters.";
KerrGeoPlungeISSOFunction::usage = "KerrGeoISSOPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*ISSO Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoISSOPlungeMino[a_, RI_, initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts, assoc,M=1, t, r, \[Theta], \[Phi], \[Epsilon], L, Q, R4, RM, RP, KZ, Z1, Z2, J,velocity},
	
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	Q =- M (RI)^(5/2) ((Sqrt[(RI-RP)(RI-RM)]-2Sqrt[RI])^2-4a^2)/(4a^2 (RI^(3/2)-Sqrt[RI]-Sqrt[(RI-RP)(RI-RM)]));
	\[Epsilon]=Sqrt[a^2 Q-2 RI^3+3 RI^4]/(Sqrt[3] RI^2);
	L=Sqrt[3 a^2 Q-a^2 RI^2-Q RI^2+3 RI^4+a^2 RI^2 \[Epsilon]^2-3 RI^4 \[Epsilon]^2]/RI;
	J = 1-\[Epsilon]^2;
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};
	consts = {\[Epsilon], L, Q};
	R4 = (a^2 Q)/(J*RI^3);
	KZ = a^2*J(Z1^2/Z2^2);
	
r[\[Lambda]_] := ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4));

z[\[Lambda]_]:= Z1*JacobiSN[Z2*\[Lambda] + EllipticK[KZ],KZ];


tr[\[Lambda]_]:= ((a^2+RI^2) (-a L+(a^2+RI^2) \[Epsilon]) \[Lambda])/((RI-RM) (RI-RP))+(2 (R4-RI)^2 \[Epsilon] \[Lambda])/(4+J (R4-RI)^2 \[Lambda]^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[Epsilon] ArcTan[(\[Lambda] (RI-R4)Sqrt[J ]  )/2 ]+(2 (a^2+RM^2) (-a L+a^2 \[Epsilon]+RM^2 \[Epsilon]) ArcTanh[(2  Sqrt[RM-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM])])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP)Sqrt[J])+(2 (a^2+RP^2) (-a L+a^2 \[Epsilon]+RP^2 \[Epsilon]) ArcTanh[(2  Sqrt[RP-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP])])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP));

\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a ((\[Lambda] (RI-R4)Sqrt[J ] (-a L+a^2 \[Epsilon]+RI^2 \[Epsilon]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a L+a^2 \[Epsilon]+RM^2 \[Epsilon]) ArcTanh[(2  Sqrt[RM-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM])])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP))+((-a L+a^2 \[Epsilon]+RP^2 \[Epsilon]) ArcTanh[(2  Sqrt[RP-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP])])/(Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP))) ;

tz[\[Lambda]_]:= 1/J \[Epsilon]  (-Z2 EllipticE[JacobiAmplitude[Z2 \[Lambda] + EllipticK[KZ],KZ],KZ]+(Z2-a^2 J/Z2) EllipticF[JacobiAmplitude[Z2 \[Lambda] + EllipticK[KZ],KZ],KZ]);
\[Phi]z[\[Lambda]_]:= L/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] + EllipticK[KZ],KZ],KZ];


	t=Function[{Global`\[Lambda]}, Evaluate[ a*L Global`\[Lambda]  tr[Global`\[Lambda]] + tz[Global`\[Lambda]]], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda]]] ], Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[-a*\[Epsilon] Global`\[Lambda] + \[Phi]r[Global`\[Lambda]] + \[Phi]z[Global`\[Lambda]]], Listable];



	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {RI,RI,RI,R4},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"rI" -> RI,
		"InitialPhases" -> initPhases
	];
	
	KerrGeoISSOPlungeFunction[a, RI, assoc]

]



(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)





KerrGeoISSOPlunge[a_,RI_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=KerrGeoISSOPlungeMino[a, RI, initPhases]
(*FIXME: add stability check but make it possible to turn it off*)



KerrGeoISSOPlungeFunction /:
 MakeBoxes[kgof:KerrGeoISSOPlungeFunction[a_, RI_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"RI: ", RI}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}],
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]};
  BoxForm`ArrangeSummaryBox[
    KerrGeoISSOPlungeFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


KerrGeoISSOPlungeFunction[a_, RI_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrGeoISSOPlungeFunction[a_,RI_, assoc_][y_?StringQ] := assoc[y]
Keys[g_KerrGeoISSOPlungeFunction]^:=Keys[g[[5]]]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
