(* ::Package:: *)

(* ::Title:: *)
(*KerrGeoOrbit subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`KerrGeoPlunge`",
	{"KerrGeodesics`ConstantsOfMotion`",
	 "KerrGeodesics`OrbitalFrequencies`",
	 "KerrGeodesics`SpecialOrbits`",
	 "KerrGeodesics`FourVelocity`"}];


KerrGeoPlunge::usage = "Takes either KerrGeoPlunge[a, {En,\[ScriptCapitalL],\[ScriptCapitalQ]}] or KerrGeoPlunge[a, ISSO , RI] for generic Plunges or ISSO plunges respectively and returns a KerrGeoPlungeFunction[..] which stores the orbital trajectory and parameters. Here the ISSO plunges are paramaterised in terms of the choice of the radius of the ISSO for a given a there are range of allowed RI's which correspond to differeing inclinations in the prograde and retrograde directions, if an RI outside this range is given the used is provided with the range of allowed values to re-run the funciton.";
KerrGeoPlungeFunction::usage = "KerrGeoPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";
Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*Homoclinic Plunges*)


KerrGeoPlunge::r0outofbounds = "Intial radius `1` is not between ISSO radius `2`, and inner radial root `3`.";
KerrGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";
KerrGeoPlunge::InnerOutter = "The choice of inner or outter trajectory must be given by a -1 (plunging) or +1 (bound), you provided `1`";


KerrGeoUSOPlunge[a_, RS_, \[ScriptCapitalQ]_, InnerOutter_, initCoords_] := Module[
	{consts, \[Theta]0, z, r0,Mino, z0, \[Phi]0, initConditions, t0, assoc,M=1,\[Chi]INC, Minoz, MinozFunc, MinoRFunc,\[CapitalUpsilon]z,MCMP,MCPP,MCPM,MCMM,DN,RISSOMIN,RISSOMAX ,roots,\[CapitalUpsilon],t,r1,r2 ,RI , r, \[Theta], \[Phi], \[ScriptCapitalE], \[ScriptCapitalL], R4, RM, RP,t\[Delta], kz, Z1, Z2, J,velocity,MinoR,tr,tz,\[Phi]r,\[Phi]z,LRoot},

	RM = 1-Sqrt[1-a^2];
	RP = 1+Sqrt[1-a^2];
	
	\[CapitalUpsilon] = M RS^5 - \[ScriptCapitalQ] (RS - 3 M )RS^3 + a^2 \[ScriptCapitalQ]^2;
	\[ScriptCapitalE] = (RS^3 (RS - 2 M ) - a (a \[ScriptCapitalQ]  - Sqrt[\[CapitalUpsilon]]))/(RS^2 Sqrt[RS^3 (RS - 3 M ) - 2a (a \[ScriptCapitalQ]  - Sqrt[\[CapitalUpsilon]])]);
	\[ScriptCapitalL]=-((2 M a RS^3+(RS^2+a^2)(a \[ScriptCapitalQ]  -Sqrt[\[CapitalUpsilon]]))/(RS^2 Sqrt[RS^3 (RS - 3 M ) - 2a (a \[ScriptCapitalQ]  - Sqrt[\[CapitalUpsilon]])]));

	
	J = 1-\[ScriptCapitalE]^2;
	
	PolyR[r_]:= ((\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2*r+a^2)(r^2+ \[ScriptCapitalQ]+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2 ));
	
	roots = Re[Values[Solve[PolyR[r]==0,r]]];
	
	If[InnerOutter!=-1 Or InnerOutter!=+1,
		Message[KerrGeoPlunge::InnerOutter,InnerOutter];
		Return[$Failed]
		];
	
	
	If[InnerOutter==-1,{r2,r1} ={roots[[1]][[1]], roots[[4]][[1]]}];
	If[InnerOutter==1,{r2,r1} ={roots[[4]][[1]], roots[[1]][[1]]}];
	
	
	
	If[a==0||Round[\[ScriptCapitalQ] ,10^-14]== 0, {Z1,Z2}={0,\[ScriptCapitalL]}, {Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]}];

	kz = a*Sqrt[J](Z1/Z2);
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[kz^2]);
	If[Arg==-(\[Pi]/2), \[ScriptCapitalL]= -\[ScriptCapitalL]];
	If[a<0,\[ScriptCapitalL]=-\[ScriptCapitalL]];
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalQ]"->\[ScriptCapitalQ]|>;
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,r2,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	z0 = Cos[\[Theta]0];
	
	Minoz[z_] :=InverseJacobiSN[z/Z1,kz^2]/Z2;
	
	If[a==0, Minoz[z_]:=0];
	If[Round[\[ScriptCapitalQ],10^-14]==0, Minoz[z_]:=0];
	
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];

	MinoR[r_] :=( 2 ArcTanh[(Sqrt[r-r2] Sqrt[r1-RS])/(Sqrt[r1-r] Sqrt[RS-r2])])/(Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[(1-\[ScriptCapitalE]^2)]);

	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r]-MinoR[r0]  ],Listable];

	

	r[\[Lambda]_] := (r2 (r1-RS)+r1 (RS-r2)   Tanh[1/2  Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]]^2)/(r1-RS+(RS-r2)  Tanh[1/2  Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]]^2);


	
	z[\[Lambda]_]:= If[Round[\[ScriptCapitalQ],10^-14]==0,0,Z1*JacobiSN[Z2*\[Lambda] ,kz^2]];
	
	
	If[Round[\[ScriptCapitalQ],10^-14]==0,Minoz[z_] :=0; z[\[Lambda]_]:=0;
	 
	];
	
		
	
	If[InnerOutter==-1,tr[\[Lambda]_]:=1/Sqrt[1-\[ScriptCapitalE]^2] (((a^2+RS^2) Sqrt[1-\[ScriptCapitalE]^2] (a^2 \[ScriptCapitalE]+RS^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) \[Lambda])/((RS-RM) (RS-RP))-(r1+r2+2 (RM+RP+RS)) \[ScriptCapitalE] ArcTan[(Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])/Sqrt[r1-RS]]+((a^2+RM^2) (a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[RM-r2] Sqrt[r1-RS]+Sqrt[r1-RM] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[RM-r2] Sqrt[r1-RS]-Sqrt[r1-RM] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(2 Sqrt[r1-RM] Sqrt[-r2+RM] (RM-RP) (RM-RS))+((a^2+RP^2) (a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[RP-r2] Sqrt[r1-RS]+Sqrt[r1-RP] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[RP-r2] Sqrt[r1-RS]-Sqrt[r1-RP] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(2 Sqrt[r1-RP] Sqrt[RP-r2] (RP-RM) (RP-RS))+((r1-r2) Sqrt[(r1-RS) (RS-r2)] \[ScriptCapitalE] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])/(r1-RS-(r2-RS) Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]]^2))];
	If[InnerOutter==+1,tr[\[Lambda]_]:=1/Sqrt[1-\[ScriptCapitalE]^2] (((a^2+RS^2) Sqrt[1-\[ScriptCapitalE]^2] (a^2 \[ScriptCapitalE]+RS^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) \[Lambda])/((RS-RM) (RS-RP))-(r1+r2+2 (RM+RP+RS)) \[ScriptCapitalE] ArcTan[(Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])/Sqrt[r1-RS]]+((a^2+RM^2) (a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[RM-r2] Sqrt[r1-RS]+Sqrt[r1-RM] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[RM-r2] Sqrt[r1-RS]-Sqrt[r1-RM] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(2 Sqrt[r1-RM] Sqrt[-r2+RM] (RM-RP) (RM-RS))+((a^2+RP^2) (a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[RP-r2] Sqrt[r1-RS]+Sqrt[r1-RP] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[RP-r2] Sqrt[r1-RS]-Sqrt[r1-RP] Sqrt[RS-r2] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(2 Sqrt[r1-RP] Sqrt[RP-r2] (RP-RM) (RP-RS))-((r1-r2) Sqrt[(r1-RS) (RS-r2)] \[ScriptCapitalE] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])/(r1-RS-(r2-RS) Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]]^2))];
	
	
	\[Phi]r[\[Lambda]_]:= 1/Sqrt[1-\[ScriptCapitalE]^2] 2 a ((Sqrt[1-\[ScriptCapitalE]^2] (a^2 \[ScriptCapitalE]+RS^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) \[Lambda])/(2 (-RM+RS) (-RP+RS))+((a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[(RM-r2) (r1-RS)]+Sqrt[(r1-RM) (RS-r2)] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[(RM-r2) (r1-RS)]-Sqrt[(r1-RM) (RS-r2)] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(4 Sqrt[r1-RM] Sqrt[RM-r2] (RM-RP) (RM-RS))+((a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) Log[(Sqrt[(RP-r2) (r1-RS)]+Sqrt[(r1-RP) (RS-r2)] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2/(Sqrt[(RP-r2) (r1-RS)]-Sqrt[(r1-RP) (RS-r2)] Tanh[1/2 Sqrt[r1-RS] Sqrt[RS-r2] Sqrt[1-\[ScriptCapitalE]^2] \[Lambda]])^2])/(4 Sqrt[r1-RP] Sqrt[RP-r2] (RP-RM) (RP-RS))) ;
	
	
	tz[\[Lambda]_]:= 1/J \[ScriptCapitalE]  (-Z2 EllipticE[JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2]+(Z2-a^2 J/Z2) Z2 \[Lambda]);
	\[Phi]z[\[Lambda]_]:= \[ScriptCapitalL]/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2];
	
	t=Function[{Global`\[Lambda]}, Evaluate[ a*\[ScriptCapitalL] Global`\[Lambda] + tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]] -  tr[MinoR[r0]] - tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[-a*\[ScriptCapitalE] Global`\[Lambda] + \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];
	MCPP = -Abs[MinoR[RP]] -MinoR[r0];
	MCMP = -Abs[MinoR[RM]] -MinoR[r0];
	MCMM = Abs[MinoR[RM]] - MinoR[r0];
	MCPM = Abs[MinoR[RP]] - MinoR[r0];
	assoc = Association[
		"a" -> a,
		"rI" -> RI,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> 0,
		"PolarFrequency"-> \[CapitalUpsilon]z,
		"Frequencies"-> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"-> 0, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]z|>,
		"RadialPeriod"-> Infinity,
		"PolarPeriod"-> 2\[Pi]/\[CapitalUpsilon]z,
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> Infinity, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]z|>,
		"RadialRoots"-> roots[[;;,1]] ,
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialMinoTime"-> MinoR,
		"AllowedISSORangeForGivenSpin"->  <|"\!\(\*SubscriptBox[\(r\), \(IMIN\)]\)"-> RISSOMIN, "\!\(\*SubscriptBox[\(r\), \(IMAX\)]\)" -> RISSOMAX|>,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(I\)], \(\)]\)"-> \[PlusMinus] Infinity(* Time to ISSO *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(inner\)], \(\)]\)"-> Abs[MinoR[R4]] - MinoR[r0](* Time to Inner*)|>,
		"HorizonCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(+\)]\)"-> MCPP(* Time ingoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(+\)]\)"-> MCMP(* Time ingoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(-\)]\)"-> MCMM(* Time outgoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(-\)]\)"-> MCPM(* Time outgoing branch crosses outer horizon *)
			|>,
	"PolarRoots"-> {Z1,Z2},
	"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];

	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]



(* ::Subsection:: *)
(*ISSO Plunges*)


KerrGeoPlunge::r0outofbounds = "Intial radius `1` is not between ISSO radius `2`, and inner radial root `3`.";
KerrGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";


KerrGeoISSOPlunge[a_, PlungeType_  ,Arg_, initCoords_] := Module[
	{consts, \[Theta]0, z, r0,Mino, z0, \[Phi]0, initConditions, t0, assoc,M=1,\[Chi]INC, Minoz, MinozFunc, MinoRFunc,\[CapitalUpsilon]z,MCMP,MCPP,MCPM,MCMM,DN,RISSOMIN,RISSOMAX ,t ,RI , r, \[Theta], \[Phi], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], R4, RM, RP,t\[Delta], kz, Z1, Z2, J,velocity,MinoR,tr,tz,\[Phi]r,\[Phi]z,LRoot},
	

	RM = 1-Sqrt[1-a^2];
	RP = 1+Sqrt[1-a^2];

	DN = (27-45 a^2+17 a^4+a^6+8 a^3 (1-a^2))^(1/3);
	RISSOMIN=  3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]-1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	RISSOMAX=3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]+1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	
	(*This finds the correct sign of the angular momentum by determining the root of where the non-analytic form of the equation touches the axis on a \[ScriptCapitalL](RI) plot*)
	If[PlungeType == "ISSORadialParam",
		RI = Arg;
		If[a==0,\[ScriptCapitalQ]=0,\[ScriptCapitalQ] = -M (RI)^(5/2) ((Sqrt[(RI-RP)(RI-RM)]-2Sqrt[RI])^2-4a^2)/(4a^2 (RI^(3/2)-Sqrt[RI]-Sqrt[(RI-RP)(RI-RM)]));
		]; 
		\[ScriptCapitalE] = Sqrt[a^2 \[ScriptCapitalQ]-2 RI^3+3 RI^4]/(Sqrt[3] RI^2);
		\[ScriptCapitalL] = Sqrt[3 a^2 \[ScriptCapitalQ]-a^2 RI^2-\[ScriptCapitalQ] RI^2+3 RI^4+a^2 RI^2 \[ScriptCapitalE]^2-3 RI^4 \[ScriptCapitalE]^2]/RI;
	RootCheckFunc[y_]:= (a^6-3 a^2 y^(5/2) (3 y^(3/2)+6 Sqrt[a^2+(-2+y) y]-2 y Sqrt[a^2+(-2+y) y])+y^(9/2) (20 Sqrt[y]-11 y^(3/2)+5 Sqrt[a^2+(-2+y) y]+3 y Sqrt[a^2+(-2+y) y])+a^4 (-4 y+3 y^2+Sqrt[y] Sqrt[a^2+(-2+y) y]+3 y^(3/2) Sqrt[a^2+(-2+y) y]));
	If[a==0,LRoot=6,y/.FindRoot[RootCheckFunc[y],{y,6}]];
	If[a!=0,If[RI>LRoot ,\[ScriptCapitalL]=-\[ScriptCapitalL]]];
	];

	If[PlungeType == "ISSOIncParam",
		\[Chi]INC = Cos[Arg - \[Pi]/2];
		RI = KerrGeoISSO[a,\[Chi]INC];
		{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]} = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[1;;3]]
		];
	J = 1-\[ScriptCapitalE]^2;
	If[a==0, {Z1,Z2}={0,\[ScriptCapitalL]}, {Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]}];
	R4 = (a^2 \[ScriptCapitalQ])/(J*RI^3);
	kz = a*Sqrt[J](Z1/Z2);
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[kz^2]);
	If[Arg==-(\[Pi]/2), \[ScriptCapitalL]= -\[ScriptCapitalL]];
	If[a<0,\[ScriptCapitalL]=-\[ScriptCapitalL]];
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalQ]"->\[ScriptCapitalQ]|>;
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,R4,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0>RI||r0<R4,
		Message[KerrGeoPlunge::r0outofbounds,r0,RI,R4];
		Return[$Failed]
		];
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[KerrGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
	z0 = Cos[\[Theta]0];
	Minoz[z_] :=InverseJacobiSN[z/Z1,kz^2]/Z2;
	If[a==0, Minoz[z_]:=0];
	If[PlungeType == "ISSOIncParam",
		If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
		MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];];
	
	If[PlungeType == "ISSORadialParam",
		If[Round[Arg - RISSOMIN,10^-14]==0, Minoz[z_] :=0];
		MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
		If[Round[Arg - RISSOMAX,10^-14]==0,Minoz[z_] :=0];
		MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];];
	MinoR[r_] :=-(2 Sqrt[r-R4])/Sqrt[J (RI-r) (R4-RI)^2];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r]-MinoR[r0]  ],Listable];

	r[\[Lambda]_] := ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4));
	z[\[Lambda]_]:= Z1*JacobiSN[Z2*\[Lambda] ,kz^2];
	If[a!=0,tr[\[Lambda]_]:= ((a^2+RI^2) (-a \[ScriptCapitalL]+(a^2+RI^2) \[ScriptCapitalE])(\[Lambda]) )/((RI-RM) (RI-RP))+(2 ((R4-RI)^2) \[ScriptCapitalE] (\[Lambda]) )/(4+J (R4-RI)^2 (\[Lambda]) ^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[ScriptCapitalE] ArcTan[((\[Lambda])  (RI-R4)Sqrt[J ]  )/2 ]+( (a^2+RM^2) (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP)Sqrt[J])+( (a^2+RP^2) (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]) Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP))];
	If[a==0,tr[\[Lambda]_]:= ((a^2+RI^2) (-a \[ScriptCapitalL]+(a^2+RI^2) \[ScriptCapitalE])(\[Lambda]) )/((RI-RM) (RI-RP))+(2 ((R4-RI)^2) \[ScriptCapitalE] (\[Lambda]) )/(4+J (R4-RI)^2 (\[Lambda]) ^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[ScriptCapitalE] ArcTan[((\[Lambda])  (RI-R4)Sqrt[J ]  )/2 ]+( (a^2+RP^2) (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]) Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP))];
	If[a!=0,\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a (((\[Lambda]) (RI-R4)Sqrt[J ] (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RI^2 \[ScriptCapitalE]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(2Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP))+((-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(2Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP))) ] ;
	If[a==0,\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a (((\[Lambda]) (RI-R4)Sqrt[J ] (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RI^2 \[ScriptCapitalE]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(2Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP)))  ]  ;
	tz[\[Lambda]_]:= 1/J \[ScriptCapitalE]  (-Z2 EllipticE[JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2]+(Z2-a^2 J/Z2) Z2 \[Lambda]);
	\[Phi]z[\[Lambda]_]:= \[ScriptCapitalL]/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2];

	t=Function[{Global`\[Lambda]}, Evaluate[ a*\[ScriptCapitalL] Global`\[Lambda] + tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]]-tr[MinoR[r0]]-tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[-a*\[ScriptCapitalE] Global`\[Lambda] + \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];
	MCPP = -Abs[MinoR[RP]] -MinoR[r0];
	MCMP = -Abs[MinoR[RM]] -MinoR[r0];
	MCMM = Abs[MinoR[RM]] - MinoR[r0];
	MCPM = Abs[MinoR[RP]] - MinoR[r0];
	assoc = Association[
		"a" -> a,
		"rI" -> RI,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> 0,
		"PolarFrequency"-> \[CapitalUpsilon]z,
		"Frequencies"-> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"-> 0, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]z|>,
		"RadialPeriod"-> Infinity,
		"PolarPeriod"-> 2\[Pi]/\[CapitalUpsilon]z,
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> Infinity, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]z|>,
		"RadialRoots"-> {RI,RI,RI,R4},
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialMinoTime"-> MinoR,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(I\)], \(\)]\)"-> \[PlusMinus] Infinity(* Time to ISSO *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(inner\)], \(\)]\)"-> Abs[MinoR[R4]] - MinoR[r0](* Time to Inner*)|>,
		"HorizonCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(+\)]\)"-> MCPP(* Time ingoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(+\)]\)"-> MCMP(* Time ingoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(-\)]\)"-> MCMM(* Time outgoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(-\)]\)"-> MCPM(* Time outgoing branch crosses outer horizon *)
			|>,
	"PolarRoots"-> {Z1,Z2},
	"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]



(* ::Subsection:: *)
(*Real Roots Generic Plunge (Mino)*)


KerrGeoPlunge::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
KerrGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";


Options[KerrGeoRealPlungeMino] = {"Roots"-> Automatic};
KerrGeoRealPlungeMino[a_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalQ]_ , initCoords_, OptionsPattern[Options[KerrGeoRealPlungeMino]]] := Module[
	{consts,assoc,M,Mino,t0,r0,\[Phi]0,kr,hR,hP,hM,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]tr,\[CapitalUpsilon]tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z,z0,\[Theta]0, \[Theta], J,\[Phi], R1,R2,R3,R4,MCMP,MCPP,MCPM,MCMM, Minoz, MinozFunc, MinoRFunc , RM, RP,ROOTS,RealRoots,ComplexRoots, kz, Z1, Z2, MinoR},

	M=1;
	J = 1-\[ScriptCapitalE]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalQ]"->\[ScriptCapitalQ]|>;
	 
	If[
		OptionValue["Roots"]===Automatic,
		ROOTS = Sort[r/.Solve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2M*r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r]],
		ROOTS = OptionValue["Roots"]
		];
	If[Length[Select[ROOTS,#>=RP&]]>=3,
	R1= ROOTS[[4]];
	R2= ROOTS[[3]];
	R3= ROOTS[[2]];
	R4= ROOTS[[1]];];
	
	If[Length[Select[ROOTS,#>=RP&]]<3,
	R1= ROOTS[[2]];
	R2= ROOTS[[1]];
	R3= ROOTS[[4]];
	R4= ROOTS[[3]];];
	

	
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,R4,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<R4||r0>R3,
		Message[KerrGeoPlunge::r0outofboundsGen,r0,R3,R4];
		Return[$Failed]
		];
		
		
	z0 = Cos[\[Theta]0];
	
	Z1=Sqrt[1/2 (1+\[ScriptCapitalL]^2/(a^2 (1-\[ScriptCapitalE]^2))+\[ScriptCapitalQ]/(a^2 (1-\[ScriptCapitalE]^2))-Sqrt[-((4 \[ScriptCapitalQ])/(a^2 (1-\[ScriptCapitalE]^2)))+(-1-\[ScriptCapitalL]^2/(a^2 (1-\[ScriptCapitalE]^2))-\[ScriptCapitalQ]/(a^2 (1-\[ScriptCapitalE]^2)))^2])];
	Z2 =Sqrt[ (a^2 (1-\[ScriptCapitalE]^2))/2 (1+\[ScriptCapitalL]^2/(a^2 (1-\[ScriptCapitalE]^2))+\[ScriptCapitalQ]/(a^2 (1-\[ScriptCapitalE]^2))+Sqrt[-((4 \[ScriptCapitalQ])/(a^2 (1-\[ScriptCapitalE]^2)))+(-1-\[ScriptCapitalL]^2/(a^2 (1-\[ScriptCapitalE]^2))-\[ScriptCapitalQ]/(a^2 (1-\[ScriptCapitalE]^2)))^2])];
	kz = a^2*J(Z1^2/Z2^2);
	kr = ((R3-R4)(R1-R2))/((R3-R1)(R4-R2));
	
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[KerrGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
		


	hR = (R3-R4)/(R3-R1);
	hP = hR (R1-RP)/(R4-RP);
	hM =  hR (R1-RM)/(R4-RM);
	
	Minoz[z_] :=InverseJacobiSN[z/Z1,kz^2]/Z2;
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	MinoR[r_]:= 2 InverseJacobiSN[Sqrt[(-R1+R3) (r-R4)]/Sqrt[(r-R1) (R3-R4)] ,kr]/Sqrt[J (R1-R3) (R2-R4)];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r]- MinoR[r0]],Listable];
	
	If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	If[\[ScriptCapitalQ]==0, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];
	
	
	\[CapitalUpsilon]r= \[Pi]/(2*EllipticK[kr]) Sqrt[J(R3-R1)(R4-R2)];
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[kz]);
	
	\[CapitalUpsilon]\[Phi]r = a/(RP-RM) ((2\[ScriptCapitalE]*RP-a*\[ScriptCapitalL])/(R1-RP) (1-(R4-R1)/(R4-RP) (EllipticPi[hP,kr]/EllipticK[kr]) )- (2\[ScriptCapitalE]*RM-a*\[ScriptCapitalL])/(R1-RM) (1-(R4-R1)/(R4-RM) (EllipticPi[hM,kr]/EllipticK[kr]))  );
	\[CapitalUpsilon]\[Phi]z  = \[ScriptCapitalL]/EllipticK[kz] EllipticPi[Z1^2,kz];
	\[CapitalUpsilon]\[Phi]= \[CapitalUpsilon]\[Phi]r+\[CapitalUpsilon]\[Phi]z;
	
	\[Phi]Tr[\[Xi]_]:= (-2*a*\[ScriptCapitalE]*(R4-R1))/((RP-RM)Sqrt[J*(R3-R1)(R4-R2)])*((2*RP-a*\[ScriptCapitalL]/\[ScriptCapitalE])/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(2*RM-a*\[ScriptCapitalL]/\[ScriptCapitalE])/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr]);
	\[Phi]Tz[\[Xi]_]:=  -\[ScriptCapitalL]/Z2 EllipticPi[Z1^2,\[Xi],kz];

	\[Phi]r[\[Lambda]_] := \[Phi]Tr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - \[Phi]Tr[\[Pi]]/(2*\[Pi]) (\[CapitalUpsilon]r*\[Lambda]);
	\[Phi]z[\[Lambda]_]:= -\[Phi]Tz[JacobiAmplitude[EllipticK[kz]*(2*(\[CapitalUpsilon]z*\[Lambda]))/\[Pi],kz]] + \[Phi]Tz[\[Pi]]/\[Pi] (\[CapitalUpsilon]z*\[Lambda]);
	

	\[CapitalUpsilon]tr = (4+a^2)\[ScriptCapitalE]+\[ScriptCapitalE](1/2 ((4+R3+R4+R1)R1-R3*R4+(R3-R1)(R4-R2) EllipticE[kr]/EllipticK[kr]+(4+R3+R4+R1+R2)(R4-R1) EllipticPi[hR,kr]/EllipticK[kr])+2/(RP-RM) (((4-a*\[ScriptCapitalL]/\[ScriptCapitalE])RP-2*a^2)/(R1-RP) (1-(R4-R1)/(R4-RP) EllipticPi[hP,kr]/EllipticK[kr])-((4-a*\[ScriptCapitalL]/\[ScriptCapitalE])RM-2*a^2)/(R1-RM) (1-(R4-R1)/(R4-RM) EllipticPi[hM,kr]/EllipticK[kr])));
	\[CapitalUpsilon]tz = -a^2*\[ScriptCapitalE]+(\[ScriptCapitalE]*\[ScriptCapitalQ])/(J(Z1^2)) (1-EllipticE[kz]/EllipticK[kz]);
	If[\[ScriptCapitalQ]==0, \[CapitalUpsilon]tz = -a^2*\[ScriptCapitalE]];
	\[CapitalUpsilon]t= \[CapitalUpsilon]tr +\[CapitalUpsilon]tz;
	
		
	tTr[\[Xi]_]:= (\[ScriptCapitalE](R4-R1))/Sqrt[J(R3-R1)(R4-R2)] ((4+R3+R4+R1+R2)EllipticPi[hR,\[Xi],kr]-4/(RP-RM) ((RP(4-a*\[ScriptCapitalL]/\[ScriptCapitalE])-2*a^2)/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(RM(4-a*\[ScriptCapitalL]/\[ScriptCapitalE])-2*a^2)/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr])+((R3-R1)(R4-R2))/(R4-R1) (EllipticE[\[Xi],kr]-(hR*Sin[\[Xi]]*Cos[\[Xi]]Sqrt[1-kr*(Sin[\[Xi]])^2])/(1-hR*(Sin[\[Xi]])^2)));
	tTz[\[Xi]_]:= -\[ScriptCapitalE]/J Z2*EllipticE[\[Xi],kz];
	
	r[\[Lambda]_]:=(R1(R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-R4(R3-R1))/((R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-(R3-R1));
	z[\[Lambda]_]:=  Z1*JacobiSN[EllipticK[kz]*2 (\[CapitalUpsilon]z*\[Lambda])/\[Pi],kz] ;
	
	
	tr[\[Lambda]_] := tTr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - tTr[\[Pi]]/(2*\[Pi]) (\[CapitalUpsilon]r*\[Lambda]);
	tz[\[Lambda]_]:= tTz[JacobiAmplitude[EllipticK[kz]*(2*(\[CapitalUpsilon]z*\[Lambda]))/\[Pi],kz]] - tTz[\[Pi]]/\[Pi] ( \[CapitalUpsilon]z*\[Lambda]);

	t=Function[{Global`\[Lambda]}, Evaluate[\[CapitalUpsilon]t*Global`\[Lambda] + tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]]-tr[MinoR[r0]]-tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[CapitalUpsilon]\[Phi]*Global`\[Lambda]+ \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];
	
	MCPP = -Abs[MinoR[RP]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMP = -Abs[MinoR[RM]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMM = Abs[MinoR[RM]] - MinoR[r0];
	MCPM = Abs[MinoR[RP]] - MinoR[r0];
	
	assoc = Association[
		"a" -> a,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> \[CapitalUpsilon]r,
		"PolarFrequency"-> \[CapitalUpsilon]z,
		"Frequencies"-> <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"->\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]z |>,
		"RadialPeriod"-> 2\[Pi]/\[CapitalUpsilon]r,
		"PolarPeriod"-> 2\[Pi]/\[CapitalUpsilon]z,
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> 2\[Pi]/\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]z|>,
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialRoots"-> {R4,R3,R2,R1},
		"RadialMinoTime"-> MinoRFunc,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R4]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R3]] - MinoR[r0](* Time outgoing to outer root*)|>,
		"HorizonCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(+\)]\)"-> MCPP(* Time ingoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(+\)]\)"-> MCMP(* Time ingoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(-\)]\)"-> MCMM(* Time outgoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(-\)]\)"-> MCPM(* Time outgoing branch crosses outer horizon *)
			|>,
		"PolarRoots"-> {Z1,Z2},
		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]



(* ::Subsection:: *)
(*Complex Roots Generic Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoPlunge::r0outofboundsGen = "Intial radius `1` is not between the outer radial root `2`, and inner radial root `3`.";
KerrGeoPlunge::\[Theta]0outofbounds = "Intial polar angle `1` must be between the ranges (`2`,`3`) or (-`2`,-`3`)";


Options[KerrGeoComplexPlungeMino]= {"Roots"->Automatic};
KerrGeoComplexPlungeMino[a_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalQ]_ , initCoords_, OptionsPattern[Options[KerrGeoComplexPlungeMino]]] := Module[{consts,MCPP, MCPM, MCMP , MCMM , assoc,MinoR,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],t0,r0,\[Phi]0,M,D1M,D1P,D2M,D2P,Minoz,MinozFunc,MinoRFunc,\[Theta]0,z0,e,b,c,d,A,B,chi,kr,p2,f, t, r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots,kz,Z1, Z2, AMR,CNR,SNR,DNR,AMZ},
	
	M=1;
	J = 1-\[ScriptCapitalE]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	ROOTS = r/.Solve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2M*r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r];
	RealRoots = Select[ROOTS,Im[#]==0&];
	ComplexRoots = Select[ROOTS,Im[#]!=0&];
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalQ]"->\[ScriptCapitalQ]|>;
	R1= RealRoots[[1]];
	R2= RealRoots[[2]];
	R3=ComplexRoots[[1]];
	R4= ComplexRoots[[2]];
	
	If[a==0, {Z1,Z2}={0,\[ScriptCapitalL]}, {Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]}];
	
	e = R2;
	b = R1;
	c = Re[R3];
	d = Abs[Im[R4]];
	A = Sqrt[(e-c)^2+d^2];
	B  = Sqrt[(b-c)^2+d^2];
	chi = (A*b+e*B)/(A+B);
	kr = Sqrt[((e-b)^2-(A-B)^2)/(4*A*B)];
	p2 = b*A^2+e*B^2-(e+b)*A*B;
	f = (4 A B)/(A-B)^2;
	kz = a*Sqrt[J](Z1/Z2);
	
	\[CapitalUpsilon]r = 2 \[Pi] Sqrt[A B (1-\[ScriptCapitalE]^2) ]/(4 EllipticK[kr^2]);
	\[CapitalUpsilon]\[Theta] = \[Pi]/2 EllipticK[kz^2]/Z2;
	
	
	D1M=Sqrt[ -f];
	D2M=Sqrt[(4 A B (RM-b) (e-RM))]/(A( RM-b)+B(e- RM));

	D1P=Sqrt[ -f];
	D2P=Sqrt[(4 A B (RP-b) (e-RP))]/(A(RP- b)+B( e- RP));
	
	
	AMR[\[Lambda]_]:=JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2];
	SNR[\[Lambda]_]:=JacobiSN[Sqrt[A B J] \[Lambda],kr^2];
	CNR[\[Lambda]_]:=JacobiCN[Sqrt[A B J] \[Lambda],kr^2];
	DNR[\[Lambda]_]:=JacobiDN[Sqrt[A B J] \[Lambda],kr^2];

	AMZ[\[Lambda]_]:=JacobiAmplitude[Z2 \[Lambda],kz^2];
	
	If[initCoords===Automatic,
					{t0,r0,\[Theta]0,\[Phi]0}={0,R1,\[Pi]/2,0},
					{t0,r0,\[Theta]0,\[Phi]0}=initCoords
					];
	If[r0<R1||r0>R2,
		Message[KerrGeoPlunge::r0outofboundsGen,r0,R2,R1];
		Return[$Failed]
		];
	If[Abs[\[Theta]0]<Abs[ArcCos[Z1]]||Abs[\[Theta]0]>\[Pi] - Abs[ArcCos[Z1]],
		Message[KerrGeoPlunge::\[Theta]0outofbounds,\[Theta]0,Abs[ArcCos[Z1]],\[Pi]-Abs[ArcCos[Z1]] ];
		Return[$Failed]
		];
		
			
	z0 = Cos[\[Theta]0];
	
	Minoz[z_] :=InverseJacobiSN[z/Z1,kz^2]/Z2;
	If[Abs[Arg]==\[Pi]/2, Minoz[z_] :=0];
	If[a==0, Minoz[z_] :=0];
	If[\[ScriptCapitalQ]==0, Minoz[z_] :=0];
	MinozFunc=Function[{Global`z}, Evaluate[Minoz[Global`z]-Minoz[z0]  ],Listable];


	MinoR[x_]:= -1/Sqrt[J A B] EllipticF[(\[Pi]/2 -ArcSin[(B(e-x)-A(x-b))/(   B(e-x)+A(x-b)  )  ]),kr^2];
	MinoRFunc=Function[{Global`r}, Evaluate[MinoR[Global`r] -MinoR[r0] ],Listable];
	
	
	r[\[Lambda]_] := ((A-B) (A b-B e) SNR[\[Lambda]]^2+2 (A B (b+e)+A B (b-e) CNR[\[Lambda]]))/(4 A B+(A-B)^2 SNR[\[Lambda]]^2);

	z[\[Lambda]_]:= Z1*JacobiSN[Z2 \[Lambda],kz^2];
	If[kz==0, z[\[Lambda]_]:= Z1*Sin[Z2*\[Lambda]]];

(*Integrals*)
	
	RINT\[Lambda][\[Lambda]_] := ((A b-B e)/(A-B) \[Lambda]-1/ Sqrt[ J] ArcTan[(e-b)/(2 Sqrt[A B])  SNR[\[Lambda]]/Sqrt[1-kr^2 (SNR[\[Lambda]])^2]]+((A+B) (e-b))/(2 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]);
	R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)+ (A^2+2 b^2-B^2-2 e^2)/(4 (e-b) Sqrt[ J])ArcTan[(f-(1+2 f kr^2) SNR[\[Lambda]]^2),2 SNR[\[Lambda]] DNR[\[Lambda]]Sqrt[f (1+f kr^2)]](*ArcTan[(2 SNR[\[Lambda]] DNR[\[Lambda]]Sqrt[f (1+f kr^2)])/(f-(1+2 f kr^2) SNR[\[Lambda]]^2)]*);
	If[a!=0,RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM)))/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) EllipticPi[1/D2M^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]-  Sqrt[(e-b)]/(Sqrt[ J] Sqrt[ (RM-b) (e-RM)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))]) 1/4 (Log[((D2M  Sqrt[1-D2M^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)/((D2M  Sqrt[1-D2M^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)])];
	If[a==0,RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))];
	RPINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP)))/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) EllipticPi[1/D2P^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]- Sqrt[(e-b)]/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))]) 1/4 (Log[((D2P  Sqrt[1-D2P^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)/((D2P  Sqrt[1-D2P^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)]);


	tr[\[Lambda]_]:=  ((2 a^2 +RM^2+RM RP+RP^2)\[ScriptCapitalE])\[Lambda]+(R2INT\[Lambda][\[Lambda]]+RINT\[Lambda][\[Lambda]](RM+RP)) \[ScriptCapitalE] + ((RM^2+a^2)(\[ScriptCapitalE](RM^2+a^2)-a*\[ScriptCapitalL]))/(RM-RP) RMINT\[Lambda][\[Lambda]]+ ((RP^2+a^2)(\[ScriptCapitalE](RP^2+a^2)-a*\[ScriptCapitalL]))/(RP-RM) RPINT\[Lambda][\[Lambda]];
	\[Phi]r[\[Lambda]_]:= a(((\[ScriptCapitalE](RM^2+a^2)-a*\[ScriptCapitalL])/(RM-RP))RMINT\[Lambda][\[Lambda]]+ (\[ScriptCapitalE](RP^2+a^2)-a*\[ScriptCapitalL])/(RP-RM) RPINT\[Lambda][\[Lambda]]);
	tz[\[Lambda]_]:= \[ScriptCapitalE]/ J ((Z2-a^2 J/Z2) Z2 \[Lambda] -Z2 EllipticE[AMZ[\[Lambda]],kz^2]);
	\[Phi]z[\[Lambda]_]:= (\[ScriptCapitalL] EllipticPi[Z1^2,AMZ[\[Lambda]],(a^2 J Z1^2)/Z2^2])/Z2;



(*Note: If including Precission Error May occur in the case of evaluating R2INT\[Lambda][0] due to a bug with mathematica in
 evaluating terms of the form EllipticE[0, num`n], num is any number and n is less thna machine precision*)
	t=Function[{Global`\[Lambda]}, Evaluate[  tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ Minoz[z0]]-tr[MinoR[r0]]-tz[Minoz[z0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]+ MinoR[r0]]], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + Minoz[z0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[ \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ Minoz[z0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[Minoz[z0]] + \[Phi]0], Listable];

	MCPP = -Abs[MinoR[RP]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMP = -Abs[MinoR[RM]] -MinoR[r0]+(2\[Pi])/\[CapitalUpsilon]r;
	MCMM = Abs[MinoR[RM]] - MinoR[r0];
	MCPM = Abs[MinoR[RP]] - MinoR[r0];
	assoc = Association[
		"a" -> a,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> \[CapitalUpsilon]r,
		"PolarFrequency"-> \[CapitalUpsilon]\[Theta],
		"Frequencies"-> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"-> \[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]\[Theta]|>,
		"RadialPeriod"-> 2\[Pi]/\[CapitalUpsilon]r,
		"PolarPeriod"-> 2\[Pi]/\[CapitalUpsilon]\[Theta],
		"Periods"-> <|"\!\(\*SubscriptBox[\(\[Lambda]\), \(r\)]\)"-> 2\[Pi]/\[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[Lambda]\), \(\[Theta]\)]\)" -> 2\[Pi]/\[CapitalUpsilon]\[Theta]|>,
		"RadialRoots"-> {R1,R2,R3,R4},
		"EllipticBasis" ->  <|"\!\(\*SubscriptBox[\(k\), \(r\)]\)"-> kr,"\!\(\*SubscriptBox[\(k\), \(z\)]\)"-> kz|>,
		"RadialMinoTime"-> MinoRFunc,
		"RadialRootCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Inner\)], \(\)]\)"-> Abs[MinoR[R1]] - MinoR[r0](* Time outgoing to inner root*),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(Outer\)], \(\)]\)"-> Abs[MinoR[R2]] - MinoR[r0](* Time outgoing to outer root*)|>,
		"HorizonCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(-\)]\)"-> MCMM(* Time outgoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(-\)]\)"-> MCPM(* Time outgoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(+\)]\)"-> MCPP(* Time ingoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(+\)]\)"-> MCMP(* Time ingoing branch crosses inner horizon *)
			|>,
		"PolarRoots"-> {Z1,Z2},
		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]


(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


KerrGeoPlunge::fourcomplexroots = "This function can only currently solve plunging geodesics in the paramater space where we have four finite real roots or two finite real roots and two finite complex roots.";
KerrGeoPlunge::highenergy = "This code does not currently support plunging solutions with \[ScriptCapitalE]>1, \[ScriptCapitalE]=`1` has been selected ";
KerrGeoPlunge::negativecarter = "This code does not currently support plunging solutions with Q<0, Q=`1` has been selected ";
KerrGeoPlunge::rIoutbounds = "Given ISSO radius `1` is not in the allowed range between `2` and `3` for spin `4`.";
KerrGeoPlunge::Inclinationoutofbounds = "Given inclination angle `1` is not between -\!\(\*FractionBox[\(\[Pi]\), \(2\)]\) and \!\(\*FractionBox[\(\[Pi]\), \(2\)]\)."
KerrGeoPlunge::aoutofbounds1 = "The plunges packages does not currently support the extremal a=1 case"
KerrGeoPlunge::aoutofbounds0 = "The plunges packages does not currently support the zero spin a=0 case"


KerrGeoPlunge[a_:0.9,{\[ScriptCapitalE]_:0.8,\[ScriptCapitalL]_:0.3,\[ScriptCapitalQ]_:3}, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=Module[
	{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
	
	
	If[a==1,
		Message[KerrGeoPlunge::aoutofbounds1];
		Return[$Failed]
		];
	
	If[\[ScriptCapitalE]>=1,
		Message[KerrGeoPlunge::highenergy,\[ScriptCapitalE]];
		Return[$Failed]
		];
		
	If[\[ScriptCapitalQ]<0,
		Message[KerrGeoPlunge::negativecarter,\[ScriptCapitalQ]];
		Return[$Failed]
		];

	ROOTS = r/.NSolve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r];

	RealRoots = Select[ROOTS,PossibleZeroQ@Im[#]&];
	ComplexRoots = Select[ROOTS,Not@PossibleZeroQ@Im[#]&];

	If[Length[ComplexRoots]== 2, Return[KerrGeoComplexPlungeMino[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], initCoords,"Roots"->ROOTS]]];
	If[Length[ComplexRoots]== 0, Return[KerrGeoRealPlungeMino[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], initCoords,"Roots"->ROOTS]]];
	];


KerrGeoPlunge[a_:0.9,"ISSORadialParam",Arg1_:0.8, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=
	Module[{RI,DN,RISSOMIN,RISSOMAX},
	If[a==1,
		Message[KerrGeoPlunge::aoutofbounds1,a];
		Return[$Failed]
		];
		RI = Arg1;
		DN = (27-45 a^2+17 a^4+a^6+8 a^3 (1-a^2))^(1/3);
		RISSOMIN = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]-1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
		RISSOMAX = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]+1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
		If[Not[RISSOMIN<=RI<=RISSOMAX],
			Message[KerrGeoPlunge::rIoutbounds,RI,RISSOMIN,RISSOMAX,a];
			Return[$Failed]
			];
		KerrGeoISSOPlunge[a,"ISSORadialParam" ,RI, initCoords]
	]


KerrGeoPlunge[a_:0.9,"ISSOIncParam",Arg1_:0.8, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=
	Module[{},
	If[a==1,
		Message[KerrGeoPlunge::aoutofbounds1,a];
		Return[$Failed]
		];
		(* Other inclinations should be fine as well, right?*)
		If[Not[-\[Pi]/2<=Arg1<=\[Pi]/2],
			Message[KerrGeoPlunge::Inclinationoutofbounds, Arg1];
			Return[$Failed]
			];
			
	KerrGeoISSOPlunge[a, "ISSOIncParam" ,Arg1 ,initCoords]
	]


KerrGeoPlunge[a_:0.9,"USOPlunge",{RS_:3.8,\[ScriptCapitalQ]_:2.0}, InnerOutter_, initCoords:{_,_,_,_}:Automatic,OptionsPattern[]]:=
	Module[{},

	KerrGeoUSOPlunge[a, RS, \[ScriptCapitalQ], InnerOutter, initCoords]
	]


KerrGeoPlungeFunction /:
 MakeBoxes[kgof:KerrGeoPlungeFunction[a_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalQ]_,assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalE]: ", \[ScriptCapitalE]}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalL]: ", \[ScriptCapitalL]}], "  ",
                  BoxForm`SummaryItem[{"\[ScriptCapitalQ]: ", \[ScriptCapitalQ]}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}],
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]};
  BoxForm`ArrangeSummaryBox[
    KerrGeoPlungeFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


KerrGeoPlungeFunction[a_, \[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalQ]_,assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrGeoPlungeFunction[a_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalQ]_,assoc_][y_?StringQ] := assoc[y]
Keys[g_KerrGeoPlungeFunction]^:=Keys[g[[5]]]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
