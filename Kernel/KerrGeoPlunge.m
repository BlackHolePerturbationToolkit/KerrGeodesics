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


KerrGeoPlunge::usage = "Takes either KerrGeoPlunge[a, Generic, En ,\[ScriptCapitalL],\[ScriptCapitalQ]] or KerrGeoPlunge[a, ISSO , RI] for generic Plunges or ISSO plunges respectively and returns a KerrGeoPlungeFunction[..] which stores the orbital trajectory and parameters. Here the ISSO plunges are paramaterised in terms of the choice of the radius of the ISSO for a given a there are range of allowed RI's which correspond to differeing inclinations in the prograde and retrograde directions, if an RI outside this range is given the used is provided with the range of allowed values to re-run the funciton.";
KerrGeoPlungeFunction::usage = "KerrGeoPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*ISSO Plunges*)


KerrGeoPlunge::r0outofbounds = "Intial radius `1` is not between ISSO radius `2`, and inner radial root `3`.";


KerrGeoISSOPlunge[a_, PlungeType_  ,Arg_, initCoords_] := Module[
	{consts, \[Theta]0, z, r0,Mino, z0, \[Phi]0, initConditions, t0, assoc,M=1,\[Chi]INC ,t,RI, r, \[Theta], \[Phi], \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], R4, RM, RP, kz, Z1, Z2, J,velocity,MinoR,tr,tz,\[Phi]r,\[Phi]z},

	RM = 1-Sqrt[1-a^2];
	RP = 1+Sqrt[1-a^2];
	
	If[PlungeType == "ISSORadialParam",
		RI = Arg;
		\[ScriptCapitalQ] = - M (RI)^(5/2) ((Sqrt[(RI-RP)(RI-RM)]-2Sqrt[RI])^2-4a^2)/(4a^2 (RI^(3/2)-Sqrt[RI]-Sqrt[(RI-RP)(RI-RM)]));
		\[ScriptCapitalE] = Sqrt[a^2 \[ScriptCapitalQ]-2 RI^3+3 RI^4]/(Sqrt[3] RI^2);
		\[ScriptCapitalL] = Sqrt[3 a^2 \[ScriptCapitalQ]-a^2 RI^2-\[ScriptCapitalQ] RI^2+3 RI^4+a^2 RI^2 \[ScriptCapitalE]^2-3 RI^4 \[ScriptCapitalE]^2]/RI;
		];
	
	If[PlungeType == "ISSOIncParam",
		\[Chi]INC = Cos[Arg];
		RI = KerrGeoISSO[a,\[Chi]INC];
		{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]} = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[1;;3]]
		];
	
	J = 1-\[ScriptCapitalE]^2;
	{Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]};
	consts = <|"\[ScriptCapitalE]"-> \[ScriptCapitalE],"\[ScriptCapitalL]"-> \[ScriptCapitalL], "\[ScriptCapitalQ]"->\[ScriptCapitalQ]|>;
	R4 = (a^2 \[ScriptCapitalQ])/(J*RI^3);
	kz = a*Sqrt[J](Z1/Z2);
	
	
	If[initCoords===Automatic,
					{t0,r0,\[Phi]0}={0,R4,0},
					{t0,r0,\[Phi]0}=initCoords
					];
	
	If[r0>RI||r0<R4,
		Message[KerrGeoPlunge::r0outofbounds,r0,RI,R4];
		Return[$Failed]
		];
		
	MinoR[r_] :=(2 Sqrt[r-R4])/Sqrt[J (RI-r) (R4-RI)^2];
	Mino=Function[{Global`r}, Evaluate[MinoR[Global`r] -MinoR[r0]  ],Listable];
	
	r[\[Lambda]_] := ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4));

	z[\[Lambda]_]:= Z1*JacobiSN[Z2*\[Lambda] ,kz^2];
	tr[\[Lambda]_]:= ((a^2+RI^2) (-a \[ScriptCapitalL]+(a^2+RI^2) \[ScriptCapitalE])(\[Lambda]) )/((RI-RM) (RI-RP))+(2 ((R4-RI)^2) \[ScriptCapitalE] (\[Lambda]) )/(4+J (R4-RI)^2 (\[Lambda]) ^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[ScriptCapitalE] ArcTan[((\[Lambda])  (RI-R4)Sqrt[J ]  )/2 ]+( (a^2+RM^2) (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP)Sqrt[J])+( (a^2+RP^2) (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE]) Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP));
	\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a (((\[Lambda]) (RI-R4)Sqrt[J ] (-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RI^2 \[ScriptCapitalE]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RM^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(2Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP))+((-a \[ScriptCapitalL]+a^2 \[ScriptCapitalE]+RP^2 \[ScriptCapitalE])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(2Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP)))    ;
	tz[\[Lambda]_]:= 1/J \[ScriptCapitalE]  (-Z2 EllipticE[JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2]+(Z2-a^2 J/Z2) Z2 \[Lambda]);
	\[Phi]z[\[Lambda]_]:= \[ScriptCapitalL]/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2];


	t=Function[{Global`\[Lambda]}, Evaluate[ a*\[ScriptCapitalL] Global`\[Lambda] + tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ MinoR[r0]]-tr[MinoR[r0]]-tz[MinoR[r0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + MinoR[r0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[-a*\[ScriptCapitalE] Global`\[Lambda] + \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ MinoR[r0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[MinoR[r0]] + \[Phi]0], Listable];

	assoc = Association[
		"a" -> a,
		"rI" -> RI,
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency"-> 0,
		"PolarFrequency"-> Missing[],
		"AzimuthalFrequency"-> 0,
		"Frequencies"-> <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"-> 0, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Missing[], "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"->0, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"-> Missing[]|>,
		"RadialRoots"-> {RI,RI,RI,R4},
		"Mino"-> Mino, (* This needs a better name. *)
		"HorizonCrossingTimeMino"-> <|
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(+\)]\)"-> Missing[](* Time ingoing branch crosses outer horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(+\)]\)"-> Missing[](* Time ingoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(-\)], \(-\)]\)"-> Mino[RM] (* Time outgoing branch crosses inner horizon *),
			"\!\(\*SubsuperscriptBox[\(\[CapitalLambda]\), SubscriptBox[\(r\), \(+\)], \(-\)]\)"-> Mino[RP] (* Time outgoing branch crosses outer horizon *)
			|>,
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {0,kz^2}, (*Is this output needed? *)
		"Trajectory" -> {t,r,\[Theta],\[Phi]}
	];
	
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]



(* ::Subsection:: *)
(*Generic Plunge (Mino)*)


(* ::Subsubsection:: *)
(*4 Real roots*)


KerrGeoRealPlunge::r0outofbounds = "Intial radius `1` is larger than outer turning radius `2`.";


Options[KerrGeoRealPlungeMino] = {"Roots"-> Automatic};
KerrGeoRealPlungeMino[a_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalQ]_ , initCoords_, OptionsPattern[Options[KerrGeoRealPlungeMino]]] := Module[
	{consts,assoc,M,Mino,t0,r0,\[Phi]0,kr,hR,hP,hM,\[CapitalUpsilon]t,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]tr,\[CapitalUpsilon]tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots, kz, Z1, Z2, MinoR},
	
	M=1;
	J = 1-\[ScriptCapitalE]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	If[
		OptionValue["Roots"]===Automatic,
		ROOTS = Sort[r/.Solve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2M*r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r]],
		ROOTS = OptionValue["Roots"]
		];
	
	If[Length[Select[ROOTS,#>=RP&]]>=3,
	R1= ROOTS[[3]];
	R2= ROOTS[[4]];
	R3= Min[Select[ROOTS,#>=RP&]];
	R4= Max[Select[ROOTS,#<=RM&]];];
	
	If[Length[Select[ROOTS,#>=RP&]]<3,
	R1= ROOTS[[2]];
	R2= ROOTS[[1]];
	R3= Min[Select[ROOTS,#>=RP&]];
	R4= Max[Select[ROOTS,#<=RM&]];];
		
	If[initCoords===Automatic,
					{t0,r0,\[Phi]0}={0,R4,0},
					{t0,r0,\[Phi]0}=initCoords
					];
	
	If[r0>R3,
		Message[KerrGeoRealPlunge::r0outofbounds,r0,R3];
		Return[$Failed]
		];
	
	{Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]};
	kz= a*Sqrt[J](Z1/Z2);
	
	kr = Sqrt[((R3-R4)(R1-R2))/((R3-R1)(R4-R2))];
	hR = (R3-R4)/(R3-R1);
	hP = hR (R1-RP)/(R4-RP);
	hM =  hR (R1-RM)/(R4-RM);
	

	MinoR[r_]:= 2 InverseJacobiSN[Sqrt[(-R1+R3) (r-R4)]/Sqrt[(r-R1) (R3-R4)] ,kr]/Sqrt[J (R1-R3) (R2-R4)];
	
	Mino=Function[{Global`r}, Evaluate[MinoR[Global`r]- MinoR[r0]],Listable];
	
	(*r=Function[{Global`\[Lambda]}, Evaluate[r[Global`\[Lambda]+ MinoR[r0]] ], Listable];*)
	
	\[CapitalUpsilon]r= \[Pi]/(2*EllipticK[kr]) Sqrt[J(R3-R1)(R4-R2)];
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[kz^2]);
	
	\[CapitalUpsilon]\[Phi]r = a/(RP-RM) ((2\[ScriptCapitalE]*RP-a*\[ScriptCapitalL])/(R1-RP) (1-(R4-R1)/(R4-RP) (EllipticPi[hP,kr]/EllipticK[kr]) )- (2\[ScriptCapitalE]*RM-a*\[ScriptCapitalL])/(R1-RM) (1-(R4-R1)/(R4-RM) (EllipticPi[hM,kr]/EllipticK[kr]))  );
	\[CapitalUpsilon]\[Phi]z  = \[ScriptCapitalL]/EllipticK[kz^2] EllipticPi[Z1^2,kz^2];
	\[CapitalUpsilon]\[Phi]= \[CapitalUpsilon]\[Phi]r+\[CapitalUpsilon]\[Phi]z ;
	
	\[Phi]Tr[\[Xi]_]:= (-2*a*\[ScriptCapitalE]*(R4-R1))/((RP-RM)Sqrt[J*(R3-R1)(R4-R2)])*((2*RP-a*\[ScriptCapitalL]/\[ScriptCapitalE])/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(2*RM-a*\[ScriptCapitalL]/\[ScriptCapitalE])/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr]);
	\[Phi]Tz[\[Xi]_]:= -\[ScriptCapitalL]/Z2 EllipticPi[Z1^2,\[Xi],kz^2];

	\[Phi]r[\[Lambda]_] := \[Phi]Tr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - \[Phi]Tr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	\[Phi]z[\[Lambda]_]:= -\[Phi]Tz[JacobiAmplitude[EllipticK[kz^2]*(2*(\[CapitalUpsilon]z*\[Lambda]))/\[Pi],kz^2]] + \[Phi]Tz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];
	

	\[CapitalUpsilon]tr = (4+a^2)\[ScriptCapitalE]+\[ScriptCapitalE](1/2 ((4+R3+R4+R1)R1-R3*R4+(R3-R1)(R4-R2) EllipticE[kr]/EllipticK[kr]+(4+R3+R4+R1+R2)(R4-R1) EllipticPi[hR,kr]/EllipticK[kr])+2/(RP-RM) (((4-a*\[ScriptCapitalL]/\[ScriptCapitalE])RP-2*a^2)/(R1-RP) (1-(R4-R1)/(R4-RP) EllipticPi[hP,kr]/EllipticK[kr])-((4-a*\[ScriptCapitalL]/\[ScriptCapitalE])RM-2*a^2)/(R1-RM) (1-(R4-R1)/(R4-RM) EllipticPi[hM,kr]/EllipticK[kr])));
	\[CapitalUpsilon]tz= -a^2*\[ScriptCapitalE]+(\[ScriptCapitalE]*\[ScriptCapitalQ])/(J(Z1^2)) (1-EllipticE[kz^2]/EllipticK[kz^2]);
	\[CapitalUpsilon]t= \[CapitalUpsilon]tr +\[CapitalUpsilon]tz;
	
	tTr[\[Xi]_]:= (\[ScriptCapitalE](R4-R1))/Sqrt[J(R3-R1)(R4-R2)] ((4+R3+R4+R1+R2)EllipticPi[hR,\[Xi],kr]-4/(RP-RM) ((RP(4-a*\[ScriptCapitalL]/\[ScriptCapitalE])-2*a^2)/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(RM(4-a*\[ScriptCapitalL]/\[ScriptCapitalE])-2*a^2)/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr])+((R3-R1)(R4-R2))/(R4-R1) (EllipticE[\[Xi],kr]-(hR*Sin[\[Xi]]*Cos[\[Xi]]Sqrt[1-kr*(Sin[\[Xi]])^2])/(1-hR*(Sin[\[Xi]])^2)));
	tTz[\[Xi]_]:= -\[ScriptCapitalE]/J Z2*EllipticE[\[Xi],kz^2];
	
	
	r[\[Lambda]_]:=(R1(R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-R4(R3-R1))/((R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-(R3-R1));
	z[\[Lambda]_]:=  Z1*JacobiSN[EllipticK[kz^2]*2 (\[CapitalUpsilon]z*\[Lambda])/\[Pi],kz^2];
	
	tr[\[Lambda]_] := tTr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - tTr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	tz[\[Lambda]_]:= tTz[JacobiAmplitude[EllipticK[kz^2]*(2*\[CapitalUpsilon]z*\[Lambda])/\[Pi],kz^2]] - tTz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];

	t=Function[{Global`\[Lambda]}, Evaluate[tr[Global`\[Lambda]+ MinoR[r0]]+tz[Global`\[Lambda]+ MinoR[r0]] + \[CapitalUpsilon]T*Global`\[Lambda]-tr[Global`\[Lambda]]-tz[Global`\[Lambda]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[r[Global`\[Lambda]+ MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda]+ MinoR[r0]]]], Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[Phi]r[Global`\[Lambda]+ MinoR[r0]]+\[Phi]z[Global`\[Lambda]+ MinoR[r0]] + \[CapitalUpsilon]\[Phi]*Global`\[Lambda] -\[Phi]r[ MinoR[r0]]-\[Phi]z[MinoR[r0]]+ \[Phi]0], Listable];

	assoc = Association[
		"a" -> a,
		"\[ScriptCapitalE]" -> \[ScriptCapitalE],
		"\[ScriptCapitalL]" -> \[ScriptCapitalL],
		"\[ScriptCapitalQ]" -> \[ScriptCapitalQ],
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"Mino"-> Mino,
		"HorizonCrossingTimeMino"-> {Mino[RM],Mino[RP]},
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {kr^2,kz^2},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"Frequencies" -> <|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" ->  \[CapitalUpsilon]r, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" ->  \[CapitalUpsilon]z, "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" ->  \[CapitalUpsilon]\[Phi] ,"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" ->  \[CapitalUpsilon]t |>
	];
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]



(* ::Subsection:: *)
(*Complex Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoComplexPlunge::r0outofbounds = "Intial radius `1` is larger than outer turning radius `2`.";


Options[KerrGeoComplexPlungeMino]= {"Roots"->Automatic};
KerrGeoComplexPlungeMino[a_, \[ScriptCapitalE]_, \[ScriptCapitalL]_, \[ScriptCapitalQ]_ , initCoords_, OptionsPattern[Options[KerrGeoComplexPlungeMino]]] := Module[{consts, assoc,MinoR,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],t0,r0,\[Phi]0,M,D1M,D1P,D2M,D2P,e,b,c,d,A,B,chi,kr,p2,f, t, r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots,kz,Z1, Z2, AMR,CNR,SNR,DNR,AMZ},
	
	M=1;
	J = 1-\[ScriptCapitalE]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	ROOTS = r/.Solve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2M*r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r];
	RealRoots = Select[ROOTS,Im[#]==0&];
	ComplexRoots = Select[ROOTS,Im[#]!=0&];

	R1= RealRoots[[1]];
	R2= RealRoots[[2]];
	R3=ComplexRoots[[1]];
	R4= ComplexRoots[[2]];
	

	If[initCoords===Automatic,
					{t0,r0,\[Phi]0}={0,R1,0},
					{t0,r0,\[Phi]0}=initCoords
					];
	
	If[r0>R2,
		Message[KerrGeoComplexPlunge::r0outofbounds,r0,R2];
		Return[$Failed]
		];
	
	{Z1,Z2}= {Sqrt[1/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)-Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J)+Sqrt[(1+(\[ScriptCapitalL]^2+\[ScriptCapitalQ])/(a^2 J))^2-(4 \[ScriptCapitalQ])/(a^2 J)])]};

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
	D2M=Sqrt[(4 A B (RM-b) (e-RM))/(A( b-RM)-B( e- RM))^2];

	D1P=Sqrt[ -f];
	D2P=Sqrt[(4 A B (RP-b) (e-RP))/(A( b-RP)-B( e- RP))^2];
	
	
	AMR[\[Lambda]_]:=JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2];
	SNR[\[Lambda]_]:=JacobiSN[Sqrt[A B J] \[Lambda],kr^2];
	CNR[\[Lambda]_]:=JacobiCN[Sqrt[A B J] \[Lambda],kr^2];
	DNR[\[Lambda]_]:=JacobiDN[Sqrt[A B J] \[Lambda],kr^2];

	AMZ[\[Lambda]_]:=JacobiAmplitude[Z2 \[Lambda],kz^2];
	
	MinoR[x_]:=1/Sqrt[J A B] EllipticF[(\[Pi]/2-ArcSin[(B(e-x)-A(x-b))/Sqrt[4 A B (e-x) (-b+x)+(B (e-x)-A (-b+x))^2]]),kr^2];
	Mino=Function[{Global`r}, Evaluate[MinoR[Global`r] -MinoR[r0] ],Listable];
	
	
	r[\[Lambda]_] := ((A-B) (A b-B e) SNR[\[Lambda]]^2+2 (A B (b+e)+A B (b-e) CNR[\[Lambda]]))/(4 A B+(A-B)^2 SNR[\[Lambda]]^2);

	z[\[Lambda]_]:= Z1*JacobiSN[Z2 \[Lambda],kz^2];

(*Integrals*)
	RINT\[Lambda][\[Lambda]_] := ((A b-B e)/(A-B) \[Lambda]-1/ Sqrt[ J] ArcTan[(e-b)/(2 Sqrt[A B])  SNR[\[Lambda]]/Sqrt[1-kr^2 (SNR[\[Lambda]])^2]]+((A+B) (e-b))/(2 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]);

	R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)+ (A^2+2 b^2-B^2-2 e^2)/(4 (e-b) Sqrt[ J])ArcTan[(f-(1+2 f kr^2) SNR[\[Lambda]]^2),2 SNR[\[Lambda]] DNR[\[Lambda]]Sqrt[f (1+f kr^2)]];
	RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM)))/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) EllipticPi[1/D2M^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]-  Sqrt[(e-b)]/(Sqrt[ J] Sqrt[ (RM-b) (e-RM)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))]) 1/4 (Log[((D2M  Sqrt[1-D2M^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)/((D2M  Sqrt[1-D2M^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)]);

	RPINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP)))/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) EllipticPi[1/D2P^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]- Sqrt[(e-b)]/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))]) 1/4 (Log[((D2P  Sqrt[1-D2P^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)/((D2P  Sqrt[1-D2P^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)]);
	tr[\[Lambda]_]:=  ((2 a^2 +RM^2+RM RP+RP^2)\[ScriptCapitalE])\[Lambda]+(R2INT\[Lambda][\[Lambda]]+RINT\[Lambda][\[Lambda]](RM+RP)) \[ScriptCapitalE] + ((RM^2+a^2)(\[ScriptCapitalE](RM^2+a^2)-a*\[ScriptCapitalL]))/(RM-RP) RMINT\[Lambda][\[Lambda]]+ ((RP^2+a^2)(\[ScriptCapitalE](RP^2+a^2)-a*\[ScriptCapitalL]))/(RP-RM) RPINT\[Lambda][\[Lambda]];
	\[Phi]r[\[Lambda]_]:= a(((\[ScriptCapitalE](RM^2+a^2)-a*\[ScriptCapitalL])/(RM-RP))RMINT\[Lambda][\[Lambda]]+ (\[ScriptCapitalE](RP^2+a^2)-a*\[ScriptCapitalL])/(RP-RM) RPINT\[Lambda][\[Lambda]]);
	tz[\[Lambda]_]:= \[ScriptCapitalE]/ J ((Z2-a^2 J/Z2) Z2 \[Lambda] -Z2 EllipticE[AMZ[\[Lambda]],kz^2]);
	\[Phi]z[\[Lambda]_]:= (\[ScriptCapitalL] EllipticPi[Z1^2,AMZ[\[Lambda]],(a^2 J Z1^2)/Z2^2])/Z2;

	t=Function[{Global`\[Lambda]}, Evaluate[  tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ MinoR[r0]]-tr[MinoR[r0]]-tz[MinoR[r0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]+ MinoR[r0]]], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + MinoR[r0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[ \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ MinoR[r0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[MinoR[r0]] + \[Phi]0], Listable];


	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[ScriptCapitalE], 
		"AngularMomentum" -> \[ScriptCapitalL], 
		"CarterConstant" -> \[ScriptCapitalQ], 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {Z1,Z2},
		"Mino"-> Mino,
		"HorizonCrossingTimeMino"-> {Mino[RM],Mino[RP]},
		"EllipticBasis"-> {kr^2,kz^2},
		"RadialFrequency"->\[CapitalUpsilon]r,
		"PolarFrequency"-> \[CapitalUpsilon]\[Theta],
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"\[ScriptCapitalE]" -> \[ScriptCapitalE],
		"\[ScriptCapitalL]" -> \[ScriptCapitalL],
		"\[ScriptCapitalQ]" -> \[ScriptCapitalQ],
		"InitialPhases" -> initPhases
	];
	KerrGeoPlungeFunction[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], assoc]
]


(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


KerrGeoPlunge::fourcomplexroots = "This function can only currently solve plunging geodesics in the paramater space where we have four finite real roots or two finite real roots and two finite complex roots.";
KerrGeoPlunge::highenergy = "This code does not currently support plunging solutions with \[ScriptCapitalE]>1, \[ScriptCapitalE]=`1` has been selected ";
KerrGeoPlunge::rIoutbounds = "Given ISSO radius `1` is not in the allowed range between `2` and `3` for spin `4`.";
KerrGeoPlunge::Inclinationoutofbounds = "Given inclination angle `1` is not between 0 and \[Pi]."


KerrGeoPlunge[a_:0.9,{\[ScriptCapitalE]_:0.8,\[ScriptCapitalL]_:0.3,\[ScriptCapitalQ]_:3}, initPhases:{_,_,_}:Automatic,OptionsPattern[]]:=Module[
	{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},

	If[\[ScriptCapitalE]>=1,
		Message[KerrGeoPlunge::highenergy,\[ScriptCapitalE]];
		Return[$Failed]
		];

	ROOTS = r/.NSolve[(\[ScriptCapitalE](r^2+a^2)-a*\[ScriptCapitalL])^2-(r^2-2r+a^2)(r^2+(a*\[ScriptCapitalE]-\[ScriptCapitalL])^2+\[ScriptCapitalQ])==0,r];

	RealRoots = Select[ROOTS,PossibleZeroQ@Im[#]&];
	ComplexRoots = Select[ROOTS,Not@PossibleZeroQ@Im[#]&];

	If[Length[ComplexRoots]== 2, Return[KerrGeoComplexPlungeMino[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], initPhases,"Roots"->ROOTS]]];
	If[Length[ComplexRoots]== 0, Return[KerrGeoRealPlungeMino[a, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], initPhases,"Roots"->ROOTS]]];
	];


KerrGeoPlunge[a_:0.9,"ISSORadialParam",Arg1_:0.8, initCoords:{_,_,_}:Automatic,OptionsPattern[]]:=
	Module[{RI,DN,RISSOMIN,RISSOMAX},
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


KerrGeoPlunge[a_:0.9,"ISSOIncParam",Arg1_:0.8, initCoords:{_,_,_}:Automatic,OptionsPattern[]]:=
	Module[{},
		(* Other inclinations should be fine as well, right?*)
		If[Not[0<=Arg1<=\[Pi]],
			Message[KerrGeoPlunge::Inclinationoutofbounds, Arg1];
			Return[$Failed]
			];
		KerrGeoISSOPlunge[a, "ISSOIncParam" ,Arg1 ,initCoords]
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
