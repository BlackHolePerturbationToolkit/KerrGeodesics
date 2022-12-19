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


KerrGeoPlunge::usage = "Takes either KerrGeoPlunge[a, Generic, En ,L,Q] or KerrGeoPlunge[a, ISSO , RI] for generic Plunges or ISSO plunges respectively and returns a KerrGeoPlungeFunction[..] which stores the orbital trajectory and parameters. Here the ISSO plunges are paramaterised in terms of the choice of the radius of the ISSO for a given a there are range of allowed RI's which correspond to differeing inclinations in the prograde and retrograde directions, if an RI outside this range is given the used is provided with the range of allowed values to re-run the funciton.";
KerrGeoPlungeFunction::usage = "KerrGeoPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*ISSO Plunges*)


KerrGeoISSOPlunge[a_, PlungeType_  ,Arg_, initCoords:{_,_,_}:{"NaN","NaN","NaN"}] := Module[{consts,Mino, \[Theta]0, z, r0, z0, \[Phi]0, initConditions, t0, assoc,M=1,\[Chi]INC ,t,RI, r, \[Theta], \[Phi], \[Epsilon], L, Q, R4, RM, RP, kz, Z1, Z2, J,velocity},
	
	
	RM = 1-Sqrt[1-a^2];
	RP = 1+Sqrt[1-a^2];
	
	If[PlungeType == "ISSORadialParam",
	RI = Arg;
	Q =- M (RI)^(5/2) ((Sqrt[(RI-RP)(RI-RM)]-2Sqrt[RI])^2-4a^2)/(4a^2 (RI^(3/2)-Sqrt[RI]-Sqrt[(RI-RP)(RI-RM)]));
	\[Epsilon]=Sqrt[a^2 Q-2 RI^3+3 RI^4]/(Sqrt[3] RI^2);
	L=Sqrt[3 a^2 Q-a^2 RI^2-Q RI^2+3 RI^4+a^2 RI^2 \[Epsilon]^2-3 RI^4 \[Epsilon]^2]/RI;];
	
	If[PlungeType == "ISSOIncParam",
	\[Chi]INC = Cos[Arg];
	RI = KerrGeoISSO[a,\[Chi]INC];
	\[Epsilon] = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[1]];
	L = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[2]];
	Q = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[3]];];
	
	J = 1-\[Epsilon]^2;
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};
	consts = {\[Epsilon], L, Q};
	R4 = (a^2 Q)/(J*RI^3);
	kz = a*Sqrt[J](Z1/Z2);
	
	If[initCoords[1]>RI,Return[Print["Picked r0 greater than ISSO radius"]]];
	
	
	
	If[initCoords=={"NaN","NaN","NaN"},
					{t0,r0,\[Phi]0}={0,R4,0};];	
	If[initCoords!={"NaN","NaN","NaN"},
					{t0,r0,\[Phi]0}=initCoords;];	

	
	MinoR[x_] :=(2 Sqrt[x-R4])/Sqrt[J (RI-x) (R4-RI)^2] - (2 Sqrt[r0-R4])/Sqrt[J (RI-x) (R4-RI)^2];
	Mino=Function[{Global`r}, Evaluate[MinoR[Global`r]],Listable];
	
	r[\[Lambda]_] := ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4));

	z[\[Lambda]_]:= Z1*JacobiSN[Z2*\[Lambda] ,kz^2];
	tr[\[Lambda]_]:= ((a^2+RI^2) (-a L+(a^2+RI^2) \[Epsilon])(\[Lambda]) )/((RI-RM) (RI-RP))+(2 ((R4-RI)^2) \[Epsilon] (\[Lambda]) )/(4+J (R4-RI)^2 (\[Lambda]) ^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[Epsilon] ArcTan[((\[Lambda])  (RI-R4)Sqrt[J ]  )/2 ]+( (a^2+RM^2) (-a L+a^2 \[Epsilon]+RM^2 \[Epsilon])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP)Sqrt[J])+( (a^2+RP^2) (-a L+a^2 \[Epsilon]+RP^2 \[Epsilon]) Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP));
	\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a (((\[Lambda]) (RI-R4)Sqrt[J ] (-a L+a^2 \[Epsilon]+RI^2 \[Epsilon]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a L+a^2 \[Epsilon]+RM^2 \[Epsilon])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]+2  Sqrt[RM-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM]-2  Sqrt[RM-R4])^2]])/(2Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP))+((-a L+a^2 \[Epsilon]+RP^2 \[Epsilon])Log[Sqrt[(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]+2  Sqrt[RP-R4])^2/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP]-2  Sqrt[RP-R4])^2]])/(2Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP)))    ;
	tz[\[Lambda]_]:= 1/J \[Epsilon]  (-Z2 EllipticE[JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2]+(Z2-a^2 J/Z2) EllipticF[JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2]);
	\[Phi]z[\[Lambda]_]:= L/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] ,kz^2],kz^2];


	t=Function[{Global`\[Lambda]}, Evaluate[ a*L Global`\[Lambda] + tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ MinoR[r0]]-tr[MinoR[r0]]-tz[MinoR[r0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda] + MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + MinoR[r0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[-a*\[Epsilon] Global`\[Lambda] + \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ MinoR[r0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[MinoR[r0]] + \[Phi]0], Listable];



	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {RI,RI,RI,R4},
		"Mino"-> Mino,
		"HorizonCrossingTimeMino"-> {Mino[RM],Mino[RP]},
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {0,kz^2},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"HorizonCrossingTimeMino"-> {Mino[RM],Mino[RP]},
		"a" -> a,
		"rI" -> RI,
		"InitialPhases" -> initPhases
	];
	
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]

]



Clear["Global`*"]
Solve[r== ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4)),\[Lambda]]


(* ::Subsection:: *)
(*Real Plunges (Mino)*)


KerrGeoRealPlungeMino[a_, \[Epsilon]_, L_, Q_ , initCoords:{_,_,_}:{0,0,0}] := Module[{consts,assoc,M,MinoR,T0,r0,\[Phi]0,kr,hR,hP,hM,\[CapitalUpsilon]T,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]Tr,\[CapitalUpsilon]Tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots, kz, Z1, Z2},
	
	M=1;
	J = 1-\[Epsilon]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	ROOTS = r/.Solve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2M*r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];

	If[Length[Select[ROOTS,#>=RP&]>=3],
	R1= ROOTS[[3]];
	R2= ROOTS[[4]];
	R3= Min[Select[ROOTS,#>=RP&]];
	R4= Max[Select[ROOTS,#<=RM&]];];
	If[Length[Select[ROOTS,#>=RP&]<3],
	R1= ROOTS[[2]];
	R2= ROOTS[[1]];
	R3= Min[Select[ROOTS,#>=RP&]];
	R4= Max[Select[ROOTS,#<=RM&]];];
	
	If[initCoords=={"NaN","NaN","NaN"},
					{T0,r0,\[Phi]0}={0,R4,0};];	
	If[initCoords!={"NaN","NaN","NaN"},
					{T0,r0,\[Phi]0}=initCoords;];	
	
	
	
	
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};
	kz= a*Sqrt[J](Z1/Z2);
	
	kr = Sqrt[((R3-R4)(R1-R2))/((R3-R1)(R4-R2))];
	hR = (R3-R4)/(R3-R1);
	hP = hR (R1-RP)/(R4-RP);
	hM =  hR (R1-RM)/(R4-RM);
	

	MinoR[x_]:= (2 (InverseJacobiSN[Sqrt[(-R1+R3) (r-R4)]/Sqrt[(r-R1) (R3-R4)] ,KR] - InverseJacobiSN[Sqrt[(-R1+R3) (r0-R4)]/Sqrt[(r0-R1) (R3-R4)] ,KR]))/Sqrt[J (R1-R3) (R2-R4)];
	Mino=Function[{Global`r}, Evaluate[MinoR[Global`r]],Listable];
	
	r=Function[{Global`\[Lambda]}, Evaluate[r[Global`\[Lambda]+ MinoR[r0]] ], Listable];
	
	\[CapitalUpsilon]r= \[Pi]/(2*EllipticK[kr]) Sqrt[J(R3-R1)(R4-R2)];
	
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[kz^2]);
	
	
	\[CapitalUpsilon]\[Phi]r = a/(RP-RM) ((2\[Epsilon]*RP-a*L)/(R1-RP) (1-(R4-R1)/(R4-RP) (EllipticPi[hP,kr]/EllipticK[kr]) )- (2\[Epsilon]*RM-a*L)/(R1-RM) (1-(R4-R1)/(R4-RM) (EllipticPi[hM,kr]/EllipticK[kr]))  );
	\[CapitalUpsilon]\[Phi]z  = L/EllipticK[kz^2] EllipticPi[Z1^2,kz^2];
	\[CapitalUpsilon]\[Phi]= \[CapitalUpsilon]\[Phi]r+\[CapitalUpsilon]\[Phi]z ;
	\[Phi]Tr[\[Xi]_]:= (-2*a*\[Epsilon]*(R4-R1))/((RP-RM)Sqrt[J*(R3-R1)(R4-R2)])*((2*RP-a*L/\[Epsilon])/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(2*RM-a*L/\[Epsilon])/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr]);
	\[Phi]Tz[\[Xi]_]:= -L/Z2 EllipticPi[Z1^2,\[Xi],kz^2];

	\[Phi]r[\[Lambda]_] := \[Phi]Tr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - \[Phi]Tr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	\[Phi]z[\[Lambda]_]:= -\[Phi]Tz[JacobiAmplitude[EllipticK[kz^2]*(2*(\[CapitalUpsilon]z*\[Lambda]))/\[Pi],kz^2]] + \[Phi]Tz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];
	

	\[CapitalUpsilon]Tr = (4+a^2)\[Epsilon]+\[Epsilon](1/2 ((4+R3+R4+R1)R1-R3*R4+(R3-R1)(R4-R2) EllipticE[kr]/EllipticK[kr]+(4+R3+R4+R1+R2)(R4-R1) EllipticPi[hR,kr]/EllipticK[kr])+2/(RP-RM) (((4-a*L/\[Epsilon])RP-2*a^2)/(R1-RP) (1-(R4-R1)/(R4-RP) EllipticPi[hP,kr]/EllipticK[kr])-((4-a*L/\[Epsilon])RM-2*a^2)/(R1-RM) (1-(R4-R1)/(R4-RM) EllipticPi[hM,kr]/EllipticK[kr])));
	\[CapitalUpsilon]Tz= -a^2*\[Epsilon]+(\[Epsilon]*Q)/(J(Z1^2)) (1-EllipticE[kz^2]/EllipticK[kz^2]);
	\[CapitalUpsilon]T= \[CapitalUpsilon]Tr +\[CapitalUpsilon]Tz;
	TTr[\[Xi]_]:= (\[Epsilon](R4-R1))/Sqrt[J(R3-R1)(R4-R2)] ((4+R3+R4+R1+R2)EllipticPi[hR,\[Xi],kr]-4/(RP-RM) ((RP(4-a*L/\[Epsilon])-2*a^2)/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],kr]-(RM(4-a*L/\[Epsilon])-2*a^2)/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],kr])+((R3-R1)(R4-R2))/(R4-R1) (EllipticE[\[Xi],kr]-(hR*Sin[\[Xi]]*Cos[\[Xi]]Sqrt[1-kr*(Sin[\[Xi]])^2])/(1-hR*(Sin[\[Xi]])^2)));
	TTz[\[Xi]_]:= -\[Epsilon]/J Z2*EllipticE[\[Xi],kz^2];
	
	
	r[\[Lambda]_]:=(R1(R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-R4(R3-R1))/((R3-R4)(JacobiSN[EllipticK[kr]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),kr])^2-(R3-R1));
	z[\[Lambda]_]:=  Z1*JacobiSN[EllipticK[kz^2]*2 (\[CapitalUpsilon]z*\[Lambda])/\[Pi],kz^2];
	
	Trr[\[Lambda]_] := TTr[JacobiAmplitude[EllipticK[kr]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],kr]] - TTr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	Tz[\[Lambda]_]:= TTz[JacobiAmplitude[EllipticK[kz^2]*(2*\[CapitalUpsilon]z*\[Lambda])/\[Pi],kz^2]] - TTz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];


	t=Function[{Global`\[Lambda]}, Evaluate[Trr[Global`\[Lambda]+ MinoR[r0]]+Tz[Global`\[Lambda]+ MinoR[r0]] + \[CapitalUpsilon]T*Global`\[Lambda]-Trr[Global`\[Lambda]]-Tz[Global`\[Lambda]] + T0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[r[Global`\[Lambda]+ MinoR[r0]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda]+ MinoR[r0]]]], Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[Phi]r[Global`\[Lambda]+ MinoR[r0]]+\[Phi]z[Global`\[Lambda]+ MinoR[r0]] + \[CapitalUpsilon]\[Phi]*Global`\[Lambda] -\[Phi]r[ MinoR[r0]]-\[Phi]z[MinoR[r0]]+ \[Phi]0], Listable];



	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"Mino"-> Mino,
		"HorizonCrossingTimeMino"-> {Mino[RM],Mino[RP]},
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {kr^2,kz^2},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"Frequencies"-> {\[CapitalUpsilon]T, \[CapitalUpsilon]r,\[CapitalUpsilon]z,\[CapitalUpsilon]\[Phi]},
		"\[Epsilon]" -> \[Epsilon],
		"L" -> L,
		"Q" -> Q,
		"InitialPhases" -> initPhases
	];
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]

]



(* ::Subsection::Closed:: *)
(*Complex Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoComplexPlungeMino[a_, \[Epsilon]_, L_, Q_ , initCoords:{_,_,_}:{"NaN","NaN","NaN"}] := Module[{consts, assoc,MinoR,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],t0,r0,\[Phi]0,M,D1M,D1P,D2M,D2P,e,b,c,d,A,B,chi,kr,p2,f, t, r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots,kz,Z1, Z2, AMR,CNR,SNR,DNR,AMZ},
	
	M=1;
	J = 1-\[Epsilon]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	ROOTS = r/.Solve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2M*r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];
	RealRoots = Select[ROOTS,Im[#]==0&];
	ComplexRoots = Select[ROOTS,Im[#]!=0&];

	R1= RealRoots[[1]];
	R2= RealRoots[[2]];
	R3=ComplexRoots[[1]];
	R4= ComplexRoots[[2]];
	
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};

	
	
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
	
	\[CapitalUpsilon]r = 2 \[Pi] EllipticK[kr^2]/Sqrt[A B (1-\[Epsilon]^2) ];
	\[CapitalUpsilon]\[Theta] = \[Pi]/2 EllipticK[kz^2]/Z2;
	
	
	If[initCoords=={"NaN","NaN","NaN"},
					{t0,r0,\[Phi]0}={0,R1,0};];	
	If[initCoords!={"NaN","NaN","NaN"},
					{t0,r0,\[Phi]0}=initCoords;];	
	
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

	R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)- (A^2+2 b^2-B^2-2 e^2)/(4 (e-b) Sqrt[ J]) ArcTan[(2( kr  Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])(Sqrt[f] kr Sqrt[1+f kr^2]))/(  kr^4 ( f +   SNR[\[Lambda]]^2)^2 -f kr^2 ( 1 +   f kr^2) +kr^2  (1-kr^2 SNR[\[Lambda]]^2)SNR[\[Lambda]]^2)];
	RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM)))/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) EllipticPi[1/D2M^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]-  Sqrt[(e-b)]/(Sqrt[ J] Sqrt[ (RM-b) (e-RM)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))]) 1/4 (Log[((D2M  Sqrt[1-D2M^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)/((D2M  Sqrt[1-D2M^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2M^2-SNR[\[Lambda]]^2))^2)]);

	RPINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP)))/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) EllipticPi[1/D2P^2,JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2],kr^2]- Sqrt[(e-b)]/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))]) 1/4 (Log[((D2P  Sqrt[1-D2P^2 kr^2]+Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)/((D2P  Sqrt[1-D2P^2 kr^2]- Sqrt[1-kr^2 SNR[\[Lambda]]^2]SNR[\[Lambda]])^2+( kr (D2P^2-SNR[\[Lambda]]^2))^2)]);
	tr[\[Lambda]_]:=  ((2 a^2 +RM^2+RM RP+RP^2)\[Epsilon])\[Lambda]+(R2INT\[Lambda][\[Lambda]]+RINT\[Lambda][\[Lambda]](RM+RP)) \[Epsilon] + ((RM^2+a^2)(\[Epsilon](RM^2+a^2)-a*L))/(RM-RP) RMINT\[Lambda][\[Lambda]]+ ((RP^2+a^2)(\[Epsilon](RP^2+a^2)-a*L))/(RP-RM) RPINT\[Lambda][\[Lambda]];
	\[Phi]r[\[Lambda]_]:= a(((\[Epsilon](RM^2+a^2)-a*L)/(RM-RP))RMINT\[Lambda][\[Lambda]]+ (\[Epsilon](RP^2+a^2)-a*L)/(RP-RM) RPINT\[Lambda][\[Lambda]]);
	tz[\[Lambda]_]:= \[Epsilon]/ J ((Z2-a^2 J/Z2) EllipticF[AMZ[\[Lambda]],kz^2]-Z2 EllipticE[AMZ[\[Lambda]],kz^2]);
	\[Phi]z[\[Lambda]_]:= (L EllipticPi[Z1^2,AMZ[\[Lambda]],(a^2 J Z1^2)/Z2^2])/Z2;

	t=Function[{Global`\[Lambda]}, Evaluate[  tr[Global`\[Lambda]+ MinoR[r0]] + tz[Global`\[Lambda]+ MinoR[r0]]-tr[MinoR[r0]]-tz[MinoR[r0]] + t0], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]+ MinoR[r0]]], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda] + MinoR[r0]]]] , Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[ \[Phi]r[Global`\[Lambda]+ MinoR[r0]] + \[Phi]z[Global`\[Lambda]+ MinoR[r0]] -  \[Phi]r[MinoR[r0]] - \[Phi]z[MinoR[r0]] + \[Phi]0], Listable];


	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
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
		"\[Epsilon]" -> \[Epsilon],
		"L" -> L,
		"Q" -> Q,
		"InitialPhases" -> initPhases
	];
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]
]



(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


KerrGeoPlunge[a_:0.9,{\[Epsilon]_:0.8,L_:0.3,Q_:3}, initPhases:{_,_,_}:{"NaN","NaN","NaN"},OptionsPattern[]]:=Module[{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
(*FIXME: add stability check but make it possible to turn it off*)
	
ROOTS = r/.NSolve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];


RealRoots = Select[ROOTS,Im[#]==0&];
ComplexRoots = Select[ROOTS,Im[#]!=0&];

	If[Length[ComplexRoots]!=0, Return[KerrGeoComplexPlungeMino[a, \[Epsilon], L, Q, initPhases]]];
	If[Length[ComplexRoots]==0, Return[KerrGeoRealPlungeMino[a, \[Epsilon], L, Q, initPhases]]];
			];



KerrGeoPlunge[a_:0.9,PlungeType_:"ISSORadialParam",Arg1_:0.8, initCoords:{_,_,_}:{"NaN","NaN","NaN"},OptionsPattern[]]:=Module[{param, \[Epsilon], method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
(*FIXME: add stability check but make it possible to turn it off*)


If[PlungeType== "ISSORadialParam" ,
	RI = Arg1;
	DN = (27-45 a^2+17 a^4+a^6+8 a^3 (1-a^2))^(1/3);
	RISSOMIN = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]-1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	RISSOMAX = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]+1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	If[Between[RI, {RISSOMIN,RISSOMAX} ], Return[KerrGeoISSOPlunge[a,"ISSORadialParam" ,RI, initCoords]]];
	Return[Print[StringForm["Please Provide valid value of rI. For a=``, the range of allowed rI's is given by ``." ,a,{RISSOMIN,RISSOMAX} ]]]];
	
If[PlungeType== "ISSOIncParam" ,
	If[Between[Arg1, {0,\[Pi]} ], Return[KerrGeoISSOPlunge[a, "ISSOIncParam" ,Arg1 ,initCoords]]];
	
	Return[Print[StringForm["Please Provide valid value of Inclination between,``." ,{0,\[Pi]} ]]]];
If[PlungeType!="ISSORadialParam",
	If[ PlungeType!= "ISSORadialParam" , Return[Print["Give a valid ISSO paramaterisation"]]]  ];

]





KerrGeoPlungeFunction /:
 MakeBoxes[kgof:KerrGeoPlungeFunction[a_,\[Epsilon]_,L_,Q_,assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"\[Epsilon]: ", \[Epsilon]}], "  ",
                  BoxForm`SummaryItem[{"L: ", L}], "  ",
                  BoxForm`SummaryItem[{"Q: ", Q}]}],
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


KerrGeoPlungeFunction[a_, \[Epsilon]_,L_,Q_,assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrGeoPlungeFunction[a_,\[Epsilon]_,L_,Q_,assoc_][y_?StringQ] := assoc[y]
Keys[g_KerrGeoPlungeFunction]^:=Keys[g[[5]]]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
