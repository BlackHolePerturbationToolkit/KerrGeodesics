(* ::Package:: *)

(* ::Title:: *)
(*KerrGeoOrbit subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`KerrGeoPlunge`",
	{"KerrGeodesics`SpecialOrbits`",
	"KerrGeodesics`ConstantsOfMotion`"}];

KerrGeoPlunge::usage = "Takes either KerrGeoPlunge[a, Generic, En ,L,Q] or KerrGeoPlunge[a, ISSO , RI] for generic Plunges or ISSO plunges respectively and returns a KerrGeoPlungeFunction[..] which stores the orbital trajectory and parameters. Here the ISSO plunges are paramaterised in terms of the choice of the radius of the ISSO for a given a there are range of allowed RI's which correspond to differeing inclinations in the prograde and retrograde directions, if an RI outside this range is given the used is provided with the range of allowed values to re-run the funciton.";
KerrGeoPlungeFunction::usage = "KerrGeoPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*ISSO Plunges (Inclination Param)*)


KerrGeoISSOPlungeInc[a_, \[Theta]INC_, initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts, assoc,M=1,\[Chi]INC ,t,RI, r, \[Theta], \[Phi], \[Epsilon], L, Q, R4, RM, RP, KZ, Z1, Z2, J,velocity},
	
	\[Chi]INC = Cos[\[Theta]INC];
	RI = KerrGeoISSO[a,\[Chi]INC];
	\[Epsilon] = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[1]];
	L = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[2]];
	Q = Values[KerrGeoConstantsOfMotion[a,RI,0,\[Chi]INC]][[3]];
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];

	J = 1-\[Epsilon]^2;
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};
	consts = {\[Epsilon], L, Q};
	R4 = (a^2 Q)/(J*RI^3);
	KZ = a^2*J(Z1^2/Z2^2);
	
	r[\[Lambda]_] := ((RI (RI-R4)^2 J*\[Lambda]^2+4*R4)/((RI-R4)^2 J*\[Lambda]^2+4));

	z[\[Lambda]_]:= Z1*JacobiSN[Z2*\[Lambda] + EllipticK[KZ],KZ];


	tr[\[Lambda]_]:= ((a^2+RI^2) (-a L+(a^2+RI^2) \[Epsilon])(\[Lambda]) )/((RI-RM) (RI-RP))+(2 ((R4-RI)^2) \[Epsilon] (\[Lambda]) )/(4+J (R4-RI)^2 (\[Lambda]) ^2)-(R4+3 RI+2 (RM+RP))/Sqrt[J] \[Epsilon] ArcTan[((\[Lambda])  (RI-R4)Sqrt[J ]  )/2 ]+(2 (a^2+RM^2) (-a L+a^2 \[Epsilon]+RM^2 \[Epsilon]) ArcTanh[(2  Sqrt[RM-R4])/((\[Lambda])  (RI-R4)Sqrt[J ] Sqrt[RI-RM])])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP)Sqrt[J])+(2 (a^2+RP^2) (-a L+a^2 \[Epsilon]+RP^2 \[Epsilon]) ArcTanh[(2  Sqrt[RP-R4])/((\[Lambda])  (RI-R4)Sqrt[J ] Sqrt[RI-RP])])/(Sqrt[J]Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP));
	\[Phi]r[\[Lambda]_]:= 1/Sqrt[J] 2 a ((\[Lambda] (RI-R4)Sqrt[J ] (-a L+a^2 \[Epsilon]+RI^2 \[Epsilon]))/( 2(-R4+RI) (RI-RM) (RI-RP))+((-a L+a^2 \[Epsilon]+RM^2 \[Epsilon]) ArcTanh[(2  Sqrt[RM-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RM])])/(Sqrt[RM-R4] (RI-RM)^(3/2) (-RM+RP))+((-a L+a^2 \[Epsilon]+RP^2 \[Epsilon]) ArcTanh[(2  Sqrt[RP-R4])/(\[Lambda] (RI-R4)Sqrt[J ] Sqrt[RI-RP])])/(Sqrt[RP-R4] (RI-RP)^(3/2) (RM-RP))) ;

	tz[\[Lambda]_]:= 1/J \[Epsilon]  (-Z2 EllipticE[JacobiAmplitude[Z2 (\[Lambda]  ) + EllipticK[KZ],KZ],KZ]+(Z2-a^2 J/Z2) EllipticF[JacobiAmplitude[Z2 (\[Lambda]  ) + EllipticK[KZ],KZ],KZ]);
	\[Phi]z[\[Lambda]_]:= L/Z2 EllipticPi[Z1^2,JacobiAmplitude[Z2 \[Lambda] + EllipticK[KZ],KZ],KZ];


	t=Function[{Global`\[Lambda]}, Evaluate[ a*L Global`\[Lambda] + tr[Global`\[Lambda]] + tz[Global`\[Lambda]]], Listable];
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
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {0,KZ},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"rI" -> RI,
		"InitialPhases" -> initPhases
	];
	
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]

]



(* ::Subsection:: *)
(*ISSO Plunges (Radial Param)*)


KerrGeoISSOPlungeRad[a_, RI_, initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts, assoc,M=1, t, r, \[Theta], \[Phi], \[Epsilon], L, Q, R4, RM, RP, KZ, Z1, Z2, J,velocity},
	
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


	t=Function[{Global`\[Lambda]}, Evaluate[ a*L Global`\[Lambda]  +tr[Global`\[Lambda]] + tz[Global`\[Lambda]]], Listable];
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
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {0,KZ},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"rI" -> RI,
		"InitialPhases" -> initPhases
	];
	
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]

]



(* ::Subsection:: *)
(*Real Plunges (Mino)*)


KerrGeoRealPlungeMino[a_, \[Epsilon]_, L_, Q_ , initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts,assoc,M,KR,hR,hP,hM,\[CapitalUpsilon]T,\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]Tr,\[CapitalUpsilon]Tz,\[CapitalUpsilon]z,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Phi]r,\[CapitalUpsilon]\[Phi]z,t,r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots, KZ, Z1, Z2},
	
	M=1;
	J = 1-\[Epsilon]^2;
	RM = M-Sqrt[M^2-a^2];
	RP = M+Sqrt[M^2-a^2];
	
	ROOTS = r/.Solve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2M*r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];

	R1= ROOTS[[3]];
	R2= ROOTS[[4]];
	R3= ROOTS[[1]];
	R4= ROOTS[[2]];
	
	
	
	
	{Z1,Z2}= {Sqrt[1/2 (1+(L^2+Q)/(a^2 J)-Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])],Sqrt[ (a^2 J)/2 (1+(L^2+Q)/(a^2 J)+Sqrt[(1+(L^2+Q)/(a^2 J))^2-(4 Q)/(a^2 J)])]};
	KZ= a^2*J(Z1^2/Z2^2);
	
	KR = ((R3-R4)(R1-R2))/((R3-R1)(R4-R2));
	hR = (R3-R4)/(R3-R1);
	hP = hR (R1-RP)/(R4-RP);
	hM =  hR (R1-RM)/(R4-RM);
	
	\[CapitalUpsilon]r= \[Pi]/(2*EllipticK[KR]) Sqrt[J(R3-R1)(R4-R2)];
	
	\[CapitalUpsilon]z= (\[Pi]*Z2)/(2*EllipticK[KZ]);
	
	
	\[CapitalUpsilon]\[Phi]r = a/(RP-RM) ((2\[Epsilon]*RP-a*L)/(R1-RP) (1-(R4-R1)/(R4-RP) (EllipticPi[hP,KR]/EllipticK[KR]) )- (2\[Epsilon]*RM-a*L)/(R1-RM) (1-(R4-R1)/(R4-RM) (EllipticPi[hM,KR]/EllipticK[KR]))  );
	\[CapitalUpsilon]\[Phi]z  = L/EllipticK[KZ] EllipticPi[Z1^2,KZ];
	\[CapitalUpsilon]\[Phi]= \[CapitalUpsilon]\[Phi]r+\[CapitalUpsilon]\[Phi]z ;
	\[Phi]Tr[\[Xi]_]:= (-2*a*\[Epsilon]*(R4-R1))/((RP-RM)Sqrt[J*(R3-R1)(R4-R2)])*((2*RP-a*L/\[Epsilon])/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],KR]-(2*RM-a*L/\[Epsilon])/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],KR]);
	\[Phi]Tz[\[Xi]_]:= -L/Z2 EllipticPi[Z1^2,\[Xi],KZ];

	\[Phi]r[\[Lambda]_] := \[Phi]Tr[JacobiAmplitude[EllipticK[KR]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],KR]] - \[Phi]Tr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	\[Phi]z[\[Lambda]_]:= -\[Phi]Tz[JacobiAmplitude[EllipticK[KZ]*(2*(\[CapitalUpsilon]z*\[Lambda]))/\[Pi],KZ]] + \[Phi]Tz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];
	

	\[CapitalUpsilon]Tr = (4+a^2)\[Epsilon]+\[Epsilon](1/2 ((4+R3+R4+R1)R1-R3*R4+(R3-R1)(R4-R2) EllipticE[KR]/EllipticK[KR]+(4+R3+R4+R1+R2)(R4-R1) EllipticPi[hR,KR]/EllipticK[KR])+2/(RP-RM) (((4-a*L/\[Epsilon])RP-2*a^2)/(R1-RP) (1-(R4-R1)/(R4-RP) EllipticPi[hP,KR]/EllipticK[KR])-((4-a*L/\[Epsilon])RM-2*a^2)/(R1-RM) (1-(R4-R1)/(R4-RM) EllipticPi[hM,KR]/EllipticK[KR])));
	\[CapitalUpsilon]Tz= -a^2*\[Epsilon]+(\[Epsilon]*Q)/(J(Z1^2)) (1-EllipticE[KZ]/EllipticK[KZ]);
	\[CapitalUpsilon]T= \[CapitalUpsilon]Tr +\[CapitalUpsilon]Tz;
	TTr[\[Xi]_]:= (\[Epsilon](R4-R1))/Sqrt[J(R3-R1)(R4-R2)] ((4+R3+R4+R1+R2)EllipticPi[hR,\[Xi],KR]-4/(RP-RM) ((RP(4-a*L/\[Epsilon])-2*a^2)/((R4-RP)(R1-RP)) EllipticPi[hP,\[Xi],KR]-(RM(4-a*L/\[Epsilon])-2*a^2)/((R4-RM)(R1-RM)) EllipticPi[hM,\[Xi],KR])+((R3-R1)(R4-R2))/(R4-R1) (EllipticE[\[Xi],KR]-(hR*Sin[\[Xi]]*Cos[\[Xi]]Sqrt[1-KR*(Sin[\[Xi]])^2])/(1-hR*(Sin[\[Xi]])^2)));
	TTz[\[Xi]_]:= -\[Epsilon]/J Z2*EllipticE[\[Xi],KZ];
	
	
	r[\[Lambda]_]:=(R1(R3-R4)(JacobiSN[EllipticK[KR]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),KR])^2-R4(R3-R1))/((R3-R4)(JacobiSN[EllipticK[KR]/\[Pi]* (\[CapitalUpsilon]r*\[Lambda]),KR])^2-(R3-R1));
	z[\[Lambda]_]:=  Z1*JacobiSN[EllipticK[KZ]*2 (\[CapitalUpsilon]z*\[Lambda])/\[Pi],KZ];
	
	Trr[\[Lambda]_] := TTr[JacobiAmplitude[EllipticK[KR]*(\[CapitalUpsilon]r*\[Lambda])/\[Pi],KR]] - TTr[\[Pi]]/(2*\[Pi]) \[CapitalUpsilon]r*\[Lambda];
	Tz[\[Lambda]_]:= TTz[JacobiAmplitude[EllipticK[KZ]*(2*\[CapitalUpsilon]z*\[Lambda])/\[Pi],KZ]] - TTz[\[Pi]]/\[Pi] \[CapitalUpsilon]z*\[Lambda];


	t=Function[{Global`\[Lambda]}, Evaluate[Trr[Global`\[Lambda]]+Tz[Global`\[Lambda]] + \[CapitalUpsilon]T*Global`\[Lambda]], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[r[Global`\[Lambda]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda]]]], Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[Phi]r[Global`\[Lambda]]+\[Phi]z[Global`\[Lambda]] + \[CapitalUpsilon]\[Phi]*Global`\[Lambda] ], Listable];



	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {KR,KZ},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"a" -> a,
		"\[Epsilon]" -> \[Epsilon],
		"L" -> L,
		"Q" -> Q,
		"InitialPhases" -> initPhases
	];
	KerrGeoPlungeFunction[a, \[Epsilon], L, Q, assoc]

]



(* ::Subsection:: *)
(*Complex Plunge (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoComplexPlungeMino[a_, \[Epsilon]_, L_, Q_ , initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts, assoc,M,D1M,D1P,D2M,D2P,e,b,c,d,A,B,chi,kr,p2,f, t, r,z, \[Theta], J,\[Phi], R1,R2,R3,R4, RM, RP,ROOTS,RealRoots,ComplexRoots, KZ, Z1, Z2, AMR,CNR,SNR,DNR,AMZ},
	
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
	KZ = a^2*J(Z1^2/Z2^2);
	
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
	
	D1M=Sqrt[ -f];
	D2M=Sqrt[(4 A B (RM-b) (e-RM))/(A( b-RM)-B( e- RM))^2];

	D1P=Sqrt[ -f];
	D2P=Sqrt[(4 A B (RP-b) (e-RP))/(A( b-RP)-B( e- RP))^2];
	
	
	AMR[\[Lambda]_]:=JacobiAmplitude[Sqrt[A B J] \[Lambda],kr^2];
	SNR[\[Lambda]_]:=JacobiSN[Sqrt[A B J] \[Lambda],kr^2];
	CNR[\[Lambda]_]:=JacobiCN[Sqrt[A B J] \[Lambda],kr^2];
	DNR[\[Lambda]_]:=JacobiDN[Sqrt[A B J] \[Lambda],kr^2];

	AMZ[\[Lambda]_]:=JacobiAmplitude[Z2 \[Lambda],KZ];
	
	r[\[Lambda]_] := ((A-B) (A b-B e) SNR[\[Lambda]]^2+2 (A B (b+e)+A B (b-e) CNR[\[Lambda]]))/(4 A B+(A-B)^2 SNR[\[Lambda]]^2);

	z[\[Lambda]_]:= Z1*JacobiSN[Z2 \[Lambda],KZ];


(*Integrals*)
RINT\[Lambda][\[Lambda]_] := ((A b-B e)/(A-B) \[Lambda]-1/ Sqrt[ J] ArcTan[(e-b)/(2 Sqrt[A B])  SNR[\[Lambda]]/Sqrt[1-kr^2 (SNR[\[Lambda]])^2]]+((A+B) (e-b))/(2 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]);

R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)+ (A^2+2 b^2-B^2-2 e^2)/(2 (e-b) Sqrt[ J]) ArcTan[( kr  DNR[\[Lambda]]SNR[\[Lambda]] -I kr^2 ( f +   SNR[\[Lambda]]^2))/(Sqrt[f] kr Sqrt[1+f kr^2])];
RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM))EllipticPi[1/D2M^2,AMR[\[Lambda]],kr^2])/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) -(Sqrt[(e-b)]ArcTanh[( DNR[\[Lambda]]SNR[\[Lambda]]-I kr (D2M^2-SNR[\[Lambda]]^2))/(D2M  Sqrt[1-D2M^2 kr^2])])/(Sqrt[ J] Sqrt[ (e-RM) (RM-b)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))]);


RPINT\[Lambda][\[Lambda]_] := ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP))EllipticPi[1/D2P^2,AMR[\[Lambda]],kr^2])/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) -( Sqrt[(e-b)]ArcTanh[( DNR[\[Lambda]]SNR[\[Lambda]]-I kr (D2P^2-SNR[\[Lambda]]^2))/(D2P  Sqrt[1-D2P^2 kr^2])])/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))]);

tr[\[Lambda]_]:=  ((2 a^2 +RM^2+RM RP+RP^2)\[Epsilon])\[Lambda]+(R2INT\[Lambda][\[Lambda]]+RINT\[Lambda][\[Lambda]](RM+RP)) \[Epsilon] + ((RM^2+a^2)(\[Epsilon](RM^2+a^2)-a*L))/(RM-RP) RMINT\[Lambda][\[Lambda]]+ ((RP^2+a^2)(\[Epsilon](RP^2+a^2)-a*L))/(RP-RM) RPINT\[Lambda][\[Lambda]];
\[Phi]r[\[Lambda]_]:= a( + ((\[Epsilon](RM^2+a^2)-a*L)/(RM-RP))RMINT\[Lambda][\[Lambda]]+ (\[Epsilon](RP^2+a^2)-a*L)/(RP-RM) RPINT\[Lambda][\[Lambda]]);
tz[\[Lambda]_]:= \[Epsilon]/ J ((Z2-a^2 J/Z2) EllipticF[AMZ[\[Lambda]],KZ]-Z2 EllipticE[AMZ[\[Lambda]],KZ]);
\[Phi]z[\[Lambda]_]:= ( L EllipticPi[Z1^2,AMZ[\[Lambda]],(a^2 J Z1^2)/Z2^2])/Z2;


	t=Function[{Global`\[Lambda]}, Evaluate[tr[Global`\[Lambda]] + tz[Global`\[Lambda]]], Listable];
	r=Function[{Global`\[Lambda]}, Evaluate[ r[Global`\[Lambda]] ], Listable];
	\[Theta]=Function[{Global`\[Lambda]}, Evaluate[ ArcCos[z[Global`\[Lambda]]] ], Listable];
	\[Phi]=Function[{Global`\[Lambda]}, Evaluate[\[Phi]r[Global`\[Lambda]] + \[Phi]z[Global`\[Lambda]]], Listable];



	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> \[Epsilon], 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialRoots"-> {R1,R2,R3,R4},
		"PolarRoots"-> {Z1,Z2},
		"EllipticBasis"-> {kr^2,KZ},
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


KerrGeoPlunge[a_:0.9,{\[Epsilon]_:0.8,L_:0.3,Q_:3}, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
(*FIXME: add stability check but make it possible to turn it off*)
	
ROOTS = r/.NSolve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];


RealRoots = Select[ROOTS,Im[#]==0&];
ComplexRoots = Select[ROOTS,Im[#]!=0&];

	If[Length[ComplexRoots]!=0, Return[KerrGeoComplexPlungeMino[a, \[Epsilon], L, Q, initPhases]]];
	If[Length[ComplexRoots]==0, Return[KerrGeoRealPlungeMino[a, \[Epsilon], L, Q, initPhases]]];
			];



KerrGeoPlunge[a_:0.9,PlungeType_:"ISSORadialParam",Arg1_:0.8, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, \[Epsilon], method, RI, DN, RISSOMIN, RISSOMAX, ROOTS, RealRoots,ComplexRoots},
(*FIXME: add stability check but make it possible to turn it off*)

If[PlungeType== "ISSORadialParam" ,
	RI = Arg1;
	DN = (27-45 a^2+17 a^4+a^6+8 a^3 (1-a^2))^(1/3);
	RISSOMIN = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]-1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	RISSOMAX = 3+Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]+1/2 \[Sqrt](72+8 (-6+a^2)-(4 (9-10 a^2+a^4))/DN-4 DN+(64 a^2)/Sqrt[3+a^2+(9-10 a^2+a^4)/DN+DN]);
	If[Between[RI, {RISSOMIN,RISSOMAX} ], Return[KerrGeoISSOPlungeRad[a, RI, initPhases]]];
	Return[Print[StringForm["Please Provide valid value of rI. For a=``, the range of allowed rI's is given by ``." ,a,{RISSOMIN,RISSOMAX} ]]]];
	
If[PlungeType== "ISSOIncParam" ,
	If[Between[Arg1, {0,\[Pi]} ], Return[KerrGeoISSOPlungeInc[a, Arg1, initPhases]]];
	
	Return[Print[StringForm["Please Provide valid value of Inclination between,``." ,{0,\[Pi]} ]]]];
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
