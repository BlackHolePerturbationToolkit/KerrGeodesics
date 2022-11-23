(* ::Package:: *)

(* ::Title:: *)
(*KerrGeoOrbit subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`KerrGeoPlunge`"];

KerrGeoPlunge::usage = "KerrGeoPlunge[a,rI] returns a KerrGeoPlungeFunction[..] which stores the orbital trajectory and parameters.";
KerrGeoPlungeFunction::usage = "KerrGeoPlungeFunction[a,rI,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
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

R2INT\[Lambda][\[Lambda]_]:=\[Lambda] /(A-B) (A b^2-B e^2)+ Sqrt[A B ]/Sqrt[ J] (EllipticE[AMR[\[Lambda]],kr^2])-((A+B) (A^2+2 b^2-B^2-2 e^2))/(4 (A-B) Sqrt[A B J]) EllipticPi[-(1/f),AMR[\[Lambda]],kr^2]-(Sqrt[A B ] (A+B-(A-B)CNR[\[Lambda]]))/((A-B) Sqrt[ J]) (SNR[\[Lambda]] DNR[\[Lambda]])/(f+(SNR[\[Lambda]])^2)+ (A^2+2 b^2-B^2-2 e^2)/(2 (e-b) Sqrt[ J]) Re[ArcTan[( kr  DNR[\[Lambda]]SNR[\[Lambda]] -I kr^2 ( f +   SNR[\[Lambda]]^2))/(Sqrt[f] kr Sqrt[1+f kr^2])]];
RMINT\[Lambda][\[Lambda]_] :=  ((A-B) \[Lambda])/(A (b-RM)-B (e-RM))+((e-b) (A (b-RM)+B (e-RM))EllipticPi[1/D2M^2,AMR[\[Lambda]],kr^2])/(2 Sqrt[A B J] (b-RM) (-e+RM) (A (b-RM)-B (e-RM))) -Re[(Sqrt[(e-b)]ArcTanh[( DNR[\[Lambda]]SNR[\[Lambda]]-I kr (D2M^2-SNR[\[Lambda]]^2))/(D2M  Sqrt[1-D2M^2 kr^2])])/(Sqrt[ J] Sqrt[ (e-RM) (RM-b)] Sqrt[ (A^2 (RM-b)-(e-RM) (b^2-B^2+e RM-b (e+RM)))])];


RPINT\[Lambda][\[Lambda]_] := ((A-B) \[Lambda])/(A (b-RP)+B (-e+RP))+((e-b) (A (b-RP)+B (e-RP))EllipticPi[1/D2P^2,AMR[\[Lambda]],kr^2])/(2 Sqrt[A B J] (b-RP) (-e+RP) (A (b-RP)+B (-e+RP))) -Re[( Sqrt[(e-b)]ArcTanh[( DNR[\[Lambda]]SNR[\[Lambda]]-I kr (D2P^2-SNR[\[Lambda]]^2))/(D2P  Sqrt[1-D2P^2 kr^2])])/(Sqrt[ J] Sqrt[(RP-b) (e-RP)] Sqrt[ (A^2 (RP-b)-(e-RP) (b^2-B^2+e RP-b (e+RP)))])];

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


KerrGeoPlunge[a_,\[Epsilon]_,L_,Q_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, method, ROOTS, RealRoots,ComplexRoots},
(*FIXME: add stability check but make it possible to turn it off*)

method = "Analytical";
param = "MINO";

ROOTS = r/.Solve[(\[Epsilon](r^2+a^2)-a*L)^2-(r^2-2r+a^2)(r^2+(a*\[Epsilon]-L)^2+Q)==0,r];

RealRoots = Select[ROOTS,Im[#]==0&];
ComplexRoots = Select[ROOTS,Im[#]!=0&];

If[Length[ComplexRoots]!=0, Return[KerrGeoComplexPlungeMino[a, \[Epsilon], L, Q, initPhases]]];
If[Length[ComplexRoots]==0, Return[KerrGeoRealPlungeMino[a, \[Epsilon], L, Q, initPhases]]];



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
