(* ::Package:: *)

(* ::Title:: *)
(*SpecialOrbits subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`SpecialOrbits`",
	{"KerrGeodesics`ConstantsOfMotion`"}];

KerrGeoPhotonSphereRadius::usage = "KerrGeoPhotonSphereRadius[a,x] returns the radius of the photon sphere."

KerrGeoISCO::usage = "KerrGeoISCO[a,x] returns the location of the innermost stable circular orbit (ISCO) for pro- and retrograde orbits."
KerrGeoISSO::usage = "KerrGeoISSO[a,x] returns the location of the innermost stable spherical orbit (ISSO)."
KerrGeoIBSO::usage = "KerrGeoISBO[a,x] returns the location of the innermost bound spherical orbit (IBSO)."

KerrGeoSeparatrix::usage = "KerrGeoSeparatrix[a,e,x] returns the value of p at the separatrix."
KerrGeoBoundOrbitQ::usage = "KerrGeoBoundOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a bound orbit."

KerrGeoPeriastron::usage = "KerrGeoPeriastron[a,p,e,x] returns the value of the periastron of a scatter orbit."
KerrGeoScatteringAngle::usage = "KerrGeoScatteringAngle[a,p,e,x] returns the value of the scattering angle (\[Phi](\!\(\*SuperscriptBox[\(\[ScriptCapitalI]\), \(+\)]\))-\[Phi](\!\(\*SuperscriptBox[\(\[ScriptCapitalI]\), \(-\)]\))) of a hyperbolic orbit."
KerrGeoDarwinBoundsChi::usage = "KerrGeoDarwinBoundsChi[e] returns the values {\[Chi](\!\(\*SuperscriptBox[\(\[ScriptCapitalI]\), \(-\)]\)),\[Chi](\!\(\*SuperscriptBox[\(\[ScriptCapitalI]\), \(+\)]\)} from the Darwin paramaterisation for a Schwarzschild scatter orbit."
KerrGeoScatterOrbitQ::usage = "KerrGeoScatterOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a scatter orbit."

KerrGeoPlungeOrbitQ::usage = "KerrGeoPlungeOrbitQ[a,p,e,x] tests if the orbital parameters correspond to a plunge orbit."

KerrGeoOrbitType::usage = "KerrGeoOrbitType[a,p,e,x] outputs whether the parameters correspond to a bound, scatter or plunge orbit."

Begin["`Private`"];


(* ::Chapter:: *)
(*Special orbits (separatrix, ISCO, ISSO etc...) *)


(* ::Section::Closed:: *)
(*Innermost stable circular orbit (ISCO)*)


(* ::Text:: *)
(*Schwarzschild ISCO is at r=6M*)


KerrGeoISCO[(0|0.),x_]:=6


(* ::Text:: *)
(*Kerr inner-most circular orbit ISCO from Bardeen, Press, Teukolsky ApJ, 178, p347 (1972), Eq. 2.21*)


KerrGeoISCO[a_,x_/;x^2==1]:=Module[{M=1,Z1,Z2},
	Z1=1+(1-a^2/M^2)^(1/3) ((1+a/M)^(1/3)+(1-a/M)^(1/3));
	Z2=(3a^2/M^2 + Z1^2)^(1/2);
	M(3+Z2-x ((3-Z1)(3+Z1+2Z2)/x^2)^(1/2))
];


(* ::Section::Closed:: *)
(*Photon Sphere*)


(* ::Text:: *)
(*The photon sphere is at 3M for all radii in Schwarzschild*)


KerrGeoPhotonSphereRadius[(0|0.),x_]:=3


(* ::Text:: *)
(*Radius of photon sphere  for equatorial orbits from Bardeen, Press, Teukolsky ApJ, 178, p347 (1972), Eq. 2.18*)


KerrGeoPhotonSphereRadius[a_,1]:=2(1+Cos[2/3 ArcCos[-a]])
KerrGeoPhotonSphereRadius[a_,-1]:=2(1+Cos[2/3 ArcCos[a]])


(* ::Text:: *)
(*For polar orbits the radius was given by E. Teo, General Relativity and Gravitation, v. 35, Issue 11, p. 1909-1926 (2003), Eq. (14)*)


KerrGeoPhotonSphereRadius[a_,(0|0.)]:=1+2Sqrt[1-1/3 a^2]Cos[1/3 ArcCos[(1-a^2)/(1-1/3 a^2)^(3/2)]]


(* ::Text:: *)
(*In the extremal limit we can find the photon sphere radius exactly*)


KerrGeoPhotonSphereRadius[1,x_]:=If[x < Sqrt[3]-1, 1+Sqrt[2] Sqrt[1-x]-x, 1];


(* ::Text:: *)
(*For all other inclinations we have to numerically find the photon sphere radius*)


KerrGeoPhotonSphereRadius[a1_?NumericQ,x0_?NumericQ/;Abs[x0]<=1]/;Precision[{a1,x0}]!=\[Infinity]:=Module[{M=1,a=a1,req,rpolar,\[CapitalPhi],Q,r,u0Sq,prec},
prec=Precision[{a1,x0}];
req=KerrGeoPhotonSphereRadius[a,Sign[x0]];
rpolar=KerrGeoPhotonSphereRadius[a,0];

\[CapitalPhi]=-((r^3-3M r^2+a^2 r+a^2 M)/(a(r-M)));
Q=-((r^3 (r^3-6M r^2+9M^2 r-4a^2 M))/(a^2 (r-M)^2));

u0Sq=((a^2-Q-\[CapitalPhi]^2)+Sqrt[(a^2-Q-\[CapitalPhi]^2)^2+4a^2  Q])/(2a^2);

r/.FindRoot[1-u0Sq-x0^2,Flatten[{r,(req+rpolar)/2,Sort[{req,rpolar}]}],WorkingPrecision->Max[MachinePrecision,prec-1]]//Quiet 
(*The final Quiet[] is there to stop FindRoot complaining about the precision of the argument. 
This seems to be fine near the equatorial plane but might not be ideal for inclincation near the polar orbit*)

]


(* ::Section::Closed:: *)
(*Innermost bound spherical orbits (IBSO)*)


KerrGeoIBSO[0,x_]:= 4


(* ::Text:: *)
(*Equatorial IBSO results from Bardeen, Press, Teukolsky 1972*)


KerrGeoIBSO[a_,1]:= 2-a+2(1-a)^(1/2)
KerrGeoIBSO[a_,-1]:= 2+a+2(1+a)^(1/2)


(* ::Text:: *)
(*At the IBSO E=1. Solve[KerrGeo[a,p,0,0]==1,p] to get the formula for the IBSO for polar orbits*)


KerrGeoIBSO[a_,0]:=Module[{\[Delta]},
	\[Delta]=27 a^4-8 a^6+3 Sqrt[3] Sqrt[27 a^8-16 a^10];
	1+Sqrt[12-4 a^2-(6 Sqrt[6] (-2+a^2))/Sqrt[6-2 a^2+(4 a^4)/\[Delta]^(1/3)+\[Delta]^(1/3)]-(4 a^4)/\[Delta]^(1/3)-\[Delta]^(1/3)]/Sqrt[6]+Sqrt[6-2 a^2+(4 a^4)/\[Delta]^(1/3)+\[Delta]^(1/3)]/Sqrt[6]
]


KerrGeoIBSO[1,(0|0.)]:=1/3 (3+(54-6 Sqrt[33])^(1/3)+(6 (9+Sqrt[33]))^(1/3))


(* ::Text:: *)
(*The below methods come from L. Stein and N. Warburton arXiv:1912.07609*)


IBSOPoly=(-4+p)^2 p^6+a^8 (-1+x^2)^2+2 a^2 p^5 (-8+2 p+4 x^2-3 p x^2)+2 a^6 p^2 (2-5 x^2+3 x^4)+a^4 p^3 (-8 (1-3 x^2+2 x^4)+p (6-14 x^2+9 x^4));


KerrGeoIBSO[a1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,x1}]!=\[Infinity])&&(1>=x1>=0)):=With[{prec=Precision[{a1,x1}]},
p/.FindRoot[IBSOPoly/.{a->a1,x->x1},{p,KerrGeoIBSO[a1,1],KerrGeoIBSO[a1,0]},WorkingPrecision->Max[MachinePrecision,prec-1]]];

KerrGeoIBSO[a1_?NumericQ,x1_?NumericQ]/;(Precision[{a1,x1}]!=\[Infinity])&&(-1<=x1<0):=With[{prec=Precision[{a1,x1}]},
p/.FindRoot[IBSOPoly/.{a->a1,x->x1},{p,KerrGeoIBSO[a1,0],KerrGeoIBSO[a1,-1]},WorkingPrecision->Max[MachinePrecision,prec-1]]];


(* ::Section::Closed:: *)
(*Separatrix*)


(* ::Text:: *)
(*Schwarzschild*)


KerrGeoSeparatrix[0,e_,x_]:= 6+2e;


(* ::Text:: *)
(*From Glampedakis and Kennefick arXiv:gr-qc/0203086, for a=M we have Subscript[p, s]=1+e*)


KerrGeoSeparatrix[1,e_,1]:= 1+e


(* ::Text:: *)
(*Polar ISSO in extremal case found from playing around with the equations (see L. Stein and N. Warburton arXiv:1912.07609)*)


KerrGeoSeparatrix[1,0,0]:=1+Sqrt[3]+Sqrt[3+2 Sqrt[3]]
KerrGeoSeparatrix[1,1,0]:=2/3 (3+(54-6Sqrt[33])^(1/3)+(6(9+Sqrt[33]))^(1/3))


(* ::Text:: *)
(*For e=1 the Subscript[p, s] is at 2 Subscript[r, ibso]*)


KerrGeoSeparatrix[a_,1,x_]:=2KerrGeoIBSO[a,x]


(* ::Text:: *)
(*The below methods come from L. Stein and N. Warburton arXiv:1912.07609*)


SepPoly=-4 (3+e) p^11+p^12+a^12 (-1+e)^4 (1+e)^8 (-1+x)^4 (1+x)^4-4 a^10 (-3+e) (-1+e)^3 (1+e)^7 p (-1+x^2)^4-4 a^8 (-1+e) (1+e)^5 p^3 (-1+x)^3 (1+x)^3 (7-7 x^2-e^2 (-13+x^2)+e^3 (-5+x^2)+7 e (-1+x^2))+8 a^6 (-1+e) (1+e)^3 p^5 (-1+x^2)^2 (3+e+12 x^2+4 e x^2+e^3 (-5+2 x^2)+e^2 (1+2 x^2))-8 a^4 (1+e)^2 p^7 (-1+x) (1+x) (-3+e+15 x^2-5 e x^2+e^3 (-5+3 x^2)+e^2 (-1+3 x^2))+4 a^2 p^9 (-7-7 e+e^3 (-5+4 x^2)+e^2 (-13+12 x^2))+2 a^8 (-1+e)^2 (1+e)^6 p^2 (-1+x^2)^3 (2 (-3+e)^2 (-1+x^2)+a^2 (e^2 (-3+x^2)-3 (1+x^2)+2 e (1+x^2)))-2 p^10 (-2 (3+e)^2+a^2 (-3+6 x^2+e^2 (-3+2 x^2)+e (-2+4 x^2)))+a^6 (1+e)^4 p^4 (-1+x^2)^2 (-16 (-1+e)^2 (-3-2 e+e^2) (-1+x^2)+a^2 (15+6 x^2+9 x^4+e^2 (26+20 x^2-2 x^4)+e^4 (15-10 x^2+x^4)+4 e^3 (-5-2 x^2+x^4)-4 e (5+2 x^2+3 x^4)))-4 a^4 (1+e)^2 p^6 (-1+x) (1+x) (-2 (11-14 e^2+3 e^4) (-1+x^2)+a^2 (5-5 x^2-9 x^4+4 e^3 x^2 (-2+x^2)+e^4 (5-5 x^2+x^4)+e^2 (6-6 x^2+4 x^4)))+a^2 p^8 (-16 (1+e)^2 (-3+2 e+e^2) (-1+x^2)+a^2 (15-36 x^2+30 x^4+e^4 (15-20 x^2+6 x^4)+4 e^3 (5-12 x^2+6 x^4)+4 e (5-12 x^2+10 x^4)+e^2 (26-72 x^2+44 x^4)));
SepEquat=a^4 (-3-2 e+e^2)^2+p^2 (-6-2 e+p)^2-2 a^2 (1+e) p (14+2 e^2+3 p-e p);
SepPolar=a^6 (-1+e)^2 (1+e)^4+p^5 (-6-2 e+p)+a^2 p^3 (-4 (-1+e) (1+e)^2+(3+e (2+3 e)) p)-a^4 (1+e)^2 p (6+2 e^3+2 e (-1+p)-3 p-3 e^2 (2+p));


pEquatPro[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepEquat/.{a->a1,e->e1},{p,1+e1,6+2e1},WorkingPrecision->Max[MachinePrecision,prec-2]]]
pEquatRet[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepEquat/.{a->a1,e->e1},{p,6+2e1,5+e1+4Sqrt[1+e1]},WorkingPrecision->Max[MachinePrecision,prec-2]]]
pPolar[a1_?NumericQ,e1_?NumericQ]/;(Precision[{a1,e1}]!=\[Infinity]):=With[{prec=Precision[{a1,e1}]},
p/.FindRoot[SepPolar/.{a->a1,e->e1},{p,1+Sqrt[3]+Sqrt[3+2Sqrt[3]],8},WorkingPrecision->Max[MachinePrecision,prec-2]]]


KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,1]/;(Precision[{a1,e1}]!=\[Infinity]):=pEquatPro[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,0]/;(Precision[{a1,e1}]!=\[Infinity]):=pPolar[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,-1]/;(Precision[{a1,e1}]!=\[Infinity]):=pEquatRet[a1,e1]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,e1,x1}]!=\[Infinity])&&(1>x1>0)):=With[{prec=Precision[{a1,e1,x1}]},
p/.FindRoot[SepPoly/.{a->a1,x->x1,e->e1},{p,pEquatPro[a1,e1],pPolar[a1,e1]},WorkingPrecision->Max[MachinePrecision,prec-2]]]
KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,x1_?NumericQ]/;((Precision[{a1,e1,x1}]!=\[Infinity])&&(-1<x1<0)):=With[{prec=Precision[{a1,e1,x1}]},
p/.FindRoot[SepPoly/.{a->a1,x->x1,e->e1},{p,pPolar[a1,e1],12},WorkingPrecision->Max[MachinePrecision,prec-2]]]


KerrGeoBoundOrbitQ[a_?NumericQ,p_?NumericQ,e_?NumericQ,x_?NumericQ]:=Module[{ps},
	ps = KerrGeoSeparatrix[a,e,x];
	If[p >= ps && 0 <= e < 1, True, False]
]


(* ::Section::Closed:: *)
(*Innermost stable spherical orbit (ISSO)*)


KerrGeoISSO[a_,x_/;Abs[x]==1]:=KerrGeoISCO[a,x]


KerrGeoISSO[a_,x_]:=KerrGeoSeparatrix[a,0,x]


(* ::Section:: *)
(*Scatter Orbits*)


(* ::Text:: *)
(*Periastron*)


KerrGeoPeriastron[a_,p_,e_/;e>=1,x_]:= p/(1+e);


(* ::Text:: *)
(*Schwarzschild hyperbolic scatter angle*)
(*Defined as \[Phi](SuperPlus[\[ScriptCapitalI]])-\[Phi](SuperMinus[\[ScriptCapitalI]])*)
(*Derived by O. Long*)


KerrGeoScatteringAngle[0,p_,e_/;e>=1,1]:= (4 Sqrt[p]EllipticF[ArcCos[-e^(-1)]/2,(4e)/(6+2e-p)])/Sqrt[-6-2e+p]


(* ::Text:: *)
(*Values of the bounds of \[Chi] from the Darwin parameterization*)


KerrGeoDarwinBoundsChi[e_]:= { -ArcCos[-1/e], ArcCos[-1/e] }


(* ::Text:: *)
(*Test whether an orbit is a scatter orbit*)
(*Currently only for equatorial Schwarzschild orbits*)
(*Lcrit derived by O. Long*)


KerrGeoScatterOrbitQ[0,p_?NumericQ,e_?NumericQ,1]:=Module[{\[CapitalEpsilon], L, Lcrit},
	If[ e == 1 || e == 1,True,
		\[CapitalEpsilon] = KerrGeoEnergy[0,p,e,1];
		L = KerrGeoAngularMomentum[0,p,e,1];
		Lcrit = (Sqrt[8-36\[CapitalEpsilon]^2+27\[CapitalEpsilon]^4-8\[CapitalEpsilon] Sqrt[-8+9\[CapitalEpsilon]^2]+9\[CapitalEpsilon]^3 Sqrt[-8+9\[CapitalEpsilon]^2]])/(Sqrt[2] Sqrt[-1+\[CapitalEpsilon]^2]);
		If[Re[L] >= Lcrit && e > 1, True, False]]
]


(* ::Section::Closed:: *)
(*Plunge Orbits*)


(* ::Text:: *)
(*Test whether an orbit is a plunge orbit*)


KerrGeoPlungeOrbitQ[0,p_?NumericQ,e_?NumericQ,1]:=
	If[KerrGeoBoundOrbitQ[0,p,e,1] == KerrGeoScatterOrbitQ[0,p,e,1] == False, True, False]


(* ::Section::Closed:: *)
(*Orbit type*)


(* ::Text:: *)
(*Output the type of orbit based on the orbital parameters*)


KerrGeoOrbitType[0,p_?NumericQ,e_?NumericQ,1]:=
	If[KerrGeoBoundOrbitQ[0,p,e,1] == True, "Bound",
		If[KerrGeoScatterOrbitQ[0,p,e,1] == True, "Scatter", "Plunge"]]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
