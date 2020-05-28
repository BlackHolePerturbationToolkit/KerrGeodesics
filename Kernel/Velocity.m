(* ::Package:: *)

(* ::Title:: *)
(*FourVelocity subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`Velocity`",
	{"KerrGeodesics`ConstantsOfMotion`"}];

KerrGeoVelocity::usage = "KerrGeoVelocity[a,p,e,x] returns the four-velocity components as parametrized functions.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*Generic (Mino)*)


KerrGeoVelocityMino[a_,p_,e_,x_,initPhases_,index_ ]:= Module[{En,L,Q,r,z,r1,r2,r3,r4,kr,zp,zm,kz, \[CapitalUpsilon]r, \[CapitalUpsilon]z, 
qr, qz, \[Lambda]local ,qr0, qz0, rprime, zprime, \[CapitalDelta], \[CapitalSigma], \[Omega], utContra,urContra,u\[Theta]Contra,uzContra,u\[Phi]Contra, utCo, urCo, u\[Theta]Co, u\[Phi]Co},

(*Constants of Motion*)
{En,L,Q}= Values[KerrGeoConstantsOfMotion[a,p,e,x]];

(*Roots*)
r1 = p/(1-e);
r2 = p/(1+e);
zm = Sqrt[1-x^2];

(*Other Roots*)
r3 = 1/(1-En^2) - (r1 + r2)/2 + Sqrt[(-(1/(1-En^2)) + (r1 + r2)/2 )^2 - (a^2 Q)/(r1 r2 (1 - En^2))];
r4  =  (a^2  Q)/(r1 r2 r3 (1-En^2));

zp = Sqrt[a^2 (1 - En^2) + (L^2)/(1 - zm^2) ]; 
kr = ((r1-r2)(r3-r4))/((r1-r3)(r2-r4));

kz = a^2 (1-En^2) zm^2/zp^2;

(*Frequencies*)
\[CapitalUpsilon]r = \[Pi]/(2 EllipticK[kr]) Sqrt[(1 - En^2)(r1 - r3)(r2 - r4)]; 
\[CapitalUpsilon]z = (\[Pi] zp)/(2EllipticK[kz] ); 

(*Action Angle Phases*)
{ qr0, qz0} = {initPhases[[1]], initPhases[[2]]};

qr = \[Lambda]local \[CapitalUpsilon]r + qr0;
qz = \[Lambda]local \[CapitalUpsilon]z + qz0  + \[Pi]/2;

(*r(qr)*)
r = (r3(r1 - r2)JacobiSN[EllipticK[kr]/\[Pi] qr, kr] ^2- r2(r1-r3) )/((r1-r2)JacobiSN[EllipticK[kr] /\[Pi] qr, kr]^2-(r1-r3)) ;  
(*r'(qr)*)
rprime = (2 (r1-r2) (r1-r3) (r2-r3) EllipticK[kr] JacobiCN[( qr EllipticK[kr])/\[Pi],kr] JacobiDN[( qr EllipticK[kr])/\[Pi],kr] JacobiSN[( qr EllipticK[kr])/\[Pi],kr])/(\[Pi] (-r1+r3+(r1-r2) JacobiSN[( qr EllipticK[kr])/\[Pi],kr]^2)^2); 

(*z(qz)*)
z = zm JacobiSN[EllipticK[kz] (2 qz)/\[Pi], kz];
(*z'(qz)*)
zprime = (2 zm EllipticK[kz] JacobiCN[(2 qz EllipticK[kz])/\[Pi],kz] JacobiDN[(2 qz EllipticK[kz])/\[Pi],kz])/\[Pi] ;

\[CapitalDelta] = r^2 + a^2 - 2r;
\[Omega] = Sqrt[r^2+ a^2];  
\[CapitalSigma] = r^2 + a^2 z^2;

If[index == "Contravariant", 

utContra[\[Lambda]_] := 1/\[CapitalSigma] (\[Omega]^2/\[CapitalDelta] ( \[Omega]^2 En  - a L) - a^2 (1-z^2)En + a L)/. \[Lambda]local-> \[Lambda];

urContra[\[Lambda]_] := ( rprime \[CapitalUpsilon]r)/\[CapitalSigma]/. \[Lambda]local-> \[Lambda];
uzContra = (\[CapitalUpsilon]z zprime)/\[CapitalSigma];  
u\[Theta]Contra[\[Lambda]_] := -(uzContra/Sqrt[1-z^2])/. \[Lambda]local-> \[Lambda];
u\[Phi]Contra[\[Lambda]_]:= 1/\[CapitalSigma] (a/\[CapitalDelta] ( \[Omega]^2 En  - a L) - a En + L/(1-z^2))/. \[Lambda]local-> \[Lambda];

<|"ut"->utContra, "ur"->urContra, "u\[Theta]"-> u\[Theta]Contra, "u\[Phi]"->   u\[Phi]Contra|>,

(*Else if Index \[Equal] Covariant*)

utCo [\[Lambda]_]:=  -En;
urCo[\[Lambda]_]:= ( rprime \[CapitalUpsilon]r)/\[CapitalDelta]/. \[Lambda]local-> \[Lambda];
u\[Theta]Co[\[Lambda]_]:=  -((\[CapitalUpsilon]z zprime)/Sqrt[1-z^2])/. \[Lambda]local-> \[Lambda];
u\[Phi]Co[\[Lambda]_]:= L;

<|"ut"->utCo, "ur"->urCo, "u\[Theta]"-> u\[Theta]Co, "u\[Phi]"->   u\[Phi]Co|>
]
]


(* ::Subsection:: *)
(*Equatorial (Darwin)*)


KerrGeoVelocityDarwin[a_,p_,e_,x_/;x^2==1,initPhases_,index_ ]:= Module[{En,L,Q,r,z,r1,r2,r3,r4,kr, \[CapitalUpsilon]r, \[CapitalLambda]r,yr,\[Lambda]0r,r01,\[CapitalLambda]r1,\[Lambda],
\[Chi]0,\[Nu], \[Chi]local ,qr0, qz0, rprime, zprime, \[CapitalDelta], \[CapitalSigma], \[Omega], ut,ur,u\[Theta],u\[Phi], MinoVelocities},

(*Constants of Motion*)
{En,L,Q}= Values[KerrGeoConstantsOfMotion[a,p,e,x]];

(*Roots*)
r1 = p/(1-e);
r2 = p/(1+e);

(*Other Roots*)
r3 = 1/(1-En^2) - (r1 + r2)/2 + Sqrt[(-(1/(1-En^2)) + (r1 + r2)/2 )^2 - (a^2 Q)/(r1 r2 (1 - En^2))];
r4  =  (a^2  Q)/(r1 r2 r3 (1-En^2));
kr = ((r1-r2)(r3-r4))/((r1-r3)(r2-r4));

(*Frequencies*)
\[CapitalUpsilon]r = \[Pi]/(2 EllipticK[kr]) Sqrt[(1 - En^2)(r1 - r3)(r2 - r4)]; 

(*Initial Phase*)
\[Chi]0 = initPhases[[1]];

(*r(\[Nu] = \[Chi] - \[Chi]0)*)
r[\[Nu]_] = p /(1-e  Cos[\[Nu]] ) ; 

(*Conversion from Mino time to Darwin Parameter*)
\[CapitalLambda]r = (2\[Pi])/\[CapitalUpsilon]r;
yr[r_]:=Sqrt[(r1-r3)/(r1-r2) (r-r2)/(r-r3)];
\[Lambda]0r[r_]:=1/Sqrt[1-En^2] 2/Sqrt[(r1-r3)(r2-r4)] EllipticF[ArcSin[yr[r]],kr];
r01=r2;
\[CapitalLambda]r1=\[Lambda]0r[r01];

\[Lambda][\[Nu]_]:=\[CapitalLambda]r Floor[\[Nu]/(2\[Pi])]+If[Mod[\[Nu],2\[Pi]]<=\[Pi], \[Lambda]0r[r[\[Nu]]]-\[CapitalLambda]r1,\[CapitalLambda]r-\[Lambda]0r[r[\[Nu]]]];


    MinoVelocities = KerrGeoVelocityMino[a,p,e,x,{0,0}, index];

ut[\[Chi]_] := MinoVelocities ["ut"][\[Lambda][\[Chi]-\[Chi]0]];
u\[Theta][\[Chi]_]:= 0;
ur[\[Chi]_] := MinoVelocities ["ur"][\[Lambda][\[Chi]-\[Chi]0]];
u\[Phi][\[Chi]_] := MinoVelocities ["u\[Phi]"][\[Lambda][\[Chi]-\[Chi]0]];

<|"ut"-> ut, "ur"-> ur, "u\[Theta]"-> u\[Theta], "u\[Phi]"-> u\[Phi] |>


]


(* ::Section:: *)
(*KerrGeoVelocity Wrapper*)


Options[KerrGeoVelocity] = {"Index" -> "Contravariant", "Parametrization"-> "Mino"}
SyntaxInformation[KerrGeoVelocity] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};


KerrGeoVelocity[a_,p_,e_,x_,initPhases:{_,_}:{0,0},OptionsPattern[] ]:= Module[{param, index},
param = OptionValue["Parametrization"];
index= OptionValue["Index"];

If[index == "Contravariant"|| index =="Covariant",

	If[param == "Darwin",

	If[ Abs[x]!=1, 
		Print["Darwin parameterization only valid for equatorial motion"];
		Return[];,
		 Return[KerrGeoVelocityDarwin[a,p,e,x,initPhases, index]]]];


	If[param == "Mino", Return[KerrGeoVelocityMino[a,p,e,x,initPhases, index]]];

	Print["Unrecognized Paramaterization: " <>param];,

	Print["Unrecognized Index: " <>index];
]
]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
