(* ::Package:: *)

(* ::Title:: *)
(*FourVelocity subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`FourVelocity`",
	{"KerrGeodesics`ConstantsOfMotion`",
	"KerrGeodesics`OrbitalFrequencies`"}];

KerrGeoFourVelocity::usage = "KerrGeoVelocity[a,p,e,x] returns the four-velocity components as parametrized functions.";

Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Error messages*)


KerrGeoFourVelocity::parametrization = "Parameterization error: `1`"


(* ::Section::Closed:: *)
(*Schwarzschild*)


(* ::Subsection::Closed:: *)
(*Circular, Equatorial*)


KerrGeoVelocityMino[a_?PossibleZeroQ,p_,e_?PossibleZeroQ,x_,initPhases_ ]:= Module[{En,L,Q,r,z,r1,r2,r3,r4,kr,zp,zm,kz, \[CapitalUpsilon]r, \[CapitalUpsilon]z, 
qr, qz, \[Lambda]local ,qr0, qz0, rprime, zprime, \[CapitalDelta], \[CapitalSigma], \[Omega], utContra,urContra,u\[Theta]Contra,uzContra,u\[Phi]Contra, utCo, urCo, u\[Theta]Co, u\[Phi]Co},

(*Constants of Motion*)
{En,L,Q}= {"\[ScriptCapitalE]","\[ScriptCapitalL]","\[ScriptCapitalQ]"}/.KerrGeoConstantsOfMotion[0,p,0,x];

\[CapitalUpsilon]z = p/Sqrt[-3+p];

{qr0,qz0} = initPhases;

qz[\[Lambda]_] := \[Lambda] \[CapitalUpsilon]z + qz0;

 

utContra= Function[{Global`\[Lambda]},Evaluate[Sqrt[p/(-3+p)] ], Listable];  
urContra:= Function[{Global`\[Lambda]},Evaluate[0],Listable];
u\[Theta]Contra = Function[{Global`\[Lambda]}, Evaluate[(Sqrt[((1-x^2)/(-3+p))] Sin[qz[Global`\[Lambda]]] )/(p Sqrt[1+(-1+x^2) Cos[qz[Global`\[Lambda]]]^2])],Listable];  
u\[Phi]Contra = Function[{Global`\[Lambda]},Evaluate[x/(Sqrt[-3+p] (p+p (-1+x^2) Cos[qz[Global`\[Lambda]]]^2))],Listable]; 

utCo =  Function[{Global`\[Lambda]},Evaluate[-En], Listable];
urCo= Function[{Global`\[Lambda]},Evaluate[0],Listable];
u\[Theta]Co=  Function[{Global`\[Lambda]},Evaluate[(p Sqrt[(1-x^2)/(-3+p)] Sin[qz[Global`\[Lambda]]])/ Sqrt[1+(-1+x^2) Cos[qz[Global`\[Lambda]]]^2]],Listable];   
u\[Phi]Co= Function[{Global`\[Lambda]},Evaluate[L],Listable];

<|"\!\(\*SuperscriptBox[\(u\), \(t\)]\)"->utContra, "\!\(\*SuperscriptBox[\(u\), \(r\)]\)"->urContra, "\!\(\*SuperscriptBox[\(u\), \(\[Theta]\)]\)"-> u\[Theta]Contra, "\!\(\*SuperscriptBox[\(u\), \(\[Phi]\)]\)"->   u\[Phi]Contra,
"\!\(\*SubscriptBox[\(u\), \(t\)]\)"->utCo, "\!\(\*SubscriptBox[\(u\), \(r\)]\)"->urCo, "\!\(\*SubscriptBox[\(u\), \(\[Theta]\)]\)"-> u\[Theta]Co, "\!\(\*SubscriptBox[\(u\), \(\[Phi]\)]\)"->   u\[Phi]Co|>



] 


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Generic (Mino)*)


KerrGeoVelocityMino[a_,p_,e_,x_,initPhases_]:= Module[{En,L,Q,r,z,r1,r2,r3,r4,kr,zp,zm,kz, \[CapitalUpsilon]r, \[CapitalUpsilon]z, 
qr, qz, \[Lambda]local ,qr0, qz0, rprime, zprime, \[CapitalDelta], \[CapitalSigma], \[Omega], utContra,urContra,u\[Theta]Contra,uzContra,u\[Phi]Contra, utCo, urCo, u\[Theta]Co, u\[Phi]Co},

(*Constants of Motion*)
{En,L,Q}= {"\[ScriptCapitalE]","\[ScriptCapitalL]","\[ScriptCapitalQ]"}/.KerrGeoConstantsOfMotion[a,p,e,x];

{r1,r2,r3,r4}=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x,En,Q];

{zp,zm}= KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a, p, e, x];

kr = ((r1-r2)(r3-r4))/((r1-r3)(r2-r4));

kz = a^2 (1-En^2) zm^2/zp^2;

(*Frequencies*)
\[CapitalUpsilon]r = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyr[a,p,e,x,{En,L,Q},{r1,r2,r3,r4}]; 
\[CapitalUpsilon]z = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequency\[Theta][a,p,e,x,{En,L,Q},{zp,zm}];

(*Action Angle Phases*)
{ qr0, qz0} = {initPhases[[1]], initPhases[[2]]};

qr[\[Lambda]_] := \[Lambda] \[CapitalUpsilon]r + qr0;
qz[\[Lambda]_] := \[Lambda] \[CapitalUpsilon]z + qz0  + \[Pi]/2;
 
(*r(qr)*)
r[qr_] := (r3(r1 - r2)JacobiSN[EllipticK[kr]/\[Pi] qr, kr] ^2- r2(r1-r3) )/((r1-r2)JacobiSN[EllipticK[kr] /\[Pi] qr, kr]^2-(r1-r3)) ;  
(*r'(qr)*)
rprime[qr_] := (2 (r1-r2) (r1-r3) (r2-r3) EllipticK[kr] JacobiCN[( qr EllipticK[kr])/\[Pi],kr] JacobiDN[( qr EllipticK[kr])/\[Pi],kr] JacobiSN[( qr EllipticK[kr])/\[Pi],kr])/(\[Pi] (-r1+r3+(r1-r2) JacobiSN[( qr EllipticK[kr])/\[Pi],kr]^2)^2); 

(*z(qz)*)
z[qz_] := zm JacobiSN[EllipticK[kz] (2 qz)/\[Pi], kz];
(*z'(qz)*)
zprime[qz_] := (2 zm EllipticK[kz] JacobiCN[(2 qz EllipticK[kz])/\[Pi],kz] JacobiDN[(2 qz EllipticK[kz])/\[Pi],kz])/\[Pi] ;

\[CapitalDelta][qr_]:= r[qr]^2 + a^2 - 2 r[qr];
\[Omega][qr_] := Sqrt[r[qr]^2+ a^2];  
\[CapitalSigma][qr_,qz_] := r[qr]^2 + a^2 z[qz]^2; 


utContra= Function[{Global`\[Lambda]},Evaluate[1/\[CapitalSigma][qr[Global`\[Lambda]],qz[Global`\[Lambda]]] (\[Omega][qr[Global`\[Lambda]]]^2/\[CapitalDelta][qr[Global`\[Lambda]]] ( \[Omega][qr[Global`\[Lambda] ]]^2 En  - a L) - a^2 (1-z[qz[Global`\[Lambda]]]^2)En + a L)], Listable];
urContra:= Function[{Global`\[Lambda]},Evaluate[( rprime[qr[Global`\[Lambda]]] \[CapitalUpsilon]r)/\[CapitalSigma][qr[Global`\[Lambda]],qz[Global`\[Lambda]]]],Listable];
u\[Theta]Contra = Function[{Global`\[Lambda]}, Evaluate[-(\[CapitalUpsilon]z zprime[qz[Global`\[Lambda]]])/(\[CapitalSigma][qr[Global`\[Lambda]],qz[Global`\[Lambda]]]Sqrt[1-z[qz[Global`\[Lambda]]]^2])],Listable];
u\[Phi]Contra = Function[{Global`\[Lambda]},Evaluate[1/\[CapitalSigma][qr[Global`\[Lambda]],qz[Global`\[Lambda]]] (a/\[CapitalDelta][qr[Global`\[Lambda]]] ( \[Omega][qr[Global`\[Lambda]]]^2 En  - a L) - a En + L/(1-z[qz[Global`\[Lambda]]]^2))],Listable];


utCo =  Function[{Global`\[Lambda]},Evaluate[-En], Listable];
urCo= Function[{Global`\[Lambda]},Evaluate[( rprime[qr[Global`\[Lambda]]] \[CapitalUpsilon]r)/\[CapitalDelta][qr[Global`\[Lambda]]]],Listable];
u\[Theta]Co=  Function[{Global`\[Lambda]},Evaluate[-((\[CapitalUpsilon]z zprime[qz[Global`\[Lambda]]])/Sqrt[1-z[qz[Global`\[Lambda]]]^2])],Listable];
u\[Phi]Co= Function[{Global`\[Lambda]},Evaluate[L],Listable];

<|"\!\(\*SuperscriptBox[\(u\), \(t\)]\)"->utContra, "\!\(\*SuperscriptBox[\(u\), \(r\)]\)"->urContra, 
"\!\(\*SuperscriptBox[\(u\), \(\[Theta]\)]\)"-> u\[Theta]Contra, "\!\(\*SuperscriptBox[\(u\), \(\[Phi]\)]\)"->   u\[Phi]Contra,
"\!\(\*SubscriptBox[\(u\), \(t\)]\)"->utCo, "\!\(\*SubscriptBox[\(u\), \(r\)]\)"->urCo, 
"\!\(\*SubscriptBox[\(u\), \(\[Theta]\)]\)"-> u\[Theta]Co, "\!\(\*SubscriptBox[\(u\), \(\[Phi]\)]\)"->  u\[Phi]Co|>
]



(* ::Subsection::Closed:: *)
(*Equatorial (Darwin)*)


(* ::Subsubsection:: *)
(*Circular Case*)


KerrGeoVelocityDarwin[a_,p_,e_?PossibleZeroQ,x_,initPhases_]:= Module[{ut,ur,u\[Theta],u\[Phi], MinoVelocities,utContra,urContra,u\[Theta]Contra,u\[Phi]Contra,
utCo,urCo,u\[Theta]Co,u\[Phi]Co,utUp,urUp,u\[Theta]Up,u\[Phi]Up, utDown,urDown,u\[Theta]Down,u\[Phi]Down},

MinoVelocities = KerrGeoVelocityMino[a,p,e,x,{0,0}];

utUp="\!\(\*SuperscriptBox[\(u\), \(t\)]\)"; urUp="\!\(\*SuperscriptBox[\(u\), \(r\)]\)"; 
u\[Theta]Up="\!\(\*SuperscriptBox[\(u\), \(\[Theta]\)]\)"; u\[Phi]Up="\!\(\*SuperscriptBox[\(u\), \(\[Phi]\)]\)";
utDown="\!\(\*SubscriptBox[\(u\), \(t\)]\)"; urDown="\!\(\*SubscriptBox[\(u\), \(r\)]\)"; 
u\[Theta]Down="\!\(\*SubscriptBox[\(u\), \(\[Theta]\)]\)"; u\[Phi]Down="\!\(\*SubscriptBox[\(u\), \(\[Phi]\)]\)";

(*All components are Constants*)
utContra = Function[{Global`\[Chi]},Evaluate[MinoVelocities [utUp][Global`\[Chi]]],Listable];
urContra = Function[{Global`\[Chi]},Evaluate[0],Listable];
u\[Theta]Contra = Function[{Global`\[Chi]}, Evaluate[0],Listable];
u\[Phi]Contra = Function[{Global`\[Chi]}, Evaluate[MinoVelocities [u\[Phi]Up][Global`\[Chi]]],Listable];

utCo = Function[{Global`\[Chi]},Evaluate[MinoVelocities [utDown][Global`\[Chi]]],Listable];
urCo = Function[{Global`\[Chi]},Evaluate[0],Listable];
u\[Theta]Co = Function[{Global`\[Chi]}, Evaluate[0],Listable];
u\[Phi]Co = Function[{Global`\[Chi]}, Evaluate[MinoVelocities [u\[Phi]Down][Global`\[Chi]]],Listable];

<|utUp ->utContra, urUp ->urContra, 
u\[Theta]Up-> u\[Theta]Contra, u\[Phi]Up->   u\[Phi]Contra,
utDown->utCo, urDown->urCo, 
u\[Theta]Down-> u\[Theta]Co, u\[Phi]Down->  u\[Phi]Co|>

]


(* ::Subsubsection:: *)
(*Eccentric Case*)


KerrGeoVelocityDarwin[a_,p_,e_,x_/;x^2==1,initPhases_]:= Module[{En,L,Q,r,z,r1,r2,r3,r4,kr, \[CapitalUpsilon]r, \[CapitalLambda]r,yr,\[Lambda]0r,r01,\[CapitalLambda]r1,\[Lambda],
\[Chi]0,\[Nu], \[Chi]local ,qr0, qz0, rprime, zprime, \[CapitalDelta], \[CapitalSigma], \[Omega], ut,ur,u\[Theta],u\[Phi], MinoVelocities,utContra,urContra,u\[Theta]Contra,u\[Phi]Contra,
utCo,urCo,u\[Theta]Co,u\[Phi]Co,utUp,urUp,u\[Theta]Up,u\[Phi]Up, utDown,urDown,u\[Theta]Down,u\[Phi]Down},

(*Constants of Motion*)
{En,L,Q}= {"\[ScriptCapitalE]","\[ScriptCapitalL]","\[ScriptCapitalQ]"}/.KerrGeoConstantsOfMotion[a,p,e,x];

(*Roots*)
{r1,r2,r3,r4}=KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x,En,Q];

kr = ((r1-r2)(r3-r4))/((r1-r3)(r2-r4));

(*Frequencies*)
\[CapitalUpsilon]r = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencyr[a,p,e,x,{En,L,Q},{r1,r2,r3,r4}]; 

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


MinoVelocities = KerrGeoVelocityMino[a,p,e,x,{0,0}];

utUp="\!\(\*SuperscriptBox[\(u\), \(t\)]\)"; urUp="\!\(\*SuperscriptBox[\(u\), \(r\)]\)"; 
u\[Theta]Up="\!\(\*SuperscriptBox[\(u\), \(\[Theta]\)]\)"; u\[Phi]Up="\!\(\*SuperscriptBox[\(u\), \(\[Phi]\)]\)";
utDown="\!\(\*SubscriptBox[\(u\), \(t\)]\)"; urDown="\!\(\*SubscriptBox[\(u\), \(r\)]\)"; 
u\[Theta]Down="\!\(\*SubscriptBox[\(u\), \(\[Theta]\)]\)"; u\[Phi]Down="\!\(\*SubscriptBox[\(u\), \(\[Phi]\)]\)";

utContra = Function[{Global`\[Chi]},Evaluate[MinoVelocities [utUp][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];
urContra = Function[{Global`\[Chi]},Evaluate[MinoVelocities [urUp][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];
u\[Theta]Contra = Function[{Global`\[Chi]}, Evaluate[0],Listable];
u\[Phi]Contra = Function[{Global`\[Chi]}, Evaluate[MinoVelocities [u\[Phi]Up][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];

utCo = Function[{Global`\[Chi]},Evaluate[MinoVelocities [utDown][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];
urCo = Function[{Global`\[Chi]},Evaluate[MinoVelocities [urDown][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];
u\[Theta]Co = Function[{Global`\[Chi]}, Evaluate[0],Listable];
u\[Phi]Co = Function[{Global`\[Chi]}, Evaluate[MinoVelocities [u\[Phi]Down][\[Lambda][Global`\[Chi]-\[Chi]0]]],Listable];

<|utUp ->utContra, urUp ->urContra, 
u\[Theta]Up-> u\[Theta]Contra, u\[Phi]Up->   u\[Phi]Contra,
utDown->utCo, urDown->urCo, 
u\[Theta]Down-> u\[Theta]Co, u\[Phi]Down->  u\[Phi]Co|>

]


(* ::Section::Closed:: *)
(*KerrGeoFourVelocity Wrapper*)


Options[KerrGeoFourVelocity] = {"Parametrization"-> "Mino"}
SyntaxInformation[KerrGeoFourVelocity] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};


KerrGeoFourVelocity[a_,p_,e_,x_,initPhases:{_,_}:{0,0}, OptionsPattern[]]:= Module[{param},
param = OptionValue["Parametrization"];
	If[param == "Darwin",

	If[ Abs[x]!=1, 
		Message[KerrGeoFourVelocity::parametrization, "Darwin parameterization only valid for equatorial motion"];
		Return[];,
		 Return[KerrGeoVelocityDarwin[a,p,e,x,initPhases]]]];


	If[param == "Mino", Return[KerrGeoVelocityMino[a,p,e,x,initPhases]]];

	Message[KerrGeoFourVelocity::parametrization, "Unrecognized Paramaterization: " <> param];

]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
