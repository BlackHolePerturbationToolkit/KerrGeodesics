(* ::Package:: *)

(* ::Title:: *)
(*OrbitalFrequencies subpackage of KerrGeodesics*)


(* ::Section::Closed:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`OrbitalFrequencies`",
	{"KerrGeodesics`ConstantsOfMotion`"}];

KerrGeoFrequencies::usage = "KerrGeoFrequencies[a, p, e, x] returns the orbital frequencies."

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Roots of the radial and polar equations*)


(* Returns the roots of the radial equation, as given by Fujita and Hikida *)
KerrGeoRadialRoots[a_, p_, e_, x_, En1_:Null,Q1_:Null] := Module[{M=1,En=En1,Q=Q1,r1,r2,r3,r4,AplusB,AB},
If[En==Null, En = KerrGeoEnergy[a, p, e, x]];
If[Q==Null,  Q = KerrGeoCarterConstant[a, p, e, x]];

r1=p/(1-e);
r2=p/(1+e);
AplusB=(2M)/(1-En^2)-(r1+r2);(*Eq. (11)*)
AB=(a^2 Q)/((1-En^2)r1 r2);(*Eq. (11)*)
r3=(AplusB+Sqrt[(AplusB)^2-4AB])/2;(*Eq. (11)*)
r4=AB/r3;

{r1,r2,r3,r4}

]


KerrGeoRadialRoots[a_, p_, e_/;e==1, x_, En1_:Null,Q1_:Null] := Module[{M=1,En=En1,Q=Q1,r1,r2,r3,r4,\[Rho]2},
r1=\[Infinity];
r2=\[Rho]2=p/(1+e);

r3=1/4 (1-2 \[Rho]2+Sqrt[(8 a^2 (-1+x^2) (-2 a x \[Rho]2+Sqrt[2] Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)])^2)/(\[Rho]2 (a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^2)+(4 \[Rho]2^2 (a^4 (x^2-x^4)+2 (-2+\[Rho]2) \[Rho]2^2+a^2 \[Rho]2 (2+x^2 \[Rho]2)-2 Sqrt[2] a x Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)])^2)/(a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^4]+(-a^4 (-1+x^2) (-1+x^2+2 \[Rho]2)-4 Sqrt[2] a x \[Rho]2 Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)]+\[Rho]2^2 (-4+4 \[Rho]2-5 \[Rho]2^2+2 \[Rho]2^3)-2 a^2 \[Rho]2 (-2+3 \[Rho]2-2 \[Rho]2^2+x^2 (2-5 \[Rho]2+\[Rho]2^2)))/(a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^2);

r4=1/4 (1-2 \[Rho]2-Sqrt[(8 a^2 (-1+x^2) (-2 a x \[Rho]2+Sqrt[2] Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)])^2)/(\[Rho]2 (a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^2)+(4 \[Rho]2^2 (a^4 (x^2-x^4)+2 (-2+\[Rho]2) \[Rho]2^2+a^2 \[Rho]2 (2+x^2 \[Rho]2)-2 Sqrt[2] a x Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)])^2)/(a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^4]+(-a^4 (-1+x^2) (-1+x^2+2 \[Rho]2)-4 Sqrt[2] a x \[Rho]2 Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (-a^2 (-1+x^2)+\[Rho]2^2)]+\[Rho]2^2 (-4+4 \[Rho]2-5 \[Rho]2^2+2 \[Rho]2^3)-2 a^2 \[Rho]2 (-2+3 \[Rho]2-2 \[Rho]2^2+x^2 (2-5 \[Rho]2+\[Rho]2^2)))/(a^2 (-1+x^2)-(-2+\[Rho]2) \[Rho]2)^2);

{r1,r2,r3,r4}

]


(* ::Text:: *)
(*This code uses the polar equation (z^2-zm^2)(a^2(1-E0^2)z^2-zp^2)==0 as the Polar equation. Hence zp is a*Sqrt[1-E0^2]*zp in other sources.*)


KerrGeoPolarRoots[a_, p_, e_, x_] := Module[{En,L,Q,zm,zp},
  {En,L,Q} = {"\[ScriptCapitalE]","\[ScriptCapitalL]","\[ScriptCapitalQ]"}/.KerrGeoConstantsOfMotion[a, p, e, x];
  zm = Sqrt[1-x^2];
  zp = (a^2 (1-En^2)+L^2/(1-zm^2))^(1/2);
  {zp,zm}
]


KerrGeoPolarRoots[a_, p_, e_, x_?PossibleZeroQ] := Module[{En,L,Q,zm,zp},
  {En,L,Q} = Values[KerrGeoConstantsOfMotion[a, p, e, x]];
  zm = Sqrt[1-x^2];
  zp = (Q)^(1/2);
  {zp,zm}
]


(* ::Section::Closed:: *)
(*Orbital Frequencies*)


(* ::Text:: *)
(*Orbital frequency calculations from Fujita and Hikida, Class. Quantum Grav .26 (2009) 135002, arXiv:0906.1420*)


(* ::Subsection::Closed:: *)
(*Schwarzschild*)


sgn[x_]=Piecewise[{{-1,Negative[x]}},1]


KerrGeoMinoFrequencies[a_?PossibleZeroQ, p_, e_?PossibleZeroQ,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> Sqrt[((-6+p) p)/(-3+p)],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> p Sqrt[1/((-3+p))] ,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> (p sgn[x])/Sqrt[(-3+p)],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"  -> Sqrt[p^5/(-3+p)] |>;


KerrGeoMinoFrequencies[a_?PossibleZeroQ, p_,e_,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> (Sqrt[-((p (-6+2 e+p))/(3+e^2-p))] \[Pi])/(2 EllipticK[(4 e)/(-6+2 e+p)]),
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> p/Sqrt[-3-e^2+p],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> (p sgn[x])/(Sqrt[-3-e^2+p]),
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> 1/2 Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))] (8+1/((-4+p)^2 EllipticK[(4 e)/(-6+2 e+p)]) (-(((-4+p) p^2 (-6+2 e+p) EllipticE[(4 e)/(-6+2 e+p)])/(-1+e^2))+(p^2 (28+4 e^2-12 p+p^2) EllipticK[(4 e)/(-6+2 e+p)])/(-1+e^2)-(2 (6+2 e-p) (3+e^2-p) p^2 EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)])/((-1+e) (1+e)^2)+(4 (-4+p) p (2 (1+e) EllipticK[(4 e)/(-6+2 e+p)]+(-6-2 e+p) EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)]))/(1+e)+2 (-4+p)^2 ((-4+p) EllipticK[(4 e)/(-6+2 e+p)]-((6+2 e-p) p EllipticPi[(16 e)/(12+8 e-4 e^2-8 p+p^2),(4 e)/(-6+2 e+p)])/(2+2 e-p)))) |>;


KerrGeoMinoFrequencies[a_?PossibleZeroQ, p_,e_/;e==1,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> (Sqrt[-((p (-6+2 e+p))/(3+e^2-p))] \[Pi])/(2 EllipticK[(4 e)/(-6+2 e+p)]),
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> p/Sqrt[-3-e^2+p],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> (p sgn[x])/(Sqrt[-3-e^2+p]),
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"  -> \[Infinity] |>;


KerrGeoBoyerLindquistFrequencies[a_?PossibleZeroQ, p_, e_?PossibleZeroQ,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> Sqrt[-6+p]/p^2,
    "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> (1)/p^(3/2),
    "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> (p Sign[x])/Sqrt[p^5] |>;


KerrGeoProperFrequencyFactor[_,p_,e_/;e==1,x_]:=\[Infinity]


KerrGeoProperFrequencyFactor[a_?PossibleZeroQ, p_, e_?PossibleZeroQ,x_]:=p^2


KerrGeoProperFrequencyFactor[a_?PossibleZeroQ ,p_,e_,x_]:=(p^2 ((1+e) (28+4 e^2+(-12+p) p)-((1+e) (-4+p) (-6+2 e+p) EllipticE[(4 e)/(-6+2 e+p)]+2 (6+2 e-p) (3+e^2-p) EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)])/EllipticK[(4 e)/(-6+2 e+p)]))/(2 (-1+e) (1+e)^2 (-4+p)^2)


(* ::Subsection::Closed:: *)
(*Kerr*)


(* ::Subsubsection::Closed:: *)
(*KerrGeoMinoFrequencyr*)


KerrGeoMinoFrequencyr[a_,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{kr},
	kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	(\[Pi] Sqrt[(1-En^2)(\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]4)])/(2 EllipticK[kr])
]


KerrGeoMinoFrequencyr[a_,p_,e_/;e==1,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{kr},
	kr=(\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	(\[Pi] Sqrt[2(\[Rho]2-\[Rho]4)])/(2 EllipticK[kr])
]


KerrGeoMinoFrequencyr[a_?PossibleZeroQ,p_,e_?PossibleZeroQ,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{},
	Sqrt[((-6+p) p)/(-3+p)]
]
		


(* ::Subsubsection::Closed:: *)
(*KerrGeoMinoFrequency\[Theta]*)


KerrGeoMinoFrequency\[Theta][a_,p_,e_,x_,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
   (\[Pi]  zp)/(2EllipticK[a^2 (1-En^2)(zm/zp)^2])
]

KerrGeoMinoFrequency\[Theta][a_,p_,e_/;e==1,x_,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
  zp
]



(* ::Subsubsection::Closed:: *)
(*KerrGeoMinoFrequency\[Phi]*)


KerrGeoMinoFrequency\[Phi][a_,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_},{zp_,zm_}]:=KerrGeoMinoFrequency\[Phi]r[a,p,e,x,{En,L,Q},{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4}]+KerrGeoMinoFrequency\[Phi]\[Theta][a,p,e,x,{En,L,Q},{zp,zm}]


KerrGeoMinoFrequency\[Phi]r[a_,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{\[Rho]in,\[Rho]out,kr,hin,hout},
	\[Rho]in= 1-Sqrt[1-a^2];\[Rho]out= 1+Sqrt[1-a^2];
	kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hout= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]out)/(\[Rho]2-\[Rho]out);
	hin= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]in)/(\[Rho]2-\[Rho]in);
	 a/(2Sqrt[1-a^2]) ((2 En \[Rho]out-a L)/(\[Rho]3-\[Rho]out) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]out) EllipticPi[hout,kr]/EllipticK[kr])-(2 En \[Rho]in-a L)/(\[Rho]3-\[Rho]in) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]in) EllipticPi[hin,kr]/EllipticK[kr]))
]

KerrGeoMinoFrequency\[Phi]r[a_/;a^2==1,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{\[Rho]in,\[Rho]out,kr,hM},
	\[Rho]in= 1-Sqrt[1-a^2];\[Rho]out= 1+Sqrt[1-a^2];
	kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hM= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-1)/(\[Rho]2-1);
	a En(2/(\[Rho]3-1) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-1) EllipticPi[hM,kr]/EllipticK[kr])+(2-a L/En)/(2(\[Rho]3-1)^2) ((2-((\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]3))/((\[Rho]1-1)(\[Rho]2-1)))+((\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]4)(\[Rho]3-1))/((\[Rho]1-1)(\[Rho]2-1)(\[Rho]4-1)) EllipticE[kr]/EllipticK[kr]+(\[Rho]2-\[Rho]3)/(\[Rho]2-1) ((\[Rho]1-\[Rho]3)/(\[Rho]1-1)+(\[Rho]2-\[Rho]3)/(\[Rho]2-1)+(\[Rho]4-\[Rho]3)/(\[Rho]4-1)-4) EllipticPi[hM,kr]/EllipticK[kr]))
]



KerrGeoMinoFrequency\[Phi]r[a_,p_,e_/;e==1,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{\[Rho]in,\[Rho]out,kr,hin,hout},
	\[Rho]in= 1-Sqrt[1-a^2];\[Rho]out= 1+Sqrt[1-a^2];
	kr= (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hout= (\[Rho]3-\[Rho]out)/(\[Rho]2-\[Rho]out);
	hin= (\[Rho]3-\[Rho]in)/(\[Rho]2-\[Rho]in);
	  a/(2Sqrt[1-a^2]) ((2  \[Rho]out-a L)/(\[Rho]3-\[Rho]out) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]out) EllipticPi[hout,kr]/EllipticK[kr])-(2  \[Rho]in-a L)/(\[Rho]3-\[Rho]in) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]in) EllipticPi[hin,kr]/EllipticK[kr]))
]


KerrGeoMinoFrequency\[Phi]r[a_/;a^2==1,p_,e_/;e==1,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{kr,hM},
	kr= (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hM= (\[Rho]3-1)/(\[Rho]2-1);
	  a(2/(\[Rho]3-1) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-1) EllipticPi[hM,kr]/EllipticK[kr])+(2-a L)/(2(\[Rho]3-1)^2) ((2-(\[Rho]2-\[Rho]3)/(\[Rho]2-1))+((\[Rho]2-\[Rho]4)(\[Rho]3-1))/((\[Rho]2-1)(\[Rho]4-1)) EllipticE[kr]/EllipticK[kr]+(\[Rho]2-\[Rho]3)/(\[Rho]2-1) (1+(\[Rho]2-\[Rho]3)/(\[Rho]2-1)+(\[Rho]4-\[Rho]3)/(\[Rho]4-1)-4) EllipticPi[hM,kr]/EllipticK[kr]))
]


KerrGeoMinoFrequency\[Phi]\[Theta][a_,p_,e_,x_,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
    ( L EllipticPi[zm^2,a^2 (1-En^2)(zm/zp)^2])/EllipticK[a^2 (1-En^2)(zm/zp)^2]
]


KerrGeoMinoFrequency\[Phi]\[Theta][a_,p_,e_,x_?PossibleZeroQ,{En_,L_,Q_},{zp_,zm_}]:=Module[{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4},
	{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4} = KerrGeoRadialRoots[a,p,e,x,En,Q];
	\[Pi] Sqrt[((\[Rho]1 \[Rho]2 (a^4+\[Rho]1^2 \[Rho]2^2+a^2 ((-2+\[Rho]1) \[Rho]1+(-2+\[Rho]2) \[Rho]2))/2)/(a^4 (2+\[Rho]1+\[Rho]2)+\[Rho]1 \[Rho]2 (\[Rho]1^2 (-2+\[Rho]2)+\[Rho]1 (-2+\[Rho]2) \[Rho]2-2 \[Rho]2^2)+a^2 (\[Rho]1^3+\[Rho]1^2 \[Rho]2+\[Rho]1 (-4+\[Rho]2) \[Rho]2+\[Rho]2^3)))]/EllipticK[(a^2 (a^4+\[Rho]1 (\[Rho]1 (-2+\[Rho]2)-2 \[Rho]2) \[Rho]2+a^2 (\[Rho]1^2+\[Rho]2^2)))/(\[Rho]1 \[Rho]2 (a^4+\[Rho]1^2 \[Rho]2^2+a^2 ((-2+\[Rho]1) \[Rho]1+(-2+\[Rho]2) \[Rho]2)))]
]


KerrGeoMinoFrequency\[Phi]\[Theta][a_,p_,e_/;e==1,x_?PossibleZeroQ,{En_,L_,Q_},{zp_,zm_}]:=Module[{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4},
	{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4} = KerrGeoRadialRoots[a,p,e,x,En,Q];
    Sqrt[2] Sqrt[(\[Rho]2 (a^2+\[Rho]2^2))/(a^2+(-2+\[Rho]2) \[Rho]2)]
]


(* ::Subsubsection::Closed:: *)
(*KerrGeoMinoFrequencyt*)


KerrGeoMinoFrequencyt[a_,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_},{zp_,zm_}]:=KerrGeoMinoFrequencytr[a,p,e,x,{En,L,Q},{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4}]+KerrGeoMinoFrequencyt\[Theta][a,p,e,x,{En,L,Q},{zp,zm}]


KerrGeoMinoFrequencytr[a_,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{\[Rho]in,\[Rho]out,kr,hin,hout,hr},
	\[Rho]in= 1-Sqrt[1-a^2];\[Rho]out= 1+Sqrt[1-a^2];
	kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hout= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]out)/(\[Rho]2-\[Rho]out);
	hin= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]in)/(\[Rho]2-\[Rho]in);
	hr=(\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3);
	(a^2+4)En+En(1/2 (\[Rho]3(\[Rho]1+\[Rho]2+\[Rho]3)-\[Rho]1 \[Rho]2+(\[Rho]1+\[Rho]2+\[Rho]3+\[Rho]4)(\[Rho]2-\[Rho]3) EllipticPi[hr,kr]/EllipticK[kr]+(\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]4) EllipticE[kr]/EllipticK[kr])+2(\[Rho]3+(\[Rho]2-\[Rho]3) EllipticPi[hr,kr]/EllipticK[kr])+1/Sqrt[1-a^2] (((4-a L/En)\[Rho]out-2a^2)/(\[Rho]3-\[Rho]out) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]out) EllipticPi[hout,kr]/EllipticK[kr])-((4-a L/En)\[Rho]in-2a^2)/(\[Rho]3-\[Rho]in) (1-( \[Rho]2-\[Rho]3)/(\[Rho]2-\[Rho]in) EllipticPi[hin,kr]/EllipticK[kr])))
]


KerrGeoMinoFrequencytr[a_,p_,e_/;e==1,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=\[Infinity]




KerrGeoMinoFrequencytr[a_/;a^2==1,p_,e_,x_,{En_,L_,Q_},{\[Rho]1_,\[Rho]2_,\[Rho]3_,\[Rho]4_}]:=Module[{\[Rho]in,\[Rho]out,kr,hM,hr},
	\[Rho]in= 1-Sqrt[1-a^2];\[Rho]out= 1+Sqrt[1-a^2];
	kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4);
	hM= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-1)/(\[Rho]2-1);
	hr=(\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3);
	5En+En(1/2 ((\[Rho]3(\[Rho]1+\[Rho]2+\[Rho]3)-\[Rho]1 \[Rho]2)+(\[Rho]1+\[Rho]2+\[Rho]3+\[Rho]4)(\[Rho]2-\[Rho]3) EllipticPi[hr,kr]/EllipticK[kr]+(\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]4)  EllipticE[kr]/EllipticK[kr])+2(\[Rho]3+(\[Rho]2-\[Rho]3) EllipticPi[hr,kr]/EllipticK[kr])+(2(4-a L/En))/(\[Rho]3-1) (1-(\[Rho]2-\[Rho]3)/(\[Rho]2-1) EllipticPi[hM,kr]/EllipticK[kr])+(2-a L/En)/(\[Rho]3-1)^2 ((2-((\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]3))/((\[Rho]1-1)(\[Rho]2-1)))+((\[Rho]1-\[Rho]3)(\[Rho]2-\[Rho]4)(\[Rho]3-1))/((\[Rho]1-1)(\[Rho]2-1)(\[Rho]4-1)) EllipticE[kr]/EllipticK[kr]+(\[Rho]2-\[Rho]3)/(\[Rho]2-1) ((\[Rho]1-\[Rho]3)/(\[Rho]1-1)+(\[Rho]2-\[Rho]3)/(\[Rho]2-1)+(\[Rho]4-\[Rho]3)/(\[Rho]4-1)-4) EllipticPi[hM,kr]/EllipticK[kr]))
]


KerrGeoMinoFrequencyt\[Theta][a_,p_,e_,x_,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
   (En Q)/((1-En^2) zm^2) (1- EllipticE[a^2 (1-En^2)(zm/zp)^2]/EllipticK[a^2 (1-En^2)(zm/zp)^2])-a^2 En
]


KerrGeoMinoFrequencyt\[Theta][a_,p_,e_/;e==1,x_,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
   -a^2+(a^2 Q)/(2 zp^2)
]


KerrGeoMinoFrequencyt\[Theta][a_,p_,e_,x_/;x^2==1,{En_,L_,Q_},{zp_,zm_}]:=Module[{},
   -a^2En
]






(* ::Subsubsection::Closed:: *)
(*KerrGeoMinoFrequencies*)


KerrGeoMinoFrequencies[a_?Negative,p_,e_,x_]:=<|"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"->1,"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"->1,"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"->-1,"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)"->1|>KerrGeoMinoFrequencies[-a,p,e,-x]


KerrGeoMinoFrequencies[a_,p_,e_?PossibleZeroQ, 1] := Module[{\[CapitalUpsilon]r, \[CapitalUpsilon]\[Phi], \[CapitalUpsilon]\[Theta], \[CapitalGamma]},

\[CapitalUpsilon]r = Sqrt[(p (-2 a^2+6 a Sqrt[p]+(-5+p) p+((a-Sqrt[p])^2 (a^2-4 a Sqrt[p]-(-4+p) p))/Abs[a^2-4 a Sqrt[p]-(-4+p) p]))/(2 a Sqrt[p]+(-3+p) p)];
\[CapitalUpsilon]\[Theta] = Abs[(p^(1/4) Sqrt[3 a^2-4 a Sqrt[p]+p^2])/Sqrt[2 a+(-3+p) Sqrt[p]]];
\[CapitalUpsilon]\[Phi] = p^(5/4)/Sqrt[2 a+(-3+p) Sqrt[p]];
\[CapitalGamma] = (p^(5/4) (a+p^(3/2)))/Sqrt[2 a+(-3+p) Sqrt[p]];

 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Abs[\[CapitalUpsilon]\[Theta]],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> \[CapitalGamma] |>

]


KerrGeoMinoFrequencies[a_,p_,e_,x_]:=Module[{En,L,Q,r1,r2,r3,r4,zm,zp,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t},
{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];

{r1,r2,r3,r4} = KerrGeoRadialRoots[a,p,e,x,En,Q];
{zp,zm}=KerrGeoPolarRoots[a,p,e,x];

\[CapitalUpsilon]r= KerrGeoMinoFrequencyr[a,p,e,x,{En,L,Q},{r1,r2,r3,r4}];
\[CapitalUpsilon]\[Theta]= KerrGeoMinoFrequency\[Theta][a,p,e,x,{En,L,Q},{zp,zm}];
\[CapitalUpsilon]\[Phi]= KerrGeoMinoFrequency\[Phi][a,p,e,x,{En,L,Q},{r1,r2,r3,r4},{zp,zm}];
\[CapitalUpsilon]t= KerrGeoMinoFrequencyt[a,p,e,x,{En,L,Q},{r1,r2,r3,r4},{zp,zm}];

 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Abs[\[CapitalUpsilon]\[Theta]],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(t\)]\)" -> \[CapitalUpsilon]t |>

]


KerrGeoBoyerLindquistFrequencies[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]} = Values[KerrGeoMinoFrequencies[a,p,e,x]];

  <| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> \[CapitalUpsilon]r,
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]\[Theta],
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi] 
   |> / \[CapitalGamma]
]


KerrGeoBoyerLindquistFrequencies[a_,p_,e_/;e>1,x_]:=Module[{},
    <| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> 0,
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> 0,
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> 0 
   |> 
]


KerrGeoProperFrequencyFactor[a_,p_,e_,x_]:=
Module[{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4,zm,zp,T},
	{\[Rho]1,\[Rho]2,\[Rho]3,\[Rho]4}=KerrGeoRadialRoots[a,p,e,x];
	{zp,zm}=KerrGeoPolarRoots[a,p,e,x];
	T=KerrGeoEnergy[a,p,e,x];
	With[{kr= (\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3) (\[Rho]3-\[Rho]4)/(\[Rho]2-\[Rho]4), k\[Theta]=a^2 (1-T^2)(zm/zp)^2, hr=(\[Rho]1-\[Rho]2)/(\[Rho]1-\[Rho]3)},

		1/2 (-((2 zp^2)/(-1+T^2))+\[Rho]1 (-\[Rho]2+\[Rho]3)+\[Rho]3 (\[Rho]2+\[Rho]3))
		+((\[Rho]1-\[Rho]3) (\[Rho]2-\[Rho]4) EllipticE[kr])/(2 EllipticK[kr])
		+(zp^2 EllipticE[k\[Theta]])/((-1+T^2) EllipticK[k\[Theta]])+((\[Rho]2-\[Rho]3) (\[Rho]1+\[Rho]2+\[Rho]3+\[Rho]4) EllipticPi[hr,kr])/(2 EllipticK[kr])
	]
]


KerrGeoProperFrequencyFactor[a_,p_,e_/;e>=1,x_]:=\[Infinity]

KerrGeoProperFrequencies[a_,p_,e_,x_]:=Module[{MinoFreqs,P},
	MinoFreqs = KerrGeoMinoFrequencies[a,p,e,x];
	P=KerrGeoProperFrequencyFactor[a,p,e,x];
	<|"\!\(\*SubscriptBox[\(\[Omega]\), \(r\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"]/P, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Theta]\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"]/P, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Phi]\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"]/P |>
]


KerrGeoProperFrequencies[a_,p_,e_/;e==1,x_]:=Module[{MinoFreqs,P},
	<|"\!\(\*SubscriptBox[\(\[Omega]\), \(r\)]\)"-> 0, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Theta]\)]\)"-> 0, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Phi]\)]\)"-> 0 |>
]


(* ::Subsection::Closed:: *)
(*Generic function for choosing between frequencies w.r.t different time coordinates*)


Options[KerrGeoFrequencies] = {"Time" -> "BoyerLindquist"}
SyntaxInformation[KerrGeoFrequencies] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
KerrGeoFrequencies[a_,p_,e_,x_,OptionsPattern[]] := Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},


If[OptionValue["Time"]=="Mino",Return[KerrGeoMinoFrequencies[a,p,e,x][[1;;4]]]];

If[OptionValue["Time"]=="BoyerLindquist", Return[KerrGeoBoyerLindquistFrequencies[a,p,e,x]]];

If[OptionValue["Time"]=="Proper",Return[KerrGeoProperFrequencies[a,p,e,x][[1;;3]]]];

]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
