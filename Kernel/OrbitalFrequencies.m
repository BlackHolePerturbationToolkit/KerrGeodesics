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


(* ::Text:: *)
(*This code uses the polar equation (z^2-zm^2)(a^2(1-E0^2)z^2-zp^2)==0 as the Polar equation. Hence zp is a*Sqrt[1-E0^2]*zp in other sources.*)


KerrGeoPolarRoots[a_, p_, e_, x_] := Module[{En,L,Q,zm,zp},
  {En,L,Q} = Values[KerrGeoConstantsOfMotion[a, p, e, x]];
  zm = Sqrt[1-x^2];
  zp = (a^2 (1-En^2)+L^2/(1-zm^2))^(1/2);
  {zp,zm}
]


(* ::Section::Closed:: *)
(*Orbital Frequencies*)


(* ::Text:: *)
(*Orbital frequency calculations from Fujita and Hikida, Class. Quantum Grav .26 (2009) 135002, arXiv:0906.1420*)


(* ::Subsection::Closed:: *)
(*Schwarzschild*)


KerrGeoMinoFrequencies[0|0., p_,0,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> Sqrt[((-6+p) p)/(-3+p)],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> p Sqrt[1/((-3+p) x^2)] x,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> (p x)/Sqrt[(-3+p) x^2],
    "\[CapitalGamma]" -> Sqrt[p^5/(-3+p)] |>;


KerrGeoMinoFrequencies[0|0., p_,e_,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> (Sqrt[-((p (-6+2 e+p))/(3+e^2-p))] \[Pi])/(2 EllipticK[(4 e)/(-6+2 e+p)]),
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> p/Sqrt[-3-e^2+p],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> (p x)/(Sqrt[-3-e^2+p] Abs[x]),
    "\[CapitalGamma]" -> 1/2 Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))] (8+1/((-4+p)^2 EllipticK[(4 e)/(-6+2 e+p)]) (-(((-4+p) p^2 (-6+2 e+p) EllipticE[(4 e)/(-6+2 e+p)])/(-1+e^2))+(p^2 (28+4 e^2-12 p+p^2) EllipticK[(4 e)/(-6+2 e+p)])/(-1+e^2)-(2 (6+2 e-p) (3+e^2-p) p^2 EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)])/((-1+e) (1+e)^2)+(4 (-4+p) p (2 (1+e) EllipticK[(4 e)/(-6+2 e+p)]+(-6-2 e+p) EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)]))/(1+e)+2 (-4+p)^2 ((-4+p) EllipticK[(4 e)/(-6+2 e+p)]-((6+2 e-p) p EllipticPi[(16 e)/(12+8 e-4 e^2-8 p+p^2),(4 e)/(-6+2 e+p)])/(2+2 e-p)))) |>;


KerrGeoBoyerLindquistFrequencies[0|0., p_,0,x_] :=
 <| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> Sqrt[-6+p]/p^2,
    "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> (Sqrt[1/x^2] x)/p^(3/2),
    "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> (p x)/Sqrt[p^5 x^2] |>;


KerrGeoProperFrequencyFactor[0|0., p_,0,x_]:=p^2


KerrGeoProperFrequencyFactor[0|0. ,p_,e_,x_]:=(p^2 ((1+e) (28+4 e^2+(-12+p) p)-((1+e) (-4+p) (-6+2 e+p) EllipticE[(4 e)/(-6+2 e+p)]+2 (6+2 e-p) (3+e^2-p) EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)])/EllipticK[(4 e)/(-6+2 e+p)]))/(2 (-1+e) (1+e)^2 (-4+p)^2)


(* ::Subsection::Closed:: *)
(*Kerr*)


KerrGeoMinoFrequencies[a_,p_,(0|0.),1] := Module[{\[CapitalUpsilon]r, \[CapitalUpsilon]\[Phi], \[CapitalUpsilon]\[Theta], \[CapitalGamma]},

\[CapitalUpsilon]r = Sqrt[(p (-2 a^2+6 a Sqrt[p]+(-5+p) p+((a-Sqrt[p])^2 (a^2-4 a Sqrt[p]-(-4+p) p))/Abs[a^2-4 a Sqrt[p]-(-4+p) p]))/(2 a Sqrt[p]+(-3+p) p)];
\[CapitalUpsilon]\[Theta] = Abs[(p^(1/4) Sqrt[3 a^2-4 a Sqrt[p]+p^2])/Sqrt[2 a+(-3+p) Sqrt[p]]];
\[CapitalUpsilon]\[Phi] = p^(5/4)/Sqrt[2 a+(-3+p) Sqrt[p]];
\[CapitalGamma] = (p^(5/4) (a+p^(3/2)))/Sqrt[2 a+(-3+p) Sqrt[p]];

 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Abs[\[CapitalUpsilon]\[Theta]],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
    "\[CapitalGamma]" -> \[CapitalGamma] |>

]


KerrGeoMinoFrequencies[a_,p_,e_,x_]:=Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},
{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];

{r1,r2,r3,r4} = KerrGeoRadialRoots[a,p,e,x,En,Q];
\[Epsilon]0=a^2 (1-En^2)/L^2;
zm=1-x^2;
a2zp=(L^2+a^2 (-1+En^2) (-1+zm))/( (-1+En^2) (-1+zm));

\[Epsilon]0zp=-((L^2+a^2 (-1+En^2) (-1+zm))/(L^2 (-1+zm)));

(*zmOverZp=If[a==0,0,zm/((L^2+a^2 (-1+En^2) (-1+zm))/(a^2 (-1+En^2) (-1+zm)))];*)
zmOverZp=zm/((L^2+a^2 (-1+En^2) (-1+zm))/(a^2 (-1+En^2) (-1+zm)));


kr=Sqrt[(r1-r2)/(r1-r3) (r3-r4)/(r2-r4)];(*Eq.(13)*)
k\[Theta]=Sqrt[zmOverZp];(*Eq.(13)*)
\[CapitalUpsilon]r=(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)])/(2EllipticK[kr^2]);(*Eq.(15)*)
\[CapitalUpsilon]\[Theta]=(Pi L Sqrt[\[Epsilon]0zp])/(2EllipticK[k\[Theta]^2]);(*Eq.(15)*)

rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));

(*Eq. (21)*)
\[CapitalUpsilon]\[Phi]=(2\[CapitalUpsilon]\[Theta])/(Pi Sqrt[\[Epsilon]0zp]) EllipticPi[zm,k\[Theta]^2]+(2a \[CapitalUpsilon]r)/(Pi(rp-rm)Sqrt[(1-En^2)(r1-r3)(r2-r4)])((2M En rp - a L)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-(2M En rm - a L)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2]));

\[CapitalGamma]=4M^2 En + (2a2zp En  \[CapitalUpsilon]\[Theta])/(Pi L Sqrt[\[Epsilon]0zp]) (EllipticK[k\[Theta]^2]- EllipticE[k\[Theta]^2]) + (2\[CapitalUpsilon]r)/(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)]) (En/2 ((r3(r1+r2+r3)-r1 r2)EllipticK[kr^2]+(r2-r3)(r1+r2+r3+r4)EllipticPi[hr,kr^2]+(r1-r3)(r2-r4)EllipticE[kr^2])+2M En(r3 EllipticK[kr^2]+(r2-r3)EllipticPi[hr,kr^2])+(2M)/(rp-rm) (((4M^2 En-a L)rp-2M a^2 En)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-((4M^2 En-a L)rm-2M a^2 En)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2])));

 <| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Abs[\[CapitalUpsilon]\[Theta]],
    "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
    "\[CapitalGamma]" -> \[CapitalGamma] |>

]


KerrGeoMinoFrequencies[(1|1.),p_,e_,x_]:=Module[{M=1,a=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hM,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];

{r1,r2,r3,r4} = KerrGeoRadialRoots[a,p,e,x,En,Q];
\[Epsilon]0=a^2 (1-En^2)/L^2;
zm=1-x^2;
a2zp=(L^2+a^2 (-1+En^2) (-1+zm))/( (-1+En^2) (-1+zm));

\[Epsilon]0zp=-((L^2+a^2 (-1+En^2) (-1+zm))/(L^2 (-1+zm)));

(*zmOverZp=If[a==0,0,zm/((L^2+a^2 (-1+En^2) (-1+zm))/(a^2 (-1+En^2) (-1+zm)))];*)
zmOverZp=zm/((L^2+a^2 (-1+En^2) (-1+zm))/(a^2 (-1+En^2) (-1+zm)));


kr=Sqrt[(r1-r2)/(r1-r3) (r3-r4)/(r2-r4)];(*Eq.(13)*)
k\[Theta]=Sqrt[zmOverZp];(*Eq.(13)*)
\[CapitalUpsilon]r=(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)])/(2EllipticK[kr^2]);(*Eq.(15)*)
\[CapitalUpsilon]\[Theta]=(Pi L Sqrt[\[Epsilon]0zp])/(2EllipticK[k\[Theta]^2]);(*Eq.(15)*)

hM = ((r1-r2)(r3-M))/((r1-r3)(r2-M));

hr=(r1-r2)/(r1-r3);

(*\[CapitalUpsilon]\[Phi] and \[CapitalGamma] from Appendix B for a=M case*)

\[CapitalUpsilon]\[Phi]= (2\[CapitalUpsilon]\[Theta])/(\[Pi] Sqrt[\[Epsilon]0zp]) EllipticPi[zm, k\[Theta]^2]+ (2 a \[CapitalUpsilon]r)/(\[Pi] Sqrt[(1-En^2)(r1-r3)(r2-r4)]) ((2M En)/(r3-M) (EllipticK[kr^2] - (r2-r3)/(r2-M) EllipticPi[hM,kr^2])+(2M^2 En-a L)/(2(r3-M)^2) ((2-((r1-r3)(r2-r3))/((r1-M)(r2-M)))EllipticK[kr^2] + ((r1-r3)(r2-r4)(r3-M))/((r1-M)(r2-M)(r4-M)) EllipticE[kr^2]+(r2-r3)/(r2-M) ((r1-r3)/(r1-M)+(r2-r3)/(r2-M)+(r4-r3)/(r4-M)-4)EllipticPi[hM, kr^2]));

\[CapitalGamma]=4M^2 En+(2a^2 En a2zp \[CapitalUpsilon]\[Theta])/(\[Pi] L Sqrt[\[Epsilon]0zp]) (EllipticK[k\[Theta]^2]-EllipticE[k\[Theta]^2]) + (2 \[CapitalUpsilon]r)/(\[Pi] Sqrt[(1-En^2)(r1-r3)(r2-r4)]) (En/2 ((r3(r1+r2+r3)-r1 r2)EllipticK[kr^2]+(r2-r3)(r1+r2+r3+r4)EllipticPi[hr,kr^2]+(r1-r3)(r2-r4)EllipticE[kr^2])+2M En(r3 EllipticK[kr^2]+(r2-r3)EllipticPi[hr,kr^2]) + (2M(4M^2 En - a L))/(r3-M) (EllipticK[kr^2]-(r2-r3)/(r2-M) EllipticPi[hM,kr^2])+(M^2 (2M^2 En-a L ))/(r3-M)^2 ((2-((r1-r3)(r2-r3))/((r1-M)(r2-M)))EllipticK[kr^2]+((r1-r3)(r2-r4)(r3-M))/((r1-M)(r2-M)(r4-M)) EllipticE[kr^2]+(r2-r3)/(r2-M) ((r1-r3)/(r1-M)+(r2-r3)/(r2-M)+(r4-r3)/(r4-M)-4)EllipticPi[hM,kr^2]));

<| "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)" -> \[CapitalUpsilon]r,
   "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> Abs[\[CapitalUpsilon]\[Theta]],
   "\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi],
   "\[CapitalGamma]" -> \[CapitalGamma] |>
]


KerrGeoBoyerLindquistFrequencies[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]} = Values[KerrGeoMinoFrequencies[a,p,e,x]];

  <| "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> \[CapitalUpsilon]r,
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> \[CapitalUpsilon]\[Theta],
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> \[CapitalUpsilon]\[Phi] 
   |> / \[CapitalGamma]
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



KerrGeoProperFrequencies[a_,p_,e_,x_]:=Module[{MinoFreqs,P},
	MinoFreqs = KerrGeoMinoFrequencies[a,p,e,x];
	P=KerrGeoProperFrequencyFactor[a,p,e,x];
	<|"\!\(\*SubscriptBox[\(\[Omega]\), \(r\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(r\)]\)"]/P, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Theta]\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)"]/P, "\!\(\*SubscriptBox[\(\[Omega]\), \(\[Phi]\)]\)"-> MinoFreqs["\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Phi]\)]\)"]/P |>
]


(* ::Subsection::Closed:: *)
(*Generic function for choosing between frequencies w.r.t different time coordinates*)


Options[KerrGeoFrequencies] = {"Time" -> "BoyerLindquist"}
SyntaxInformation[KerrGeoFrequencies] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
KerrGeoFrequencies[a_,p_,e_,x_,OptionsPattern[]] := Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},


If[OptionValue["Time"]=="Mino",Return[KerrGeoMinoFrequencies[a,p,e,x][[1;;3]]]];

If[OptionValue["Time"]=="BoyerLindquist", Return[KerrGeoBoyerLindquistFrequencies[a,p,e,x]]];

If[OptionValue["Time"]=="Proper",Return[KerrGeoProperFrequencies[a,p,e,x][[1;;3]]]];

]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
