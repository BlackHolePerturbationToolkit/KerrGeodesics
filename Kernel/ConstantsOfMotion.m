(* ::Package:: *)

(* ::Title:: *)
(*ConstantsOfMotion subpackage of KerrGeodesics*)


(* ::Section::Closed:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`ConstantsOfMotion`"];

KerrGeoEnergy::usage = "KerrGeoEnergy[a, p, e, x] returns the orbital energy."
KerrGeoAngularMomentum::usage = "KerrGeoAngularMomentum[a, p, e, x] returns the orbital angular momentum about the symmetry axis."
KerrGeoCarterConstant::usage = "KerrGeoCarterConstant[a, p, e, x] returns the Carter constant of the orbit."
KerrGeoConstantsOfMotion::usage = "KerrGeoConstantsOfMotion[a, p, e, x] returns the three constants of motion."

(*KerrGeoVelocityAtInfinity::usage = "KerrGeoVelocityAtInfinity[a, p, e, x] returns the magnitude of the velocity at infinity of a scatter orbit."
KerrGeoImpactParameter::usage = "KerrGeoImpactParameter[a, p, e, x] returns the impact parameter of a scatter orbit."*)

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Schwarzschild (a=0)*)


(* ::Subsection::Closed:: *)
(*Circular (e=0)*)


KerrGeoEnergy[0,p_,0,x_]:=(-2+p)/Sqrt[(-3+p) p]


KerrGeoAngularMomentum[0,p_,0,x_]:=(p x)/Sqrt[-3+p]


KerrGeoCarterConstant[0,p_,0,x_]:=-((p^2 (-1+x^2))/(-3+p))


(* ::Subsection::Closed:: *)
(*Eccentric*)


KerrGeoEnergy[0,p_,e_,x_]:=Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))]


KerrGeoAngularMomentum[0,p_,e_,x_]:=(p x)/Sqrt[-3-e^2+p]


KerrGeoCarterConstant[0,p_,e_,x_]:=(p^2 (-1+x^2))/(3+e^2-p)


(* ::Subsection::Closed:: *)
(*Scatter*)


KerrGeoVelocityAtInfinity[0,p_,e_/;e>1,x_]:=Module[{En},
	En = KerrGeoEnergy[0,p,e,x];
	Sqrt[En^2-1]/En
]


KerrGeoImpactParameter[0,p_,e_/;e>1,x_]:=KerrGeoAngularMomentum[0,p,e,x]/Sqrt[KerrGeoEnergy[0,p,e,x]^2-1]


(* ::Text:: *)
(*Schwarzschild hyperbolic scatter angle*)
(*Defined as \[Phi](SuperPlus[\[ScriptCapitalI]])-\[Phi](SuperMinus[\[ScriptCapitalI]]) - \[Pi] e.g. Eq. (29) of arXiv:2209.03740 *)


KerrGeoScatteringAngle[0,p_,e_/;e>=1,1]:= -\[Pi]+(4 Sqrt[p]EllipticF[ArcCos[-e^(-1)]/2,(4e)/(6+2e-p)])/Sqrt[-6-2e+p]


(* ::Subsection::Closed:: *)
(*Convenience function to compute all three constants of motion*)


KerrGeoConstantsOfMotion[0,p_,e_,x_]:= 
 <|"\[ScriptCapitalE]" -> KerrGeoEnergy[0,p,e,x],
   "\[ScriptCapitalL]" -> KerrGeoAngularMomentum[0,p,e,x],
   "\[ScriptCapitalQ]" -> KerrGeoCarterConstant[0,p,e,x] |>


KerrGeoConstantsOfMotion[0,p_,e_/;e>=1,x_]:= 
 <|"\[ScriptCapitalE]" -> KerrGeoEnergy[0,p,e,x],
   "\[ScriptCapitalL]" -> KerrGeoAngularMomentum[0,p,e,x],
   "\[ScriptCapitalQ]" -> KerrGeoCarterConstant[0,p,e,x], 
   "\!\(\*SubscriptBox[\(v\), \(\[Infinity]\)]\)" -> KerrGeoVelocityAtInfinity[0,p,e,x],
   "b" -> KerrGeoImpactParameter[0,p,e,x],
   "\[Psi]" -> KerrGeoScatteringAngle[0,p,e,x]
   |>


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Negative a*)


KerrGeoEnergy[a_?Negative,p_,e_,x_]:=KerrGeoEnergy[-a,p,e,-x]
KerrGeoAngularMomentum[a_?Negative,p_,e_,x_]:=-KerrGeoAngularMomentum[-a,p,e,-x]


(* ::Subsection::Closed:: *)
(*Equatorial orbits (x^2 = 1)*)


(* ::Text:: *)
(*The Carter constant is zero for all equatorial orbits*)


KerrGeoCarterConstant[a_,p_,e_,x_/;x^2==1]:=0


(* ::Subsubsection::Closed:: *)
(*Circular (e=0)*)


KerrGeoEnergy[a_,p_,0,x_/;x^2==1]:=((-2+p) Sqrt[p]+a/x)/Sqrt[2 a/x p^(3/2)+(-3+p) p^2]


KerrGeoAngularMomentum[a_,p_,0,x_/;x^2==1]:=((a^2+p^2)x-2 a Sqrt[p]) /(p^(3/4) Sqrt[x^2 (-3+p) Sqrt[p]+2 a x])


(* ::Subsubsection::Closed:: *)
(*Eccentric*)


(* ::Text:: *)
(*Simplified from Glampedakis and Kennefick, Phys. Rev. D66 (2002) 044002, arXiv:gr-qc/0203086, Eq. 7 and appendix A*)


KerrGeoEnergy[a_,p_,e_,x_/;x^2==1]:= Sqrt[1-((1-e^2) (1+((-1+e^2) (a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(p^3 x^2)])))/(-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)))/p];


KerrGeoAngularMomentum[a_,p_,e_,x_/;x^2==1]:= p x Sqrt[(a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(x^2 p^3)]))/((-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)x^2)]+a Sqrt[1-((1-e^2) (1+((-1+e^2) (a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(p^3 x^2)])))/(-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)))/p];


(* ::Subsubsection::Closed:: *)
(*Convenience function to compute all three constants of motion*)


KerrGeoConstantsOfMotion[a_,p_,e_,x:(1|-1)]:=
 <|"\[ScriptCapitalE]" -> KerrGeoEnergy[a,p,e,x],
   "\[ScriptCapitalL]" -> KerrGeoAngularMomentum[a,p,e,x],
   "\[ScriptCapitalQ]" -> KerrGeoCarterConstant[a,p,e,x] |>


(* ::Subsection::Closed:: *)
(*Polar orbits (x=0)*)


(* ::Text:: *)
(*The angular momentum is zero for all polar orbits*)


KerrGeoAngularMomentum[a_,p_,e_,(0|0.)]:=0


(* ::Subsubsection::Closed:: *)
(*Spherical (e=0)*)


(* ::Text:: *)
(*Simplified formula starting from Stoghianidis & Tsoubelis, Gen. Rel, Grav., vol. 19, No. 12, p. 1235 (1987), Eqs. (17)-(19)*)


KerrGeoEnergy[a_,p_,(0|0.),(0|0.)]:=Sqrt[(p (a^2-2 p+p^2)^2)/((a^2+p^2) (a^2+a^2 p-3 p^2+p^3))]


KerrGeoCarterConstant[a_,p_,(0|0.),(0|0.)]:=(p^2 (a^4+2 a^2 (-2+p) p+p^4))/((a^2+p^2) ((-3+p) p^2+a^2 (1+p)))


(* ::Subsubsection::Closed:: *)
(*Eccentric*)


(* ::Text:: *)
(*These equations were worked out by N. Warburton starting with Schmidt's formula*)


KerrGeoEnergy[a_,p_,e_,(0|0.)]:=Sqrt[-((p (a^4 (-1+e^2)^2+(-4 e^2+(-2+p)^2) p^2+2 a^2 p (-2+p+e^2 (2+p))))/(a^4 (-1+e^2)^2 (-1+e^2-p)+(3+e^2-p) p^4-2 a^2 p^2 (-1-e^4+p+e^2 (2+p))))]


KerrGeoCarterConstant[a_,p_,e_,(0|0.)]:= -((p^2 (a^4 (-1+e^2)^2+p^4+2 a^2 p (-2+p+e^2 (2+p))))/(a^4 (-1+e^2)^2 (-1+e^2-p)+(3+e^2-p) p^4-2 a^2 p^2 (-1-e^4+p+e^2 (2+p))))


(* ::Subsubsection::Closed:: *)
(*Convenience function to compute all three constants of motion*)


KerrGeoConstantsOfMotion[a_,p_,e_,(0|0.)] :=
 <|"\[ScriptCapitalE]" -> KerrGeoEnergy[a,p,e,0],
   "\[ScriptCapitalL]" -> KerrGeoAngularMomentum[a,p,e,0],
   "\[ScriptCapitalQ]" -> KerrGeoCarterConstant[a,p,e,0] |>


(* ::Subsection::Closed:: *)
(*Spherical orbits (e=0)*)


KerrGeoEnergy[a_,p_,(0|0.),x_]:=\[Sqrt](((-3+p) (-2+p)^2 p^5-2 a^5 x (-1+x^2) Sqrt[p^3+a^2 p (-1+x^2)]+a^4 p^2 (-1+x^2) (4-5 p (-1+x^2)+3 p^2 (-1+x^2))-a^6 (-1+x^2)^2 (x^2+p^2 (-1+x^2)-p (1+2 x^2))+a^2 p^3 (4-4 x^2+p (12-7 x^2)-3 p^3 (-1+x^2)+p^2 (-13+10 x^2))+a (-2 p^(9/2) x Sqrt[p^2+a^2 (-1+x^2)]+4 p^3 x Sqrt[p^3+a^2 p (-1+x^2)])+2 a^3 (2 p x (-1+x^2) Sqrt[p^3+a^2 p (-1+x^2)]-x^3 Sqrt[p^7+a^2 p^5 (-1+x^2)]))/((p^2-a^2 (-1+x^2)) ((-3+p)^2 p^4-2 a^2 p^2 (3+2 p-3 x^2+p^2 (-1+x^2))+a^4 (-1+x^2) (-1+x^2+p^2 (-1+x^2)-2 p (1+x^2)))))


KerrGeoAngularMomentum[a_,p_,(0|0.),x_,En1_:Null]:=Block[{En=En1,g,d,h,f},
If[En==Null,En=KerrGeoEnergy[a,p,0,x]];

g=2 a p;
d=(a^2+(-2+p) p) (p^2-a^2 (-1+x^2));
h=((-2+p) p-a^2 (-1+x^2))/x^2;
f=p^4+a^2 (p (2+p)-(a^2+(-2+p) p) (-1+x^2));

(-En g + x Sqrt[(-d h + En^2 (g^2+ f h))/x^2])/h

]


(* ::Text:: *)
(*CarterConstant and ConstantsOfMotion calculations are covered by the generic case*)


(* ::Subsection::Closed:: *)
(*Marginally bound orbits (e = 1)*)


KerrGeoEnergy[a_,p_,e_/;e==1,x_]:=1


KerrGeoAngularMomentum[a_,p_,e_/;e==1,x_,En1_:Null]:= Module[{En=En1,\[Rho]2},
	If[En==Null,En=KerrGeoEnergy[a,p,e,x]];
	\[Rho]2=p/(1+e);
	((x^2) (-2 a \[Rho]2+Sqrt[2] /x Sqrt[\[Rho]2 (a^2+(-2+\[Rho]2) \[Rho]2) (a^2 (1-x^2)+\[Rho]2^2)]) )/(a^2 (1-x^2)+(-2+\[Rho]2) \[Rho]2)
]



(* ::Subsection::Closed:: *)
(*Generic orbits*)


KerrGeoEnergy[a_,p_,e_,x_]:= Module[{r1,r2,zm,\[CapitalDelta],f,g,h,d,\[Kappa],\[Rho],\[Epsilon],\[Sigma],\[Eta],r},
	
	r1 = p/(1-e);
	r2 = p/(1+e);

	zm = Sqrt[1-x^2];

    \[CapitalDelta][r_] = r^2 - 2 r + a^2;

    f[r_] = r^4 + a^2 (r (r + 2) + zm^2 \[CapitalDelta][r]);
    g[r_] = 2 a r;
    h[r_] = r (r - 2) + zm^2/(1 - zm^2) \[CapitalDelta][r];
    d[r_] = (r^2 + a^2 zm^2) \[CapitalDelta][r];
    
    \[Kappa] = d[r1] h[r2] - h[r1] d[r2];
    \[Epsilon] = d[r1] g[r2] - g[r1] d[r2];
    \[Rho] = f[r1] h[r2] - h[r1] f[r2];
    \[Eta] = f[r1] g[r2] - g[r1] f[r2];
    \[Sigma] = g[r1] h[r2] - h[r1] g[r2];

	Sqrt[(\[Kappa] \[Rho] + 2 \[Epsilon] \[Sigma] - x 2 Sqrt[\[Sigma] (\[Sigma] \[Epsilon]^2 + \[Rho] \[Epsilon] \[Kappa] - \[Eta] \[Kappa]^2)/x^2])/(\[Rho]^2 + 4 \[Eta] \[Sigma])]
]


KerrGeoAngularMomentum[a_,p_,e_,x_,En1_:Null]:= Module[{En=En1,r1,zm,\[CapitalDelta],f,g,h,d,r},
	If[En==Null,En=KerrGeoEnergy[a,p,e,x]];
	
	r1 = p/(1-e);

	zm = Sqrt[1-x^2];

    \[CapitalDelta][r_] = r^2 - 2 r + a^2;

    f[r_] = r^4 + a^2 (r (r + 2) + zm^2 \[CapitalDelta][r]);
    g[r_] = 2 a r;
    h[r_] = r (r - 2) + zm^2/(1 - zm^2) \[CapitalDelta][r];
    d[r_] = (r^2 + a^2 zm^2) \[CapitalDelta][r];
    
    (-En g[r1] + x Sqrt[(-d[r1] h[r1] + En^2 (g[r1]^2+ f[r1] h[r1]))/x^2])/h[r1]
]


KerrGeoCarterConstant[a_,p_,e_,x_,En1_:Null,L1_:Null]:= Module[{En=En1,L=L1,zm},
	If[En==Null,En=KerrGeoEnergy[a,p,e,x]];
	If[L==Null,L= KerrGeoAngularMomentum[a,p,e,x,En]];
    zm = Sqrt[1-x^2];
	
	zm^2 (a^2 (1 - En^2) + L^2/(1 - zm^2))
]


KerrGeoConstantsOfMotion[a_,p_,e_,x_] :=Module[{En,L},
  En = KerrGeoEnergy[a,p,e,x];
  L = KerrGeoAngularMomentum[a,p,e,x];

 <|"\[ScriptCapitalE]" -> En, "\[ScriptCapitalL]" -> L, "\[ScriptCapitalQ]" -> KerrGeoCarterConstant[a,p,e,x,En,L] |>
 ]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
