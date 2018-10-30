(* ::Package:: *)

(* ::Title:: *)
(*Package for the calculation of bound time-like geodesics and their properties in Kerr spacetime*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`"];

KerrGeoEnergy::usage = "KerrGeoEnergy[a, p, e, x] returns the orbital energy."
KerrGeoAngularMomentum::usage = "KerrGeoAngularMomentum[a, p, e, x] returns the orbital angular momentum about the symmetry axis."
KerrGeoCarterConstant::usage = "KerrGeoCarterConstant[a, p, e, x] returns the Carter constant of the orbit."
KerrGeoConstantsOfMotion::usage = "KerrGeoConstantsOfMotion[a, p, e, x] returns the three constants of motion {E,L,Q}."

KerrGeoFrequencies::usage = "KerrGeoFrequencies[a, p, e, x] returns the orbital frequencies"

KerrGeoOrbit::usage = "KerrGeoOrbit[a,p,e,x] returns a KerrGeoOrbitFunction[..] which stores the orbital trajectory and parameters.";
KerrGeoOrbitFunction::usage = "KerrGeoOrbitFunction[a,p,e,x,assoc] an object for storing the trajectory and orbital parameters in the assoc Association."

KerrGeoPhotonSphereRadius::usage = "KerrGeoPhotonSphereRadius[a,x] returns the radius of the photon sphere."

KerrGeoISCO::usage = "KerrGeoISCO[a,x] returns the location of the innermost stable circular orbit (ISCO) for pro- and retrograde orbits"
KerrGeoISSO::usage = "KerrGeoISCO[a,x] returns the location of the innermost stable spherical orbit (ISSO)"
KerrGeoIBSO::usage = "KerrGeoISCO[a,x] returns the location of the innermost bound spherical orbit (IBSO)"

KerrGeoSeparatrix::usage = "KerrGeoSeparatrix[a,e,x] returns the value of \!\(\*SubscriptBox[\(p\), \(s\)]\) at the separatrix"

Begin["`Private`"];


(* ::Chapter::Closed:: *)
(*Constants of Motion*)


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
(*Convenience function to compute all three constants of motion*)


KerrGeoConstantsOfMotion[0,p_,e_,x_]:= {KerrGeoEnergy[0,p,e,x],KerrGeoAngularMomentum[0,p,e,x],KerrGeoCarterConstant[0,p,e,x]}


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Equatorial orbits (x^2 = 1)*)


(* ::Text:: *)
(*The Carter constant is zero for all equatorial orbits*)


KerrGeoCarterConstant[a_,p_,e_,x_/;x^2==1]:=0


(* ::Subsubsection:: *)
(*Circular (e=0)*)


KerrGeoEnergy[a_,p_,0,x_/;x^2==1]:=((-2+p) Sqrt[p]+a/x)/Sqrt[2 a/x p^(3/2)+(-3+p) p^2]


KerrGeoAngularMomentum[a_,p_,0,x_/;x^2==1]:=(a^2-2 a/x Sqrt[p]+p^2)/(Sqrt[2 a/x+(-3+p) Sqrt[p]] p^(3/4))


(* ::Subsubsection:: *)
(*Eccentric*)


(* ::Text:: *)
(*Simplified from Glampedakis and Kennefick, Phys. Rev. D66 (2002) 044002, arXiv:gr-qc/0203086, Eq. 7 and appendix A*)


KerrGeoEnergy[a_,p_,e_,x_/;x^2==1]:= Sqrt[1-((1-e^2) (1+((-1+e^2) (a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(p^3 x^2)])))/(-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)))/p];


KerrGeoAngularMomentum[a_,p_,e_,x_/;x^2==1]:= p x Sqrt[(a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(x^2 p^3)]))/((-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)x^2)]+a Sqrt[1-((1-e^2) (1+((-1+e^2) (a^2 (1+3 e^2+p)+p (-3-e^2+p-2x Sqrt[(a^6 (-1+e^2)^2+a^2 (-4 e^2+(-2+p)^2) p^2+2 a^4 p (-2+p+e^2 (2+p)))/(p^3 x^2)])))/(-4 a^2 (-1+e^2)^2+(3+e^2-p)^2 p)))/p];


(* ::Subsubsection:: *)
(*Convenience function to compute all three constants of motion*)


KerrGeoConstantsOfMotion[a_,p_,e_,x_/;x^2==1]:= {KerrGeoEnergy[a,p,e,x],KerrGeoAngularMomentum[a,p,e,x],KerrGeoCarterConstant[a,p,e,x]}


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


KerrGeoConstantsOfMotion[a_,p_,e_,(0|0.)]:={KerrGeoEnergy[a,p,e,0],KerrGeoAngularMomentum[a,p,e,0],KerrGeoCarterConstant[a,p,e,0]}


(* ::Subsection::Closed:: *)
(*Spherical orbits (e=0)*)


KerrGeoEnergy[a_,p_,0,x_]:=\[Sqrt](((-3+p) (-2+p)^2 p^5-2 a^5 x (-1+x^2) Sqrt[p^3+a^2 p (-1+x^2)]+a^4 p^2 (-1+x^2) (4-5 p (-1+x^2)+3 p^2 (-1+x^2))-a^6 (-1+x^2)^2 (x^2+p^2 (-1+x^2)-p (1+2 x^2))+a^2 p^3 (4-4 x^2+p (12-7 x^2)-3 p^3 (-1+x^2)+p^2 (-13+10 x^2))+a (-2 p^(9/2) x Sqrt[p^2+a^2 (-1+x^2)]+4 p^3 x Sqrt[p^3+a^2 p (-1+x^2)])+2 a^3 (2 p x (-1+x^2) Sqrt[p^3+a^2 p (-1+x^2)]-x^3 Sqrt[p^7+a^2 p^5 (-1+x^2)]))/((p^2-a^2 (-1+x^2)) ((-3+p)^2 p^4-2 a^2 p^2 (3+2 p-3 x^2+p^2 (-1+x^2))+a^4 (-1+x^2) (-1+x^2+p^2 (-1+x^2)-2 p (1+x^2)))))


KerrGeoAngularMomentum[a_,p_,0,x_,En1_:Null]:=Block[{En=En1,g,d,h,f},
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


KerrGeoConstantsOfMotion[a_,p_,e_,x_]:=Module[{En,L,Q},
	En=KerrGeoEnergy[a,p,e,x];
	L=KerrGeoAngularMomentum[a,p,e,x,En];
	Q=KerrGeoCarterConstant[a,p,e,x,En,L];
	{En,L,Q}
]


(* ::Chapter::Closed:: *)
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


KerrGeoPolarRoots[a_, p_, e_, x_] := Module[{En,L,Q,zm,zp},
  {En,L,Q} = KerrGeoConstantsOfMotion[a, p, e, x];
  zm = Sqrt[1-x^2];
  zp = (a^2 (1-En^2)+L^2/(1-zm^2))^(1/2);
  {zp,zm}
]


(* ::Chapter::Closed:: *)
(*Orbital Frequencies*)


(* ::Text:: *)
(*Orbital frequency calculations from Fujita and Hikida, Class. Quantum Grav .26 (2009) 135002, arXiv:0906.1420*)


(* ::Section::Closed:: *)
(*Schwarzschild*)


KerrGeoMinoFrequencies[0,p_,0,x_]:={Sqrt[((-6+p) p)/(-3+p)],p Sqrt[1/((-3+p) x^2)] x,(p x)/Sqrt[(-3+p) x^2],Sqrt[p^5/(-3+p)]}


KerrGeoMinoFrequencies[0,p_,e_,x_]:={(Sqrt[-((p (-6+2 e+p))/(3+e^2-p))] \[Pi])/(2 EllipticK[(4 e)/(-6+2 e+p)]),p Sqrt[1/((-3-e^2+p) x^2)] x,(p x)/Sqrt[(-3-e^2+p) x^2],1/2 Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))] (8+1/((-4+p)^2 EllipticK[(4 e)/(-6+2 e+p)]) (-(((-4+p) p^2 (-6+2 e+p) EllipticE[(4 e)/(-6+2 e+p)])/(-1+e^2))+(p^2 (28+4 e^2-12 p+p^2) EllipticK[(4 e)/(-6+2 e+p)])/(-1+e^2)-(2 (6+2 e-p) (3+e^2-p) p^2 EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)])/((-1+e) (1+e)^2)+(4 (-4+p) p (2 (1+e) EllipticK[(4 e)/(-6+2 e+p)]+(-6-2 e+p) EllipticPi[(2 e (-4+p))/((1+e) (-6+2 e+p)),(4 e)/(-6+2 e+p)]))/(1+e)+2 (-4+p)^2 ((-4+p) EllipticK[(4 e)/(-6+2 e+p)]-((6+2 e-p) p EllipticPi[(16 e)/(12+8 e-4 e^2-8 p+p^2),(4 e)/(-6+2 e+p)])/(2+2 e-p))))}


KerrGeoBoyerLindquistFrequencies[0,p_,0,x_]:={Sqrt[-6+p]/p^2,(Sqrt[1/x^2] x)/p^(3/2),(p x)/Sqrt[p^5 x^2]}


(* ::Section::Closed:: *)
(*Kerr*)


KerrGeoMinoFrequencies[a_,p_,e_,x_]:=Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},
{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];

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


{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]}

]


KerrGeoBoyerLindquistFrequencies[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]},

{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma]} = KerrGeoMinoFrequencies[a,p,e,x];

{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]}/\[CapitalGamma]

]


(* ::Section:: *)
(*Generic function for choosing between frequencies w.r.t different time coordinates*)


Options[KerrGeoFrequencies] = {Time -> "BoyerLindquist"}
SyntaxInformation[KerrGeoFrequencies] = {"ArgumentsPattern"->{_,_,_,_,OptionsPattern[]}};
KerrGeoFrequencies[a_,p_,e_,x_,OptionsPattern[]]:=Module[{M=1,En,L,Q,r1,r2,r3,r4,\[Epsilon]0,zm,a2zp,\[Epsilon]0zp,zmOverZp,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],rp,rm,hr,hp,hm,\[CapitalUpsilon]\[Phi],\[CapitalGamma]},


If[OptionValue["Time"]=="Mino",Return[KerrGeoMinoFrequencies[a,p,e,x][[1;;3]]]];

If[OptionValue["Time"]=="BoyerLindquist", Return[KerrGeoBoyerLindquistFrequencies[a,p,e,x]]];

If[OptionValue["Time"]=="Proper",Print["Propertime frequencies not implemented yet"]];

]


(* ::Chapter::Closed:: *)
(*Orbital Trajectory*)


(* ::Section::Closed:: *)
(*Schwarzschild*)


(* ::Text:: *)
(*The analytic equations below are taken from Appendix B of "Fast Self-forced Inspirals" by M. van de Meent and N. Warburton, Class. Quant. Grav. 35:144003 (2018), arXiv:1802.05281*)


(*t and \[Phi] accumulated over one orbit*)
\[CapitalPhi]SchwarzDarwin[p_,e_]:=4 Sqrt[p/(p-6+2e)] EllipticK[(4 e)/(p-6+2e)]
TSchwarzDarwin[p_,e_]:=(2p Sqrt[(p-6+2e)((p-2)^2-4e^2)])/((1-e^2)(p-4)) EllipticE[(4e)/(p-6+2e)]-2p Sqrt[(p-2)^2-4e^2]/((1-e^2)Sqrt[p-6+2e]) EllipticK[(4e)/(p-6+2e)]-(4(8(1-e^2)+p(1+3e^2-p))Sqrt[(p-2)^2-4e^2])/((1-e)(1-e^2)(p-4)Sqrt[p-6+2e]) EllipticPi[-((2e)/(1-e)),(4e)/(p-6+2e)]+(16Sqrt[(p-2)^2-4e^2])/((p-2+2e)Sqrt[p-6+2e]) EllipticPi[(4e)/(p-2+2e),(4e)/(p-6+2e)]


tSchwarzDarwin[p_,e_,\[Xi]_]:=TSchwarzDarwin[p,e]/2+((p Sqrt[(p-6+2e)((p-2)^2-4e^2)])/((1-e^2)(p-4)) EllipticE[\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-p Sqrt[(p-2)^2-4e^2]/((1-e^2)Sqrt[p-6+2e]) EllipticF[\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-(2(8(1-e^2)+p(1+3e^2-p))Sqrt[(p-2)^2-4e^2])/((1-e)(1-e^2)(p-4)Sqrt[p-6+2e]) EllipticPi[-((2e)/(1-e)),\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]+(8Sqrt[(p-2)^2-4e^2])/((p-2+2e)Sqrt[p-6+2e]) EllipticPi[(4e)/(p-2+2e),\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-e p Sqrt[((p-2)^2-4e^2)(p-6-2e Cos[\[Xi]])]/((1-e^2)(p-4)(1+e Cos[\[Xi]])) Sin[\[Xi]])
rSchwarzDarwin[p_,e_,\[Chi]_]:=p/(1 + e Cos[\[Chi]])
\[Theta]SchwarzDarwin[p_,e_,\[Chi]_]:= \[Pi]/2 
\[Phi]SchwarzDarwin[p_,e_,\[Xi]_]:=\[CapitalPhi]SchwarzDarwin[p,e]/2+2Sqrt[p/(p-6+2e)]EllipticF[\[Xi]/2-\[Pi]/2,(4 e)/(p-6+2e)]


(* ::Text:: *)
(*FIXME: make the below work for inclined orbits and accept initial phases*)


KerrGeoOrbitSchwarzDarwin[p_, e_]:=Module[{t, r, \[Theta], \[Phi], assoc},

t[\[Chi]_] := tSchwarzDarwin[p,e,\[Chi]];
r[\[Chi]_] := rSchwarzDarwin[p,e,\[Chi]];
\[Theta][\[Chi]_] := \[Theta]SchwarzDarwin[p,e,\[Chi]];
\[Phi][\[Chi]_] := \[Phi]SchwarzDarwin[p,e,\[Chi]];

assoc = Association["Trajectory" -> {t,r,\[Theta],\[Phi]}, "Parametrization"->"Darwin"];

KerrGeoOrbitFunction[0, p, e, 0, assoc]

]


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Equatorial (Darwin)*)


(* ::Text:: *)
(*Compute the orbit using Mino time and then convert to Darwin time using \[Lambda][r[\[Chi]]] where \[Lambda][r] is found in Fujita and Hikida (2009).*)


KerrGeoOrbitEquatorialDarwin[a_,p_,e_,x_/;x^2==1]:=Module[{orbitMino,freqs,r1,r2,r3,r4,\[CapitalLambda]r,yr,kr,\[Lambda]0r,r,r01,\[CapitalLambda]r1,\[Lambda],En,L,Q,tMino,rMino,\[Theta]Mino,\[Phi]Mino,tDarwin,rDarwin,\[Theta]Darwin,\[Phi]Darwin,assoc},

orbitMino = KerrGeoOrbit[a,p,e,x];

{r1,r2,r3,r4} = orbitMino["RadialRoots"];
freqs = orbitMino["Frequencies"];
{En,L,Q} = orbitMino["ConstantsOfMotion"];
\[CapitalLambda]r = (2\[Pi])/freqs[[1]];

yr[r_]:=Sqrt[(r1-r3)/(r1-r2) (r-r2)/(r-r3)];
kr=(r1-r2)/(r1-r3) (r3-r4)/(r2-r4);
\[Lambda]0r[r_]:=1/Sqrt[1-En^2] 2/Sqrt[(r1-r3)(r2-r4)] EllipticF[ArcSin[yr[r]],kr];


r[\[Chi]_]:=p/(1+e Cos[\[Chi]]);

r01=r2;
\[CapitalLambda]r1=\[Lambda]0r[r01];


\[Lambda][\[Chi]_]:=\[CapitalLambda]r Floor[\[Chi]/(2\[Pi])]+If[Mod[\[Chi],2\[Pi]]<=\[Pi], \[Lambda]0r[r[\[Chi]]]-\[CapitalLambda]r1,\[CapitalLambda]r-\[Lambda]0r[r[\[Chi]]]];
{tMino, rMino, \[Theta]Mino, \[Phi]Mino} = orbitMino["Trajectory"];

tDarwin[\[Chi]_]:= tMino[\[Lambda][\[Chi]]];
rDarwin[\[Chi]_]:= rMino[\[Lambda][\[Chi]]];
\[Theta]Darwin[\[Chi]_]:= \[Theta]Mino[\[Lambda][\[Chi]]];
\[Phi]Darwin[\[Chi]_]:= \[Phi]Mino[\[Lambda][\[Chi]]];

assoc = Association[
			"Trajectory" -> {tDarwin,rDarwin,\[Theta]Darwin,\[Phi]Darwin}, 
			"Parametrization" -> "Darwin", 
			"ConstantsOfMotion"->{En,L,Q}, 
			"RadialRoots"->{r1,r2,r3,r4},
			"Energy" -> En,
			"AngularMomentum" -> L,
			"CarterConstant" -> Q
		];

KerrGeoOrbitFunction[a, p, e, 0, assoc]

]


(* ::Subsection::Closed:: *)
(*Equatorial (Fast Spec - Darwin)*)


(* Hopper, Forseth, Osburn, and Evans, PRD 92 (2015)*)


(* ::Subsubsection::Closed:: *)
(*Subroutines that checks for the number of samples necessary for spectral integration*)


Clear[DarwinFastSpecIntegrateAndConvergenceCheck];
DarwinFastSpecIntegrateAndConvergenceCheck[func_]:=
Module[{test,compare,res,NInit,iter=1,sampledFunc,phaseList,pg,eps,coeffs,
	coeffsList,coeffsN,\[CapitalDelta]integratedFunc,growthRate,fn,nIter,sampleMax},
	(*DarwinFastSpecIntegrateAndConvergenceCheck takes a function 'func' and integrates
	'func' with respect to the Darwin parameter \[Chi] using spectral methods:
				-- func: a function that takes an integer 'N' as an argument and returns
						 function values at N points. These points are sampled at evenly
						 spaced values of the Darwin parameter \[Chi]  *)
	
	(* Memoize function that we are integrating with respect to the Darwin parameter *)
	sampledFunc[NN_]:=sampledFunc[NN]=func[NN];
	(* Determine precision of sampled points. 
	Use precision to check for convergence of spectral methods *)
	pg=Precision[sampledFunc[16][[2]]];

	(* Use Mathematica's built in DCT solver to determine DCT coefficients *)
	coeffs[Nr_]:=coeffs[Nr]=FourierDCT[sampledFunc[Nr],1]Sqrt[2/(Nr-1)]/.var_?NumericQ/;var==0:>0;
	(* Find the relative accuracy of the DCT coefficients by comparing the last two 
	DCT coefficients to the n=0 coefficient *)
	res[Nr_]:=Min[Abs@RealExponent[coeffs[Nr][[-1]]/coeffs[Nr][[1]]],Abs@RealExponent[coeffs[Nr][[-2]]/coeffs[Nr][[1]]]];

	(* Treat machine precision calculations differently from arbitrary precision *)
	If[pg==$MachinePrecision,
		(* eps sets the precision goal/numerical tolerance for our solutions *)
		eps=15;
		(* Set some initial value of sample points *)
		NInit=2^4;
		(* Find number of sample points necessary to match precision goal 'eps' *)
		While[res[NInit]<eps&&iter<20,NInit=2*NInit;iter++]
		,
		(* Set initial value of sample points at slightly larger value for extended precision calculations *)
		NInit=2^5;
		(* Also test for convergence by increasing the sample size until the n=0 coefficient is unchanged *)
		compare=coeffs[NInit/2][[1]]; (*n=0 coefficient *)
		While[((compare =!= (compare = coeffs[NInit][[1]]))||res[NInit]<pg+1)&&iter<20,NInit=2*NInit; iter++]
	];
	(* After determining the number of sampled points necessary for convergence, store DCT coefficients *)
	coeffsList=coeffs[NInit];
	fn[n_]:=fn[n]=coeffsList[[n+1]];
	growthRate=coeffsList[[1]]/2;
	(* Evaluate new weighted coefficients due to integrating the interpolated inverse DCT transformation *)
	If[pg==MachinePrecision,
		coeffsN[NInit-1]=fn[NInit-1]/(NInit-1);
		nIter=2^2+2;
		Do[coeffsN[n]=coeffsList[[n+1]]/n,{n,1,nIter-2}];
		eps=15-RealExponent[Sum[Abs[fn[n]],{n,0,nIter-2}]];
		While[(-RealExponent[fn[nIter]]<=eps||-RealExponent[fn[nIter-1]]<=eps)&&nIter<NInit-2,coeffsN[nIter]=fn[nIter]/nIter;coeffsN[nIter-1]=fn[nIter-1]/(nIter-1);nIter+=2];
		coeffsN[nIter]=fn[nIter]/nIter;
		coeffsN[nIter-1]=fn[nIter-1]/(nIter-1);
		sampleMax=Min[nIter,NInit-2],
		Do[coeffsN[n]=coeffsList[[n+1]]/n,{n,1,NInit-1}];
		sampleMax=NInit-2;
	];
	(* Construct integrated series solution *)
	\[CapitalDelta]integratedFunc[\[Chi]_?NumericQ]:=coeffsN[NInit-1]/2 Sin[(NInit-1)*\[Chi]]+Sum[coeffsN[nIter] Sin[nIter \[Chi]],{nIter,1,sampleMax}];
	(* Allow function to evaluate lists by threading over them *)
	\[CapitalDelta]integratedFunc[\[Chi]List_List]:=\[CapitalDelta]integratedFunc[#]&/@\[Chi]List;
	(* Return the linear rate of growth and the oscillatory function \[CapitalDelta]integratedFunc *)
	{growthRate,\[CapitalDelta]integratedFunc}
];


(* ::Subsubsection::Closed:: *)
(*Main file that calculates geodesics using spectral integration*)


Clear[KerrGeoOrbitFastSpecDarwin];
KerrGeoOrbitFastSpecDarwin[a_,p_,e_,x_/;x^2==1,initPhases:{_,_,_,_}:{0,0,0,0}]:=
Module[{M=1,En,L,Q,r1,r2,r3,r4,p3,p4,assoc,var,t0, \[Chi]0, \[Phi]0,r0,\[Theta]0,t,r,\[Theta],\[Phi],\[Chi],growthRateT,growthRatePh,
		\[Chi]r,NrMax,pg,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[Phi]C,tC,Pr,r0Sample,PrSample,dtd\[Chi],d\[Phi]d\[Chi],TVr,PVr},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, x, En, Q];
	p3=(1-e)r3/M;
	p4=(1+e)r4/M;
	
	(*Precision of sampling depends on precision of arguments*)
	pg=Min[{Precision[{a,p,e,x}],Precision[initPhases]}];
	
	(*Parameterize r in terms of Darwin parameter \[Chi]*)
	r0[chi_]:=p M/(1+e Cos[chi]);
	\[Theta]0[chi_?NumericQ]:=N[Pi/2,pg];
	\[Theta]0[chi_List]:=\[Theta]0[#]&/@chi;
	
	(* Expressions for dt/d\[Lambda] = TVr and d\[Phi]/d\[Lambda] = PVr *)
	TVr[rp_]:=(En*(a^2 + rp^2)^2)/(a^2 - 2*M*rp + rp^2) + a*L*(1 - (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2))-a^2*En;
	PVr[rp_]:=-((a^2*L)/(a^2 - 2*M*rp + rp^2)) + a*En*(-1 + (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2))+L;
	
	(* Sampling of radial position using evenly spaced values of Darwin parameter \[Chi]*)
	If[pg==$MachinePrecision,
		\[Chi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}]],
		\[Chi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}],1.5pg]
	];
	r0Sample[Nr_]:=r0[\[Chi]r[Nr]];
	
	(* Expression for d\[Lambda]/d\[Chi]*)
	Pr[chi_]:=(1-e^2)/Sqrt[(p-p4)+e(p-p4 Cos[chi])]/(M (1-En^2)^(1/2)Sqrt[(p-p3)-e(p+p3 Cos[chi])]);
	PrSample[Nr_]:=Pr[\[Chi]r[Nr]];
	
	(* Sampling of expressions for dt/d\[Chi] and d\[Phi]/d\[Chi] *)
	dtd\[Chi][Nr_]:=TVr[r0Sample[Nr]]PrSample[Nr];
	d\[Phi]d\[Chi][Nr_]:=PVr[r0Sample[Nr]]PrSample[Nr];
	
	(*Spectral integration of t and \[Phi] as functions of \[Chi]*)
	{growthRateT,\[CapitalDelta]tr}=DarwinFastSpecIntegrateAndConvergenceCheck[dtd\[Chi]];
	{growthRatePh,\[CapitalDelta]\[Phi]r}=DarwinFastSpecIntegrateAndConvergenceCheck[d\[Phi]d\[Chi]];
	
	(*Collect initial phases*)
	{t0, \[Chi]0, \[Phi]0} = {initPhases[[1]], initPhases[[2]], initPhases[[4]]};
	(* Find integration constants for t and \[Phi], so that t(\[Chi]=0)=t0 and \[Phi](\[Chi]=0)=\[Phi]0 *)
	\[Phi]C=\[CapitalDelta]\[Phi]r[\[Chi]0];
	tC=\[CapitalDelta]tr[\[Chi]0];
	
	t[\[Chi]_]:=\[CapitalDelta]tr[\[Chi]+\[Chi]0]+growthRateT \[Chi]+t0-tC;
	r[\[Chi]_]:=r0[\[Chi]+\[Chi]0];
	\[Theta][\[Chi]_]:=\[Theta]0[\[Chi]];
	\[Phi][\[Chi]_]:=\[CapitalDelta]\[Phi]r[\[Chi]+\[Chi]0]+growthRatePh \[Chi]+\[Phi]0-\[Phi]C;
	
	assoc = Association[
		"Parametrization"->"Darwin",
		"Energy" -> En, 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> {En,L,Q},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"RadialRoots" -> {r1,r2,r3,r4}
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Subsection::Closed:: *)
(*Circular (Fast Spec - Darwin)*)


(* Hopper, Forseth, Osburn, and Evans, PRD 92 (2015)*)


(* ::Subsubsection::Closed:: *)
(*Main file that calculates geodesics using spectral integration*)


KerrGeoOrbitFastSpecDarwin[a_,p_,e_/;e==0,x_,initPhases:{_,_,_,_}:{0,0,0,0}]:=
Module[{M=1,En,L,Q,zp,zm,assoc,var,t0, \[Chi]0, \[Phi]0,r0,\[Theta]0,t,r,\[Theta],\[Phi],\[Chi],freqT,freqPh,
		\[Chi]\[Theta],pg,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]\[Theta],\[Phi]C,tC,P\[Theta],\[Theta]0Sample,P\[Theta]Sample,dtd\[Chi],d\[Phi]d\[Chi],TV\[Theta],PV\[Theta],\[Beta],\[Alpha],zRoots},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	
	(* Useful constants for \[Theta]-dependent calculations *)
	\[Beta]=a^2(1-En^2);
	\[Alpha]=L^2+Q+\[Beta];
	zp=Sqrt[(\[Alpha]+Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zm=Sqrt[(\[Alpha]-Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zRoots={ArcCos[zm],Pi-ArcCos[zm]};
	
	(*Precision of sampling depends on precision of arguments*)
	pg=Min[{Precision[{a,p,e,x}],Precision[initPhases]}];
	
	(*Parameterize \[Theta] in terms of Darwin-like parameter \[Chi]*)
	r0[chi_?NumericQ]:=N[p M,pg];
	r0[chi_List]:=r0[#]&/@chi;
	\[Theta]0[chi_]:=ArcCos[zm Cos[chi]];

	(* Expressions for dt/d\[Lambda] = TV\[Theta] and d\[Phi]/d\[Lambda] = PV\[Theta] *)
	TV\[Theta][\[Theta]p_]:=(En*(a^2 + M^2 p^2)^2)/(a^2 - 2*M^2*p + M^2p^2) + a*L*(1 - (a^2 + p^2)/(a^2 - 2*M^2*p + M^2 p^2))-a^2*En Sin[\[Theta]p]^2;
	PV\[Theta][\[Theta]p_]:=-((a^2*L)/(a^2 - 2*M^2*p + M^2p^2)) + a*En*(-1 + (a^2 + M^2p^2)/(a^2 - 2*M^2*p + M^2p^2))+L Csc[\[Theta]p]^2;
	
	(* Sampling of radial position using evenly spaced values of Darwin parameter \[Chi]*)
	If[pg==$MachinePrecision,
		\[Chi]\[Theta][Nth_]:=N[Table[i Pi/(Nth-1),{i,0,Nth-1}]],
		\[Chi]\[Theta][Nth_]:=N[Table[i Pi/(Nth-1),{i,0,Nth-1}],1.5pg]
	];
	\[Theta]0Sample[Nth_]:=\[Theta]0[\[Chi]\[Theta][Nth]];
	
	(* Expression for d\[Lambda]/d\[Chi]*)
	P\[Theta][chi_]:=(\[Beta](zp^2-zm^2 Cos[chi]^2))^(-1/2);
	P\[Theta]Sample[Nth_]:=P\[Theta][\[Chi]\[Theta][Nth]];
	
	(* Sampling of expressions for dt/d\[Chi] and d\[Phi]/d\[Chi] *)
	dtd\[Chi][Nr_]:=TV\[Theta][\[Theta]0Sample[Nr]]P\[Theta]Sample[Nr];
	d\[Phi]d\[Chi][Nr_]:=PV\[Theta][\[Theta]0Sample[Nr]]P\[Theta]Sample[Nr];

	(*Spectral integration of t and \[Phi] as functions of \[Chi]*)
	{freqT,\[CapitalDelta]t\[Theta]}=DarwinFastSpecIntegrateAndConvergenceCheck[dtd\[Chi]];
	{freqPh,\[CapitalDelta]\[Phi]\[Theta]}=DarwinFastSpecIntegrateAndConvergenceCheck[d\[Phi]d\[Chi]];
	
	(*Collect initial phases*)
	{t0, \[Chi]0, \[Phi]0} = {initPhases[[1]], initPhases[[3]], initPhases[[4]]};
	(* Find integration constants for t and \[Phi], so that t(\[Chi]=0)=t0 and \[Phi](\[Chi]=0)=\[Phi]0 *)
	\[Phi]C=\[CapitalDelta]\[Phi]\[Theta][\[Chi]0];
	tC=\[CapitalDelta]t\[Theta][\[Chi]0];
	
	t[\[Chi]_]:=\[CapitalDelta]t\[Theta][\[Chi]+\[Chi]0]+freqT \[Chi]+t0-tC;
	r[\[Chi]_]:=r0[\[Chi]];
	\[Theta][\[Chi]_]:=\[Theta]0[\[Chi]+\[Chi]0];
	\[Phi][\[Chi]_]:=\[CapitalDelta]\[Phi]\[Theta][\[Chi]+\[Chi]0]+freqPh \[Chi]+\[Phi]0-\[Phi]C;
	
	assoc = Association[
		"Parametrization"->"Darwin",
		"Energy" -> En, 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> {En,L,Q},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"PolarRoots" -> zRoots
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Subsection::Closed:: *)
(*Generic (Mino)*)


KerrGeoOrbitMino[a_,p_,e_,x_,initPhases:{_,_,_,_}:{0,0,0,0}]:=Module[{M=1,En,L,Q,assoc,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,zp,zm,kr,k\[Theta],rp,rm,hr,hp,hm,rq,zq,\[Psi]r,tr,\[Phi]f,\[Psi]z,tz,\[Phi]z,qt0,qr0,qz0,q\[Phi]0,t,r,\[Theta],\[Phi],\[Phi]t,\[Phi]r,Ct,C\[Phi]},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = KerrGeoMinoFrequencies[a,p,e,x];
	{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, x, En, Q];
	{zp,zm} = KerrGeoPolarRoots[a, p, e, x];
	
	kr = (r1-r2)/(r1-r3) (r3-r4)/(r2-r4);
	k\[Theta] = a^2 (1-En^2)(zm/zp)^2;


rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));

rq = Function[{qr},(r3(r1 - r2)JacobiSN[EllipticK[kr]/\[Pi] qr,kr]^2-r2(r1-r3))/((r1-r2)JacobiSN[EllipticK[kr]/\[Pi] qr,kr]^2-(r1-r3))];

zq = Function[{qz}, -zm JacobiSN[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]]];

\[Psi]r[qr_]:=\[Psi]r[qr]= JacobiAmplitude[EllipticK[kr]/\[Pi] qr,kr];

tr[qr_]:= -En/Sqrt[(1-En^2) (r1-r3) (r2-r4)] (
4(r2-r3) (EllipticPi[hr,kr] qr/\[Pi]-EllipticPi[hr,\[Psi]r[qr],kr])
-4 (r2-r3)/(rp-rm) (
-(1/((-rm+r2) (-rm+r3)))(-2 a^2+rm (4-(a L)/En)) (EllipticPi[hm,kr] qr/\[Pi]-EllipticPi[hm,\[Psi]r[qr],kr] )
+1/((-rp+r2) (-rp+r3)) (-2 a^2+rp (4-(a L)/En)) (EllipticPi[hp,kr] qr/\[Pi]-EllipticPi[hp,\[Psi]r[qr],kr])
)
+(r2-r3) (r1+r2+r3+r4) (EllipticPi[hr,kr] qr/\[Pi]-EllipticPi[hr,\[Psi]r[qr],kr] )
+(r1-r3) (r2-r4) (EllipticE[kr] qr/\[Pi]-EllipticE[\[Psi]r[qr],kr]+hr((Sin[\[Psi]r[qr]]Cos[\[Psi]r[qr]] Sqrt[1-kr Sin[\[Psi]r[qr]]^2])/(1-hr Sin[\[Psi]r[qr]]^2))) );

\[Phi]r[qr_]:= (2 a En (-1/((-rm+r2) (-rm+r3))(2 rm-(a L)/En) (r2-r3) (EllipticPi[hm,kr] qr/\[Pi]-EllipticPi[hm,\[Psi]r[qr],kr])+1/((-rp+r2) (-rp+r3))(2 rp-(a L)/En) (r2-r3) (EllipticPi[hp,kr] qr/\[Pi]-EllipticPi[hp,\[Psi]r[qr],kr] )))/((-rm+rp) Sqrt[(1-En^2) (r1-r3) (r2-r4)]);

\[Psi]z[qz_]:=\[Psi]z[rq]= JacobiAmplitude[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]];
tz[qz_]:= 1/(1-En^2) En zp ( EllipticE[k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticE[\[Psi]z[qz],k\[Theta]]);
\[Phi]z[qz_]:= -1/zp L ( EllipticPi[zm^2,k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticPi[zm^2,\[Psi]z[qz],k\[Theta]]);

{qt0, qr0, qz0, q\[Phi]0} = {initPhases[[1]], initPhases[[2]], initPhases[[3]], initPhases[[4]]};
(*Calculate normalization constants so that t=0 and \[Phi]=0 at \[Lambda]=0 when qt0=0 and q\[Phi]0=0 *)
Ct=tr[qr0]+tz[qz0];
C\[Phi]=\[Phi]r[qr0]+\[Phi]z[qz0];

t[\[Lambda]_]:= qt0 + \[CapitalUpsilon]t \[Lambda] + tr[\[CapitalUpsilon]r \[Lambda] + qr0] + tz[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]-Ct;
r[\[Lambda]_]:= rq[\[CapitalUpsilon]r \[Lambda]+ qr0];
\[Theta][\[Lambda]_]:= ArcCos[zq[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]];
\[Phi][\[Lambda]_]:= q\[Phi]0 + \[CapitalUpsilon]\[Phi] \[Lambda] + \[Phi]r[\[CapitalUpsilon]r \[Lambda]+ qr0] + \[Phi]z[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]-C\[Phi];

	assoc = Association[
	"Parametrization"->"Mino", 
	"Energy" -> En, 
	"AngularMomentum" -> L, 
	"CarterConstant" -> Q, 
	"ConstantsOfMotion" -> {En,L,Q},
	"RadialFrequency" -> \[CapitalUpsilon]r,
	"PolarFrequency" ->  \[CapitalUpsilon]\[Theta],
	"AzimuthalFrequency" -> \[CapitalUpsilon]\[Phi],
	"Frequencies" -> {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
	"Trajectory" -> {t,r,\[Theta],\[Phi]},
	"RadialRoots" -> {r1,r2,r3,r4}
	];

	KerrGeoOrbitFunction[a,p,e,x,assoc]

]


(* ::Subsection::Closed:: *)
(*Generic (Fast Spec - Mino)*)


(* Hopper, Forseth, Osburn, and Evans, PRD 92 (2015)*)


(* ::Subsubsection::Closed:: *)
(*Subroutines for calculating \[Lambda](\[Psi]) and \[Lambda](\[Chi])*)


Clear[MinoRFastSpec];
MinoRFastSpec[a_,p_,e_,x_]:=
Module[{M=1,En,L,Q,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,\[Psi]r,r0,p3,p4,Pr,PrSample,NrInit=2^4,
		\[ScriptCapitalP]rSample,\[ScriptCapitalP]rn,\[ScriptCapitalP]rList,\[CapitalDelta]\[Lambda]r,pg,iter=1,compare,res,rate},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, x, En, Q];
	p3=(1-e)r3/M;
	p4=(1+e)r4/M;
	pg=Precision[{a,p,e,x}];
	
	If[pg==$MachinePrecision,
		\[Psi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}]],
		\[Psi]r[Nr_]:=N[Table[i Pi/(Nr-1),{i,0,Nr-1}],1.5pg]
	];
	
	Pr[psi_]:=Pr[psi]=(1-e^2)((p-p4)+e(p-p4 Cos[psi]))^(-1/2)/(M (1-En^2)^(1/2)((p-p3)-e(p+p3 Cos[psi]))^(1/2));
	PrSample[Nr_]:=Pr[\[Psi]r[Nr]];
	
	{rate,\[CapitalDelta]\[Lambda]r}=DarwinFastSpecIntegrateAndConvergenceCheck[PrSample];
	\[CapitalDelta]\[Lambda]r
];


Clear[MinoThetaFastSpec];
MinoThetaFastSpec[a_,p_,e_,x_]:=
Module[{M=1,En,L,Q,zp,zm,\[Chi]\[Theta],\[Beta],P\[Theta],P\[Theta]Sample,NthInit=2^4,\[ScriptCapitalP]\[Theta]List,\[ScriptCapitalP]\[Theta]Sample,
	\[ScriptCapitalP]\[Theta]k,\[CapitalDelta]\[Lambda]\[Theta],pg,iter=1,compare,res,\[Alpha],rate},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	\[Beta]=a^2(1-En^2);
	\[Alpha]=L^2+Q+\[Beta];
	zp=Sqrt[(\[Alpha]+Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zm=Sqrt[(\[Alpha]-Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	pg=Precision[{a,p,e,x}];
	
	If[pg==$MachinePrecision,
		\[Chi]\[Theta][Nth_]:=N[Table[i Pi/(Nth-1),{i,0,Nth-1}]],
		\[Chi]\[Theta][Nth_]:=N[Table[i Pi/(Nth-1),{i,0,Nth-1}],1.5pg]
	];
	P\[Theta][chi_]:=P\[Theta][chi]=(\[Beta](zp^2-zm^2 Cos[chi]^2))^(-1/2);
	P\[Theta]Sample[Nth_]:=P\[Theta][\[Chi]\[Theta][Nth]];
	
	{rate,\[CapitalDelta]\[Lambda]\[Theta]}=DarwinFastSpecIntegrateAndConvergenceCheck[P\[Theta]Sample];
	\[CapitalDelta]\[Lambda]\[Theta]
];


(* ::Subsubsection::Closed:: *)
(*Subroutine for calculating \[Psi](\[Lambda]) and \[Chi](\[Lambda])*)


PhaseOfMinoFastSpec[\[CapitalUpsilon]_,sampledMino_]:=
Module[{sampledFunc,NInit,phase},
	NInit=Length[sampledMino];
	sampledFunc[NN_]:=ConstantArray[1,NN];
	phase=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFunc,{\[CapitalUpsilon],sampledMino}];
	phase
];


(* ::Subsubsection::Closed:: *)
(*Subroutines for calculating \[CapitalDelta]\[Phi]r(\[Lambda]), \[CapitalDelta]\[Phi]\[Theta](\[Lambda]), \[CapitalDelta]tr(\[Lambda]), \[CapitalDelta]t\[Theta](\[Lambda])*)


PhiOfMinoFastSpecR[a_,p_,e_,x_,{\[CapitalUpsilon]r_,minoSampleR_}]:=
Module[{M=1,En,L,Q,sampledFuncR,sampledMinoR,PVr,\[CapitalDelta]\[Phi]r},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	PVr[rp_]:=-((a^2*L)/(a^2 - 2*M*rp + rp^2)) + a*En*(-1 + (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2));

	sampledFuncR=LambdaToPsiRTransform[a,p,e,x,PVr];
	\[CapitalDelta]\[Phi]r=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncR,{\[CapitalUpsilon]r,minoSampleR}];
	\[CapitalDelta]\[Phi]r
];

PhiOfMinoFastSpecR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[Lambda]rSample,pg,\[Psi]r},
	\[CapitalUpsilon]r = Re[KerrGeoMinoFrequencies[a,p,e,x]][[1]];
	pg=Precision[{a,p,e,x}];

	If[pg==$MachinePrecision,
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}]],
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}],1.5pg]
	];
	\[CapitalDelta]\[Lambda]r=MinoRFastSpec[a,p,e,x];
	\[Lambda]r[psi_]:=\[Lambda]r[psi]=\[CapitalDelta]\[Lambda]r[psi];
	\[Lambda]rSample[Nr_]:=\[Lambda]r[\[Psi]r[Nr]];
	PhiOfMinoFastSpecR[a,p,e,x,{\[CapitalUpsilon]r,\[Lambda]rSample}]
];


PhiOfMinoFastSpecTheta[a_,p_,e_,x_,{\[CapitalUpsilon]\[Theta]_,minoSampleTh_}]:=
Module[{M=1,En,L,Q,sampledFuncTheta,sampledMinoTheta,PV\[Theta],\[CapitalDelta]\[Phi]\[Theta]},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	PV\[Theta][\[Theta]p_]:=L*Csc[\[Theta]p]^2;
	
	sampledFuncTheta=LambdaToChiThetaTransform[a,p,e,x,PV\[Theta]];
	\[CapitalDelta]\[Phi]\[Theta]=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncTheta,{\[CapitalUpsilon]\[Theta],minoSampleTh}];
	\[CapitalDelta]\[Phi]\[Theta]
];

PhiOfMinoFastSpecTheta[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]\[Theta],\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],\[Lambda]\[Theta]Sample,pg,\[Chi]\[Theta]},
	\[CapitalUpsilon]\[Theta] = Re[KerrGeoMinoFrequencies[a,p,e,x]][[2]];
	pg=Precision[{a,p,e,x}];

	If[pg==$MachinePrecision,
		\[Chi]\[Theta][Nth_]:=N[Table[i 2 Pi/Nth,{i,0,Nth-1}]],
		\[Chi]\[Theta][Nth_]:=N[Table[i 2 Pi/Nth,{i,0,Nth-1}],1.5pg]
	];
	\[CapitalDelta]\[Lambda]\[Theta]=MinoRFastSpec[a,p,e,x];
	\[Lambda]\[Theta][psi_]:=\[Lambda]\[Theta][psi]=\[CapitalDelta]\[Lambda]\[Theta][psi];
	\[Lambda]\[Theta]Sample[Nr_]:=\[Lambda]\[Theta][\[Chi]\[Theta][Nr]];
	PhiOfMinoFastSpecTheta[a,p,e,x,{\[CapitalUpsilon]\[Theta],\[Lambda]\[Theta]Sample}]
];


TimeOfMinoFastSpecR[a_,p_,e_,x_,{\[CapitalUpsilon]r_,minoSampleR_}]:=
Module[{M=1,En,L,Q,sampledFuncR,sampledMinoR,TVr,\[CapitalDelta]tr},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	TVr[rp_]:=(En*(a^2 + rp^2)^2)/(a^2 - 2*M*rp + rp^2) + a*L*(1 - (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2));

	sampledFuncR=LambdaToPsiRTransform[a,p,e,x,TVr];
	\[CapitalDelta]tr=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncR,{\[CapitalUpsilon]r,minoSampleR}];
	\[CapitalDelta]tr
];

TimeOfMinoFastSpecR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[Lambda]rSample,pg,\[Psi]r},
	\[CapitalUpsilon]r = Re[KerrGeoMinoFrequencies[a,p,e,x]][[1]];
	pg=Precision[{a,p,e,x}];

	If[pg==$MachinePrecision,
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}]],
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}],1.5pg]
	];
	\[CapitalDelta]\[Lambda]r=MinoRFastSpec[a,p,e,x];
	\[Lambda]r[psi_]:=\[Lambda]r[psi]=\[CapitalDelta]\[Lambda]r[psi];
	\[Lambda]rSample[Nr_]:=\[Lambda]r[\[Psi]r[Nr]];
	TimeOfMinoFastSpecR[a,p,e,x,{\[CapitalUpsilon]r,\[Lambda]rSample}]
];


TimeOfMinoFastSpecTheta[a_,p_,e_,x_,{\[CapitalUpsilon]\[Theta]_,minoSampleTh_}]:=
Module[{M=1,En,L,Q,sampledFuncTheta,sampledMinoTheta,TV\[Theta],\[CapitalDelta]t\[Theta]},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	TV\[Theta][\[Theta]p_]:=-(a^2*En*Sin[\[Theta]p]^2);

	sampledFuncTheta=LambdaToChiThetaTransform[a,p,e,x,TV\[Theta]];
	\[CapitalDelta]t\[Theta]=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncTheta,{\[CapitalUpsilon]\[Theta],minoSampleTh}];
	\[CapitalDelta]t\[Theta]
];

TimeOfMinoFastSpecTheta[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]\[Theta],\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],\[Lambda]\[Theta]Sample,pg,\[Chi]\[Theta]},
	\[CapitalUpsilon]\[Theta] = Re[KerrGeoMinoFrequencies[a,p,e,x]][[2]];
	pg=Precision[{a,p,e,x}];

	If[pg==$MachinePrecision,
		\[Chi]\[Theta][Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}]],
		\[Chi]\[Theta][Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}],1.5pg]
	];
	\[CapitalDelta]\[Lambda]\[Theta]=MinoRFastSpec[a,p,e,x];
	\[Lambda]\[Theta][psi_]:=\[Lambda]\[Theta][psi]=\[CapitalDelta]\[Lambda]\[Theta][psi];
	\[Lambda]\[Theta]Sample[Nr_]:=\[Lambda]\[Theta][\[Chi]\[Theta][Nr]];
	TimeOfMinoFastSpecTheta[a,p,e,x,{\[CapitalUpsilon]\[Theta],\[Lambda]\[Theta]Sample}]
];


(* ::Subsubsection::Closed:: *)
(*Generic subroutines that transform functions from *)
(*\[Lambda] dependence to \[Psi] or \[Chi] dependence*)


LambdaToPsiRTransform[a_,p_,e_,x_,rFunc_]:=
Module[{M=1,En,L,Q,r1,r2,r3,r4,p3,p4,\[Psi]r,r0,r0Sample,Pr,PrSample,rFuncSample,pg},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, x, En, Q];
	p3=(1-e)r3/M;
	p4=(1+e)r4/M;
	pg=Precision[{a,p,e,x}];

	If[pg==$MachinePrecision,
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}]],
		\[Psi]r[Nr_]:=N[Table[i 2 Pi/Nr,{i,0,Nr-1}],1.5pg]
	];
	r0[psi_]:=p M/(1+e Cos[psi]);
	r0Sample[Nr_]:=r0[\[Psi]r[Nr]];
	
	Pr[psi_]:=(1-e^2)/Sqrt[(p-p4)+e(p-p4 Cos[psi])]/(M (1-En^2)^(1/2)Sqrt[(p-p3)-e(p+p3 Cos[psi])]);
	PrSample[Nr_]:=Pr[\[Psi]r[Nr]];
	
	rFuncSample[Nr_]:=rFunc[r0Sample[Nr]]PrSample[Nr];
	rFuncSample
];


LambdaToChiThetaTransform[a_,p_,e_,x_,thFunc_]:=
Module[{M=1,En,L,Q,zp,zm,\[Beta],\[Chi]\[Theta],\[Theta]0,\[Theta]0Sample,P\[Theta],P\[Theta]Sample,thFuncSample,pg,\[Alpha]},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	\[Beta]=a^2(1-En^2);
	\[Alpha]=L^2+Q+\[Beta];
	zp=Sqrt[(\[Alpha]+Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zm=Sqrt[(\[Alpha]-Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	pg=Precision[{a,p,e,x}];
	
	If[pg==$MachinePrecision,
		\[Chi]\[Theta][Nth_]:=N[Table[i 2 Pi/Nth,{i,0,Nth-1}]],
		\[Chi]\[Theta][Nth_]:=N[Table[i 2 Pi/Nth,{i,0,Nth-1}],1.5pg]
	];
	\[Theta]0[chi_]:=ArcCos[zm Cos[chi]];
	\[Theta]0Sample[Nth_]:=\[Theta]0[\[Chi]\[Theta][Nth]];
	
	P\[Theta][chi_]:=(\[Beta](zp^2-zm^2 Cos[chi]^2))^(-1/2);
	P\[Theta]Sample[Nth_]:=P\[Theta][\[Chi]\[Theta][Nth]];
	
	thFuncSample[Nth_]:=thFunc[\[Theta]0Sample[Nth]]P\[Theta]Sample[Nth];
	thFuncSample
];


(* ::Subsubsection::Closed:: *)
(*Subroutines that checks for the number of samples necessary for spectral integration of an even function*)


Clear[DarwinFastSpecMinoIntegrateAndConvergenceCheck];
DarwinFastSpecMinoIntegrateAndConvergenceCheck[func_,{freq_,mino_}]:=
Module[{test,compare,res,NInit,iter=1,fn,sampledFunc,phaseList,pg,eps,nTest},
	(*
	DarwinFastSpecMinoIntegrateAndConvergenceCheck takes an even function 'func' 
	and determines the number of samples necessary to integrate 'func' with
	respect to Mino time \[Lambda] using spectral sampling in the Darwin-like parameters 
	\[Psi] and \[Chi]:
				-- func: a function that includes function values 
				that result from sampling the function at evenly spaced values 
				of the Darwin-like parameters \[Psi] or \[Chi]
				-- freq: Mino time frequency with respect to r or \[Theta] (\[CapitalUpsilon]r or \[CapitalUpsilon]\[Theta])
				-- mino: a list of Mino time values that results from 
				sampling the function at evenly spaced values of the Darwin-like 
				parameters \[Psi] or \[Chi] 
	*)
	
	(* Memoize function that we are integrating with respect to the Darwin parameter *)
	sampledFunc[NN_]:=sampledFunc[NN]=func[NN];
	(* Determine precision of arguments. 
	Use precision to check for convergence of spectral methods *)
	pg=Precision[freq];

	(* Treat machine precision calculations differently from arbitrary precision *)
	If[pg==MachinePrecision,
		(* eps sets the precision goal/numerical tolerance for our solutions *)
		eps=15;
		(* Set some initial value of sample points *)
		NInit=2^3;
		(* Sampled points need to also be weighted by the Mino time *)
		phaseList[NN_]:=phaseList[NN]=freq*mino[NN]+N[Table[2Pi i/NN,{i,0,NN-1}]];
		(* Create functions for discrete cosine series coefficient fn *)
		fn[NN_,nn_]:=Total[sampledFunc[NN]Cos[nn phaseList[NN]]]/NN;
		(* Test convergence by comparing coefficients to the n=0 coefficient*)
		nTest=0;
		test[NN_]:=fn[NN,nTest];
		(* If the n=0 coefficient vanishes or is unity, compare against n=1 coefficient *)
		If[test[NInit]==0||test[NInit]==1,nTest=1];
		(* Find the relative accuracy of the DFT coefficients by comparing the last 
		DFT coefficient to the test coefficient set above *)
		res[NN_]:=Abs@RealExponent[(fn[NN,NN/2]/.{y_/;y==0:>0})/fn[NN,0]];
		(* Find number of sample points necessary to match precision goal 'eps' *)
		While[res[NInit]<eps&&iter<8,NInit=2*NInit; iter++]
		,
		(* Set some initial value of sample points *)
		NInit=2^6;
		(* Same process as above but with higher precision tolerances *)
		phaseList[NN_]:=phaseList[NN]=freq*mino[NN]+N[Table[2Pi i/NN,{i,0,NN-1}],1.5pg];
		fn[NN_,nn_]:=Total[sampledFunc[NN]Cos[nn phaseList[NN]]]/NN;
		nTest=0;
		test[NN_]:=fn[NN,nTest];
		If[test[NInit]==0||test[NInit]==1,nTest=1];
		res[NN_]:=Abs@RealExponent[(fn[NN,NN/2]/.{y_/;y==0:>0})/fn[NN,0]];
		(* Also test for convergence by increasing the sample size until the test coefficient is unchanged *)
		compare=test[NInit/2];
		While[((compare =!= (compare = test[NInit]))||res[NInit/2]<pg+1)&&iter<10,NInit=2*NInit; iter++];
		(* Compare residuals by comparing to a different Fourier coefficient to ensure proper convergence *)
		(* This part might be overkill *)
		nTest++;
		compare=test[NInit/2];
		iter=1;
		While[((compare =!= (compare = test[NInit]))||res[NInit/2]<pg+1)&&iter<10,NInit=2*NInit; iter++];
		(* If the Fourier coefficients were unaffected after doubling the sampling size, 
		then we can sample using the previous number of sample points *)
		If[res[NInit]==0&&res[NInit/2]!=0,NInit,NInit=NInit/2];
	];
	(* After determining the necessary number of sample points, we sample the function
	we want to integrate and the Mino time, and pass both lists of sampled points to our
	spectral integrator *)
	DarwinFastSpecMinoIntegrateEven[sampledFunc[NInit],{freq,mino[NInit]}]
];


(* ::Subsubsection::Closed:: *)
(*Subroutine that performs spectral integration on even functions*)


DarwinFastSpecMinoIntegrateEven[sampledFunctionTemp_,{freq_,minoTimeTemp_}]:=
Module[{sampledF,fn,fList,f,sampleN,\[Lambda],integratedF,phaseList,pg,nn,samplePhase,nList,
		sampleMax,eps,sampleHalf},
	(*
	DarwinFastSpecMinoIntegrateEven takes an even function 'sampledFunctionTemp' 
	and integrates 'sampledFunctionTemp' with respect to Mino time \[Lambda] using spectral 
	sampling in the Darwin-like parameters \[Psi] and \[Chi]:
				-- sampledFunctionTemp: a list that includes function values 
				that result from sampling the function at evenly spaced values 
				of the Darwin-like parameters \[Psi] or \[Chi]
				-- freq: Mino time frequency with respect to r or \[Theta] (\[CapitalUpsilon]r or \[CapitalUpsilon]\[Theta])
				-- minoTimeTemp: a list of Mino time values that results from 
				sampling the function at evenly spaced values of the Darwin-like 
				parameters \[Psi] or \[Chi] 
	*)
	
	(* Represent values equal to zero with infinite precision *)
	sampledF=sampledFunctionTemp/.x_?NumericQ/;x==0:>0;
	\[Lambda]=minoTimeTemp/.x_?NumericQ/;x==0:>0;
	
	(* Determine the number of sampled points and precision of arguments *)
	sampleN=Length[sampledF];
	pg=Precision[freq];
	
	(* Use precision and number of sampled points to generate list of
	evenly spaced values of \[Psi] or \[Chi] *)
	If[pg==$MachinePrecision,
		phaseList=N[Table[2Pi i/sampleN,{i,0,sampleN-1}]],
		phaseList=N[Table[2Pi i/sampleN,{i,0,sampleN-1}],1.5pg]
	];
	
	(* Create functions for discrete cosine series coefficient fn *)
	samplePhase=(freq \[Lambda]+phaseList);
	fn[n_]:=fn[n]=(sampledF.Cos[nn*samplePhase])/.nn->n;
	
	(* Calculate series coefficients until they equal 0 (with respect 
	to the precision being used) *)
	If[pg==MachinePrecision,
		fList=Block[{halfSample,nIter},
			nIter=2^2+2;
			eps=15-RealExponent[Sum[Abs@fn[n],{n,0,nIter-2}]];
			While[(-RealExponent[fn[nIter]]<=eps||-RealExponent[fn[nIter-1]]<=eps)&&nIter<sampleN/2,nIter+=2];
			sampleMax=Min[nIter,sampleN/2];
			1/sampleN fn[Table[n,{n,0,sampleMax}]]/.x_?NumericQ/;x==0:>0
		],
		fList=Block[{halfSample,nIter},
			nIter=2^5+2;
			While[(fn[nIter]!=0||fn[nIter-1]!=0)&&nIter<sampleN/2,nIter+=2];
			sampleMax=Min[nIter,sampleN/2];
			1/sampleN fn[Table[n,{n,0,sampleMax}]]/.x_?NumericQ/;x==0:>0
		]
	];
	f[n_]:=fList[[n+1]];
	
	(* Construct integrated series solution *)
	integratedF[mino_?NumericQ]:=2*Sum[f[n]/n Sin[n freq mino],{n,1,sampleMax}];
	(* Allow function to evaluate lists by threading over them *)
	integratedF[minoList_List]:=integratedF[#]&/@minoList;
	integratedF
];


(* ::Subsubsection::Closed:: *)
(*Main file that calculates geodesics using spectral integration*)


Clear[KerrGeoOrbitFastSpec];
Options[KerrGeoOrbitFastSpec]={InitialPosition->{0,0,0,0}};
KerrGeoOrbitFastSpec[a_,p_,e_,x_,initPhases:{_,_,_,_}:{0,0,0,0},opts:OptionsPattern[]]:=
Module[{M=1,En,L,Q,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,p3,p4,\[Alpha],\[Beta],zp,zm,assoc,var,\[Chi]0,\[Psi]0,
		r0,\[Theta]0,qt0,qr0,q\[Theta]0,q\[Phi]0,\[Lambda]t0,\[Lambda]r0,\[Lambda]\[Theta]0,\[Lambda]\[Phi]0,t,r,\[Theta],\[Phi],\[Psi],\[Chi],\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],rC,\[Theta]C,\[CapitalDelta]r,\[CapitalDelta]\[Theta],
		\[Psi]r,\[Chi]\[Theta],NrMax,NthMax,pg,\[Lambda]rSample,\[Lambda]\[Theta]Sample,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]\[Theta],\[Phi]C,tC,zRoots,tInit,rInit,\[Theta]Init,\[Phi]Init},
	{En,L,Q} = KerrGeoConstantsOfMotion[a,p,e,x];
	
	(* Useful constants for \[Theta]-dependent calculations *)
	\[Beta]=a^2(1-En^2);
	\[Alpha]=L^2+Q+\[Beta];
	zp=Sqrt[(\[Alpha]+Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zm=Sqrt[(\[Alpha]-Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zRoots={ArcCos[zm],Pi-ArcCos[zm]};
	
	(* Useful constants for r-dependent calculations *)
	{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, x, En, Q];
	
	(* Mino frequencies of orbit. I call Re, because some frequencies 
	are given as imaginary for restricted orbits. 
	Maybe something to fix with KerrGeoMinoFrequencies*)
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = KerrGeoMinoFrequencies[a,p,e,x];
	If[e>0&&(Im[\[CapitalUpsilon]r]!=0||Re[\[CapitalUpsilon]r]==0),Print["Unstable orbit. Aborting."]; Abort[]];
	If[r2<1+Sqrt[1-a^2],Print["Unstable orbit. Aborting."]; Abort[]];
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Re[{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t}];
	
	(*Parameterize r and \[Theta] in terms of Darwin-like parameters \[Psi] and \[Chi]*)
	r0[psi_]:=p M/(1+e Cos[psi]);
	\[Theta]0[chi_]:=ArcCos[zm Cos[chi]];
	
	(*Precision of sampling depends on precision of arguments*)
	pg=Min[{Precision[{a,p,e,x}],Precision[initPhases]}];
	If[pg==$MachinePrecision,
		\[Psi]r[Nr_]:=N[Table[2Pi i/Nr,{i,0,Nr-1}]];
		\[Chi]\[Theta][Nth_]:=N[Table[2Pi i/Nth,{i,0,Nth-1}]],
		\[Psi]r[Nr_]:=N[Table[2Pi i/Nr,{i,0,Nr-1}],1.3pg];
		\[Chi]\[Theta][Nth_]:=N[Table[2Pi i/Nth,{i,0,Nth-1}],1.3pg]
	];
	
	(*Solve for Mino time as a function of \[Psi] and \[Chi]*)
	\[CapitalDelta]\[Lambda]r=MinoRFastSpec[a,p,e,x];
	\[Lambda]r[psi_]:=\[Lambda]r[psi]=\[CapitalDelta]\[Lambda]r[psi];
	\[Lambda]rSample[Nr_]:=\[Lambda]r[\[Psi]r[Nr]];
	\[Lambda]\[Theta]=MinoThetaFastSpec[a,p,e,x];
	\[Lambda]\[Theta][chi_]:=\[Lambda]\[Theta][chi]=\[CapitalDelta]\[Lambda]\[Theta][chi];
	\[Lambda]\[Theta]Sample[Nth_]:=\[Lambda]\[Theta][\[Chi]\[Theta][Nth]];
	
	(*Find the inverse transformation of \[Psi] and \[Chi] as functions of \[Lambda] 
	using spectral integration*)
	If[e>0,
		\[Psi]=PhaseOfMinoFastSpec[\[CapitalUpsilon]r,\[Lambda]rSample],
		\[Psi][\[Lambda]_]:=0
	];
	If[x^2<1,
		\[Chi]=PhaseOfMinoFastSpec[\[CapitalUpsilon]\[Theta],\[Lambda]\[Theta]Sample],
		\[Chi][\[Lambda]_]:=0
	];

	(*Spectral integration of t and \[Phi] as functions of \[Lambda]*)
	If[e>0,
		\[CapitalDelta]tr=TimeOfMinoFastSpecR[a,p,e,x,{\[CapitalUpsilon]r,\[Lambda]rSample}];
		\[CapitalDelta]\[Phi]r=PhiOfMinoFastSpecR[a,p,e,x,{\[CapitalUpsilon]r,\[Lambda]rSample}],
		\[CapitalDelta]tr[\[Lambda]_?NumericQ]:=0;
		\[CapitalDelta]tr[\[Lambda]_List]:=\[CapitalDelta]tr[#]&/@\[Lambda];
		\[CapitalDelta]\[Phi]r[\[Lambda]_?NumericQ]:=0;
		\[CapitalDelta]\[Phi]r[\[Lambda]_List]:=\[CapitalDelta]\[Phi]r[#]&/@\[Lambda];
	];
	If[x^2<1,
		(* Calculate theta dependence for non-equatorial orbits *)
		\[CapitalDelta]t\[Theta]=TimeOfMinoFastSpecTheta[a,p,e,x,{\[CapitalUpsilon]\[Theta],\[Lambda]\[Theta]Sample}];
		\[CapitalDelta]\[Phi]\[Theta]=PhiOfMinoFastSpecTheta[a,p,e,x,{\[CapitalUpsilon]\[Theta],\[Lambda]\[Theta]Sample}],
		(* No theta dependence for equatorial orbits *)
		\[CapitalDelta]t\[Theta][\[Lambda]_?NumericQ]:=0;
		\[CapitalDelta]t\[Theta][\[Lambda]_List]:=\[CapitalDelta]t\[Theta][#]&/@\[Lambda];
		\[CapitalDelta]\[Phi]\[Theta][\[Lambda]_?NumericQ]:=0;
		\[CapitalDelta]\[Phi]\[Theta][\[Lambda]_List]:=\[CapitalDelta]\[Phi]\[Theta] [#]&/@\[Lambda];
	];

	(*Collect initial Mino time phases*)
	{qt0, qr0, q\[Theta]0, q\[Phi]0} = {initPhases[[1]], initPhases[[2]], initPhases[[3]], initPhases[[4]]};
	
	(* If the user specifies a valid set of initial coordinate positions, 
	find phases that give these initial positions *)
	{tInit,rInit,\[Theta]Init,\[Phi]Init}=OptionValue[InitialPosition];
	If[{tInit,rInit,\[Theta]Init,\[Phi]Init}!={0,0,0,0},
		If[rInit<=r1&&rInit>=r2&&\[Theta]Init>=zRoots[[1]]&&\[Theta]Init<=zRoots[[2]],
			If[e==0,\[Psi]0=0,\[Psi]0=ArcCos[p M/(e rInit)-1/e]];
			If[zm==0,\[Chi]0=0,\[Chi]0=ArcCos[Cos[\[Theta]Init]/zm]];
			qt0=tInit;
			qr0=\[CapitalUpsilon]r*\[Lambda]r[\[Psi]0]+\[Psi]0;
			q\[Theta]0=\[CapitalUpsilon]\[Theta]*\[Lambda]\[Theta][\[Chi]0]+\[Chi]0;
			q\[Phi]0=\[Phi]Init,
			Print["{t,r,\[Theta],\[Phi]} = "<>ToString[{tInit,rInit,\[Theta]Init,\[Phi]Init}]<>" is not a valid set of initial coordinates"]
		]
	];
	If[\[CapitalUpsilon]r==0,\[Lambda]r0=0,\[Lambda]r0=qr0/\[CapitalUpsilon]r];
	If[\[CapitalUpsilon]\[Theta]==0,\[Lambda]\[Theta]0=0,\[Lambda]\[Theta]0=q\[Theta]0/\[CapitalUpsilon]\[Theta]];
	
	(* Find integration constants for t and \[Phi], so that t(\[Lambda]=0)=qt0 and \[Phi](\[Lambda]=0)=q\[Phi]0 *)
	\[Phi]C=\[CapitalDelta]\[Phi]r[\[Lambda]r0]+\[CapitalDelta]\[Phi]\[Theta][\[Lambda]\[Theta]0];
	tC=\[CapitalDelta]tr[\[Lambda]r0]+\[CapitalDelta]t\[Theta][\[Lambda]\[Theta]0];

	t[\[Lambda]_]:=\[CapitalDelta]tr[\[Lambda]+\[Lambda]r0]+\[CapitalDelta]t\[Theta][\[Lambda]+\[Lambda]\[Theta]0]+\[CapitalUpsilon]t \[Lambda]+qt0-tC;
	r[\[Lambda]_]:=r0[\[Psi][\[Lambda]+\[Lambda]r0]+\[CapitalUpsilon]r \[Lambda]+qr0];
	\[Theta][\[Lambda]_]:=\[Theta]0[\[Chi][\[Lambda]+\[Lambda]\[Theta]0]+\[CapitalUpsilon]\[Theta] \[Lambda]+q\[Theta]0];
	\[Phi][\[Lambda]_]:=\[CapitalDelta]\[Phi]r[\[Lambda]+\[Lambda]r0]+\[CapitalDelta]\[Phi]\[Theta][\[Lambda]+\[Lambda]\[Theta]0]+\[CapitalUpsilon]\[Phi] \[Lambda]+q\[Phi]0-\[Phi]C;

	assoc = Association[
		"Parametrization"->"Mino",
		"Energy" -> En, 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> {En,L,Q},
		"RadialFrequency" -> \[CapitalUpsilon]r,
		"PolarFrequency" ->  \[CapitalUpsilon]\[Theta],
		"AzimuthalFrequency" -> \[CapitalUpsilon]\[Phi],
		"Frequencies" -> {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"RadialRoots" -> {r1,r2,r3,r4},
		"PolarRoots" -> zRoots
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


Options[KerrGeoOrbit] = {Parametrization -> "Mino", Method -> "FastSpec"}
SyntaxInformation[KerrGeoOrbit] = {"ArgumentsPattern"->{_,_,OptionsPattern[]}};


KerrGeoOrbit[a_,p_,e_,x_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, method},
(*FIXME: add stability check but make it possible to turn it off*)

method = OptionValue["Method"];
param = OptionValue["Parametrization"];

If[param == "Darwin" && Abs[x]!=1, Print["Darwin parameterization only valid for equatorial motion"];Return[];];

If[Precision[{a,p,e,x}] > 30, method = "Analytic"];

If[method == "FastSpec",

	If[param == "Mino", Return[KerrGeoOrbitFastSpec[a, p, e, x, initPhases]]];
	If[param == "Darwin", 
		If[a==0,Return[KerrGeoOrbitSchwarzDarwin[p, e]], Return[KerrGeoOrbitFastSpecDarwin[a,p,e,x,initPhases]]]
	];
	Print["Unrecognized parametrization: " <> OptionValue["Parametrization"]];
	
];

If[method == "Analytic",

	If[param == "Mino", Return[KerrGeoOrbitMino[a, p, e, x, initPhases]]];
	If[param == "Darwin", 
		If[a==0,Return[KerrGeoOrbitSchwarzDarwin[p, e]], Return[KerrGeoOrbitDarwin[a,p,e,x,initPhases]]]
	];
	Print["Unrecognized parametrization: " <> OptionValue["Parametrization"]];

];

Print["Unrecognized method: " <> method];

]


Format[KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_]] := "KerrGeoOrbitFunction["<>ToString[a]<>","<>ToString[p]<>","<>ToString[e]<>","<>ToString[N[x]]<>",<<>>]";
KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_][y_?StringQ] := assoc[y]


(* ::Chapter::Closed:: *)
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


KerrGeoPhotonSphereRadius[0,x_]:=3


(* ::Text:: *)
(*Radius of photon sphere  for equatorial orbits from Bardeen, Press, Teukolsky ApJ, 178, p347 (1972), Eq. 2.18*)


KerrGeoPhotonSphereRadius[a_,1]:=2(1+Cos[2/3 ArcCos[-a]])
KerrGeoPhotonSphereRadius[a_,-1]:=2(1+Cos[2/3 ArcCos[a]])


(* ::Text:: *)
(*For polar orbits the radius was given by E. Teo, General Relativity and Gravitation, v. 35, Issue 11, p. 1909-1926 (2003), Eq. (14)*)


KerrGeoPhotonSphereRadius[a_,0]:=1+2Sqrt[1-1/3 a^2]Cos[1/3 ArcCos[(1-a^2)/(1-1/3 a^2)^(3/2)]]


(* ::Text:: *)
(*In the extremal limit we can find the photon sphere radius exactly*)


KerrGeoPhotonSphereRadius[1,x_]:=If[x < Sqrt[3]-1,1+Sqrt[2] Sqrt[1-x]-x,1];;


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


KerrGeoIBSO[a1_?NumericQ,x1_?NumericQ]/;Precision[{a1,x1}]!=\[Infinity]:=Block[{a=a1,x=x1,rph,prec,n=1,ru,E0},
prec=Precision[{a1,x}];
rph=KerrGeoPhotonSphereRadius[a,x];

E0=KerrGeoEnergy[a,ru,0,x];

While[(E0/.ru->rph+10^-n)<1,n++];
ru/.FindRoot[E0==1,{ru,rph+10^-n,10},WorkingPrecision->Max[MachinePrecision,prec-1]]
]


(* ::Section::Closed:: *)
(*Separatrix*)


(* ::Text:: *)
(*Schwarzschild*)


KerrGeoSeparatrix[0,e_,x_]:= 6+2e;


(* ::Text:: *)
(*From Glampedakis and Kennefick arXiv:gr-qc/0203086, for a=M we have Subscript[p, s]=1+e*)


KerrGeoSeparatrix[1,e_,1]:= 1+e


(* ::Text:: *)
(*Separatrix for equatorial Kerr from Levin and Periz-Giz arXiv:1108.1819*)


KerrGeoSeparatrix[a1_,e_,x_/;Abs[x]==1]:= Module[{ru,a=a1},
	If[x==-1, a = -a];
	ru=ru/.Solve[e==(-ru^2+6 ru-8a ru^(1/2)+3a^2)/(ru^2-2ru+a^2),ru][[-1]];
	(4 ru (ru^(1/2)-a)^2)/(ru^2-2ru+a^2)
]


(* ::Text:: *)
(*Polar ISSO in extremal case found from playing around with the equations*)


KerrGeoSeparatrix[1,0,0]:=1+Sqrt[3]+Sqrt[3+2 Sqrt[3]]


(* ::Text:: *)
(*For e=1 the \!\(TraditionalForm\`*)
(*\*SubscriptBox[\(p\), \(s\)]\)is at 2\!\(TraditionalForm\`*)
(*\*SubscriptBox[\(r\), \(ibso\)]\)*)


KerrGeoSeparatrix[a_,1,x_]:=2KerrGeoIBSO[a,x]


KerrGeoSeparatrix[a1_?NumericQ,e1_?NumericQ,x1_?NumericQ]/;Precision[{a1,e1,x1}]!=\[Infinity]:=Block[{a=a1,ra2,\[Beta],E0,L0,Q0,e2,ru,x=x1,prec,r1,\[Beta]2,p,ru0},

{E0,L0,Q0}=KerrGeoConstantsOfMotion[a,ru,0,x];
\[Beta]=(-1+E0^2);

ra2=2 (-a E0+L0)^2+2 Q0+ru^2 (-1-ru \[Beta]+Sqrt[1+\[Beta] (L0^2+Q0-a^2 \[Beta]-2 ru (1+ru \[Beta]))]);

\[Beta]2=ru (2+ru \[Beta]-2 Sqrt[1+L0^2 \[Beta]+Q0 \[Beta]-2 ru \[Beta]-a^2 \[Beta]^2-2 ru^2 \[Beta]^2]);

e2=(ra2-ru \[Beta]2)/(ra2+ru \[Beta]2);


prec=Precision[{a,e1,x}];
r1=KerrGeoIBSO[a,x];
ru0=ru/.FindRoot[e2==e1,{ru,(r1+10)/2,r1,10},WorkingPrecision->Max[MachinePrecision,prec-1]];

p=(2ra2 ru)/(ra2+ru \[Beta]2)/.ru->ru0
]


(* ::Section::Closed:: *)
(*Innermost stable spherical orbit (ISSO)*)


KerrGeoISSO[a_,x_/;Abs[x]==1]:=KerrGeoISCO[a,x]


KerrGeoISSO[a_,x_]:=KerrGeoSeparatrix[a,0,x]


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
