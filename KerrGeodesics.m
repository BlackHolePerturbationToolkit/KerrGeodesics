(* ::Package:: *)

BeginPackage["KerrGeodesics`"];

KerrGeoELQ::usage = "KerrGeoELQ[a, p, e, \[Theta]inc] returns the energy, z-component of the angular momentum and the Carter constant";
KerrGeoFreqs::usage = "KerrGeoFreqs[a, p, e, \[Theta]inc] returns the radial, polar and azimuthal frequencies and the conversion factor between Boyer-Lindquist and Mino time frequencies";
KerrGeoStableOrbitQ::usage = "KerrGeoStableOrbitQ[a,p,e,\[Theta]inc] checks if given parameters correspond to a stable orbit";
KerrGeoISCO::usage = "KerrGeoISCO[a,\[Theta]inc] computes the location of the inner-most stable circular orbit (ISCO)"
KerrGeoPhotonSphereRadius::usage = "KerrGeoPhotonSphereRadius[a,\[Theta]inc] computes the radius of the photon sphere"
KerrGeoIBSO::usage = "KerrGeoIBSO[a,\[Theta]inc] computes the radius of the inner-most bound spherical orbit (IBSO)"
KerrGeoISSO::usage = "KerrGeoISSO[a,\[Theta]inc] computes the radius of the inner-most stable spherical orbit (ISSO)"
KerrGeoSeparatrix::usage = "KerrGeoSepatrix[a,e,\[Theta]inc] calculates the value of p at the sepatrix between stable and plunging/scattered orbits*)"
KerrGeoOrbit::usage = "KerrGeoOrbit[a,p,e,\[Theta]inc] calculates the orbital trajectory in Boyer-Lindquist coordinates"
KerrGeoOrbitFunction::usage = "KerrGeoOrbitFunction[a,p,e,\[Theta]inc,data] function containing information relating to a specific Kerr orbit"
KerrGeoOrbitFunction2::usage = "KerrGeoOrbitFunction[a,p,e,\[Theta]inc,data] function containing information relating to a specific Kerr orbit"


Begin["`Private`"];

(*ELQ calculation for Schwarzschild spacetime*)
(*Cutler, Kennefick and Poisson, Phys. Rev. D, 50, 6, p3816, (1994)*)
(*Eqs. 2.5 and 2.6*)
KerrGeoELQ[(0|0.0),p_,e_,\[Theta]inc_]:=Module[{M=1,E0,L0,Q},
 E0=Sqrt[((p-2-2e)(p-2+2e))/(p(p-3-e^2))];
 L0=p/Sqrt[p/M-e^2-3];
 {E0,Cos[\[Theta]inc]L0,Sin[\[Theta]inc]^2 L0^2}
]

(*ELQ calculation for Kerr circular equatorial - prograde*)
(*Bardeen, Press, Teukolsky ApJ, 178, p347 (1972)*)
(*Eqs. 2.12 and 2.13*)
KerrGeoELQ[a_,p_,0,0]:=Module[{M=1,E0,L0,Q},
E0 = (p^(3/2)-2M p^(1/2)+a M^(1/2))/(p^(3/4) (p^(3/2)-3M p^(1/2)+2a M^(1/2))^(1/2));
L0 = (M^(1/2) (p^2-2a M^(1/2) p^(1/2)+a^2))/(p^(3/4) (p^(3/2)-3M p^(1/2)+2a M^(1/2))^(1/2));
Q = 0;
{E0,L0,Q}
]

(*ELQ calculation for Kerr circular equatorial - retrograde*)
(*Bardeen, Press, Teukolsky ApJ, 178, p347 (1972)*)
(*Eqs. 2.12 and 2.13*)
KerrGeoELQ[a_,p_,0,\[Pi]]:=Module[{M=1,E0,L0,Q},
E0 = (p^(3/2)-2M p^(1/2)-a M^(1/2))/(p^(3/4) (p^(3/2)-3M p^(1/2)-2a M^(1/2))^(1/2));
L0 = (-M^(1/2)(p^2+2a M^(1/2) p^(1/2)+a^2))/(p^(3/4) (p^(3/2)-3M p^(1/2)-2a M^(1/2))^(1/2));
Q = 0;
{E0,L0,Q}
]

(*ELQ calculation for Kerr eccentric equatorial*)
(*Glampedakis and Kennefick, Phys. Rev. D66 (2002) 044002, arXiv:gr-qc/0203086*)
(*Eq. 7 and appendix A*)
KerrGeoELQ[a1_,p_,e_,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]==0]:=Module[{M=1,a=a1,F,G,B,\[CapitalDelta]1,x,E0,L0,Q=0,sign=1},

If[Mod[\[Theta]inc,2\[Pi]]==\[Pi], a=-a;sign=-1];

F = 1/p^3 (p^3-2M(3+e^2)p^2+M^2 (3+e^2)^2 p-4M a^2 (1-e^2)^2);
G=2/p (-M p^2+(M^2 (3+e^2)-a^2)p-M a^2 (1+3e^2));
B=(a^2-M p)^2;
\[CapitalDelta]1=G^2-4F B;
x=Sqrt[(-G-sign Sqrt[\[CapitalDelta]1])/(2 F)];

E0=(1-M/p (1-e^2)(1-x^2/p^2 (1-e^2)))^(1/2);
L0=a E0 + x;

{E0,sign L0,Q}

]


(*ELQ calculation for spherical polar orbits*)
(*Stoghianidis & Tsoubelis, Gen. Rel, Gravitation, vol. 19, No. 12, p. 1235 (1987)*)
(*Eqs. 17-19*)
KerrGeoELQ[a_,p_,0,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]==\[Pi]/2]:=Module[{M=1,r0,\[CapitalDelta]0,Z,En,Q},
		r0 = p;
		\[CapitalDelta]0 = r0^2-2M r0+a^2;
		Z = r0^3 - 3M r0^2 + a^2 r0+M a^2;
		En =Sqrt[(r0 \[CapitalDelta]0^2)/((r0^2+a^2)Z)];
		Q = r0 (M r0^3+a^2 r0^2-3M a^2 r0 +a^4)/Z - a^2 En^2;
		{En,0,Q}
]

(*ELQ calculation for generic orbits in Kerr*)
(*W. Schmidt, Class. Quant. Grav. 19 (2002) 2743, arXiv:gr-qc/0202090*)
(*Appendix B*)
KerrGeoELQ[a_(*/;Abs[a]<=1*), p_, e_, \[Theta]inc1_?NumericQ] := Module[{M=1,f, g, h, d, fp, gp, hp, dp, r, rp, ra, zm, \[CapitalDelta], \[Rho], \[Kappa], \[Epsilon], \[Eta], \[Sigma], En, L, Q, E1, Em1, f1, g1, h1, d1, f2, g2, h2, d2, L1, L2,r0,\[CapitalDelta]0,Z,\[Theta]min,\[Theta]inc=\[Theta]inc1},

\[Theta]inc=Mod[\[Theta]inc,2\[Pi]];
If[\[Theta]inc>\[Pi], \[Theta]inc=2\[Pi]-\[Theta]inc];
If[\[Theta]inc <= \[Pi]/2 , \[Theta]min = \[Pi]/2-\[Theta]inc, \[Theta]min=-\[Pi]/2+\[Theta]inc];

If[Mod[\[Theta]inc,\[Pi]]==\[Pi]/2 && e!=0,Print["Polar, non-spherical orbits not yet implemented"]; Return[]];

 rp = p/(1 + e);
 ra = p/(1 - e);
 
 zm = Cos[\[Theta]min];
 
 
 \[CapitalDelta][r_] = r^2 - 2 r + a^2;
 
 
 f[r_] = r^4 + a^2 (r (r + 2) + zm^2 \[CapitalDelta][r]);
 g[r_] = 2 a r;
 h[r_] = r (r - 2) + zm^2/(1 - zm^2) \[CapitalDelta][r];
 d[r_] = (r^2 + a^2 zm^2) \[CapitalDelta][r];
 
 fp[r_] = 4 r^3 + 2 a^2 ((1 + zm^2) r + (1 - zm^2));
 gp[r_] = 2 a;
 hp[r_] = 2 (r - 1)/(1 - zm^2);
 dp[r_] = 2 (2 r - 3) r^2 + 2 a^2 ((1 + zm^2) r - zm^2);
 
 {f1, g1, h1, d1} = {f[rp], g[rp], h[rp], d[rp]};
 
 If[e != 0 || !NumericQ[e],
  {f2, g2, h2, d2} = {f[ra], g[ra], h[ra], d[ra]}
  ,
  {f2, g2, h2, d2} = {fp[ra], gp[ra], hp[ra], dp[ra]}
  ];
 
 \[Kappa] = d1 h2 - h1 d2;
 \[Epsilon] = d1 g2 - g1 d2;
 \[Rho] = f1 h2 - h1 f2;
 \[Eta] = f1 g2 - g1 f2;
 \[Sigma] = g1 h2 - h1 g2;
 
 En[pm_] := Sqrt[(\[Kappa] \[Rho] + 2 \[Epsilon] \[Sigma] + pm 2 Sqrt[\[Sigma] (\[Sigma] \[Epsilon]^2 + \[Rho] \[Epsilon] \[Kappa] - \[Eta] \[Kappa]^2)])/(\[Rho]^2 + 4 \[Eta] \[Sigma])];
 E1 = En[1];
 Em1 = En[-1];
 
 L1[En_] := -((En g1 + Sqrt[-d1 h1 + En^2 (g1^2 + f1 h1)])/h1);
 L2[En_] := (-En g1 + Sqrt[-d1 h1 + En^2 (g1^2 + f1 h1)])/h1;
 Q[En_, L_] := zm^2 (a^2 (1 - En^2) + L^2/(1 - zm^2));
 

 If[ \[Theta]inc <= \[Pi]/2,
    Return[{Em1, L2[Em1], Q[Em1, L2[Em1]]}],
	Return[{E1, L1[E1], Q[E1, L1[E1]]}]
  ]
]


(*Orbital frequencies for circular orbits in Schwarzschild spacetime*)
KerrGeoFreqs[(0|0.),p_,(0|0.),\[Theta]inc_]:=Module[{M=1},
	{Sqrt[((p-6M)M)/p^4], (M/p^3)^(1/2), (M/p^3)^(1/2), Sqrt[p^5/(p-3)]}
]

(*Orbital frequencies for eccentric orbits in Schwarzschild spacetime*)
KerrGeoFreqs[(0|0.),p_,e_,\[Theta]inc_]:=Module[{En,L,Q,r1,r2,r3,kr,\[CapitalUpsilon]r,\[CapitalGamma],hp,hr,M=1},
{En,L,Q}=KerrGeoELQ[0,p,e,0];

r1=p/(1-e);
r2=p/(1+e);
r3=(2M)/(1-En^2)-(r1+r2);(*Eq. (11)*)


kr=(r1-r2)/(r1-r3) r3/r2;(*Eq.(13)*)
\[CapitalUpsilon]r=(Pi Sqrt[(1-En^2)(r1-r3)r2])/(2EllipticK[kr]);(*Eq.(15)*)
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-2M))/((r1-r3)(r2-2M));
(*Convert to frequencies w.r.t BL time Subscript[\[CapitalOmega], k]=Subscript[\[Gamma], k]/\[CapitalGamma] using Fujita and Hikida's formula Eq. (21)*)
\[CapitalGamma]=4M^2 En +En/2 ((r3(r1+r2+r3)-r1 r2)+(r2-r3)(r1+r2+r3) EllipticPi[hr,kr]/EllipticK[kr]+(r1-r3)r2 EllipticE[kr]/EllipticK[kr])+2M En(r3 +(r2-r3) EllipticPi[hr,kr]/EllipticK[kr])+((8M^3 En)/(r3-2M) (1-(r2-r3)/(r2-2M) EllipticPi[hp,kr]/EllipticK[kr]));

{\[CapitalUpsilon]r/\[CapitalGamma],L/\[CapitalGamma],L/\[CapitalGamma],\[CapitalGamma]}

]


(*Calculate the orbital frequencies, Subscript[\[CapitalOmega], \[Alpha]], w.r.t Boyer-Lindquist time and the conversion factor, \[CapitalGamma], to frequencies w.r.t. Mino time, Subscript[\[CapitalUpsilon], \[Alpha]]=\[CapitalGamma] Subscript[\[CapitalOmega], \[Alpha]]*)
(*Fujita and Hikida, Class. Quantum Grav.26 (2009) 135002, arXiv:0906.1420*)
KerrGeoFreqs[a_/;Abs[a]<1,p_,e_,\[Theta]inc1_?NumericQ]:=Module[{M=1,En,L,Q,r1,r2,AplusB,AB,r3,r4,\[Epsilon]0,zm,kr,k\[Theta],\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalGamma],rp,rm,hp,hm,hr,EnLQ,a2zp,\[Epsilon]0zp,zmOverZp,\[Theta]min,\[Theta]inc=\[Theta]inc1},

\[Theta]inc=Mod[\[Theta]inc,2\[Pi]];
If[\[Theta]inc>\[Pi], \[Theta]inc=2\[Pi]-\[Theta]inc];

If[\[Theta]inc==\[Pi]/2, Print["Equations for polar orbits not implemented yet"];Return[];];

{En,L,Q}=KerrGeoELQ[a,p,e,\[Theta]inc];
\[Theta]min=(\[Pi]/2-\[Theta]inc)/Sign[L];

r1=p/(1-e);
r2=p/(1+e);
AplusB=(2M)/(1-En^2)-(r1+r2);(*Eq. (11)*)
AB=(a^2 Q)/((1-En^2)r1 r2);(*Eq. (11)*)
r3=(AplusB+Sqrt[(AplusB)^2-4AB])/2;(*Eq. (11)*)
r4=AB/r3;(*Eq. (11)*)
\[Epsilon]0=a^2 (1-En^2)/L^2;
zm=Cos[\[Theta]min]^2;
a2zp=(L^2+a^2 (-1+En^2) (-1+zm))/( (-1+En^2) (-1+zm));

\[Epsilon]0zp=-((L^2+a^2 (-1+En^2) (-1+zm))/(L^2 (-1+zm)));
zmOverZp=If[a==0,0,zm/((L^2+a^2 (-1+En^2) (-1+zm))/(a^2 (-1+En^2) (-1+zm)))];
kr=Sqrt[(r1-r2)/(r1-r3) (r3-r4)/(r2-r4)];(*Eq.(13)*)
k\[Theta]=Sqrt[zmOverZp];(*Eq.(13)*)
\[CapitalUpsilon]r=(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)])/(2EllipticK[kr^2]);(*Eq.(15)*)
\[CapitalUpsilon]\[Theta]=(Pi L Sqrt[\[Epsilon]0zp])/(2EllipticK[k\[Theta]^2]);(*Eq.(15)*)

rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));

\[CapitalUpsilon]\[Phi]=(2\[CapitalUpsilon]\[Theta])/(Pi Sqrt[\[Epsilon]0zp]) EllipticPi[zm,k\[Theta]^2]+(2a \[CapitalUpsilon]r)/(Pi(rp-rm)Sqrt[(1-En^2)(r1-r3)(r2-r4)])(*Eq. (21)*)((2M En rp - a L)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-(2M En rm - a L)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2]));

(*Convert to frequencies w.r.t BL time Subscript[\[CapitalOmega], k]=Subscript[\[Gamma], k]/\[CapitalGamma] using Fujita and Hikida's formula Eq. (21)*)
\[CapitalGamma]=4M^2 En + (2a2zp En  \[CapitalUpsilon]\[Theta])/(Pi L Sqrt[\[Epsilon]0zp]) (EllipticK[k\[Theta]^2]- EllipticE[k\[Theta]^2]) + (2\[CapitalUpsilon]r)/(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)]) (En/2 ((r3(r1+r2+r3)-r1 r2)EllipticK[kr^2]+(r2-r3)(r1+r2+r3+r4)EllipticPi[hr,kr^2]+(r1-r3)(r2-r4)EllipticE[kr^2])+2M En(r3 EllipticK[kr^2]+(r2-r3)EllipticPi[hr,kr^2])+(2M)/(rp-rm) (((4M^2 En-a L)rp-2M a^2 En)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-((4M^2 En-a L)rm-2M a^2 En)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2])));

(*Output the BL frequencies by dividing the Mino time frequencies by the conversion factor \[CapitalGamma]*)
{\[CapitalUpsilon]r/\[CapitalGamma],Abs[\[CapitalUpsilon]\[Theta]/\[CapitalGamma]],\[CapitalUpsilon]\[Phi]/\[CapitalGamma],\[CapitalGamma]}
]

KerrGeoFreqs[a_/;a==1,p_,e_,\[Theta]inc1_?NumericQ]:=Module[{},
	Print["Frequency calculation not yet implemented for a=M (but the equations are in the appendix of Fujita and Hikida so please implement them!)"];
]



(*Orbit stability check, Schwarzschild case*)
KerrGeoStableOrbitQ[(0|0.0),p_,e_,\[Theta]inc_]:=Module[{},
  If[p>6+2e,True,False](*> not \[GreaterEqual] as orbits along the separatrix are marginally stable*)
];

(* The orbit is stable only if Subscript[\[CapitalOmega], r] is real*)
KerrGeoStableOrbitQ[a_?NumericQ/;a!=0,p_?NumericQ,e_?NumericQ,\[Theta]inc_?NumericQ]:=Module[{ps},
  ps=KerrGeoSeparatrix[a,e,\[Theta]inc];
  If[p<ps,False,True]
];


(*Kerr inner-most circular orbit ISCO*)
(*Bardeen, Press, Teukolsky ApJ, 178, p347 (1972)*)
(*Eq. 2.21*)
KerrGeoISCO[a_,\[Theta]inc1_/;Mod[\[Theta]inc1,\[Pi]]==0]:=Module[{M=1,\[Theta]inc=\[Theta]inc1,Z1,Z2},
    \[Theta]inc=Mod[\[Theta]inc,2\[Pi]];
	Z1=1+(1-a^2/M^2)^(1/3) ((1+a/M)^(1/3)+(1-a/M)^(1/3));
	Z2=(3a^2/M^2 + Z1^2)^(1/2);
	If[\[Theta]inc==0,
		Return[M(3+Z2-((3-Z1)(3+Z1+2Z2))^(1/2))],
		Return[M(3+Z2+((3-Z1)(3+Z1+2Z2))^(1/2))]
	];
];


(*Photon sphere becomes light ring at r=3M in Schwarzschild*)
KerrGeoPhotonSphereRadius[(0|0.0),\[Theta]inc_]:=3;

(*Radius of photon sphere  for equatorial orbits*)
(*Bardeen, Press, Teukolsky ApJ, 178, p347 (1972)*)
(*Eq. 2.18*)
KerrGeoPhotonSphereRadius[a_,(0|0.)]:=2(1+Cos[2/3 ArcCos[-a]])
KerrGeoPhotonSphereRadius[a_,\[Pi]]:=2(1+Cos[2/3 ArcCos[a]])

(*Photon sphere radius is where the energy for a timelike orbit diverges*)
KerrGeoPhotonSphereRadius[a_,\[Theta]inc_]:= Module[{res},
	res=Sort[Re[p/.Solve[Simplify[Denominator[KerrGeoELQ[a,p,0,\[Theta]inc][[1]]^2]==0,Assumptions->{a>=0,p>0}],p]],Greater];
	If[\[Theta]inc>=\[Pi]/2,Return[res[[1]]],Return[res[[2]]]]
]

(* The IBSO is where E=1. Use the fact that rph < r_ibso*)
KerrGeoIBSO[a_?NumericQ,\[Theta]inc_?NumericQ]:= Module[{rph},
	rph=KerrGeoPhotonSphereRadius[a,\[Theta]inc];
	p/.FindRoot[KerrGeoELQ[a,p,0,\[Theta]inc][[1]]-1,{p,rph+10^-10,9}, Method->"Brent"]
]

(*Equatorial IBSO results from Bardeen, Press, Teukolsky 1972*)
KerrGeoIBSO[a_,0]:= 2-a+2(1-a)^(1/2)
KerrGeoIBSO[a_,\[Pi]]:= 2+a+2(1+a)^(1/2)

(*Returns the roots of the radial equation and r2-r3 *)
KerrGeoRadialEqRoots[a_,p_,e_,\[Theta]inc_]:=Module[{M=1,En,L,Q,AplusB,AB,r1,r2,r3,r4},
	{En,L,Q}=KerrGeoELQ[a,p,e,\[Theta]inc];

	r1=p/(1-e);
	r2=p/(1+e);
	AplusB=(2M)/(1-En^2)-(r1+r2);(*Eq. (11)*)
	AB=(a^2 Q)/((1-En^2)r1 r2);(*Eq. (11)*)
	r3=(AplusB+Sqrt[(AplusB)^2-4AB])/2;(*Eq. (11)*)
	r4=AB/r3;(*Eq. (11)*)

	{r1,r2,r3,r4,r2-r3}
]

(*The ISSO occurs when the r2-r3=0. Also use r_ibso < r_isso*)
KerrGeoISSO[a_?NumericQ,\[Theta]inc_?NumericQ]:=Module[{rmb},
	rmb=KerrGeoIBSO[a,\[Theta]inc];
	p/.FindRoot[KerrGeoRadialEqRoots[a,p,0,\[Theta]inc][[5]],{p,rmb+10^-10,9},Method->"Brent"]
]

(*For equatorial orbits the ISSO is the ISCO*)
KerrGeoISSO[a_,0]:= KerrGeoISCO[a,0]
KerrGeoISSO[a_,\[Pi]]:= KerrGeoISCO[a,\[Pi]]



(*Separatrix for Schwarzschild*)
KerrGeoSeparatrix[0,e_,\[Theta]inc_]:= 6+2e;

(*Separatrix for equatorial Kerr from Levin and Periz-Giz arXiv:1108.1819*)
KerrGeoSeparatrix[a1_,e_,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]==0]:= Module[{ru,a=a1},
	If[Mod[\[Theta]inc,2\[Pi]] == \[Pi], a = -a];
	ru=ru/.Solve[e==(-ru^2+6 ru-8a ru^(1/2)+3a^2)/(ru^2-2ru+a^2),ru][[-1]];
	(4 ru (ru^(1/2)-a)^2)/(ru^2-2ru+a^2)
]

(* The separatrix occurs when Subscript[\[CapitalOmega], r]=0 which implies Subscript[r, 2]=Subscript[r, 3]. Search for the value of p at which this occurs which must be outside the p of the ISSO*)
(*FIXME: values near the equatorial plane seem to be problematic*)
KerrGeoSeparatrix[a_?NumericQ,e_?NumericQ,\[Theta]inc_?NumericQ]:=Module[{},
	p/.FindRoot[KerrGeoRadialEqRoots[a,p,e,\[Theta]inc][[5]],{p,KerrGeoISSO[a,\[Theta]inc],20},Method->"Brent"]
]

(*From Glampedakis and Kennefick arXiv:gr-qc/0203086*)
KerrGeoSeparatrix[1,e_,0]:=1+e




(*Calculate the trajectory and properties of an eccentric, equatorial Kerr orbit*)
(*Differential form of trajectory equations from Glampedakis and Kennefick, Phys. Rev. D66 (2002) 044002, arXiv:gr-qc/0203086*)
(*Integration of trajectory equations using method from Hopper et al., Phys. Rev. D 92, 044048 (2015), arXiv:1506.04742*)
Options[KerrGeoOrbit] = {"MaxIterations"->10}
KerrGeoOrbit[a1_?NumericQ,p_?NumericQ,e_?NumericQ,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]==0,OptionsPattern[]]:=Module[{M=1,a=a1,ELQ,freqs,F,G,B,\[CapitalDelta]1,x,\[ScriptCapitalE]0,\[ScriptCapitalL]0,Vr,Vt,V\[Phi],J,dtd\[Chi],d\[Phi]d\[Chi],\[ScriptCapitalN],\[Chi]k,dtd\[Chi]k,d\[Phi]d\[Chi]k,\[ScriptCapitalG]tn,\[ScriptCapitalG]\[Phi]n,t,\[Phi],\[Chi],\[Chi]k2,tk,rk,\[Phi]k,r\[Chi],\[CapitalDelta]\[Phi],Tr,tI,\[Phi]I,rI,\[CapitalDelta]\[ScriptCapitalN],estPrec,jmax,d\[Tau]d\[Chi],u,ut,ur,u\[Theta],u\[Phi],drd\[Chi],xp},
If[Precision[{a1,p,e,\[Theta]inc}]==Infinity, 
	Print["Cannot get infinite precision orbit trajectory, will work to MachinePrecision. Specify specific precision values for input if that precision is sought"];
	a=N[a];
];

ELQ=KerrGeoELQ[a,p,e,\[Theta]inc];
{\[ScriptCapitalE]0,\[ScriptCapitalL]0}=ELQ[[1;;2]];
If[Mod[\[Theta]inc,2\[Pi]] == \[Pi], a = -a];

F = 1/p^3 (p^3-2M(3+e^2)p^2+M^2 (3+e^2)^2 p-4M a^2 (1-e^2)^2);
G=2/p (-M p^2+(M^2 (3+e^2)-a^2)p-M a^2 (1+3e^2));
B=(a^2-M p)^2;
\[CapitalDelta]1=G^2-4F B;
x=Sqrt[(-G-Sign[\[ScriptCapitalL]0]Sqrt[\[CapitalDelta]1])/(2 F)];

Vr[\[Chi]_]:= x^2+a^2+2a x \[ScriptCapitalE]0-2M x^2/p (3+e Cos[\[Chi]]);
Vt[\[Chi]_]:= a^2 \[ScriptCapitalE]0-(2a M x)/p (1+e Cos[\[Chi]])+(\[ScriptCapitalE]0 p^2)/(1+e Cos[\[Chi]])^2;
V\[Phi][\[Chi]_]:= x+a \[ScriptCapitalE]0-(2M x)/p (1+e Cos[\[Chi]]);
J[\[Chi]_]:= 1-(2M)/p (1+e Cos[\[Chi]])+a^2/p^2 (1+e Cos[\[Chi]])^2;

dtd\[Chi][\[Chi]_]:=Vt[\[Chi]]/(J[\[Chi]] Vr[\[Chi]]^(1/2));
d\[Phi]d\[Chi][\[Chi]_]:=V\[Phi][\[Chi]]/(J[\[Chi]] Vr[\[Chi]]^(1/2));



d\[Tau]d\[Chi][\[Chi]_]:=p^2/((1+e Cos[\[Chi]])^2 Sqrt[a^2+2 a \[ScriptCapitalE]0 x+x^2-(2 M x^2 (3+e Cos[\[Chi]]))/p]);


(*Initial guesses that should give MachinePrecision. As e\[Rule]1, \[ScriptCapitalN] diverges.*)
Which[e<=0.2, \[ScriptCapitalN]=10, e<=0.5, \[ScriptCapitalN]=20, e<=0.7, \[ScriptCapitalN]=30, e<=0.9, \[ScriptCapitalN]=45, e<=0.95,\[ScriptCapitalN]=70, e>0.95, \[ScriptCapitalN]=100];

jmax=OptionValue["MaxIterations"];
freqs=KerrGeoFreqs[a,p,e,\[Theta]inc];
Table[
\[Chi]k=Table[( \[Pi] k)/(\[ScriptCapitalN]-1),{k, 0, \[ScriptCapitalN]-1}];
{dtd\[Chi]k,d\[Phi]d\[Chi]k}=Transpose[Table[{dtd\[Chi][\[Chi]k[[i]]],d\[Phi]d\[Chi][\[Chi]k[[i]]]},{i,1,Length[\[Chi]k]}]];

\[ScriptCapitalG]tn=FourierDCT[N[dtd\[Chi]k,Precision[{a,p,e}]],1];
\[ScriptCapitalG]\[Phi]n=FourierDCT[N[d\[Phi]d\[Chi]k,Precision[{a,p,e}]],1];

t[\[Chi]_]:=Sqrt[2/(\[ScriptCapitalN]-1)](1/2 \[ScriptCapitalG]tn[[0+1]]\[Chi]+1/2 \[ScriptCapitalG]tn[[\[ScriptCapitalN]-1+1]] Sin[(\[ScriptCapitalN]-1)\[Chi]]/(\[ScriptCapitalN]-1)+Sum[1/n \[ScriptCapitalG]tn[[n+1]]Sin[n \[Chi]],{n,1,\[ScriptCapitalN]-2}]);
\[Phi][\[Chi]_]:=Sqrt[2/(\[ScriptCapitalN]-1)](1/2 \[ScriptCapitalG]\[Phi]n[[0+1]]\[Chi]+1/2 \[ScriptCapitalG]\[Phi]n[[\[ScriptCapitalN]-1+1]] Sin[(\[ScriptCapitalN]-1)\[Chi]]/(\[ScriptCapitalN]-1)+Sum[1/n \[ScriptCapitalG]\[Phi]n[[n+1]]Sin[n \[Chi]],{n,1,\[ScriptCapitalN]-2}]);

\[CapitalDelta]\[Phi] = \[Phi][2\[Pi]];
Tr = t[2\[Pi]];

estPrec=Abs[MantissaExponent[Abs[1-2\[Pi]/freqs[[1]]/Tr]]][[2]];

\[CapitalDelta]\[ScriptCapitalN]=3 Abs[Precision[freqs[[1]]] - estPrec];(*FIXME, make a better e-dependent estimate here*)

(*Print["Freq precision: ", Precision[freqs[[1]]]];
Print["Estimated precision in result: ", estPrec];*)
If[Precision[{a,p,e,\[Theta]inc}] == MachinePrecision || \[CapitalDelta]\[ScriptCapitalN]<3,Return[Null,Table]];
\[ScriptCapitalN]=Ceiling[\[ScriptCapitalN]+\[CapitalDelta]\[ScriptCapitalN]];
If[j==jmax,Print["Failed to reach desired precision within ", jmax, " iterations"]];
,
{j,1,jmax}
];


r\[Chi][\[Chi]_]=p/(1+e Cos[\[Chi]]);
drd\[Chi][\[Chi]_]:=(e p Sin[\[Chi]])/(1+e Cos[\[Chi]])^2;

xp={t,r\[Chi],\[Pi]/2&,\[Phi]};

ut[\[Chi]_]:=dtd\[Chi][\[Chi]]/d\[Tau]d\[Chi][\[Chi]];
ur[\[Chi]_]:=drd\[Chi][\[Chi]]/d\[Tau]d\[Chi][\[Chi]];
u\[Theta][\[Chi]_]:=0;
u\[Phi][\[Chi]_]:=d\[Phi]d\[Chi][\[Chi]]/d\[Tau]d\[Chi][\[Chi]];

u={ut,ur,u\[Theta],u\[Phi]};


KerrGeoOrbitFunction[a,p,e,\[Theta]inc,{\[ScriptCapitalN],ELQ,freqs,xp,u}]

]

(*KerrGeoOrbit[a1_,p_,e_,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]!=0,OptionsPattern[]]:=Module[{},
	Print["Orbit trajectory not yet implemented for non-equatorial orbits"];
]*)


Format[KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]]:=KerrGeoOrbitFunction[a,p,e,\[Theta]inc,"<<>>"];

KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}][\[Chi]_?NumericQ] := Through[x[\[Chi]]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["4Velocity"][\[Chi]_?NumericQ] := Through[u[\[Chi]]]

KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["Energy"]:=ELQ[[1]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["AngMom"]:=ELQ[[2]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["CarterConst"]:=ELQ[[3]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["ELQ"]:=ELQ

KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["Freqs"]:=freqs
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalOmega]r"]:=freqs[[1]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalOmega]\[Theta]"]:=freqs[[2]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalOmega]\[Phi]"]:=freqs[[3]]

KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalUpsilon]r"]:=freqs[[1]]freqs[[4]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalUpsilon]\[Theta]"]:=freqs[[2]]freqs[[4]]
KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["\[CapitalUpsilon]\[Phi]"]:=freqs[[3]]freqs[[4]]

KerrGeoOrbitFunction[a_,p_,e_,\[Theta]inc_, {\[ScriptCapitalN]_,ELQ_,freqs_,x_,u_}]["Interpolation"]:=Module[{\[Chi],\[Chi]k,tk,rk,\[Phi]k,tI,rI,\[Phi]I,Tr,\[CapitalDelta]\[Phi]},

Tr = 2\[Pi]/freqs[[1]];
\[CapitalDelta]\[Phi] = freqs[[3]] Tr;

\[Chi]k=Table[(2\[Pi] k)/\[ScriptCapitalN],{k,0,\[ScriptCapitalN]}];
{tk,rk,\[Phi]k}=Transpose[Table[{x[[1]][\[Chi]]-Tr/(2\[Pi]) \[Chi], x[[2]][\[Chi]], x[[4]][\[Chi]]-\[CapitalDelta]\[Phi]/(2\[Pi]) \[Chi]}/.\[Chi]-> \[Chi]k[[i]],{i,1,Length[\[Chi]k]-1}]];
AppendTo[tk,tk[[1]]]; (*Fix to stop PeriodicInterpolation complaining about the end points not being equal*)
AppendTo[rk,rk[[1]]];
AppendTo[\[Phi]k,\[Phi]k[[1]]];



tI = Interpolation[Transpose[{\[Chi]k,tk}],PeriodicInterpolation->True];
rI = Interpolation[Transpose[{\[Chi]k,rk}],PeriodicInterpolation->True];
\[Phi]I = Interpolation[Transpose[{\[Chi]k,\[Phi]k}],PeriodicInterpolation->True];

t[\[Chi]_]:=Tr/(2\[Pi])\[Chi]+tI[\[Chi]];
r[\[Chi]_]:=rI[\[Chi]];
\[Phi][\[Chi]_]:=\[CapitalDelta]\[Phi]/(2\[Pi])\[Chi] +\[Phi]I[\[Chi]];

{t, r , \[Theta], \[Phi]}

]


(* ::Subsection:: *)
(*Generic orbit implementation*)


(* FIXME this function is a copy of one above, merge the two together *)
KerrGeoRadialRoots[a_, p_, e_, \[Theta]inc_] := Module[{M=1,En,L,Q,r1,r2,r3,r4,AplusB,AB},
{En,L,Q}=KerrGeoELQ[a, p, e, \[Theta]inc];

r1=p/(1-e);
r2=p/(1+e);
AplusB=(2M)/(1-En^2)-(r1+r2);(*Eq. (11)*)
AB=(a^2 Q)/((1-En^2)r1 r2);(*Eq. (11)*)
r3=(AplusB+Sqrt[(AplusB)^2-4AB])/2;(*Eq. (11)*)
r4=AB/r3;

{r1,r2,r3,r4}

]


KerrGeoPolarRoots[a_, p_, e_, \[Theta]inc_] := Module[{En,L,Q,\[Theta]min,zm,zp},
  {En,L,Q} = KerrGeoELQ[a, p, e, \[Theta]inc];
  \[Theta]min=(\[Pi]/2-\[Theta]inc)/Sign[L];
  zm = Cos[\[Theta]min];
  zp = (a^2 (1-En^2)+L^2/(1-zm^2))^(1/2);
  {zp,zm}
]


(*
Generic bound orbits. Based on code from Maarten van de Meent which is an 
implementation of Fujita and Hikida's analytic formula 
Class.Quantum Grav.26 (2009) 135002,arXiv:0906.1420
*)

KerrGeoOrbit[a_, p_, e_, \[Theta]inc_]:=Module[{M=1, r1, r2, r3, r4, kr, k\[Theta], En, L, Q, zp, zm, \[Psi]r, \[Psi]z, rq, zq, rp, rm, hr, hp, hm, \[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], \[CapitalGamma], \[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalUpsilon]\[Phi], \[CapitalUpsilon]t, qr0, qz0, qt0, q\[Phi]0, tr,tz,\[Phi]r,\[Phi]z,t,r,\[Theta],\[Phi]},

{En,L,Q} = KerrGeoELQ[a, p, e, \[Theta]inc];
{r1,r2,r3,r4} = KerrGeoRadialRoots[a, p, e, \[Theta]inc];
{zp,zm} = KerrGeoPolarRoots[a, p, e, \[Theta]inc];
{\[CapitalOmega]r, \[CapitalOmega]\[Theta], \[CapitalOmega]\[Phi], \[CapitalGamma]} = KerrGeoFreqs[a, p, e, \[Theta]inc];
{\[CapitalUpsilon]r, \[CapitalUpsilon]\[Theta], \[CapitalUpsilon]\[Phi], \[CapitalUpsilon]t} = {\[CapitalOmega]r \[CapitalGamma], \[CapitalOmega]\[Theta] \[CapitalGamma], \[CapitalOmega]\[Phi] \[CapitalGamma], \[CapitalGamma]};

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
-2 (r2-r3)/(rp-rm) (
-(1/((-rm+r2) (-rm+r3)))(-2 a^2+rm (4-(a L)/En)) (EllipticPi[hm,kr] qr/\[Pi]-EllipticPi[hm,\[Psi]r[qr],kr] )
+1/((-rp+r2) (-rp+r3)) (-2 a^2+rp (4-(a L)/En)) (EllipticPi[hp,kr] qr/\[Pi]-EllipticPi[hp,\[Psi]r[qr],kr])
)
+(r2-r3) (r1+r2+r3+r4) (EllipticPi[hr,kr] qr/\[Pi]-EllipticPi[hr,\[Psi]r[qr],kr] )
+(r1-r3) (r2-r4) (EllipticE[kr] qr/\[Pi]-EllipticE[\[Psi]r[qr],kr]+hr((Sin[\[Psi]r[qr]]Cos[\[Psi]r[qr]] Sqrt[1-kr Sin[\[Psi]r[qr]]^2])/(1-hr Sin[\[Psi]r[qr]]^2))) );

\[Phi]r[qr_]:= (2 a En (-1/((-rm+r2) (-rm+r3))(2 rm-(a L)/En) (r2-r3) (EllipticPi[hm,kr] qr/\[Pi]-EllipticPi[hm,\[Psi]r[qr],kr])+1/((-rp+r2) (-rp+r3))(2 rp-(a L)/En) (r2-r3) (EllipticPi[hp,kr] qr/\[Pi]-EllipticPi[hp,\[Psi]r[qr],kr] )))/((-rm+rp) Sqrt[(1-En^2) (r1-r3) (r2-r4)]);

\[Psi]z[qz_]:=\[Psi]z[rq]= JacobiAmplitude[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]];
tz[qz_]:= 1/(1-En^2) En zp ( EllipticE[k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticE[\[Psi]z[qz],k\[Theta]]);
\[Phi]z[qz_]:= -1/zp L ( EllipticPi[zm^2,k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticPi[zm^2,\[Psi]z[qz],k\[Theta]]);

(*FIXME; need to be able to give these as input*)
{qt0, qr0, qz0, q\[Phi]0}={0,0,0,0};

t[\[Lambda]_]:= qt0 + \[CapitalUpsilon]t \[Lambda] + tr[\[CapitalUpsilon]r \[Lambda] + qr0] + tz[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0];
r[\[Lambda]_]:= rq[\[CapitalUpsilon]r \[Lambda]+ qr0];
\[Theta][\[Lambda]_]:= ArcCos[zq[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]];
\[Phi][\[Lambda]_]:= q\[Phi]0 + \[CapitalUpsilon]\[Phi] \[Lambda] + \[Phi]r[\[CapitalUpsilon]r \[Lambda]+ qr0] + \[Phi]z[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0];

KerrGeoOrbitFunction2[a, p, e, \[Theta]inc, {t,r,\[Theta],\[Phi]}]


]


Format[KerrGeoOrbitFunction2[a_,p_,e_,\[Theta]inc_, {t_,r_,\[Theta]_,\[Phi]_}]]:=KerrGeoOrbitFunction2[a,p,e,\[Theta]inc,"<<>>"];

KerrGeoOrbitFunction2[a_,p_,e_,\[Theta]inc_, {t_,r_,\[Theta]_,\[Phi]_}][\[Lambda]_]:={t[\[Lambda]], r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]}



End[];

EndPackage[];
