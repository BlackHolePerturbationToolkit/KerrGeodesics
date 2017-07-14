(* ::Package:: *)

BeginPackage["KerrGeodesics`"];

KerrGeoELQ::usage = "KerrGeoELQ[a, p, e, \[Theta]inc] returns the energy, z-component of the angular momentum and the Carter constant";
KerrGeoFreqs::usage = "KerrGeoFreqs[a, p, e, \[Theta]inc] returns the radial, polar and azimuthal frequencies and the conversion factor between Boyer-Lindquist and Mino time frequencies";
KerrGeoStableOrbitQ::usage = "KerrGeoStableOrbitQ[a,p,e,\[Theta]inc] checks if given parameters correspond to a stable orbit";
KerrGeoISCO::usage = "KerrGeoISCO[a,opts] computes the location of the inner-most stable circular orbit (ISCO)"
KerrGeoPhotonSphereRadius::usage = "KerrGeoPhotonSphereRadius[a,\[Theta]inc] computes the radius of the photon sphere"
KerrGeoIBSO::usage = "KerrGeoIBSO[a,\[Theta]inc] computes the radius of the inner-most bound spherical orbit (IBSO)"
KerrGeoISSO::usage = "KerrGeoISSO[a,\[Theta]inc] computes the radius of the inner-most stable spherical orbit (ISSO)"
KerrGeoSeparatrix::usage = "KerrGeoSepatrix[a,e,\[Theta]inc] calculates the value of p at the sepatrix between stable and plunging/scattered orbits*)"


Begin["`Private`"];

(*ELQ calculation for Schwarzschild spacetime*)
KerrGeoELQ[(0|0.0),p_,e_,\[Theta]inc_]:=Module[{E0,L0,Q},
 E0=Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))];
 L0=p/Sqrt[-3-e^2+p];
 {E0,Cos[\[Theta]inc]L0,Sin[\[Theta]inc]^2 L0^2}
]


(*This function computes the orbital constants of motion from W. Schmidt, arXiv:gr-qc/0202090 *)
KerrGeoELQ[a_(*/;Abs[a]<=1*), p_, e_, \[Theta]inc1_?NumericQ] := Module[{M=1,f, g, h, d, fp, gp, hp, dp, r, rp, ra, zm, \[CapitalDelta], \[Rho], \[Kappa], \[Epsilon], \[Eta], \[Sigma], En, L, Q, E1, Em1, f1, g1, h1, d1, f2, g2, h2, d2, L1, L2,r0,\[CapitalDelta]0,Z,\[Theta]min,\[Theta]inc=\[Theta]inc1},

\[Theta]inc=Mod[\[Theta]inc,2\[Pi]];
If[\[Theta]inc>\[Pi], \[Theta]inc=2\[Pi]-\[Theta]inc];
If[\[Theta]inc <= \[Pi]/2 , \[Theta]min = \[Pi]/2-\[Theta]inc, \[Theta]min=-\[Pi]/2+\[Theta]inc];

(*Equations for polar orbits from Stoghianidis & Tsoubelis GRG, vol. 19, No. 12, p. 1235 (1987)*)
If[Mod[\[Theta]inc,\[Pi]]==\[Pi]/2,
	If[e!=0,Print["Polar, non-spherical orbits not yet implemented"];Return[];,
		r0 = p;
		\[CapitalDelta]0 = r0^2-2M r0+a^2;
		Z = r0^3-3M r0^2+a^2 r0+M a^2;
		En =Sqrt[(r0 \[CapitalDelta]0^2)/((r0^2+a^2)Z)];
		Q = r0 (M r0^3+a^2 r0^2-3M a^2 r0 +a^4)/Z - a^2 En^2;
		Return[{En,0,Q}]
	];
];

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

(*Calculate the frequencies [\[Gamma]r,\[Gamma]\[Theta],\[Gamma]\[Phi]] with respect to Mino time as per Fujita and Hikida [Class.Quantum Grav.26 (2009) 135002]*)	  
KerrGeoFreqs[a_/;Abs[a]<1,p_,e_,\[Theta]inc1_?NumericQ]:=Module[{M=1,En,L,Q,r1,r2,AplusB,AB,r3,r4,\[Epsilon]0,zm,kr,k\[Theta],\[Gamma]r,\[Gamma]\[Theta],\[Gamma]\[Phi],\[CapitalGamma],rp,rm,hp,hm,hr,EnLQ,a2zp,\[Epsilon]0zp,zmOverZp,\[Theta]min,\[Theta]inc=\[Theta]inc1},

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
\[Gamma]r=(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)])/(2EllipticK[kr^2]);(*Eq.(15)*)
\[Gamma]\[Theta]=(Pi L Sqrt[\[Epsilon]0zp])/(2EllipticK[k\[Theta]^2]);(*Eq.(15)*)

rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));

\[Gamma]\[Phi]=(2\[Gamma]\[Theta])/(Pi Sqrt[\[Epsilon]0zp]) EllipticPi[zm,k\[Theta]^2]+(2a \[Gamma]r)/(Pi(rp-rm)Sqrt[(1-En^2)(r1-r3)(r2-r4)])(*Eq. (21)*)((2M En rp - a L)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-(2M En rm - a L)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2]));

(*Convert to frequencies w.r.t BL time Subscript[\[CapitalOmega], k]=Subscript[\[Gamma], k]/\[CapitalGamma] using Fujita and Hikida's formula Eq. (21)*)
\[CapitalGamma]=4M^2 En + (2a2zp En  \[Gamma]\[Theta])/(Pi L Sqrt[\[Epsilon]0zp]) (EllipticK[k\[Theta]^2]- EllipticE[k\[Theta]^2]) + (2\[Gamma]r)/(Pi Sqrt[(1-En^2)(r1-r3)(r2-r4)]) (En/2 ((r3(r1+r2+r3)-r1 r2)EllipticK[kr^2]+(r2-r3)(r1+r2+r3+r4)EllipticPi[hr,kr^2]+(r1-r3)(r2-r4)EllipticE[kr^2])+2M En(r3 EllipticK[kr^2]+(r2-r3)EllipticPi[hr,kr^2])+(2M)/(rp-rm) (((4M^2 En-a L)rp-2M a^2 En)/(r3-rp) (EllipticK[kr^2]-(r2-r3)/(r2-rp) EllipticPi[hp,kr^2])-((4M^2 En-a L)rm-2M a^2 En)/(r3-rm) (EllipticK[kr^2]-(r2-r3)/(r2-rm) EllipticPi[hm,kr^2])));

(*Output the BL frequencies by dividing the Mino time frequencies by the conversion factor \[CapitalGamma]*)
{\[Gamma]r/\[CapitalGamma],Abs[\[Gamma]\[Theta]/\[CapitalGamma]],\[Gamma]\[Phi]/\[CapitalGamma],\[CapitalGamma]}
]

KerrGeoFreqs[a_/;a==1,p_,e_,\[Theta]inc1_?NumericQ]:=Module[{},
	Print["Frequency calculation not yet implemented for a=M (but the equations are in the appendix of Fujita and Hikida so please implement them!)"];
]

(*Orbit stability check, Schwarzschild case*)
KerrGeoStableOrbitQ[(0|0.0),p_,e_,\[Theta]inc_]:=Module[{},
  If[p>6+2e,True,False](*> not \[GreaterEqual] as orbits along the separatrix are marginally stable*)
];

(* The orbit is stable only if Subscript[\[CapitalOmega], r] is real*)
KerrGeoStableOrbitQ[a_?NumericQ/;a!=0,p_?NumericQ,e_?NumericQ,\[Theta]inc_?NumericQ]:=Module[{freqs},
  freqs=KerrGeoFreqs[a,p,e,\[Theta]inc];
  Element[freqs[[1]],Reals]
];


Options[KerrGeoISCO] = {"orbit" -> "Prograde"}
KerrGeoISCO[a_,OptionsPattern[]]:=Module[{Z1,Z2},
	Z1=1+(1-a^2)^(1/3) ((1+a)^(1/3)+(1-a)^(1/3));
	Z2=(3a^2+Z1^2)^(1/2);
	If[OptionValue["orbit"]=="Prograde",
		Return[3+Z2-((3-Z1)(3+Z1+2Z2))^(1/2)],
		Return[3+Z2+((3-Z1)(3+Z1+2Z2))^(1/2)]
	];
];

(*Photon sphere radius is where the energy for a timelike orbit diverges*)
KerrGeoPhotonSphereRadius[a_,\[Theta]inc_]:= Module[{res},
	res=Sort[Re[p/.Solve[Simplify[Denominator[KerrGeoELQ[a,p,0,\[Theta]inc][[1]]^2]==0,Assumptions->{a>=0,p>0}],p]],Greater];
	If[\[Theta]inc>=\[Pi]/2,Return[res[[1]]],Return[res[[2]]]]
]

(*Photon sphere becomes light ring at r=3M in Schwarzschild*)
KerrGeoPhotonSphereRadius[(0|0.0),\[Theta]inc_]:=3;

(*Equatorial photon sphere results from Bardeen, Press, Teukolsky 1972*)
KerrGeoPhotonSphereRadius[a_,0]:=2(1+Cos[2/3 ArcCos[-a]])
KerrGeoPhotonSphereRadius[a_,\[Pi]]:=2(1+Cos[2/3 ArcCos[a]])

(* The IBSO is where E=1. Use the fact that rph < r_ibso*)
KerrGeoIBSO[a_,\[Theta]inc_]:= Module[{rph},
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
KerrGeoISSO[a_,\[Theta]inc_]:=Module[{rmb},
	rmb=KerrGeoIBSO[a,\[Theta]inc];
	p/.FindRoot[KerrGeoRadialEqRoots[a,p,0,\[Theta]inc][[5]],{p,rmb+10^-10,9},Method->"Brent"]
]

(*For equatorial orbits the ISSO is the ISCO*)
KerrGeoISSO[a_,0]:= KerrGeoISCO[a]
KerrGeoISSO[a_,\[Pi]]:= KerrGeoISCO[a,orbit->"Retrograde"]

(* The separatrix occurs when Subscript[\[CapitalOmega], r]=0 which implies Subscript[r, 2]=Subscript[r, 3]. Search for the value of p at which this occurs which must be outside the p of the ISSO*)
KerrGeoSeparatrix[a_,e_,\[Theta]inc_]:=Module[{},
	p/.FindRoot[KerrGeoRadialEqRoots[a,p,e,\[Theta]inc][[5]],{p,KerrGeoISSO[a,\[Theta]inc],20},Method->"Brent"]
]

(*Separatrix for Schwarzschild*)
KerrGeoSeparatrix[0,e_,\[Theta]inc_]:= 6+2e;

(*Separatrix for equatorial Kerr from Levin and Periz-Giz arXiv:1108.1819*)
KerrGeoSeparatrix[a1_,e_,\[Theta]inc_/;Mod[\[Theta]inc,\[Pi]]==0]:= Module[{ru,a=a1},
If[\[Theta]inc==\[Pi], a = -a];
ru=ru/.Solve[e==(-ru^2+6 ru-8a ru^(1/2)+3a^2)/(ru^2-2ru+a^2),ru][[-1]];
(4 ru (ru^(1/2)-a)^2)/(ru^2-2ru+a^2)
]

(*From Glampedakis and Kennefick arXiv:gr-qc/0203086*)
KerrGeoSeparatrix[1,e_,0]:=1+e
	
End[];

EndPackage[];
