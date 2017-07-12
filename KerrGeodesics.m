(* ::Package:: *)

BeginPackage["KerrGeodesics`"];

KerrGeoELQ::usage = "KerrGeoELQ[a, p, e, \[Theta]min] returns the energy, z-component of the angular momentum and the Carter constant";
KerrGeoFreqs::usage = "KerrGeoFreqs[a, p, e, \[Theta]min] returns the radial, polar and azimuthal frequencies and the conversion factor between Boyer-Lindquist and Mino time frequencies";


Begin["`Private`"];

KerrGeoELQSchwarz[p_,e_,\[Theta]inc_]:=Module[{E0,L0,Q},
E0=Sqrt[(-4 e^2+(-2+p)^2)/(p (-3-e^2+p))];
L0=p/Sqrt[-3-e^2+p];
{E0,Cos[\[Theta]inc]L0,Sin[\[Theta]inc]^2 L0^2}
]


(*This function computes the orbital constants of motion from W. Schmidt, arXiv:0202090 [gr-qc] *)
KerrGeoELQGeneric[a_, p_, e_, \[Theta]inc1_?NumericQ] := Module[{M=1,f, g, h, d, fp, gp, hp, dp, r, rp, ra, zm, \[CapitalDelta], \[Rho], \[Kappa], \[Epsilon], \[Eta], \[Sigma], En, L, Q, E1, Em1, f1, g1, h1, d1, f2, g2, h2, d2, L1, L2,r0,\[CapitalDelta]0,Z,\[Theta]min,\[Theta]inc=\[Theta]inc1},

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

KerrGeoELQ[a_, p_, e_, \[Theta]inc_]:=Module[{},

If[! (0 <= a <= 1), Print["Domain error: 0 <= q <= 1 reqired"]; Return[];];

If[a==0,
	Return[KerrGeoELQSchwarz[p,e,\[Theta]inc]],
	Return[KerrGeoELQGeneric[a,p,e,\[Theta]inc]]
];

];

	  
KerrGeoFreqs[a_,p_,e_,\[Theta]inc1_?NumericQ]:=Module[{M=1,En,L,Q,r1,r2,AplusB,AB,r3,r4,\[Epsilon]0,zm,kr,k\[Theta],\[Gamma]r,\[Gamma]\[Theta],\[Gamma]\[Phi],\[CapitalGamma],rp,rm,hp,hm,hr,EnLQ,a2zp,\[Epsilon]0zp,zmOverZp,\[Theta]min,\[Theta]inc=\[Theta]inc1},

\[Theta]inc=Mod[\[Theta]inc,2\[Pi]];
If[\[Theta]inc>\[Pi], \[Theta]inc=2\[Pi]-\[Theta]inc];

If[\[Theta]inc==\[Pi]/2, Print["Equations for polar orbits not implemented yet"];Return[];];

{En,L,Q}=KerrGeoELQ[a,p,e,\[Theta]inc];
\[Theta]min=(\[Pi]/2-\[Theta]inc)/Sign[L];

(*Calculate the frequencies [\[Gamma]r,\[Gamma]\[Theta],\[Gamma]\[Phi]] with respect to Mino time as per Fujita and Hikida [Class.Quantum Grav.26 (2009) 135002]*)
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
	
End[];

EndPackage[];
