(* ::Package:: *)

(* ::Title:: *)
(*KerrGeoOrbit subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`KerrGeoOrbit`",
	{"KerrGeodesics`ConstantsOfMotion`",
	 "KerrGeodesics`OrbitalFrequencies`",
	 "KerrGeodesics`Velocity`"}];

KerrGeoOrbit::usage = "KerrGeoOrbit[a,p,e,x] returns a KerrGeoOrbitFunction[..] which stores the orbital trajectory and parameters.";
KerrGeoOrbitFunction::usage = "KerrGeoOrbitFunction[a,p,e,x,assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Section::Closed:: *)
(*Schwarzschild*)


(* ::Text:: *)
(*The analytic equations below are taken from Appendix B of "Fast Self-forced Inspirals" by M. van de Meent and N. Warburton, Class. Quant. Grav. 35:144003 (2018), arXiv:1802.05281*)


(* ::Text:: *)
(**)


(*t and \[Phi] accumulated over one orbit*)
\[CapitalPhi]SchwarzDarwin[p_,e_]:=4 Sqrt[p/(p-6+2e)] EllipticK[(4 e)/(p-6+2e)]
TSchwarzDarwin[p_,e_]:=(2p Sqrt[(p-6+2e)((p-2)^2-4e^2)])/((1-e^2)(p-4)) EllipticE[(4e)/(p-6+2e)]-2p Sqrt[(p-2)^2-4e^2]/((1-e^2)Sqrt[p-6+2e]) EllipticK[(4e)/(p-6+2e)]-(4(8(1-e^2)+p(1+3e^2-p))Sqrt[(p-2)^2-4e^2])/((1-e)(1-e^2)(p-4)Sqrt[p-6+2e]) EllipticPi[-((2e)/(1-e)),(4e)/(p-6+2e)]+(16Sqrt[(p-2)^2-4e^2])/((p-2+2e)Sqrt[p-6+2e]) EllipticPi[(4e)/(p-2+2e),(4e)/(p-6+2e)]


tSchwarzDarwin[p_,e_,\[Xi]_]:=TSchwarzDarwin[p,e]/2+((p Sqrt[(p-6+2e)((p-2)^2-4e^2)])/((1-e^2)(p-4)) EllipticE[\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-p Sqrt[(p-2)^2-4e^2]/((1-e^2)Sqrt[p-6+2e]) EllipticF[\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-(2(8(1-e^2)+p(1+3e^2-p))Sqrt[(p-2)^2-4e^2])/((1-e)(1-e^2)(p-4)Sqrt[p-6+2e]) EllipticPi[-((2e)/(1-e)),\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]+(8Sqrt[(p-2)^2-4e^2])/((p-2+2e)Sqrt[p-6+2e]) EllipticPi[(4e)/(p-2+2e),\[Xi]/2-\[Pi]/2,(4e)/(p-6+2e)]-e p Sqrt[((p-2)^2-4e^2)(p-6-2e Cos[\[Xi]])]/((1-e^2)(p-4)(1+e Cos[\[Xi]])) Sin[\[Xi]])
rSchwarzDarwin[p_,e_,\[Chi]_]:=p/(1 + e Cos[\[Chi]])
\[Theta]SchwarzDarwin[p_,e_,\[Chi]_]:= \[Pi]/2 
\[Phi]SchwarzDarwin[p_,e_,\[Xi]_]:=\[CapitalPhi]SchwarzDarwin[p,e]/2+2Sqrt[p/(p-6+2e)]EllipticF[\[Xi]/2-\[Pi]/2,(4 e)/(p-6+2e)]


(* ::Text:: *)
(*Specialization to circular (equatorial) orbits*)


tSchwarzDarwin[p_/;p>6, 0, \[Xi]_] := ((p^2) \[Xi] )/Sqrt[-6+p] 
rSchwarzDarwin[p_/;p>6, 0, \[Xi]_] := 0;
\[Phi]SchwarzDarwin[p_/;p>6, 0, \[Xi]_] := Sqrt[p/(-6+p)] \[Xi]


(* ::Text:: *)
(*FIXME: make the below work for inclined orbits and accept initial phases*)


KerrGeoOrbitSchwarzDarwin[p_, e_] := Module[{t, r, \[Theta], \[Phi], assoc, consts, En,L,Q},

t[\[Chi]_] := tSchwarzDarwin[p,e,\[Chi]];
r[\[Chi]_] := rSchwarzDarwin[p,e,\[Chi]];
\[Theta][\[Chi]_] := \[Theta]SchwarzDarwin[p,e,\[Chi]];
\[Phi][\[Chi]_] := \[Phi]SchwarzDarwin[p,e,\[Chi]];

consts = KerrGeoConstantsOfMotion[0,p,e,1];
{En,L,Q} = Values[consts];

assoc = Association[
			"Trajectory" -> {t,r,\[Theta],\[Phi]},
			"Parametrization" -> "Darwin", 
			"ConstantsOfMotion"-> consts, 
			"a" -> 0,
			"p" -> p,
			"e" -> e,
			"Inclination" -> 1,
			"Energy" -> En,
			"AngularMomentum" -> L,
			"CarterConstant" -> Q
			];

KerrGeoOrbitFunction[0, p, e, 1, assoc]

]


(* ::Section:: *)
(*Kerr*)


(* ::Subsection:: *)
(*Equatorial (Darwin)*)


(* ::Text:: *)
(*Compute the orbit using Mino time and then convert to Darwin time using \[Lambda][r[\[Chi]]] where \[Lambda][r] is found in Fujita and Hikida (2009).*)
(*ToFix: Currently doesn't take in initial phases in any way*)


KerrGeoOrbitEquatorialDarwin[a_,p_,e_,x_/;x^2==1, initPhases_] := Module[{orbitMino,freqs,r1,r2,r3,r4,\[CapitalLambda]r,yr,kr,\[Lambda]0r,r,r01,\[CapitalLambda]r1,\[Lambda],consts,En,L,Q,tMino,rMino,\[Theta]Mino,\[Phi]Mino,tDarwin,rDarwin,\[Theta]Darwin,\[Phi]Darwin,assoc,velocity},

orbitMino = KerrGeoOrbit[a,p,e,x,initPhases];

{r1,r2,r3,r4} = orbitMino["RadialRoots"];
freqs = orbitMino["Frequencies"];
consts = orbitMino["ConstantsOfMotion"];
{En,L,Q} = Values[consts];
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

velocity = Values[KerrGeoVelocity[a,p,e,x,"Parametrization"->"Darwin"]];

assoc = Association[
			"Trajectory" -> {tDarwin,rDarwin,\[Theta]Darwin,\[Phi]Darwin}, 
			"Velocity"-> velocity,
			"Parametrization" -> "Darwin", 
			"ConstantsOfMotion"-> consts, 
			"RadialRoots"->{r1,r2,r3,r4},
			"a" -> a,
			"p" -> p,
			"e" -> e,
			"Inclination" -> x,
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
Module[{M=1,consts,En,L,Q,r1,r2,r3,r4,p3,p4,assoc,var,t0, \[Chi]0, \[Phi]0,r0,\[Theta]0,t,r,\[Theta],\[Phi],\[Chi],growthRateT,growthRatePh,
		\[Chi]r,NrMax,pg,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[Phi]C,tC,Pr,r0Sample,PrSample,dtd\[Chi],d\[Phi]d\[Chi],TVr,PVr},
	consts = KerrGeoConstantsOfMotion[a,p,e,x];
	{En,L,Q} = Values[consts];
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
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
		"ConstantsOfMotion" -> consts,
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"RadialRoots" -> {r1,r2,r3,r4},
		"a" -> a,
		"p" -> p,
		"e" -> e,
		"Inclination" -> x
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Subsection::Closed:: *)
(*Circular (Fast Spec - Darwin)*)


(* Hopper, Forseth, Osburn, and Evans, PRD 92 (2015)*)


(* ::Subsubsection::Closed:: *)
(*Main file that calculates geodesics using spectral integration*)


KerrGeoOrbitFastSpecDarwin[a_,p_,e_/;e==0,x_,initPhases:{_,_,_,_}:{0,0,0,0}]:=
Module[{M=1,consts,En,L,Q,zp,zm,assoc,var,t0, \[Chi]0, \[Phi]0,r0,\[Theta]0,t,r,\[Theta],\[Phi],\[Chi],freqT,freqPh,
		\[Chi]\[Theta],pg,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]\[Theta],\[Phi]C,tC,P\[Theta],\[Theta]0Sample,P\[Theta]Sample,dtd\[Chi],d\[Phi]d\[Chi],TV\[Theta],PV\[Theta],\[Beta],\[Alpha],zRoots},
	consts = KerrGeoConstantsOfMotion[a,p,e,x];
	{En,L,Q} = Values[consts];
	
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
		"ConstantsOfMotion" -> consts,
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"PolarRoots" -> zRoots,
		"a" -> a,
		"p" -> p,
		"e" -> e,
		"Inclination" -> x
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Subsection::Closed:: *)
(*Circular (Mino)*)


(* ::Text:: *)
(*FIXME: make the initial phases work in this case*)


KerrGeoOrbitMino[a_, p_, (0|0.), (1|1.), initPhases:{_,_,_,_}:{0,0,0,0}] := Module[{consts, assoc, t, r, \[Theta], \[Phi], En, L, Q, \[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t, e=0, x=1, velocity},

	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]];
	consts = KerrGeoConstantsOfMotion[a,p,e,x];
	{En,L,Q} = Values[consts];

	t[\[Lambda]_] := ((a^3 Sqrt[2 a+(-3+p) Sqrt[p]] p^2+a Sqrt[2 a+(-3+p) Sqrt[p]] (-2+p) p^3+a^2 Sqrt[(2 a+(-3+p) Sqrt[p]) p^7]-2 Sqrt[(2 a+(-3+p) Sqrt[p]) p^9]+Sqrt[(2 a+(-3+p) Sqrt[p]) p^11]) \[Lambda])/((2 a+(-3+p) Sqrt[p]) p^(3/4) (a^2+(-2+p) p));
	r[\[Lambda]_] := p;
	\[Theta][\[Lambda]_] := \[Pi]/2;
	\[Phi][\[Lambda]_] := (p^(5/4) \[Lambda])/Sqrt[2 a+(-3+p) Sqrt[p]];
	
	velocity = Values[KerrGeoVelocity[a,p,e,x,{initPhases[[2]],initPhases[[3]]}]];
	
	assoc = Association[
		"Parametrization"->"Mino", 
		"Energy" -> En, 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency" -> \[CapitalUpsilon]r,
		"PolarFrequency" ->  \[CapitalUpsilon]\[Theta],
		"AzimuthalFrequency" -> \[CapitalUpsilon]\[Phi],
		"Frequencies" -> {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"Velocity" -> velocity,
		"a" -> a,
		"p" -> p,
		"e" -> e,
		"Inclination" -> x
	];
	
	KerrGeoOrbitFunction[a, p, e, x, assoc]

]


(* ::Subsection::Closed:: *)
(*Generic (Mino)*)


KerrGeoOrbitMino[a_,p_,e_,x_,initPhases:{_,_,_,_}:{0,0,0,0}]:=Module[{M=1,consts,En,L,Q,assoc,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,zp,zm,kr,k\[Theta],rp,rm,hr,hp,hm,rq,zq,\[Psi]r,tr,\[Phi]f,\[Psi]z,tz,\[Phi]z,qt0,qr0,qz0,q\[Phi]0,t,r,\[Theta],\[Phi],\[Phi]t,\[Phi]r,Ct,C\[Phi],velocity},
	consts = KerrGeoConstantsOfMotion[a,p,e,x];
	{En,L,Q} = Values[consts];
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]];
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
	{zp,zm} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoPolarRoots[a, p, e, x];
	
	kr = (r1-r2)/(r1-r3) (r3-r4)/(r2-r4);
	k\[Theta] = a^2 (1-En^2)(zm/zp)^2;

rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
hr=(r1-r2)/(r1-r3);
hp=((r1-r2)(r3-rp))/((r1-r3)(r2-rp));
hm=((r1-r2)(r3-rm))/((r1-r3)(r2-rm));

rq = Function[{qr},(r3(r1 - r2)JacobiSN[EllipticK[kr]/\[Pi] qr,kr]^2-r2(r1-r3))/((r1-r2)JacobiSN[EllipticK[kr]/\[Pi] qr,kr]^2-(r1-r3))];

zq = Function[{qz}, zm JacobiSN[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]]];

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

\[Psi]z[qz_]:= \[Psi]z[qz] = JacobiAmplitude[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]];
tz[qz_]:= 1/(1-En^2) En zp ( EllipticE[k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticE[\[Psi]z[qz],k\[Theta]]);
\[Phi]z[qz_]:= -1/zp L ( EllipticPi[zm^2,k\[Theta]]2((qz+\[Pi]/2)/\[Pi])-EllipticPi[zm^2,\[Psi]z[qz],k\[Theta]]);

{qt0, qr0, qz0, q\[Phi]0} = {initPhases[[1]], initPhases[[2]], initPhases[[3]], initPhases[[4]]};
(*Calculate normalization constants so that t=0 and \[Phi]=0 at \[Lambda]=0 when qt0=0 and q\[Phi]0=0 *)
Ct=tr[qr0]+tz[qz0]/.i_/;i==0:>0;
C\[Phi]=\[Phi]r[qr0]+\[Phi]z[qz0]/.i_/;i==0:>0;

t[\[Lambda]_]:= qt0 + \[CapitalUpsilon]t \[Lambda] + tr[\[CapitalUpsilon]r \[Lambda] + qr0] + tz[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]-Ct;
r[\[Lambda]_]:= rq[\[CapitalUpsilon]r \[Lambda]+ qr0];
\[Theta][\[Lambda]_]:= ArcCos[zq[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]];
\[Phi][\[Lambda]_]:= q\[Phi]0 + \[CapitalUpsilon]\[Phi] \[Lambda] + \[Phi]r[\[CapitalUpsilon]r \[Lambda]+ qr0] + \[Phi]z[\[CapitalUpsilon]\[Theta] \[Lambda] + qz0]-C\[Phi];

velocity = Values[KerrGeoVelocity[a,p,e,x,{initPhases[[2]],initPhases[[3]]}]];


	assoc = Association[
	"Parametrization"->"Mino", 
	"Energy" -> En, 
	"AngularMomentum" -> L, 
	"CarterConstant" -> Q, 
	"ConstantsOfMotion" -> consts,
	"RadialFrequency" -> \[CapitalUpsilon]r,
	"PolarFrequency" ->  \[CapitalUpsilon]\[Theta],
	"AzimuthalFrequency" -> \[CapitalUpsilon]\[Phi],
	"Frequencies" -> {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
	"Trajectory" -> {t,r,\[Theta],\[Phi]},
	"Velocity"-> velocity,
	"RadialRoots" -> {r1,r2,r3,r4},
	"a" -> a,
	"p" -> p,
	"e" -> e,
	"Inclination" -> x
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	PVr[rp_]:=-((a^2*L)/(a^2 - 2*M*rp + rp^2)) + a*En*(-1 + (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2));

	sampledFuncR=LambdaToPsiRTransform[a,p,e,x,PVr];
	\[CapitalDelta]\[Phi]r=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncR,{\[CapitalUpsilon]r,minoSampleR}];
	\[CapitalDelta]\[Phi]r
];

PhiOfMinoFastSpecR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[Lambda]rSample,pg,\[Psi]r},
	\[CapitalUpsilon]r = Re[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]][[1]];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	PV\[Theta][\[Theta]p_]:=L*Csc[\[Theta]p]^2;
	
	sampledFuncTheta=LambdaToChiThetaTransform[a,p,e,x,PV\[Theta]];
	\[CapitalDelta]\[Phi]\[Theta]=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncTheta,{\[CapitalUpsilon]\[Theta],minoSampleTh}];
	\[CapitalDelta]\[Phi]\[Theta]
];

PhiOfMinoFastSpecTheta[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]\[Theta],\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],\[Lambda]\[Theta]Sample,pg,\[Chi]\[Theta]},
	\[CapitalUpsilon]\[Theta] = Re[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]][[2]];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	TVr[rp_]:=(En*(a^2 + rp^2)^2)/(a^2 - 2*M*rp + rp^2) + a*L*(1 - (a^2 + rp^2)/(a^2 - 2*M*rp + rp^2));

	sampledFuncR=LambdaToPsiRTransform[a,p,e,x,TVr];
	\[CapitalDelta]tr=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncR,{\[CapitalUpsilon]r,minoSampleR}];
	\[CapitalDelta]tr
];

TimeOfMinoFastSpecR[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]r,\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[Lambda]rSample,pg,\[Psi]r},
	\[CapitalUpsilon]r = Re[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]][[1]];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	TV\[Theta][\[Theta]p_]:=-(a^2*En*Sin[\[Theta]p]^2);

	sampledFuncTheta=LambdaToChiThetaTransform[a,p,e,x,TV\[Theta]];
	\[CapitalDelta]t\[Theta]=DarwinFastSpecMinoIntegrateAndConvergenceCheck[sampledFuncTheta,{\[CapitalUpsilon]\[Theta],minoSampleTh}];
	\[CapitalDelta]t\[Theta]
];

TimeOfMinoFastSpecTheta[a_,p_,e_,x_]:=Module[{\[CapitalUpsilon]\[Theta],\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],\[Lambda]\[Theta]Sample,pg,\[Chi]\[Theta]},
	\[CapitalUpsilon]\[Theta] = Re[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]][[2]];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
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
	{En,L,Q} = Values[KerrGeoConstantsOfMotion[a,p,e,x]];
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


(* ::Subsubsection:: *)
(*Main file that calculates geodesics using spectral integration*)


Clear[KerrGeoOrbitFastSpec];
Options[KerrGeoOrbitFastSpec]={InitialPosition->{0,0,0,0}};
KerrGeoOrbitFastSpec[a_,p_,e_,x_,initPhases:{_,_,_,_}:{0,0,0,0},opts:OptionsPattern[]]:=
Module[{M=1,consts,En,L,Q,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,p3,p4,\[Alpha],\[Beta],zp,zm,assoc,var,\[Chi]0,\[Psi]0,
		r0,\[Theta]0,qt0,qr0,q\[Theta]0,q\[Phi]0,\[Lambda]t0,\[Lambda]r0,\[Lambda]\[Theta]0,\[Lambda]\[Phi]0,t,r,\[Theta],\[Phi],\[Psi],\[Chi],\[CapitalDelta]\[Lambda]r,\[Lambda]r,\[CapitalDelta]\[Lambda]\[Theta],\[Lambda]\[Theta],rC,\[Theta]C,\[CapitalDelta]r,\[CapitalDelta]\[Theta],
		\[Psi]r,\[Chi]\[Theta],NrMax,NthMax,pg,\[Lambda]rSample,\[Lambda]\[Theta]Sample,\[CapitalDelta]tr,\[CapitalDelta]\[Phi]r,\[CapitalDelta]t\[Theta],\[CapitalDelta]\[Phi]\[Theta],\[Phi]C,tC,zRoots,tInit,rInit,\[Theta]Init,\[Phi]Init, velocity},
	consts = KerrGeoConstantsOfMotion[a,p,e,x];
	{En,L,Q} = Values[consts];
	
	(* Useful constants for \[Theta]-dependent calculations *)
	\[Beta]=a^2(1-En^2);
	\[Alpha]=L^2+Q+\[Beta];
	zp=Sqrt[(\[Alpha]+Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zm=Sqrt[(\[Alpha]-Sqrt[\[Alpha]^2-4 Q \[Beta]])/(2\[Beta])];
	zRoots={ArcCos[zm],Pi-ArcCos[zm]};
	
	(* Useful constants for r-dependent calculations *)
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
	
	(* Mino frequencies of orbit. I call Re, because some frequencies 
	are given as imaginary for restricted orbits. 
	Maybe something to fix with KerrGeoMinoFrequencies*)
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]];
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
	
	velocity = Values[KerrGeoVelocity[a,p,e,x,{initPhases[[2]],initPhases[[3]]}]];
	
	
	assoc = Association[
		"Parametrization"->"Mino",
		"Energy" -> En, 
		"AngularMomentum" -> L, 
		"CarterConstant" -> Q, 
		"ConstantsOfMotion" -> consts,
		"RadialFrequency" -> \[CapitalUpsilon]r,
		"PolarFrequency" ->  \[CapitalUpsilon]\[Theta],
		"AzimuthalFrequency" -> \[CapitalUpsilon]\[Phi],
		"Frequencies" -> {\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi]},
		"Trajectory" -> {t,r,\[Theta],\[Phi]},
		"Velocity"-> velocity,
		"RadialRoots" -> {r1,r2,r3,r4},
		"PolarRoots" -> zRoots,
		"a" -> a,
		"p" -> p,
		"e" -> e,
		"Inclination" -> x
		];
	
	KerrGeoOrbitFunction[a,p,e,x,assoc]
]


(* ::Section:: *)
(*KerrGeoOrbit and KerrGeoOrbitFuction*)


Options[KerrGeoOrbit] = {"Parametrization" -> "Mino", "Method" -> "FastSpec"}
SyntaxInformation[KerrGeoOrbit] = {"ArgumentsPattern"->{_,_,OptionsPattern[]}};


KerrGeoOrbit[a_,p_,e_,x_, initPhases:{_,_,_,_}:{0,0,0,0},OptionsPattern[]]:=Module[{param, method},
(*FIXME: add stability check but make it possible to turn it off*)

method = OptionValue["Method"];
param = OptionValue["Parametrization"];

If[param == "Darwin" && Abs[x]!=1, Print["Darwin parameterization only valid for equatorial motion"];Return[];];

If[Precision[{a,p,e,x}] > 30, method = "Analytic"];

If[method == "FastSpec",

	If[param == "Mino",  If[PossibleZeroQ[a] || PossibleZeroQ[e], Return[KerrGeoOrbitMino[a, p, e, x, initPhases]], Return[KerrGeoOrbitFastSpec[a, p, e, x, initPhases]]]];
	If[param == "Darwin", 
		If[PossibleZeroQ[a], Return[KerrGeoOrbitSchwarzDarwin[p, e]], Return[KerrGeoOrbitFastSpecDarwin[a,p,e,x,initPhases]]]
	];
	Print["Unrecognized parametrization: " <> OptionValue["Parametrization"]];
	
];

If[method == "Analytic",

	If[param == "Mino", Return[KerrGeoOrbitMino[a, p, e, x, initPhases]]];
	If[param == "Darwin", 
		If[PossibleZeroQ[a], Return[KerrGeoOrbitSchwarzDarwin[p, e]], Return[KerrGeoOrbitEquatorialDarwin[a,p,e,x,initPhases]]]
	];
	Print["Unrecognized parametrization: " <> OptionValue["Parametrization"]];

];

Print["Unrecognized method: " <> method];

]


KerrGeoOrbitFunction /:
 MakeBoxes[kgof:KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"p: ", p}], "  ",
                  BoxForm`SummaryItem[{"e: ", e}], "  ",
                  BoxForm`SummaryItem[{"x: ", x}]}],
             BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}]};
  extended = {BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}],
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["AngularMomentum"]}],
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["CarterConstant"]}]};
  BoxForm`ArrangeSummaryBox[
    KerrGeoOrbitFunction,
    kgof,
    Module[{t, r, \[Theta], \[CurlyPhi]},
      {t, r, \[Theta], \[CurlyPhi]} = kgof["Trajectory"];
      Quiet[Show[ParametricPlot3D[{r[\[Lambda]] Sin[\[Theta][\[Lambda]]] Cos[\[CurlyPhi][\[Lambda]]],r[\[Lambda]] Sin[\[Theta][\[Lambda]]] Sin[\[CurlyPhi][\[Lambda]]],r[\[Lambda]] Cos[\[Theta][\[Lambda]]]}, {\[Lambda],0,2\[Pi]}, Boxed->False, Axes->False, PlotRange->All, ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"]/AbsoluteCurrentValue[Magnification]}]],
           Graphics3D[{Black,Sphere[{0,0,0},1+Sqrt[1-a^2]]}]]]
      ],
    summary,
    extended,
    form]
];


KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := Through[assoc["Trajectory"][\[Lambda]]]
KerrGeoOrbitFunction[a_, p_, e_, x_, assoc_][y_?StringQ] := assoc[y]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
