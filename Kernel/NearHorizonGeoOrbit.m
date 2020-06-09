(* ::Package:: *)

(* ::Title:: *)
(*NearHorizonGeoOrbit subpackage for KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`NearHorizonGeoOrbit`"];

NearHorizonGeoOrbitClass::usage = "NearHorizonGeoOrbitClass[spacetime, radialClass, polarClass] returns a NearHorizonGeoOrbitFunction[..] which stores the (symbolic) trajectory and orbital parameters of 
a generic near-horizon geodesic of given radial and polar class. The classification refers to arXiv:2001.03478.";
NearHorizonGeoOrbit::usage = "NearHorizonGeoOrbit[spacetime, \[ScriptCapitalE], \[ScriptCapitalL], \[ScriptCapitalQ], \[Mu]] returns a NearHorizonGeoOrbitFunction[..] which stores the (numerical) trajectory and orbital parameters of a 
near-horizon geodesic.";
NearHorizonGeoOrbitFunction::usage = "NearHorizonGeoOrbitFunction[assoc] an object for storing the trajectory and orbital parameters in the assoc Association.";


(* ::Section::Closed:: *)
(*Conventions for the output syntax*)


(* Rule to write the output with well-chosen notations *)
BHPTStyle := {0 -> 0, \[Lambda]f -> Subscript[\[Lambda], f], Rf -> Subscript[R, f], RPlus -> SubPlus[R], ti -> Subscript[t, i], \[Phi]i -> Subscript[\[Phi], i], t0 -> Subscript[t, 0], 
\[Phi]0 -> Subscript[\[Phi], 0], e -> \[ScriptE], \[Lambda]0 -> Subscript[\[Lambda], 0], T0 -> Subscript[T, 0], \[Lambda]Plus -> SubPlus[\[Lambda]], \[Lambda]Minus -> SubMinus[\[Lambda]], Rm -> Subscript[R, m], 
\[Lambda]m -> Subscript[\[Lambda], m], R0 -> Subscript[R, 0], Ti -> Subscript[T, i], \[Lambda]i -> Subscript[\[Lambda], i], EE -> \[ScriptCapitalE], CC -> \[ScriptCapitalC],  Q -> \[ScriptCapitalQ], lStar -> SubStar[\[ScriptCapitalL]], 
lNot -> "\!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(\[SmallCircle]\)]\)", l -> \[ScriptCapitalL], CNot -> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\[SmallCircle]\)]\)", zP->SubPlus[z], zM->SubMinus[z], 
si\[Theta]->"\!\(\*SubsuperscriptBox[\(s\), \(\[Theta]\), \(i\)]\)", \[Theta]i->"\!\(\*SubscriptBox[\(\[Theta]\), \(i\)]\)", \[CapitalPhi]0->"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(i\)]\)", 
\[Lambda]Minus->"\!\(\*SubscriptBox[\(\[Lambda]\), \(-\)]\)", \[Lambda]Plus->"\!\(\*SubscriptBox[\(\[Lambda]\), \(+\)]\)", M -> "M", \[Mu] -> "\[Mu]", Ri -> "\!\(\*SubscriptBox[\(R\), \(i\)]\)", Subscript[\[CapitalPhi], \[Theta]]->"\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)", \[Kappa]->"\[Kappa]"}

Begin["`Private`"];


(* ::Section::Closed:: *)
(*NHEK radial parametrizations of generating classes*)


(* ::Subsection:: *)
(*Deflecting*)


NHEKDeflectingRadial[R_]:={\[Lambda]Plus-I/Sqrt[-CC]ArcCos[(EE l-CC R)/(EE Sqrt[CC+l^2])],-Sqrt[EE^2+2EE l R-CC R^2]/(EE R),\[CapitalPhi]0+Log[(EE+l R+Sqrt[EE^2+2EE l R-CC R^2])/R]-3 l/(4Sqrt[-CC])Log[EE l-CC R+Sqrt[-CC]Sqrt[EE^2+2EE l R-CC R^2]]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NHEKDeflectingMino[\[Lambda]_]:={Sqrt[-CC]/EE Sinh[Sqrt[-CC](\[Lambda]-\[Lambda]Plus)]/(l/Sqrt[CC+l^2]-Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Plus)]),EE/CC(l-Sqrt[CC+l^2]Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Plus)]), \[CapitalPhi]0-3l/4(\[Lambda]-\[Lambda]Plus)+2 ArcTanh[(l-Sqrt[CC+l^2])/Sqrt[-CC]Tanh[Sqrt[-CC](\[Lambda]-\[Lambda]Plus)]]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*Plunging*)


NHEKPlungingRadial[R_]:={\[Lambda]Minus + 1/Sqrt[-CC]Log[EE Sqrt[CC+l^2]/(EE l-CC R+Sqrt[-CC]Sqrt[EE^2+2EE l R-CC R^2])], Sqrt[EE^2+2EE l R-CC R^2]/(EE R), \[CapitalPhi]0-Log[(EE+l R+Sqrt[EE^2+2EE l R-CC R^2])/R]+3 l/(4Sqrt[-CC])Log[EE l-CC R+Sqrt[-CC]Sqrt[EE^2+2EE l R-CC R^2]]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NHEKPlungingMino[\[Lambda]_]:={Sqrt[-CC]/EE Abs[Sinh[Sqrt[-CC](\[Lambda]-\[Lambda]Minus)]]/(l/Sqrt[l^2+CC]-Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Minus)]), EE/CC (l-Sqrt[CC+l^2]Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Minus)]), \[CapitalPhi]0-3 l/4(\[Lambda]-\[Lambda]Minus)+2ArcTanh[(l+Sqrt[CC+l^2])/Sqrt[-CC] Tanh[Sqrt[-CC]/2 (\[Lambda]-\[Lambda]Minus)]]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*SphericalStar*)


NHEKSphericalStarMino[\[Lambda]_]:={T0+lStar/R0(\[Lambda]-\[Lambda]0), R0, \[CapitalPhi]0-3/4 lStar(\[Lambda]-\[Lambda]0)+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*Marginal*)


NHEKMarginalRadial[R_]:={\[Lambda]m+1/Sqrt[-CC]Log[Rm/R], T0+l/(Sqrt[-CC] R)-l/(Sqrt[-CC]Rm), \[CapitalPhi]0+3/4 l/Sqrt[-CC]Log[R/Rm]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NHEKMarginalMino[\[Lambda]_]:={T0+l/(Rm Sqrt[-CC])Exp[Sqrt[-CC](\[Lambda]-\[Lambda]m)], Rm Exp[-Sqrt[-CC](\[Lambda]-\[Lambda]m)], \[CapitalPhi]0-3/4 l(\[Lambda]-\[Lambda]m)+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*PlungingStar*)


NHEKPlungingStarRadial[R_]:={\[Lambda]0-1/lStar Sqrt[1+2lStar R/EE],T0+1/R Sqrt[1+2 lStar R/EE],\[CapitalPhi]0+3/4 Sqrt[1+2 lStar R/EE]-2 ArcTanh[Sqrt[1+2lStar R/EE]]+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NHEKPlungingStarMino[\[Lambda]_]:={T0 + 2 lStar (\[Lambda]-\[Lambda]0)/(EE (1-lStar^2(\[Lambda]-\[Lambda]0)^2)), EE/(2 lStar) (lStar^2(\[Lambda]-\[Lambda]0)^2-1), \[CapitalPhi]0 - 3/4 lStar(\[Lambda]-\[Lambda]0) + 2ArcTanh[lStar(\[Lambda]-\[Lambda]0)]+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*BoundedSubcritical*)


NHEKBoundedSubcriticalRadial[R_]:={\[Lambda]Plus + 1/Sqrt[CC] ArcCos[(CC R-EE l)/(EE Sqrt[CC+l^2])], Sqrt[EE^2+2EE l R-CC R^2]/(EE R), \[CapitalPhi]0-Log[(EE+l R+Sqrt[EE^2+2EE l R-CC R^2])/R]+3 l/(4Sqrt[CC])ArcTan[Sqrt[CC]Sqrt[EE^2+2EE l R-CC R^2]/(EE l-CC R)]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NHEKBoundedSubcriticalMino[\[Lambda]_]:={Sqrt[CC]/EE Sin[Sqrt[CC](\[Lambda]-\[Lambda]Plus)]/(l/Sqrt[CC+l^2]+Cos[Sqrt[CC](\[Lambda]-\[Lambda]Plus)]), EE/CC(l+Sqrt[CC+l^2]Cos[Sqrt[CC](\[Lambda]-\[Lambda]Plus)]), \[CapitalPhi]0 - 3l/4(\[Lambda]-\[Lambda]Plus)+2 ArcTanh[(l-Sqrt[CC+l^2])/Sqrt[CC]Tan[Sqrt[CC]/2(\[Lambda]-\[Lambda]Plus)]]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Section::Closed:: *)
(*Near-NHEK radial parametrizations of generating classes*)


v\[Kappa][R_]:=(e+\[Kappa] l)^2+2e l(R-\[Kappa])-CC(R^2-\[Kappa]^2)


FPlus[R_]:=(e R+\[Kappa](\[Kappa] l+Sqrt[v\[Kappa][R]]))^2
FMinus[R_]:=(e R+\[Kappa](\[Kappa] l-Sqrt[v\[Kappa][R]]))^2
GPlus[R_]:=(e+l R+Sqrt[v\[Kappa][R]])^2
GMinus[R_]:=(e+l R-Sqrt[v\[Kappa][R]])^2


(* ::Subsection:: *)
(*Spherical*)


NearNHEKSphericalMino[\[Lambda]_]:={t0+l/R0 (\[Lambda]-\[Lambda]0), R0, \[Phi]0-3/4 l (\[Lambda]-\[Lambda]0)+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*PlungingStarNull*)


NearNHEKPlungingStarNullRadial[R_]:={(R-Ri)/(\[Kappa] lStar),ti-1/(2\[Kappa])Log[(R^2-\[Kappa]^2)/(Ri^2-\[Kappa]^2)],\[Phi]i+3/(4\[Kappa])R+1/2 Log[(R-\[Kappa])/(R+\[Kappa])] +lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NearNHEKPlungingStarNullMino[\[Lambda]_]:={ti-1/(2\[Kappa]) Log[1+\[Kappa] lStar (\[Lambda]-\[Lambda]i) (\[Kappa] lStar(\[Lambda]-\[Lambda]i)-2Ri)/(Ri^2-\[Kappa]^2)], Ri-\[Kappa] lStar (\[Lambda]-\[Lambda]i), \[Phi]i-3/4 lStar (\[Lambda]-\[Lambda]i)+1/2 Log[(1-(\[Kappa] lStar (\[Lambda]-\[Lambda]i))/(Ri-\[Kappa]))/(1-(\[Kappa] lStar(\[Lambda]-\[Lambda]i))/Ri+\[Kappa])]+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*PlungingStar & BoundedStar*)


NearNHEKPlungingBoundedStarRadial[R_]:={-Sqrt[v\[Kappa][R]]/(e lStar),1/\[Kappa] ArcCosh[Abs[R+\[Kappa]^2 lStar/e]/Sqrt[R^2-\[Kappa]^2]],-3/(4e)Sqrt[v\[Kappa][R]]+ArcTanh[Sqrt[v\[Kappa][R]]/(e+lStar R)]+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NearNHEKPlungingBoundedStarMino[\[Lambda]_]:={-1/\[Kappa] ArcTanh[2\[Kappa] e lStar^2(\[Lambda]-\[Lambda]0)/(\[Kappa]^2lStar^2+e^2(lStar^2(\[Lambda]-\[Lambda]0)^2-1))], R0 + e lStar/2(\[Lambda]-\[Lambda]0)^2, -3/4 lStar(\[Lambda]-\[Lambda]0)+Sign[e]ArcTanh[2\[Kappa] e lStar^2(\[Lambda]-\[Lambda]0)/(\[Kappa]^2lStar^2+e^2(lStar^2(\[Lambda]-\[Lambda]0)^2-1))]+lStar Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


(* ::Subsection:: *)
(*BoundedSupercritical (Retrograde)*)


BoundedSupercriticalR[\[Lambda]_]:=1/CC(e l-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Minus)])


NearNHEKBoundedSupercriticalRadial[R_]:={\[Lambda]Minus+1/Sqrt[-CC]ArcCosh[(e l-R CC)/Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]],1/(4\[Kappa]) Log[FPlus[R]/FMinus[R]],3l/(4 Sqrt[-CC])Log[e l-CC R+Sqrt[-CC]Sqrt[v\[Kappa][R]]]-1/2Log[GPlus[R]/((CC+l^2)(R^2-\[Kappa]^2))]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NearNHEKBoundedSupercriticalMino[\[Lambda]_]:={NearNHEKBoundedSupercriticalRadial[BoundedSupercriticalR[\[Lambda]]][[2]], BoundedSupercriticalR[\[Lambda]], NearNHEKBoundedSupercriticalRadial[BoundedSupercriticalR[\[Lambda]]][[3]]}


(* ::Subsection:: *)
(*BoundedSubcritical*)


BoundedSubcriticalR[\[Lambda]_]:=1/CC(e l+Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]Cos[Sqrt[CC](\[Lambda]-\[Lambda]Plus)])


NearNHEKBoundedSubcriticalRadial[R_]:={\[Lambda]Plus+1/Sqrt[CC]ArcCos[(-e l+R CC)/Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]],1/(4\[Kappa]) Log[FPlus[R]/FMinus[R]],3l/(4 Sqrt[CC])ArcTan[Sqrt[CC]Sqrt[v\[Kappa][R]]/(e l-CC R)]-1/2Log[GPlus[R]/((CC+l^2)(R^2-\[Kappa]^2))]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NearNHEKBoundedSubcriticalMino[\[Lambda]_]:={NearNHEKBoundedSubcriticalRadial[BoundedSubcriticalR[\[Lambda]]][[2]], BoundedSubcriticalR[\[Lambda]], NearNHEKBoundedSubcriticalRadial[BoundedSubcriticalR[\[Lambda]]][[3]]}


(* ::Subsection:: *)
(*Plunging*)


PlungingR[\[Lambda]_]:=1/CC(e l-Sqrt[-(CC+l^2)(e^2+\[Kappa]^2 CC)](Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]f)]-Sqrt[2]Sinh[Sqrt[-CC](\[Lambda]-\[Lambda]f)]))


NearNHEKPlungingRadial[R_]:={\[Lambda]f+(-Log[(2  Sqrt[-(CC+l^2) (e^2+CC \[Kappa]^2)])/Sqrt[-CC]+2 Sqrt[2] Sqrt[((CC+l^2) (e^2+CC \[Kappa]^2))/CC]]+Log[(2  (e l-CC R))/Sqrt[-CC]+2 Sqrt[e^2+2 e l R-CC R^2+(CC+l^2) \[Kappa]^2]])/Sqrt[-CC],1/(4\[Kappa]) Log[FPlus[R]/FMinus[R]],3l/(4 Sqrt[-CC])Log[e l-CC R+Sqrt[-CC]Sqrt[v\[Kappa][R]]]-1/2Log[GPlus[R]/((CC+l^2)(R^2-\[Kappa]^2))]+l Subscript[\[CapitalPhi], \[Theta]][\[Lambda]]}


NearNHEKPlungingMino[\[Lambda]_]:={NearNHEKPlungingRadial[PlungingR[\[Lambda]]][[2]], PlungingR[\[Lambda]], NearNHEKPlungingRadial[PlungingR[\[Lambda]]][[3]]}


(* ::Subsection:: *)
(*Deflecting *)


DeflectingR[\[Lambda]_]:=1/CC(e l-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]Cosh[Sqrt[-CC](\[Lambda]-\[Lambda]Minus)])


NearNHEKDeflectingRadial[R_]:={\[Lambda]Minus+1/Sqrt[-CC]ArcCosh[(e l-R CC)/Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]],-1/(4\[Kappa]) Log[FPlus[R]/FMinus[R]],-3l/(4 Sqrt[-CC])Log[e l-CC R+Sqrt[-CC]Sqrt[v\[Kappa][R]]]+1/2Log[GPlus[R]/((CC+l^2)(R^2-\[Kappa]^2))]+l Subscript[\[CapitalPhi], \[Theta]]}


NearNHEKDeflectingMino[\[Lambda]_]:={NearNHEKDeflectingRadial[DeflectingR[\[Lambda]]][[2]], DeflectingR[\[Lambda]], NearNHEKDeflectingRadial[DeflectingR[\[Lambda]]][[3]]}


(* ::Section::Closed:: *)
(*Polar parametrizations*)


(* ::Subsection:: *)
(*Kerr*)


(* ::Subsubsection:: *)
(*Potential and roots*)


\[CapitalDelta]\[Theta]=1/2(1-(Q+l^2)/\[Epsilon]0);
zPlus=\[CapitalDelta]\[Theta]+Sign[\[Epsilon]0]Sqrt[\[CapitalDelta]\[Theta]^2+Q/\[Epsilon]0];
zMinus=\[CapitalDelta]\[Theta]-Sign[\[Epsilon]0]Sqrt[\[CapitalDelta]\[Theta]^2+Q/\[Epsilon]0];
zNull=Q/(Q+l^2);


\[CapitalTheta][z_]:=-l^2 z +(Q+\[Epsilon]0 z)(1-z)


mFunction[\[Lambda]_,\[CapitalLambda]_,\[Lambda]\[Theta]_]:=Floor[2/\[CapitalLambda] (\[Lambda]-\[Lambda]\[Theta])+1/2]


(* ::Subsubsection:: *)
(*Pendular*)


\[CapitalPsi]Plus[x_]:=ArcSin[x/Sqrt[zP]]
EllipticEPrime[x_,m_]:=1/(2m)(EllipticE[x,m]-EllipticF[x,m])
\[Lambda]i\[Theta]Pendular=\[Lambda]i-si\[Theta]/Sqrt[-\[Epsilon]0 zM]EllipticF[\[CapitalPsi]Plus[Cos[\[Theta]i]],zP/zM];
\[CapitalLambda]Pendular=4 /Sqrt[-\[Epsilon]0 zM]EllipticK[zP/zM];
cos\[Theta]Pendular[\[Lambda]_]:=si\[Theta] Sqrt[zP]JacobiSN[Sqrt[-\[Epsilon]0 zM](\[Lambda]-\[Lambda]i\[Theta]Pendular),zP/zM]
\[Theta]Pendular[\[Lambda]_]:=ArcCos[cos\[Theta]Pendular[\[Lambda]]]
\[CapitalPhi]hatPendular[\[Lambda]_,m_]:=1/Sqrt[-\[Epsilon]0 zM](2m EllipticPi[zP,zP/zM]+si\[Theta] (-1)^m EllipticPi[zP,\[CapitalPsi]Plus[cos\[Theta]Pendular[\[Lambda]]],zP/zM]-si\[Theta] EllipticPi[zP,\[CapitalPsi]Plus[Cos[\[Theta]i]],zP/zM])-\[Lambda]+\[Lambda]i
ThatPendular[\[Lambda]_,m_]:=-2 zP/Sqrt[-\[Epsilon]0 zM](2m EllipticEPrime[\[Pi]/2,zP/zM]+si\[Theta] (-1)^m EllipticEPrime[\[CapitalPsi]Plus[cos\[Theta]Pendular[\[Lambda]]],zP/zM]-si\[Theta] EllipticEPrime[\[CapitalPsi]Plus[Cos[\[Theta]i]],zP/zM])


(* ::Subsubsection::Closed:: *)
(*PendularNought*)


\[Lambda]i\[Theta]PendularNought=\[Lambda]i-Sqrt[z0/Q]ArcSin[Cos[\[Theta]i]/Sqrt[z0]];
cos\[Theta]PendularNought[\[Lambda]_]:=si\[Theta] Sqrt[z0]Sin[Sqrt[Q/z0](\[Lambda]-\[Lambda]i\[Theta]PendularNought)]
\[CapitalLambda]PendularNought=2\[Pi] Sqrt[z0/Q];
\[Theta]PendularNought[\[Lambda]_]:=ArcCos[cos\[Theta]PendularNought[\[Lambda]]]
\[CapitalPhi]hatPendularNought[\[Lambda]_,m_]:=Sqrt[z0/(Q(1-z0))](m \[Pi]+(-1)^m si\[Theta] ArcSin[Sqrt[(1-z0)/z0]Cot[\[Theta]PendularNought[\[Lambda]]]]-si\[Theta] ArcSin[Sqrt[(1-z0)/z0]Cot[\[Theta]i]])-\[Lambda]+\[Lambda]i
ThatPendularNought[\[Lambda]_,m_]:=1/2(z0(\[Lambda]-\[Lambda]i)-Sqrt[z0/Q]((-1)^m si\[Theta] cos\[Theta]PendularNought[\[Lambda]] Sqrt[z0-cos\[Theta]PendularNought[\[Lambda]]^2]-si\[Theta] Cos[\[Theta]i]Sqrt[z0-Cos[\[Theta]i]^2]))


(* ::Subsection:: *)
(*Near-horizon*)


(* ::Subsubsection:: *)
(*Dictionary*)


v\[Theta][z_]:=(Q+CNot z)(1-z)-l^2 z


nearHorizonPolarDictionary:={a->M,EE->l/(2M),\[Epsilon]0->CNot,\[CapitalTheta]->v\[Theta]}


\[CapitalDelta]\[Theta]=\[CapitalDelta]\[Theta]/.nearHorizonPolarDictionary;
zPlus=zPlus/.nearHorizonPolarDictionary;
zMinus=zMinus/.nearHorizonPolarDictionary;
zNull=zNull/.nearHorizonPolarDictionary;


(* ::Subsubsection:: *)
(*Pendular*)


cos\[Theta]PendularNH[\[Lambda]_]:=cos\[Theta]Pendular[\[Lambda]]/.nearHorizonPolarDictionary
\[CapitalLambda]PendularNH:=\[CapitalLambda]Pendular/.nearHorizonPolarDictionary
\[Lambda]i\[Theta]PendularNH:=\[Lambda]i\[Theta]Pendular/.nearHorizonPolarDictionary
\[Theta]PendularNH[\[Lambda]_]:=ArcCos[cos\[Theta]PendularNH[\[Lambda]]]
\[CapitalPhi]PendularNH[\[Lambda]_,m_]:=\[CapitalPhi]hatPendular[\[Lambda],m]-1/4 ThatPendular[\[Lambda],m]/.nearHorizonPolarDictionary


(* ::Subsubsection::Closed:: *)
(*PendularNought*)


cos\[Theta]PendularNoughtNH[\[Lambda]_]:=cos\[Theta]PendularNought[\[Lambda]]/.nearHorizonPolarDictionary
\[CapitalLambda]PendularNoughtNH:=\[CapitalLambda]PendularNought/.nearHorizonPolarDictionary
\[Lambda]i\[Theta]PendularNoughtNH:=\[Lambda]i\[Theta]PendularNought/.nearHorizonPolarDictionary
\[Theta]PendularNoughtNH[\[Lambda]_]:=ArcCos[cos\[Theta]PendularNoughtNH[\[Lambda]]]
\[CapitalPhi]PendularNoughtNH[\[Lambda]_,m_]:=\[CapitalPhi]hatPendularNought[\[Lambda],m]-1/4ThatPendularNought[\[Lambda],m]/.nearHorizonPolarDictionary


(* ::Subsection:: *)
(**)


(* ::Section::Closed:: *)
(*Flips*)


LeftRightFlip:={T->-T, T0->-T0, Ti->-Ti, \[CapitalPhi]->-\[CapitalPhi], \[CapitalPhi]0->-\[CapitalPhi]0, \[CapitalPhi]i->-\[CapitalPhi]i, \[Lambda]->-\[Lambda], \[Lambda]i->-\[Lambda]i, \[Lambda]Minus->-\[Lambda]Minus, \[Lambda]0->-\[Lambda]0, \[Lambda]Plus->-\[Lambda]Plus, si\[Theta]->-si\[Theta], siR->-siR}


DiagFlip:={R->-R, Ri->-Ri, R0->-R0, RPlus->-RPlus, RMinus->-RMinus, \[CapitalPhi]->-\[CapitalPhi], \[CapitalPhi]0->-\[CapitalPhi]0, \[CapitalPhi]i->-\[CapitalPhi]i, l->-l, lStar->-lStar, lNot->-lNot, siR->-siR}


LeftRightFlipMino[assos_]:=Module[{},
	Return[{-assos["Trajectory"][[1]], assos["Trajectory"][[2]], assos["Trajectory"][[3]], -assos["Trajectory"][[4]]}]
];

DiagFlipMino[assos_]:=Module[{},
	Return[{assos["Trajectory"][[1]], -assos["Trajectory"][[2]], assos["Trajectory"][[3]], -assos["Trajectory"][[4]]}]
];

LeftRightFlipRadial[assos_]:=Module[{},
	Return[{-assos["Trajectory"][[1]], -assos["Trajectory"][[2]], assos["Trajectory"][[3]], -assos["Trajectory"][[4]]}]
];

DiagFlipRadial[assos_]:=Module[{},
	Return[{assos["Trajectory"][[1]], assos["Trajectory"][[2]], assos["Trajectory"][[3]], -assos["Trajectory"][[4]]}]
];


(* ::Section:: *)
(*Get the symbolic form of the orbits: GetOrbit method*)


Simplification[rule_, ass_, expr_]:=Module[{},
	Switch[rule,
	"Simplify", Return[Assuming[ass, Simplify[expr]]];,
	"FullSimplify", Return[Assuming[ass, FullSimplify[expr]]];,
	_, Return[expr];
	];
];


Options[GetOrbit] = {"OutStyle" -> "BHPT", "SimplificationRule" -> "Simplify", "Retrograde" -> False, "ReplaceLStar" -> False, "ReplaceC" -> False, "ReplaceLNot" -> False, "ReplaceCNot" -> False, "ReplaceRoots" -> False, "ReplaceTurningPoints" -> False, "CosTheta" -> True, "ReplacePhiTheta" -> False};

GetOrbit[spacetime_, parametrization_, radial_, polar_, type_, OptionsPattern[]]:=Module[{style, toFlip, rule, cosTheta, polarAssumptions, radialTrajectory, assos, replaceLStarRule, replaceCRule, radialAssumptions,assumptions, mu, trajectory, parameters, energy, momentum, carter, 
initialData, criticalMomentum, LNot, casimir,  kappa, criticalRadius, radialPotential, radialRoots, polarPotential, polarRoots, CCNot, PhiTheta, frequencies, turningPoints, rootDictionary, paramForm},

Clear[T, R, \[CapitalPhi], \[Lambda], EE, l, Q, CC, lStar, lNot, CNot, \[Lambda]i, zP, zM, z0, m, \[Lambda]Plus, \[Lambda]Minus, si\[Theta], \[Lambda]0, T0, R0, Rm, \[Lambda]m, \[CapitalPhi]0, t, \[Phi], e, \[Lambda]f, \[Theta], t, r];

(* Define the type of geodesic *)
Switch[type,
	"Null", mu = 0,
	"Timelike", mu = \[Mu], 
	_, Print["Unknown geodesic type: " <> type]; Return[$Failed];
];

parameters = {energy, momentum, carter, mu};
criticalMomentum = 2/Sqrt[3]Sqrt[M^2 mu^2+Q];
casimir = 3/4(lStar^2-l^2);
LNot = 2 M mu;
CCNot = (l^2-lNot^2)/4;

(* Replacement rules *)
If[OptionValue["ReplaceLStar"], lStar = criticalMomentum;];
If[OptionValue["ReplaceC"], CC = casimir;];
If[OptionValue["ReplaceLNot"], lNot = LNot;];
If[OptionValue["ReplaceCNot"], CNot = CCNot;];
If[OptionValue["ReplaceTurningPoints"], m = turningPoints;];
If[OptionValue["ReplaceRoots"], rootDictionary = {zP->zPlus, zM->zMinus, z0->zNull};, rootDictionary = {};];

(* Radial Part *)
Switch[spacetime,
	"NHEK",
	(* NHEK geodesics *)
	kappa = 0;
	criticalRadius = -energy/momentum;
	
	Switch[radial,
	"Deflecting",
		(* Parameters *)
		energy = EE;
		momentum = l;
		carter = Q;
		radialAssumptions = {EE<0 && l>lStar && CC<0 && R0>0};
		
		initialData = {si\[Theta], 0, SubPlus[R], \[Theta]i, \[CapitalPhi]0, \[Lambda]Plus};
		\[Lambda]i = \[Lambda]Plus;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions,Simplify[{SubMinus[R] -> EE l/CC - Abs[EE]/Abs[CC]Sqrt[l^2+CC], SubPlus[R] -> EE l/CC + Abs[EE]/Abs[CC]Sqrt[l^2+CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NHEKDeflectingRadial[R];
		,"Mino",
			radialTrajectory = NHEKDeflectingMino[\[Lambda]];
		];
		
	,"Plunging",
		(* Parameters *)
		energy = EE;
		momentum = l;
		carter = Q;
		radialAssumptions = {EE>0 && l>lStar && CC<0};
		
		initialData = {si\[Theta], 0, SubMinus[R], \[Theta]i, \[CapitalPhi]0, \[Lambda]Minus};
		\[Lambda]i = \[Lambda]Minus;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions,Simplify[{SubMinus[R] -> EE l/CC - Abs[EE]/Abs[CC]Sqrt[l^2+CC], SubPlus[R] -> EE l/CC + Abs[EE]/Abs[CC]Sqrt[l^2+CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NHEKPlungingRadial[R];
		,"Mino",
			radialTrajectory = NHEKPlungingMino[\[Lambda]];
		];
	
	,"Outward",
		(* FLIP *)
		toFlip = GetOrbit["NHEK", parametrization, "Plunging", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Radial Class"] = "Outward";
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = LeftRightFlipMino[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = LeftRightFlipRadial[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];
		]
		
	,"SphericalStar",
		(* Parameters *)
		energy = 0;
		momentum = lStar;
		carter = Q;
		casimir = 0;
		radialAssumptions = {EE== 0 && l==lStar && CC==0};
		
		initialData = {si\[Theta], Ti, R0, \[Theta]i, \[CapitalPhi]0, \[Lambda]i};
		T0 = Ti;
		\[CapitalPhi]0 = \[CapitalPhi]i;
		\[Lambda]0 = \[Lambda]i;
		
		(* Roots and potential *)
		radialRoots = {Undefined};
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			Print["Radial parametrization impossible for spherical geodesics !"];
			Return[$Failed];
		,"Mino",
			radialTrajectory = NHEKSphericalStarMino[\[Lambda]];
		];
	
	,"Marginal",
		(* Parameters *)
		energy = 0;
		momentum = l;
		carter = Q;
		radialAssumptions = {EE==0 && Rm>0 && CC<0 && l>lStar};
		
		initialData = {si\[Theta], Ti, Rm, \[Theta]i, \[CapitalPhi]0, \[Lambda]m};
		T0 = Ti;
		\[CapitalPhi]0 = \[CapitalPhi]i;
		\[Lambda]i = \[Lambda]m;
		
		(* Roots and potential *)
		radialRoots = {Undefined};
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NHEKMarginalRadial[R];
		,"Mino",
			radialTrajectory = NHEKMarginalMino[\[Lambda]];
		];
	
	,"PlungingStar",
		(* Parameters *)
		energy = EE;
		momentum = lStar;
		l = lStar;
		carter = Q;
		casimir = 0;
		radialAssumptions = {EE>0 && CC==0 && l==lStar};
		
		initialData = {si\[Theta], T0, R0, \[Theta]i, \[CapitalPhi]0, \[Lambda]0};
		\[Lambda]i = \[Lambda]0;
		
		(* Roots and potential *)
		radialRoots = {R0->-EE/(2 lStar)};
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NHEKPlungingStarRadial[R];
		,"Mino",
			radialTrajectory = NHEKPlungingStarMino[\[Lambda]];
		];
	
	,"OutwardStar",
		(* FLIP *)
		toFlip = GetOrbit["NHEK", parametrization, "PlungingStar", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Radial Class"] = "OutwardStar";
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = LeftRightFlipMino[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = LeftRightFlipRadial[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];
		]
	
	,"BoundedSubcritical",
		(* Parameters *)
		energy = EE;
		momentum = l;
		carter = Q;
		radialAssumptions = {EE>0 && l^2<lStar^2 && CC>0};
		
		initialData = {si\[Theta], 0, SubPlus[R], \[Theta]i, \[CapitalPhi]0, \[Lambda]Plus};
		\[Lambda]i = \[Lambda]Plus;
		\[Lambda]0 = \[Lambda]Plus;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions,Simplify[{SubMinus[R] -> EE l/CC - Abs[EE]/Abs[CC]Sqrt[l^2+CC], SubPlus[R] -> EE l/CC + Abs[EE]/Abs[CC]Sqrt[l^2+CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[EE^2+2EE l R-CC R^2]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NHEKBoundedSubcriticalRadial[R];
		,"Mino",
			radialTrajectory = NHEKBoundedSubcriticalMino[\[Lambda]];
		];
	
	,"BoundedStarMinus",
		(* FLIP *)
		toFlip = GetOrbit["NHEK", parametrization, "PlungingStar", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> True, "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.DiagFlip;
		toFlip["Radial Class"] = "BoundedStarMinus";
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Geodesic Parameters"] = toFlip["Geodesic Parameters"]/.lStar->-lStar;
		toFlip["Assumptions"] = toFlip["Assumptions"]/.lStar->-lStar;
		toFlip["Radial Roots"] = toFlip["Radial Roots"]/.lStar->-lStar;
		toFlip["Radial Potential"] = toFlip["Radial Potential"]/.lStar->-lStar;
		toFlip["Angular Momentum"] = toFlip["Angular Momentum"]/.lStar->-lStar;
		toFlip["Critical Radius"] = toFlip["Critical Radius"]/.lStar->-lStar;
		toFlip["Frequencies"] = toFlip["Frequencies"]/.lStar->-lStar;
		toFlip["Turning Points"] = toFlip["Turning Points"]/.lStar->-lStar;
		
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = DiagFlipMino[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = DiagFlipRadial[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];
		]
		
	,"BoundedSupercritical",
		(* FLIP *)
		toFlip = GetOrbit["NHEK", parametrization, "Plunging", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.DiagFlip;
		toFlip["Radial Class"] = "BoundedSupercritical";
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Geodesic Parameters"] = toFlip["Geodesic Parameters"]/.l->-l;
		toFlip["Assumptions"] = toFlip["Assumptions"]/.l->-l;
		toFlip["Radial Roots"] = toFlip["Radial Roots"]/.l->-l;
		toFlip["Radial Potential"] = toFlip["Radial Potential"]/.l->-l;
		toFlip["Angular Momentum"] = toFlip["Angular Momentum"]/.l->-l;
		toFlip["Critical Radius"] = toFlip["Critical Radius"]/.l->-l;
		toFlip["Frequencies"] = toFlip["Frequencies"]/.l->-l;
		toFlip["Turning Points"] = toFlip["Turning Points"]/.l->-l;
		
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = DiagFlipMino[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = DiagFlipRadial[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];
		]
	
	,_ ,Print["Unknown radial class: " <> radial]; Return[$Failed];	
	];
	
	(****************)
	
	,"Near-NHEK", 
	(* Near-NHEK geodesics *)
	kappa = \[Kappa];
	criticalRadius = -energy/momentum;
	
	Switch[radial,
	"Spherical",
	(* Parameters *)
		energy = -\[Kappa] Sqrt[-CC];
		momentum = l;
		carter = Q;
		casimir = CC;
		radialAssumptions = { e==-\[Kappa] Sqrt[-CC] && l>lStar && CC<0};
		
		initialData = {si\[Theta], t0, R0, \[Theta]i, \[Phi]0, \[Lambda]0}/.R0 ->  \[Kappa] l/Sqrt[-CC];
		\[Lambda]i = \[Lambda]0;
		
		(* Roots and potential *)
		radialRoots = {Undefined};
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][\[Kappa] l/Sqrt[-CC]]/.{e->energy}]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			Print["Radial parametrization impossible for spherical geodesics !"];
			Return[$Failed];
		,"Mino",
			radialTrajectory = NearNHEKSphericalMino[\[Lambda]]/.R0 ->  \[Kappa] l/Sqrt[-CC];
		];
	
	,"PlungingStarNull",
		(* Parameters *)
		energy = 0;
		momentum = lStar;
		l = lStar;
		carter = Q;
		casimir = 0;
		radialAssumptions = { e==0 && l==lStar && CC==0 && Ri>0};
		
		initialData = {si\[Theta], ti, Ri, \[Theta]i, \[Phi]i, \[Lambda]i};
		
		(* Roots and potential *)
		radialRoots = {Undefined};
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKPlungingStarNullRadial[R];
		,"Mino",
			radialTrajectory = NearNHEKPlungingStarNullMino[\[Lambda]];
		];
	
	,"OutwardStarNull",
	(* FLIP *)
		toFlip = GetOrbit["Near-NHEK", parametrization, "PlungingStarNull", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Radial Class"] = "OutwardStarNull";
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = LeftRightFlipMino[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = LeftRightFlipRadial[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];
		]
	
	,"PlungingStar",
		(* Parameters *)
		energy = e;
		momentum = lStar;
		l = lStar;
		carter = Q;
		casimir = 0;
		radialAssumptions = { e>0 && l==lStar && CC==0};
		
		initialData = {si\[Theta], 0, R0, \[Theta]i, 0, \[Lambda]i}/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		
		(* Roots and potential *)
		radialRoots = {Subscript[R, 0]->-(e^2+\[Kappa]^2 lStar^2)/(2e lStar)};
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKPlungingBoundedStarRadial[R]/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		,"Mino",
			radialTrajectory = NearNHEKPlungingBoundedStarMino[\[Lambda]]/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		];
	
	,"OutwardStar",
		(* FLIP *)
		toFlip = GetOrbit["Near-NHEK", parametrization, "PlungingStar", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Radial Class"] = "OutwardStar";
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = LeftRightFlipMino[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = LeftRightFlipRadial[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];s
		]
	
	,"BoundedStar",
		(* Parameters *)
		energy = e;
		momentum = lStar;
		carter = Q;
		casimir = 0;
		radialAssumptions = { e<0 && l==lStar && CC==0};
		
		initialData = {si\[Theta], 0, R0, \[Theta]i, 0, \[Lambda]i}/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		
		(* Roots and potential *)
		radialRoots = {Subscript[R, 0]->-(e^2+\[Kappa]^2 lStar^2)/(2e lStar)};
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKPlungingBoundedStarRadial[R]/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		,"Mino",
			radialTrajectory = NearNHEKPlungingBoundedStarMino[\[Lambda]]/.R0 -> -(e^2+\[Kappa]^2 lStar^2)/(2e lStar);
		];
		
	,"BoundedStarMinus",
		(* FLIP *)
		toFlip = GetOrbit["Near-NHEK", parametrization, "PlungingStar", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> True, "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.DiagFlip;
		toFlip["Radial Class"] = "BoundedStarMinus";
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Geodesic Parameters"] = toFlip["Geodesic Parameters"]/.lStar->-lStar;
		toFlip["Assumptions"] = {toFlip["Assumptions"][[1]]/.lStar->-lStar && e>\[Kappa] lStar};
		toFlip["Radial Roots"] = toFlip["Radial Roots"]/.lStar->-lStar;
		toFlip["Radial Potential"] = toFlip["Radial Potential"]/.lStar->-lStar;
		toFlip["Angular Momentum"] = toFlip["Angular Momentum"]/.lStar->-lStar;
		toFlip["Critical Radius"] = toFlip["Critical Radius"]/.lStar->-lStar;
		toFlip["Frequencies"] = toFlip["Frequencies"]/.lStar->-lStar;
		toFlip["Turning Points"] = toFlip["Turning Points"]/.lStar->-lStar;
		
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = DiagFlipMino[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = DiagFlipRadial[toFlip]/.DiagFlip;
			Return[toFlip/.BHPTStyle];
		]
	
	,"BoundedSupercritical",
		(* Parameters *)
		energy = e;
		momentum = l;
		carter = Q;
		casimir = CC;
		radialAssumptions = {e>-\[Kappa] l && l<-lStar && CC<0};
		
		initialData = {si\[Theta], 0, RMinus, \[Theta]i, 0, \[Lambda]Minus};
		\[Lambda]i = \[Lambda]Minus;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions, Simplify[{SubMinus[R]-> e l/CC-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC], SubPlus[R]-> e l/CC+Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKBoundedSupercriticalRadial[R];
		,"Mino",
			radialTrajectory = NearNHEKBoundedSupercriticalMino[\[Lambda]];
		];
	
	,"BoundedSubcritical",
		(* Parameters *)
		energy = e;
		momentum = l;
		carter = Q;
		casimir = CC;
		radialAssumptions = {e>-\[Kappa] l && l^2<lStar^2 && CC>0};
		
		initialData = {si\[Theta], 0, RPlus, \[Theta]i, 0, \[Lambda]Plus};
		\[Lambda]i = \[Lambda]Plus;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions, Simplify[{SubMinus[R]-> e l/CC-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC], SubPlus[R]-> e l/CC+Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKBoundedSubcriticalRadial[R];
		,"Mino",
			radialTrajectory = NearNHEKBoundedSubcriticalMino[\[Lambda]];
		];

	,"Plunging",
		(* Parameters *)
		energy = e;
		momentum = l;
		carter = Q;
		casimir = CC;
		radialAssumptions = {e>-\[Kappa] Sqrt[-CC] && l^2>lStar^2 && CC<0 && Rf==1/CC(e l-Sqrt[-(CC+l^2)(e^2+\[Kappa]^2 CC)])};
		
		initialData = {si\[Theta], 0, Rf, \[Theta]i, 0, \[Lambda]f};
		\[Lambda]i = \[Lambda]f;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions, Simplify[{SubMinus[R]-> e l/CC-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC], SubPlus[R]-> e l/CC+Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKPlungingRadial[R];
		,"Mino",
			radialTrajectory = NearNHEKPlungingMino[\[Lambda]];
		];
		
	,"Outward",
		(* FLIP *)
		toFlip = GetOrbit["Near-NHEK", parametrization, "Plunging", polar, type, "OutStyle"-> "None", "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], "ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNot"], "ReplaceCNot" -> OptionValue["ReplaceCNot"], "ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]; 
		toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"] = toFlip["\!\(\*SubscriptBox[\(\[CapitalPhi]\), \(\[Theta]\)]\)"]/.LeftRightFlip;
		toFlip["Radial Class"] = "Outward";
		Switch[parametrization,
			"Mino", 
			toFlip["Trajectory"] = LeftRightFlipMino[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];,
			"Radial",
			toFlip["Trajectory"] = LeftRightFlipRadial[toFlip]/.LeftRightFlip;
			Return[toFlip/.BHPTStyle];
		]
	
	,"Deflecting",
		(* Parameters *)
		energy = e;
		momentum = l;
		carter = Q;
		casimir = CC;
		radialAssumptions = {e<-\[Kappa] Sqrt[-CC] && l^2>lStar^2 && CC<0};
		
		initialData = {si\[Theta], 0, RMinus, \[Theta]i, 0, \[Lambda]Minus};
		\[Lambda]i = \[Lambda]f;
		
		(* Roots and potential *)
		radialRoots = Assuming[radialAssumptions, Simplify[{SubMinus[R]-> e l/CC-Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC], SubPlus[R]-> e l/CC+Sqrt[(CC+l^2)(e^2+\[Kappa]^2 CC)]/Abs[CC]}]];
		radialPotential = Assuming[radialAssumptions, Simplify[v\[Kappa][R]]];
		
		(* Trajectory *)
		Switch[parametrization,
		"Radial",
			radialTrajectory = NearNHEKDeflectingRadial[R];
		,"Mino",
			radialTrajectory = NearNHEKDeflectingMino[\[Lambda]];
		];
	
	,_ ,Print["Unknown radial class: " <> radial]; Return[$Failed];
	];
	
	(*********************)
	,_ ,Print["Unknown spacetime: " <> spacetime]; Return[$Failed];
];

(* Polar part *)
Switch[polar,
"Pendular",
	(*** Pendular ***)
	polarAssumptions = {Q>0 && l^2!=lNot^2 && CNot!=0 && si\[Theta]^2==1};
	polarPotential = Assuming[polarAssumptions, Simplify[v\[Theta][z]]];
	polarRoots = {SubMinus[z] -> zMinus, SubPlus[z] -> zPlus}/.nearHorizonPolarDictionary;
	frequencies = Assuming[polarAssumptions,Simplify[{"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalLambda]PendularNH}/.rootDictionary]];
	PhiTheta = \[CapitalPhi]PendularNH[\[Lambda], m]/.rootDictionary;
	turningPoints = mFunction[\[Lambda], \[CapitalLambda]PendularNH, \[Lambda]i\[Theta]PendularNH]/.rootDictionary;
	cosTheta = cos\[Theta]PendularNH[\[Lambda]]/.rootDictionary;
	
,"PendularNought",
	(*** PendularNot ***)
	CCNot = 0;
	If[OptionValue["Retrograde"], 
		polarAssumptions = {Q>0 && l==-lNot  && si\[Theta]^2==1};
		l = -lNot;,
		polarAssumptions = {Q>0 && l==lNot && si\[Theta]^2==1};
		l = lNot; 
	];
	polarPotential = Assuming[polarAssumptions, Simplify[v\[Theta][z]]];
	polarRoots = {"\!\(\*SubscriptBox[\(z\), \(0\)]\)" -> zNull}/.nearHorizonPolarDictionary;
	frequencies = Assuming[polarAssumptions,Simplify[{"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> \[CapitalLambda]PendularNoughtNH}/.rootDictionary]];
	PhiTheta = \[CapitalPhi]PendularNoughtNH[\[Lambda], m]/.rootDictionary;
	turningPoints = mFunction[\[Lambda], \[CapitalLambda]PendularNoughtNH, \[Lambda]i\[Theta]PendularNoughtNH]/.rootDictionary;
	cosTheta = cos\[Theta]PendularNoughtNH[\[Lambda]]/.rootDictionary;
	
,"Equatorial",
	(*** Equatorial ***)
	carter = 0;
	Q = 0;
	If[OptionValue["Retrograde"], 
		polarAssumptions = {Q==0 && 0>=l>=-lNot && si\[Theta]^2==1 && z==0};,
		polarAssumptions = {Q==0 && 0<=l<=lNot && si\[Theta]^2==1 && z==0};
	];
	polarPotential = Assuming[polarAssumptions, Simplify[v\[Theta][z]]];
	polarRoots = Undefined;
	frequencies = {"\!\(\*SubscriptBox[\(\[CapitalUpsilon]\), \(\[Theta]\)]\)" -> 0};
	PhiTheta = 0;
	turningPoints = Undefined;
	cosTheta = 0;
	
,_ ,Print["Unknown polar class: " <> polar]; Return[$Failed];
];

(* Putting all pieces together *)
assumptions = Join[radialAssumptions, polarAssumptions, {lStar>0}];

Switch[parametrization,
"Mino",
	trajectory = {radialTrajectory[[1]] , radialTrajectory[[2]], cosTheta, radialTrajectory[[3]]};
,"Radial",
	trajectory = {radialTrajectory[[1]], radialTrajectory[[2]], cosTheta, radialTrajectory[[3]]};
];

 If[Not[OptionValue["CosTheta"]], trajectory[[3]] = ArcCos[trajectory[[3]]];];
 If[OptionValue["ReplacePhiTheta"], trajectory = trajectory/.Subscript[\[CapitalPhi], \[Theta]][\[Lambda]] -> PhiTheta;];
 If[polar=="Equatorial", trajectory = trajectory/.Subscript[\[CapitalPhi], \[Theta]][\[Lambda]] -> 0;];
 If[parametrization == "Radial", trajectory = trajectory/. \[Lambda] -> radialTrajectory[[1]];];

(* Parametrization form *)
paramForm = Switch[parametrization,
	"Mino",
		If[OptionValue["CosTheta"], "{T(\[Lambda]), R(\[Lambda]), cos \[Theta](\[Lambda]), \[CapitalPhi](\[Lambda])}", "{T(\[Lambda]), R(\[Lambda]), \[Theta](\[Lambda]), \[CapitalPhi](\[Lambda])}"]
	,"Radial",
		If[OptionValue["CosTheta"], "{\[Lambda](R), T(R), cos \[Theta](R), \[CapitalPhi](R)}", "{\[Lambda](R), T(R), \[Theta](R), \[CapitalPhi](R)}"]
];

(* Output *)
rule = OptionValue["SimplificationRule"];

If[OptionValue["OutStyle"]=="BHPT", style = BHPTStyle;, style = {}; ];

assos = Association[
	"Spacetime" -> spacetime,
	"Radial Class" -> radial,
	"Polar Class" -> polar,
	"Type" -> type,
	"Parametrization" -> parametrization,
	Switch[parametrization,
		"Mino", "Trajectory" -> ToFunction[trajectory, \[Lambda], style, rule, assumptions],
		"Radial", "Trajectory" -> ToFunction[trajectory, R, style, rule, assumptions]], 
	"Geodesic Parameters" -> parameters,
	"Initial Data" -> initialData,
	"Assumptions" -> assumptions,
	"Radial Roots" -> radialRoots,
	"Radial Potential" -> ToFunction[radialPotential, R, style, rule, assumptions],
	"LStar" -> criticalMomentum,
	"Casimir" -> casimir,
	"Energy" -> energy,
	"Angular Momentum" -> momentum,
	"Carter Constant" -> carter,
	"LNought" -> LNot,
	"Kappa" -> kappa,
	"Critical Radius" -> criticalRadius,
	"Polar Potential" -> ToFunction[polarPotential, z, style, rule, assumptions],
	"Polar Roots" -> Simplification[rule, assumptions, polarRoots],
	"PhiTheta" -> Simplification[rule, assumptions, PhiTheta],
	"Frequencies" -> Simplification[rule, assumptions, frequencies],
	"CNought" -> CCNot,
	"Turning Points" -> ToFunction[turningPoints, \[Lambda], style, rule, assumptions],
	"Parametrization Form" -> paramForm
];

Return[assos/.style];

];

ToFunction[expr_, var_, style_, rule_, assumptions_]:=Function[var, Simplification[rule, assumptions, expr]/.style]


(* ::Section:: *)
(*Output: NearHorizonGeoOrbit & NearHorizonGeoOrbitFunction*)


(* SYMBOLIC TRAJECTORY *)

Options[NearHorizonGeoOrbitClass] = {"Parametrization" -> "Mino", "Type" -> "Timelike", "SimplificationRule" -> "Simplify", "Retrograde" -> False, "ReplaceLStar" -> False, "ReplaceC" -> False, "ReplaceLNought" -> False, 
"ReplaceCNought" -> False, "ReplaceRoots" -> False, "ReplaceTurningPoints" -> False, "CosTheta" -> True, "ReplacePhiTheta" -> False, "Style"->"BHPT"};

NearHorizonGeoOrbitClass[spacetime_String, radial_String, polar_String, OptionsPattern[]]:=Module[{assoc, type, parametrization},
	type = OptionValue["Type"];
	parametrization = OptionValue["Parametrization"];
	
	assoc = GetOrbit[spacetime, parametrization, radial, polar, type, "OutStyle"-> OptionValue["Style"], "SimplificationRule" -> OptionValue["SimplificationRule"], "Retrograde" -> OptionValue["Retrograde"], 
	"ReplaceLStar" -> OptionValue["ReplaceLStar"], "ReplaceC" -> OptionValue["ReplaceC"], "ReplaceLNot" -> OptionValue["ReplaceLNought"], "ReplaceCNot" -> OptionValue["ReplaceCNought"], 
	"ReplaceRoots" -> OptionValue["ReplaceRoots"], "ReplaceTurningPoints" -> OptionValue["ReplaceTurningPoints"], "CosTheta" -> OptionValue["CosTheta"], 
	"ReplacePhiTheta" -> OptionValue["ReplacePhiTheta"]]//Quiet;
	
	If[assoc["Spacetime"]==$Failed["Spacetime"], Return["An error occured in the process !"]; Return[$Failed];];
	
	Return[NearHorizonGeoOrbitFunction[assoc]];
];


(* NUMERICAL TRAJECTORY *)

Options[NearHorizonGeoOrbit] = {"RadialMotion" -> "Ingoing", "ExplicitMass" -> 1, "Parametrization" -> "Mino", "SimplificationRule" -> "Simplify", "CosTheta" -> False, "Numerical" -> True};

NearHorizonGeoOrbit[st_String, EE_, l_, Q_, mu_, \[Kappa]_:\[Kappa], initData:{si\[Theta]_, Ti_, Ri_, \[Theta]i_, \[CapitalPhi]i_, \[Lambda]i_}:{1, 0, Ri, \[Pi]/2, 0, 0}, OptionsPattern[]]:=Module[{CC, paramRule, mass, parametrization, assoc, lStar, isRetro, lNot, radial, polar, type},
	
	mass = OptionValue["ExplicitMass"];
	parametrization = OptionValue["Parametrization"];
	
	(* Determine the spacetime *)
	If[\[Kappa]<0, Print["Uncorrect input parameters: \[Kappa] must be positive !"]; Return[$Failed];];
	
	If[Q<0, Print["Uncorrect input parameters: Q must be positive !"]; Return[$Failed];];
	
	If[st=="Near-NHEK" && EE<=-\[Kappa] l, Print["Uncorrect input parameters: \[Epsilon]>-\[Kappa] \[ScriptCapitalL] must be satisfied !"]; Return[$Failed];];
	
	(* Determine the type *)
	Which[
		mu>0, type = "Timelike";,
		mu==0, type = "Null";,
		_, Print["Uncorrect input parameters: \[Mu] must be positive or null !"]; Return[$Failed];
	];
	
	(* Determine the radial class *)
	lStar := 2/Sqrt[3]Sqrt[mass^2 mu^2+Q];
	CC := 3/4(lStar^2-l^2);
	
	Switch[st,
		"NHEK",
		Which[
			l>lStar, Which[
				EE<0, radial = "Deflecting";,
				EE==0, radial = "Marginal";,
				EE>0, Switch[OptionValue["RadialMotion"],
								"Ingoing", radial = "Plunging";,
								"Outgoing", radial = "Outward";
								,_, Print["Unable to choose between Plunging and Outward classes. Please provide 'Radial Motion' option !"]; Return[$Failed];
							];
				,_, Print["Uncorrect input parameters: unable to determine the radial class !"]; Return[$Failed];
			];
			
			,l==lStar, Which[
				EE==0, radial = "SphericalStar";,
				EE>0, Switch[OptionValue["RadialMotion"],
								"Ingoing", radial = "PlungingStar";,
								"Outgoing", radial = "OutwardStar";
								,_, Print["Unable to choose between Plunging and Outward classes. Please provide 'Radial Motion' option !"]; Return[$Failed];
							];
				,_, Print["Uncorrect input parameters: unable to determine the radial class !"]; Return[$Failed];
			];
			
			,l^2<lStar^2 && EE>0, radial = "BoundedSubcritical";
			
			,l==-lStar && EE>0, radial = "BoundedStarMinus";
			
			,l<-lStar && EE>0, radial = "BoundedSupercritical";
			
			,_, Print["Uncorrect input parameters: unable to determine the radial class !"]; Return[$Failed];
		];
		,"Near-NHEK",
		Which[
		l>lStar && EE>-\[Kappa] l, 
			Which[
				EE<-\[Kappa] Sqrt[-CC], 
					radial = "Deflecting";				
				,EE==-\[Kappa] Sqrt[-CC],
					radial = "Spherical";
				,EE>-\[Kappa] Sqrt[-CC],
					Switch[OptionValue["RadialMotion"],
								"Ingoing", radial = "Plunging";,
								"Outgoing", radial = "Outward";
								,_, Print["Unable to choose between Plunging and Outward classes. Please provide 'Radial Motion' option !"]; Return[$Failed];
							];
				,_, Print["Uncorrect input parameters: unable to determine the radial class !"]; Return[$Failed];
			];
		
		,l==lStar && EE>-\[Kappa] l,
			Which[
				EE<0, radial = "BoundedStar";
				
				,EE==0, Switch[OptionValue["RadialMotion"],
								"Ingoing", radial = "PlungingStarNull";,
								"Outgoing", radial = "OutwardStarNull";
								,_, Print["Unable to choose between Plunging and Outward classes. Please provide 'Radial Motion' option !"]; Return[$Failed];
							];
				
				,EE>0, Switch[OptionValue["RadialMotion"],
								"Ingoing", radial = "PlungingStar";,
 								"Outgoing", radial = "OutwardStar";
								,_, Print["Unable to choose between Plunging and Outward classes. Please provide 'Radial Motion' option !"]; Return[$Failed];
							];
			];
		
		,l^2<lStar^2 && EE>-\[Kappa] l, radial = "BoundedSubcritical";
		
		,l==-lStar && EE>-\[Kappa] l, radial = "BoundedStarMinus";
		
		,l<-lStar && EE>-\[Kappa] l, radial = "BoundedSupercritical";
		
		,_, Print["Uncorrect input parameters: unable to determine the radial class !"]; Return[$Failed];		
		];
	];
	
	
	(* Determine the polar class *)
	lNot = 2 mass mu;
	
	Which[
		l==lNot && Q>0, polar = "PendularNought"; isRetro = False;,
		l==-lNot && Q>0, polar = "PendularNought"; isRetro = True;,
		l^2!=lNot^2 && l>0 && Q>0, polar = "Pendular"; isRetro = False;,
		l^2!=lNot^2 && l<0 && Q>0, polar = "Pendular"; isRetro = True;,
		l^2<=lNot^2 && l>0 && Q==0, polar = "Equatorial"; isRetro = False;,
		l^2<=lNot^2 && l<0 && Q==0, polar = "Equatorial"; isRetro = True;,
		_, Print["Uncorrect input parameters: unable to determine the polar class !"]; Return[$Failed];	
	];
	
	(* Produce the output *)
	assoc = GetOrbit[st, parametrization, radial, polar, type, "OutStyle"-> "None", "SimplificationRule" -> "None", 
	"Retrograde" -> isRetro, "ReplaceLStar" -> True, "ReplaceC" ->True, "ReplaceLNot" -> True, "ReplaceCNot" -> True, "ReplaceRoots" -> True, 
	"ReplaceTurningPoints" -> True, "CosTheta" -> OptionValue["CosTheta"], "ReplacePhiTheta" -> True]//Quiet;
	
	(* Rule to replace initial data and parameters*)
	
	paramRule = {assoc["Kappa"] -> \[Kappa], assoc["Energy"] -> EE, assoc["Angular Momentum"] -> l, assoc["Carter Constant"] -> Q, \[Mu] -> mu, M -> mass, assoc["Initial Data"][[1]] -> initData[[1]], 
	assoc["Initial Data"][[2]] -> initData[[2]], assoc["Initial Data"][[3]] -> initData[[3]], assoc["Initial Data"][[4]] -> initData[[4]], assoc["Initial Data"][[5]] -> initData[[5]], 
	assoc["Initial Data"][[6]] -> initData[[6]]};
	
	If[OptionValue["Numerical"], assoc["Trajectory"] =  N[assoc["Trajectory"]];];
	
	If[assoc["Spacetime"]==$Failed["Spacetime"], Return["An error occured in the process !"]; Return[$Failed];];
	
	Return[NearHorizonGeoOrbitFunction[assoc/.paramRule]];
	
];


NearHorizonGeoOrbitFunction /:
 MakeBoxes[kgof:NearHorizonGeoOrbitFunction[assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended, Q, EE, l},

  summary = {Row[{BoxForm`SummaryItem[{"Spacetime: ", assoc["Spacetime"]}], "\t",
                  BoxForm`SummaryItem[{"Radial class: ", ClassName[assoc["Radial Class"]]}], "\t",
                  BoxForm`SummaryItem[{"Polar class: ", ClassName[assoc["Polar Class"]]}]}],
             Row[{BoxForm`SummaryItem[{"Parametrization: ", assoc["Parametrization"]}], "\t",
                  BoxForm`SummaryItem[{"Parametrization form: ", assoc["Parametrization Form"]}]}]};
  extended = {Row[{BoxForm`SummaryItem[{"Energy: ", assoc["Energy"]}], "\t",
              BoxForm`SummaryItem[{"Angular Momentum: ", assoc["Angular Momentum"]}], "\t",
              BoxForm`SummaryItem[{"Carter Constant: ", assoc["Carter Constant"]}]}],
              Row[{BoxForm`SummaryItem[{"\[ScriptCapitalC]: ", assoc["Casimir"]}], " \t", BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(*\)]\): ", assoc["LStar"]}]}],
              Row[{BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(\[SmallCircle]\)]\): ", assoc["CNought"]}], "\t", BoxForm`SummaryItem[{"\!\(\*SubscriptBox[\(\[ScriptCapitalL]\), \(\[SmallCircle]\)]\): ", assoc["LNought"]}]}]
              };
  BoxForm`ArrangeSummaryBox[
    NearHorizonGeoOrbitFunction,
    kgof,
    None,
    summary,
    extended,
    form]
];


NearHorizonGeoOrbitFunction[assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := assoc["Trajectory"][\[Lambda]]
NearHorizonGeoOrbitFunction[assoc_][y_?StringQ] := assoc[y]


ClassName[class_]:=Switch[class
	,"PendularNought", Return["\!\(\*SubscriptBox[\(Pendular\), \(\[SmallCircle]\)]\)"];
	,"SphericalStar", Return["\!\(\*SubscriptBox[\(Spherical\), \(*\)]\)(ISSO)"];
	,"PlungingStar", Return["\!\(\*SubscriptBox[\(Plunging\), \(*\)]\)"];
	,"OutwardStar", Return["\!\(\*SubscriptBox[\(Outward\), \(*\)]\)"];
	,"BoundedSubcritical", Return["\!\(\*SubscriptBox[\(Bounded\), \(<\)]\)"];
	,"BoundedStarMinus", Return["\!\(\*SubsuperscriptBox[\(Bounded\), \(*\), \(-\)]\)"];
	,"BoundedStar", Return["\!\(\*SubscriptBox[\(Bounded\), \(*\)]\)"];
	,"BoundedSupercritical", Return["\!\(\*SubscriptBox[\(Bounded\), \(>\)]\)"];
	,"PlungingStarNull", Return["\!\(\*SubscriptBox[\(Plunging\), \(*\)]\)(\[ScriptE]=0)"];
	,"OutwardStarNull", Return["\!\(\*SubscriptBox[\(Outward\), \(*\)]\)(\[ScriptE]=0)"];
	,_, Return[class];
];


(* ::Section:: *)
(*Close the package*)


End[];

EndPackage[];
