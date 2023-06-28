(* ::Package:: *)

(* ::Title:: *)
(*ParallelTransport subpackage of KerrGeodesics*)


(* ::Text:: *)
(*This subpackage implements the  analytic equations for parallel transport along a geodesic in Kerr spacetime given by M. van de Meent in https://arxiv.org/abs/1906.05090*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`ParallelTransport`",
	{
	 "KerrGeodesics`KerrGeoOrbit`"}];

KerrParallelTransportFrame::usage = "KerrParallelTransportFrame[a,p,e,x] returns a KerrParallelTransportFrameFunction[..] which stores the parallely transported frame, orbital trajectory, and parameters.";
KerrParallelTransportFrameFunction::usage = "KerrParallelTransportFrameFunction[a,p,e,x,assoc] an object for storing the parallely transported frame, trajectory and orbital parameters in the assoc Association.";

Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Error messages*)


KerrParallelTransportFrame::general = "`1`"
KerrParallelTransportFrame::parametrization = "Parametrization error: `1`"


(* ::Section::Closed:: *)
(*Kerr*)


(* ::Subsection::Closed:: *)
(*Generic (Mino)*)


(* ::Subsubsection::Closed:: *)
(*Frequencies*)


(* ::Input::Initialization:: *)
MinoPrecessionFrequency[orbit_KerrGeoOrbitFunction]:=MinoPrecessionFrequencyr[orbit]+MinoPrecessionFrequency\[Theta][orbit]


(* ::Input::Initialization:: *)
MinoPrecessionFrequencyr[orbit_KerrGeoOrbitFunction]:=
With[
{
a=orbit["a"],
\[ScriptCapitalE]=orbit["Energy"],
\[ScriptCapitalL]=orbit["AngularMomentum"],
\[ScriptCapitalQ]=orbit["CarterConstant"],
r1=orbit["RadialRoots"][[1]],
r2=orbit["RadialRoots"][[2]],
r3=orbit["RadialRoots"][[3]],
r4=orbit["RadialRoots"][[4]],
\[CapitalUpsilon]r=orbit["RadialFrequency"]
},
With[
{
\[ScriptCapitalK]=\[ScriptCapitalQ]+(a \[ScriptCapitalE]-\[ScriptCapitalL])^2,
kr= (r1-r2)/(r1-r3) (r3-r4)/(r2-r4)
},
(\[CapitalUpsilon]r Sqrt[\[ScriptCapitalK]])/\[Pi] ((2 (a^2 \[ScriptCapitalE]+r3^2 \[ScriptCapitalE]-a \[ScriptCapitalL]) EllipticK[kr])/(Sqrt[-(r1-r3) (r2-r4) (-1+\[ScriptCapitalE]^2)] (r3^2+\[ScriptCapitalK]))+( (r2-r3) (-a^2 \[ScriptCapitalE]+\[ScriptCapitalE] \[ScriptCapitalK]+a \[ScriptCapitalL]) (2Im[(r2+I Sqrt[\[ScriptCapitalK]]) (r3+I Sqrt[\[ScriptCapitalK]]) EllipticPi[((r1-r2) (r3-I Sqrt[\[ScriptCapitalK]]))/((r1-r3) (r2-I Sqrt[\[ScriptCapitalK]])),kr]]))/(Sqrt[(r1-r3) (r2-r4) (1-\[ScriptCapitalE]^2)] Sqrt[\[ScriptCapitalK]] (r2^2+\[ScriptCapitalK]) (r3^2+\[ScriptCapitalK])))
]
];


(* ::Input::Initialization:: *)
MinoPrecessionFrequency\[Theta][orbit_KerrGeoOrbitFunction]:=
With[
{
a=orbit["a"],
x=orbit["Inclination"],
\[ScriptCapitalE]=orbit["Energy"],
\[ScriptCapitalL]=orbit["AngularMomentum"],
\[ScriptCapitalQ]=orbit["CarterConstant"],
\[CapitalUpsilon]z=orbit["PolarFrequency"]
},
With[
{
zm=Sqrt[1-x^2],
zp=(a^2 (1-\[ScriptCapitalE]^2)+\[ScriptCapitalL]^2/x^2)^(1/2)
},
	With[
	{
	\[ScriptCapitalK]=\[ScriptCapitalQ]+(a \[ScriptCapitalE]-\[ScriptCapitalL])^2,
	k\[Theta]= a^2 (1-\[ScriptCapitalE]^2)(zm/zp)^2
	},
	-((2\[CapitalUpsilon]z Sqrt[\[ScriptCapitalK]])/(\[Pi] zp))(\[ScriptCapitalE] EllipticK[k\[Theta]]+((a^2 \[ScriptCapitalE]-\[ScriptCapitalE] \[ScriptCapitalK]-a \[ScriptCapitalL])/\[ScriptCapitalK]) EllipticPi[(a^2 zm^2)/\[ScriptCapitalK],k\[Theta]] )
	]
]
];


(* ::Subsubsection::Closed:: *)
(*PrecessionPhase*)


PrecessionPhaser[orbit_KerrGeoOrbitFunction]:=
	Function[{qr},Evaluate@With[
		{
			a=orbit["a"],
			\[ScriptCapitalE]=orbit["Energy"],
			\[ScriptCapitalL]=orbit["AngularMomentum"],
			\[ScriptCapitalQ]=orbit["CarterConstant"],
			r1=orbit["RadialRoots"][[1]],
			r2=orbit["RadialRoots"][[2]],
			r3=orbit["RadialRoots"][[3]],
			r4=orbit["RadialRoots"][[4]],
			\[CapitalUpsilon]r=orbit["RadialFrequency"],
			\[CapitalUpsilon]\[Psi]r=MinoPrecessionFrequencyr[orbit]
		},
			With[
				{
					\[ScriptCapitalK]=\[ScriptCapitalQ]+(a \[ScriptCapitalE]-\[ScriptCapitalL])^2,
					kr= (r1-r2)/(r1-r3) (r3-r4)/(r2-r4)
				},
					With[
						{
							\[Gamma]=JacobiAmplitude[EllipticK[kr] qr/\[Pi],kr]
						},
							Sqrt[\[ScriptCapitalK]](( (r2-r3) (-a^2 \[ScriptCapitalE]+\[ScriptCapitalE] \[ScriptCapitalK]+a \[ScriptCapitalL]) (2Im[(r2+I Sqrt[\[ScriptCapitalK]]) (r3+I Sqrt[\[ScriptCapitalK]])( EllipticPi[((r1-r2) (r3-I Sqrt[\[ScriptCapitalK]]))/((r1-r3) (r2-I Sqrt[\[ScriptCapitalK]])),\[Gamma],kr]-EllipticPi[((r1-r2) (r3-I Sqrt[\[ScriptCapitalK]]))/((r1-r3) (r2-I Sqrt[\[ScriptCapitalK]])),kr] qr/\[Pi])]))/(Sqrt[(r1-r3) (r2-r4) (1-\[ScriptCapitalE]^2)] Sqrt[\[ScriptCapitalK]] (r2^2+\[ScriptCapitalK]) (r3^2+\[ScriptCapitalK])))
					]
			]
	]
]


PrecessionPhasez[orbit_KerrGeoOrbitFunction]:=
Function[{qz},Evaluate@With[
	{
		a=orbit["a"],
		x=orbit["Inclination"],
		\[ScriptCapitalE]=orbit["Energy"],
		\[ScriptCapitalL]=orbit["AngularMomentum"],
		\[ScriptCapitalQ]=orbit["CarterConstant"],
		\[CapitalUpsilon]z=orbit["PolarFrequency"]
	},
		With[
			{
				zm=Sqrt[1-x^2],
				zp=(a^2 (1-\[ScriptCapitalE]^2)+\[ScriptCapitalL]^2/x^2)^(1/2)
			},
				With[
					{
						\[ScriptCapitalK]=\[ScriptCapitalQ]+(a \[ScriptCapitalE]-\[ScriptCapitalL])^2,
						k\[Theta]= a^2 (1-\[ScriptCapitalE]^2)(zm/zp)^2
					},
					With[
						{
							\[Gamma]=JacobiAmplitude[EllipticK[k\[Theta]] 2/\[Pi] (qz+\[Pi]/2),k\[Theta]]
						},
						-(Sqrt[\[ScriptCapitalK]]/zp)(((a^2 \[ScriptCapitalE]-\[ScriptCapitalE] \[ScriptCapitalK]-a \[ScriptCapitalL])/\[ScriptCapitalK]) (EllipticPi[(a^2 zm^2)/\[ScriptCapitalK],\[Gamma],k\[Theta]] -EllipticPi[(a^2 zm^2)/\[ScriptCapitalK],k\[Theta]] 2/\[Pi] (qz+\[Pi]/2)))
					]
				]
			]
		]
];


PrecessionPhase[orbit_KerrGeoOrbitFunction,\[CapitalUpsilon]\[Psi]_,qr0_,qz0_,q\[Psi]0_:0]:=
Function[{\[Lambda]},
	Evaluate@With[
		{
			\[CapitalUpsilon]r=orbit["RadialFrequency"],
			\[CapitalUpsilon]z=orbit["PolarFrequency"]
		},
		q\[Psi]0+\[CapitalUpsilon]\[Psi] \[Lambda] +PrecessionPhaser[orbit][\[CapitalUpsilon]r \[Lambda] + qr0]+PrecessionPhasez[orbit][\[CapitalUpsilon]z \[Lambda] + qz0]
]
]





(* ::Subsubsection::Closed:: *)
(*Tetrads*)


MarckCarterFrame[orbit_KerrGeoOrbitFunction]:=
	Module[
		{
			a,\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ],\[ScriptCapitalK],rf,\[Theta]f
		},
				a=orbit["a"];
				{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]}=Values@orbit["ConstantsOfMotion"];
				\[ScriptCapitalK]=\[ScriptCapitalQ]+(a \[ScriptCapitalE]-\[ScriptCapitalL])^2;
				rf=Function[\[Lambda],orbit["Trajectory"][[2]][\[Lambda]]];
				\[Theta]f=Function[\[Lambda],orbit["Trajectory"][[3]][\[Lambda]]];
				Function[{\[Lambda]},
					Evaluate@With[{
						r=rf[\[Lambda]],\[Theta]=\[Theta]f[\[Lambda]],rd=rf'[\[Lambda]],\[Theta]d=\[Theta]f'[\[Lambda]]
					},
					{(*(Subscript[e, i])^(a)Subscript[\[Eta], (a)(b)]Subscript[(\[Omega]^(b)), \[Mu]]*)
						{-\[ScriptCapitalE], rd/\[CapitalDelta][r],\[Theta]d, \[ScriptCapitalL]},
						{(-r rd \[Alpha]-a^2 \[Beta] \[Theta]d Cos[\[Theta]] Sin[\[Theta]])/(Sqrt[\[ScriptCapitalK]] \[CapitalSigma][r,\[Theta]]),(r ((a^2+r^2) \[ScriptCapitalE]-a \[ScriptCapitalL]) \[Alpha])/(Sqrt[\[ScriptCapitalK]] \[CapitalDelta][r]),(a \[Beta] Cos[\[Theta]] (-\[ScriptCapitalL] Csc[\[Theta]]+a \[ScriptCapitalE] Sin[\[Theta]]))/Sqrt[\[ScriptCapitalK]],(a Sin[\[Theta]] ((a^2+r^2) \[Beta] \[Theta]d Cos[\[Theta]]+r rd \[Alpha] Sin[\[Theta]]))/(Sqrt[\[ScriptCapitalK]] \[CapitalSigma][r,\[Theta]])},
						{(a^2 \[ScriptCapitalE] \[Alpha]+r^2 \[ScriptCapitalE] \[Alpha]+a \[ScriptCapitalL] (-\[Alpha]+\[Beta])-a^2 \[ScriptCapitalE] \[Beta] Sin[\[Theta]]^2)/\[CapitalSigma][r,\[Theta]],-((rd \[Alpha])/\[CapitalDelta][r]),-\[Beta] \[Theta]d,(-(a^2+r^2) \[ScriptCapitalL] \[Beta]+a (a \[ScriptCapitalL] \[Alpha]+a^2 \[ScriptCapitalE] (-\[Alpha]+\[Beta])+r^2 \[ScriptCapitalE] (-\[Alpha]+\[Beta])) Sin[\[Theta]]^2)/\[CapitalSigma][r,\[Theta]]},
						{(a (-rd Cos[\[Theta]]+r \[Theta]d Sin[\[Theta]]))/(Sqrt[\[ScriptCapitalK]] \[CapitalSigma][r,\[Theta]]),(a ((a^2+r^2) \[ScriptCapitalE]-a \[ScriptCapitalL]) Cos[\[Theta]])/(Sqrt[\[ScriptCapitalK]] \[CapitalDelta][r]),(r \[ScriptCapitalL] Csc[\[Theta]]-a r \[ScriptCapitalE] Sin[\[Theta]])/Sqrt[\[ScriptCapitalK]],(Sin[\[Theta]] (-2 r (a^2+r^2) \[Theta]d+a^2 rd Sin[2 \[Theta]]))/(2 Sqrt[\[ScriptCapitalK]] \[CapitalSigma][r,\[Theta]])}
					}/.{
						\[CapitalSigma][r,\[Theta]]-> r^2+a^2 Cos[\[Theta]]^2,
						\[CapitalDelta][r]-> r^2-2r+a^2,
						\[Alpha]-> Sqrt[(\[ScriptCapitalK]-a^2 Cos[\[Theta]]^2)/(r^2+\[ScriptCapitalK])],
						\[Beta]-> 1/Sqrt[((\[ScriptCapitalK]-a^2 Cos[\[Theta]]^2)/(r^2+\[ScriptCapitalK]))]
					}
					]
				]
	]


KerrParallelTransportFrameMino[a_,p_,e_,x_,initPhases:{_,_,_,_,_}:{0,0,0,0,0}]:=Module[
	{orbit,mcf,\[CapitalUpsilon]\[Psi],pt\[Psi],ptf,assoc},
		orbit=KerrGeoOrbit[a,p,e,x,initPhases[[1;;4]],"Method"->"Analytic","Parametrization"-> "Mino"];
		mcf=MarckCarterFrame[orbit];
		\[CapitalUpsilon]\[Psi] = MinoPrecessionFrequency[orbit];
		pt\[Psi]=PrecessionPhase[orbit,\[CapitalUpsilon]\[Psi],initPhases[[2]],initPhases[[3]],initPhases[[5]]];
		
		ptf=Function[{\[Lambda]},
			Evaluate[
				{
					mcf[\[Lambda]][[1]],
					mcf[\[Lambda]][[2]]Cos[pt\[Psi][\[Lambda]]]+mcf[\[Lambda]][[3]]Sin[pt\[Psi][\[Lambda]]],
					-mcf[\[Lambda]][[2]]Sin[pt\[Psi][\[Lambda]]]+mcf[\[Lambda]][[3]]Cos[pt\[Psi][\[Lambda]]],
					mcf[\[Lambda]][[4]]
				}
			]
		];
		
		assoc = Last@orbit;
		assoc["PrecessionFrequency"]=\[CapitalUpsilon]\[Psi];
		assoc["ParallelTransportedFrame"]=ptf;

		KerrParallelTransportFrameFunction[a,p,e,x,assoc]

]


(* ::Section::Closed:: *)
(*KerrParallelTransportFrame and KerrParallelTransportFrameFuction*)


Options[KerrParallelTransportFrame] = {"Parametrization" -> "Mino", "Method" -> "Analytic"}
SyntaxInformation[KerrParallelTransportFrame] = {"ArgumentsPattern"->{_,_,OptionsPattern[]}};


KerrParallelTransportFrame[a_,p_,e_,x_, initPhases:{_,_,_,_,_}:{0,0,0,0,0},OptionsPattern[]]:=Module[{param, method},
(*FIXME: add stability check but make it possible to turn it off*)

method = OptionValue["Method"];
param = OptionValue["Parametrization"];

If[param =!= "Mino", Message[KerrParallelTransportFrame::parametrization, "Only Mino time parametrization has been implemented for parallel transport."]; Return[];];
If[method =!= "Analytic", Message[KerrParallelTransportFrame::general, "Only Analytic method has been implemented for parallel transport."]; Return[];];

If[method == "Analytic",
	If[param == "Mino", Return[KerrParallelTransportFrameMino[a, p, e, x, initPhases]]];
];

Message[KerrParallelTransportFrame::general, "Method: " <> method <> "is not one of {Analytic}."];

]


KerrParallelTransportFrameFunction /:
 MakeBoxes[kptff:KerrParallelTransportFrameFunction[a_, p_, e_, x_, assoc_], form:(StandardForm|TraditionalForm)] :=
 Module[{summary, extended},
  summary = {Row[{BoxForm`SummaryItem[{"a: ", a}], "  ",
                  BoxForm`SummaryItem[{"p: ", p}], "  ",
                  BoxForm`SummaryItem[{"e: ", e}], "  ",
                  BoxForm`SummaryItem[{"x: ", x}]}]
             };
  extended = {BoxForm`SummaryItem[{"Precession Frequency: ", assoc["PrecessionFrequency"]}]};
  BoxForm`ArrangeSummaryBox[
    KerrParallelTransportFrameFunction,
    kptff,
    None,
    summary,
    extended,
    form]
];


KerrParallelTransportFrameFunction[a_, p_, e_, x_, assoc_][\[Lambda]_/;StringQ[\[Lambda]] == False] := assoc["ParallelTransportedFrame"][\[Lambda]]
KerrParallelTransportFrameFunction[a_, p_, e_, x_, assoc_][y_?StringQ] := assoc[y]


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
