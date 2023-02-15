BeginTestSection["PlungesTest"]

VerificationTest[(* 1 *)
	Block[{t, r, \[Theta], \[Phi]}, orbit=KerrGeoPlunge[0.05, "ISSORadialParam", 5.9 ];
{t, r, \[Theta], \[Phi]}  =orbit["Trajectory"];
{t[\[Lambda]], r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]}/.\[Lambda]->5]
	,
	{166.93605140000858,  5.669341902321057, 2.4860551451282125, 17.222045750314397}	
	,
	TestID -> "KerrGeoPlunge - ISSORadialParam"
]

VerificationTest[(* 2 *)
	orbit=KerrGeoPlunge[0.05, "ISSORadialParam", 200.9 ]
	,
	$Failed
	,
	{KerrGeoPlunge::rIoutbounds}
	,
	TestID -> "KerrGeoPlunge -ISSORadialParam - ROutOfBoundsUpper"
]

VerificationTest[(* 3 *)
	orbit=KerrGeoPlunge[0.05, "ISSORadialParam", 0 ]
	,
	$Failed
	,
	{KerrGeoPlunge::rIoutbounds}
	,
	TestID -> "KerrGeoPlunge -ISSORadialParam - ROutOfBoundsUnder"
]

VerificationTest[(* 4 *)
	Block[{t, r, \[Theta], \[Phi]}, orbit=KerrGeoPlunge[0.08, "ISSORadialParam", 6.06, {63 , 1.1,2.4, -140}];
{t, r, \[Theta], \[Phi]}  =orbit["Trajectory"];
{t[\[Lambda]], r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]}/.\[Lambda]->8.9]
	,
	{456.9800421782238,  5.99235372581383, 2.771355583007127, -170.597372084777}	
	,
	TestID -> "KerrGeoPlunge - ISSORadialParamICs"
]

VerificationTest[(* 5 *)
	KerrGeoPlunge[0.08, "ISSORadialParam", 6.06, {63 , 78.1,2.4, -140}]
	,
	$Failed
	,
	{KerrGeoPlunge::r0outofbounds}
	,
	TestID -> "KerrGeoPlunge - ISSORadialParamICs - RadialOutOfBounds"
]

VerificationTest[(* 6 *)
	KerrGeoPlunge[0.08, "ISSORadialParam", 6.06, {63 , 2.1,12.4, -140}]
	,
	$Failed
	,
	{KerrGeoPlunge::\[Theta]0outofbounds}
	,
	TestID -> "KerrGeoPlunge - ISSORadialParamICs - PolarOutOfBounds"
]

VerificationTest[(* 7 *)
	Block[{t, r, \[Theta], \[Phi]}, orbit=KerrGeoPlunge[0.08, "ISSORadialParam", 6.06, {63` ,(11/10),(12/5), -140`}];
{t, r, \[Theta], \[Phi]}  =orbit["Trajectory"];
Round[{t[\[Lambda]], r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]}/.\[Lambda]->0, 10^(-8)]]
	,
	{63,  11/10, 12/5, -140}	
	,
	TestID -> "KerrGeoPlunge - ISSORadialParamICs - TimeZero"
]

VerificationTest[(* 8 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.99,  "ISSORadialParam", 8.717352279606489]; 
   {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 8.99]
	,
	{709.2968761641665,  8.643214446769269, 1.6120407986455219, -34.4819996184413}	
	,
	TestID -> "KerrGeoPlunge - ISSORadialParam MaximalRI"
]

VerificationTest[(* 9 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.7,  "ISSOIncParam", Pi/2.5]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> -7.8]
	,
	{-136.21547228952554,  3.389631843209988, 1.8758140096570342, -24.6802217009732}	
	,
	TestID -> "KerrGeoPlunge - ISSOIncParam"
]

VerificationTest[(* 10 *)
	KerrGeoPlunge[0.7,  "ISSOIncParam", 15]
	,
	$Failed
	,
	{KerrGeoPlunge::Inclinationoutofbounds}
	,
	TestID -> "KerrGeoPlunge - ISSOIncParam MaxIncOutOfBounds"
]

VerificationTest[(* 11 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.77,  "ISSOIncParam", -(Pi/4.2), {123,  4.3, Pi/2.2, 12.8}]; 
   {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> -2.09]
	,
	{95.35960150614794,  4.282702895995593, 2.391128629110148, 21.014943131587707}	
	,
	TestID -> "KerrGeoPlunge - ISSOIncParam MaxIncICs"
]

VerificationTest[(* 12 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.77,  "ISSOIncParam", -(Pi/2), {123,  4.3, Pi/2, 12.8}]; 
   {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> -2.09]
	,
	{87.79211403209042,  5.614689072605751, 1.5707963267948966, 21.38125127100541}	
	,
	TestID -> "KerrGeoPlunge - ISSOIncParamICs - ExtremeInc Negative"
]

VerificationTest[(* 13 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.977,  "ISSOIncParam", Pi/2, {123,  1.3, Pi/2, 12.8}]; 
   {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> -2.09]
	,
	{105.04498996779,  0.6649464562686522, 1.5707963267948966, 4.598376768609648}	
	,
	TestID -> "KerrGeoPlunge - ISSOIncParamICs - ExtremeInc Positive"
]

VerificationTest[(* 14 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.92,  {0.8,  10, 6}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 3.66]
	,
	{20.661343459257957,  0.15028185959212795, 1.5630908171743572, 38.368616042198404}	
	,
	TestID -> "KerrGeoPlunge - Generic Plunge - Complex Roots"
]

VerificationTest[(* 15 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.987,  {0.101,  0.39, 0}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 3.66]
	,
	{2.4976443942888236,  1.0611065184188126, 1.5707963267948966, -0.259629783302759}
	,
	{Power::infy, Power::infy, Infinity::indet, Power::infy, General::stop, Infinity::indet, Infinity::indet, General::stop}
	,
	TestID -> "KerrGeoPlunge - Generic Plunge - Real Roots 3 Inside rm"
]

VerificationTest[(* 16 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.1,  {0.956,  2.55, 6.54}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 3.66]
	,
	{48.654713642077695,  0.8239191438889363, 1.126915405459793, 13.12077316222189}	
	,
	TestID -> "KerrGeoPlunge - Generic Plunge - Real Roots 3 Outside rp"
]

VerificationTest[(* 17 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.92,  {0.8,  10, 0}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 3.66]
	,
	{18.729898545961788,  0.9316108659154707, 1.5707963267948966, 37.204816126657136}
	,
	{Power::infy, Power::infy, Infinity::indet}
	,
	TestID -> " KerrGeoPlunge - Generic - Zero Carter"
]

VerificationTest[(* 18 *)
	KerrGeoPlunge[0.92,  {10,  10, 2}]
	,
	$Failed
	,
	{KerrGeoPlunge::highenergy}
	,
	TestID -> "KerrGeoPlunge - Generic - High Energy"
]

VerificationTest[(* 19 *)
	KerrGeoPlunge[0.92,  {0.6,  2.4, -5}]
	,
	$Failed
	,
	{KerrGeoPlunge::negativecarter}
	,
	TestID -> "KerrGeoPlunge - Generic - Negative Carter"
]

VerificationTest[(* 20 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.9,  {0.3,  4, 10}, {-221,  1, Pi/2.9, 42}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 8.66]
	,
	{-186.53361138627804,  1.430951644567102, 0.9346228824493619, 91.0672074418539}	
	,
	TestID -> "KerrGeoPlunge - Generic - ICs"
]

VerificationTest[(* 21 *)
	Block[{t,  r, \[Theta], \[Phi]},  orbit = KerrGeoPlunge[0.987,  {0.101,  0.39, 0}]; {t,  r, \[Theta], \[Phi]} = orbit["Trajectory"]; 
   {t[\[Lambda]],  r[\[Lambda]], \[Theta][\[Lambda]], \[Phi][\[Lambda]]} /. \[Lambda] -> 3.66]
	,
	{2.4976443942888236,  1.0611065184188126, 1.5707963267948966, -0.259629783302759}
	,
	{Power::infy, Power::infy, Infinity::indet, Power::infy, General::stop, Infinity::indet, Infinity::indet, General::stop}
	,
	TestID -> "KerrGeoPlunge - RealRoot QZero"
]

VerificationTest[(* 22 *)
	KerrGeoPlunge[0.9,  {0.3,  4, 10}, {-221,  1101, Pi/2.9, 42}]
	,
	$Failed
	,
	{KerrGeoPlunge::r0outofboundsGen}
	,
	TestID -> "KerrGeoPlunge - Generic - IC radius out of bounds"
]

VerificationTest[(* 23 *)
	KerrGeoPlunge[0.9,  {0.3,  4, 10}, {-221,  1, Pi/0.9, 42}]
	,
	$Failed
	,
	{KerrGeoPlunge::\[Theta]0outofbounds}
	,
	TestID -> "KerrGeoPlunge - Generic - IC polar out of bounds"
]

EndTestSection[]
