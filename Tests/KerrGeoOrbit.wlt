(****************************************************************)
(* KerrGeoOrbit Schwarzschild/Circular/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0., p=10., e=0. x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Schwarzschild/Circular/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Schwarzschild/Circular/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0., p=10., e=0. x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Schwarzschild/Circular/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Schwarzschild/eccentric/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0., p=10., e=0.5 x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Schwarzschild/eccentric/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Schwarzschild/eccentric/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0., p=10., e=0.5 x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Schwarzschild/eccentric/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Circular/Prograde/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0. x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Circular/Prograde/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Circular/Prograde/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0. x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Circular/Prograde/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Circular/Retrograde/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0. x=-1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Circular/Retrograde/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Circular/Retrograde/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0. x=-1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Circular/Retrograde/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Eccentric/Prograde/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0.5 x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Eccentric/Prograde/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Eccentric/Prograde/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0.5 x=1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Eccentric/Prograde/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Eccentric/Retrograde/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0.5 x=-1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Eccentric/Retrograde/Mino Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Eccentric/Retrograde/Darwin Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0.5 x=-1},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Darwin"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Darwin"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Eccentric/Retrograde/Darwin Analytic vs. FastSpec"
]

(****************************************************************)
(* KerrGeoOrbit Kerr/Generic/Prograde/Mino Analytic vs. FastSpec *)
(****************************************************************)
VerificationTest[
    Block[{Internal`$EqualTolerance = 4, a=0.9, p=10., e=0.5 x=0.5},
		KerrGeoOrbit[a, p, e, x, Method -> "Analytic", "Parametrization" -> "Mino"][42.]
		==
		KerrGeoOrbit[a, p, e, x, Method -> "FastSpec", "Parametrization" -> "Mino"][42.]
 ]
    ,
    TestID->"KerrGeoOrbit Kerr/Generic/Prograde/Mino Analytic vs. FastSpec"
]

