(****************************************************************)
(* KerrGeoAngularMomentum formula for polar in Schw             *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[0, p, e, 1]
    ,
    p/Sqrt[-3 - e^2 + p]
    ,
    TestID -> "KerrGeoEnergy for polar orbit in Schwarzschild formula"
]

(****************************************************************)
(* KerrGeoAngularMomentum formula for equatorial circular orbit *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[a, p, 0, 0]
    ,
    0
    ,
    TestID -> "KerrGeoEnergy for equatorial circular orbit formula"
]

(****************************************************************)
(* KerrGeoAngularMomentum numerics                              *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[0.3771034125400172, 0.20945249623845186, 0.49967208040572664, 0.014887423705454594]
    ,
    -0.0006984251378633733` - 0.01742473413668152` I
    ,
    TestID -> "KerrGeoAngularMomentum with random numerical input"
]
