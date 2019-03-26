(****************************************************************)
(* KerrGeoAngularMomentum documentation example 1               *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[0.9`20, 10, 0.1`20, Cos[Pi/3]]
    ,
    1.79647169810973255845128705122101046531
    ,
    TestID->"KerrGeoAngularMomentum documentation example 1"
]

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
