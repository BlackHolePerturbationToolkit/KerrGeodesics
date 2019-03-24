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
(* KerrGeoAngularMomentum documentation example 2               *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[0, p, e, 1]
    ,
    p/Sqrt[-3 - e^2 + p]
    ,
    TestID->"KerrGeoAngularMomentum documentation example 2"
]

(****************************************************************)
(* KerrGeoAngularMomentum documentation example 3               *)
(****************************************************************)
VerificationTest[
    KerrGeoAngularMomentum[a, p, 0, 0]
    ,
    0
    ,
    TestID->"KerrGeoAngularMomentum documentation example 3"
]
