(****************************************************************)
(* KerrGeoEnergy documentation example 1                        *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0.9`20, 10, 0.1`20, Cos[Pi/3]]
    ,
    0.95426997475311240549349843416842480591
    ,
    TestID->"KerrGeoEnergy generic numerical test"
]

(****************************************************************)
(* KerrGeoEnergy documentation example 2                        *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0, p, e, 1]
    ,
    Sqrt[(-4 e^2 + (-2 + p)^2)/(p (-3 - e^2 + p))]
    ,
    TestID->"KerrGeoEnergy for Schwarzschild equatorial circular orbit formula"
]



