(****************************************************************)
(* KerrGeoEnergy documentation example 1                        *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0.9`20, 10, 0.1`20, Cos[Pi/3]]
    ,
    0.95426997475311240549349843416842480591
    ,
    TestID->"KerrGeoEnergy documentation example 1"
]

(****************************************************************)
(* KerrGeoEnergy documentation example 2                        *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0, p, e, 1]
    ,
    Sqrt[(-4 e^2 + (-2 + p)^2)/(p (-3 - e^2 + p))]
    ,
    TestID->"KerrGeoEnergy documentation example 2"
]

(****************************************************************)
(* KerrGeoEnergy documentation example 3                        *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[a, p, 0, 0]
    ,
    Sqrt[(p (a^2 - 2 p + p^2)^2)/((a^2 + p^2) (a^2 + a^2 p - 3 p^2 + p^3))]
    ,
    TestID->"KerrGeoEnergy documentation example 3"
]
