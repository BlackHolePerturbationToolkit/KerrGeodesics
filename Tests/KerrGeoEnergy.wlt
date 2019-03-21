(****************************************************************)
(* KerrGeoEnergy formula for polar in Schw                      *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0, r0, e, 1]
    ,
    Sqrt[(-4 e^2 + (-2 + p)^2)/(p (-3 - e^2 + p))]
    ,
    TestID -> "KerrGeoEnergy for polar orbit in Schwarzschild formula"
]

(****************************************************************)
(* KerrGeoEnergy formula for equatorial circular orbit          *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[a, p, 0, 0]
    ,
    Sqrt[(p (a^2 - 2 p + p^2)^2)/((a^2 + p^2) (a^2 + a^2 p - 3 p^2 + p^3))]
    ,
    TestID -> "KerrGeoEnergy for equatorial circular orbit formula"
]

(****************************************************************)
(* KerrGeoEnergy numerics                                       *)
(****************************************************************)
VerificationTest[
    KerrGeoEnergy[0.3771034125400172, 0.20945249623845186, 0.49967208040572664, 0.014887423705454594]
    ,
    2.327726198784594` - 0.17607223312831283` I
    ,
    TestID -> "KerrGeoEnergy with random numerical input"
]
