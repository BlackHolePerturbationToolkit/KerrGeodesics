(****************************************************************)
(* KerrGeoPhotonSphereRadius polar formula                      *)
(****************************************************************)
VerificationTest[
    KerrGeoPhotonSphereRadius[a, 0]
    ,
    1 + 2 Sqrt[1 - a^2/3] Cos[1/3 ArcCos[(1 - a^2)/(1 - a^2/3)^(3/2)]]
    ,
    TestID -> "KerrGeoPhotonSphereRadius polar formula"
]

(****************************************************************)
(* KerrGeoPhotonSphereRadius numerics                           *)
(****************************************************************)
VerificationTest[
    KerrGeoPhotonSphereRadius[0.974953488795457, 0.49967208040572664]
    ,
    1.6613320588698588`
    ,
    TestID -> "KerrGeoPhotonSphereRadius Numerics"
    ,
    SameTest -> withinRoundoff
]