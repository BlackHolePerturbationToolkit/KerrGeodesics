(****************************************************************)
(* KerrGeoPhotonSphereRadius documentation example 2            *)
(****************************************************************)
VerificationTest[
    KerrGeoPhotonSphereRadius[a, 0]
    ,
    1 + 2 Sqrt[1 - a^2/3] Cos[1/3 ArcCos[(1 - a^2)/(1 - a^2/3)^(3/2)]]
    ,
    TestID->"KerrGeoPhotonSphereRadius documentation example 2"
]
