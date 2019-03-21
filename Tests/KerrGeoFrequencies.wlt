(****************************************************************)
(* KerrGeoFrequencies formula for circular polar in Schw        *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0, r0, 0, 1]
    ,
    {Sqrt[-6 + r0]/r0^2, 1/r0^(3/2), r0/Sqrt[r0^5]}
    ,
    TestID -> "KerrGeoFrequencies for circular polar orbit in Schwarzschild formula"
]

(****************************************************************)
(* KerrGeoFrequencies numerics                                  *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0.9, 6.0, 0.3, 1]
    ,
    {0.037645648405229756`, 0.05262560724893245`, 0.05816122061091471`}
    ,
    TestID -> "KerrGeoFrequencies numerics"
]

