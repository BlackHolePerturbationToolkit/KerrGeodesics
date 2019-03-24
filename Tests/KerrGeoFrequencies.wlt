(****************************************************************)
(* KerrGeoFrequencies generic                                   *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0.9, 6.0, 0.3, 1]
    ,
    <|"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> 0.037645648405229756, 
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> 0.05262560724893245, 
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> 0.05816122061091471|>
    ,
    TestID->"KerrGeoFrequencies generic"
]

(****************************************************************)
(* KerrGeoFrequencies documentation example 1                   *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0.9`20, 5, 0.7`20, Cos[\[Pi]/4]]
    ,
    {
        0.02147956794769575998110000119989314899,
        0.0401121096645215624490993245979412002,
        0.04729735868384957938674379795253462165
    }
    ,
    TestID->"KerrGeoFrequencies documentation example 1"
]

(****************************************************************)
(* KerrGeoFrequencies documentation example 2                   *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0.9`20, 5, 0.7`20, Cos[\[Pi]/4], Time -> "Mino"]
    ,
    {
        1.56234759143678127542701361075015211086,
        2.9176125923210750323363531819102210566,
        3.44024212223323910871580068419286366985
    }
    ,
    TestID->"KerrGeoFrequencies documentation example 2"
]

(****************************************************************)
(* KerrGeoFrequencies documentation example 3                   *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0, r0, 0, 1]
    ,
    {
        Sqrt[-6 + r0]/r0^2,
        1/r0^(3/2),
        r0/Sqrt[r0^5]
    }
    ,
    TestID->"KerrGeoFrequencies documentation example 3"
]
