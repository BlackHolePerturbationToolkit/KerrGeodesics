(****************************************************************)
(* SpinWeightedSpheroidalEigenvalue                             *)
(****************************************************************)
VerificationTest[
    KerrGeoFrequencies[0.9, 6.0, 0.3, 1]
    ,
    <|"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(r\)]\)" -> 0.037645648405229756, 
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Theta]\)]\)" -> 0.05262560724893245, 
     "\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Phi]\)]\)" -> 0.05816122061091471|>
    ,
    TestID->"KerrGeoFrequencies generic",
    SameTest -> withinRoundoff
]
