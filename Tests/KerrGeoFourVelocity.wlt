(****************************************************************************************************)
(* KerrGeoFourVelocity correctness test for the contravariant four-velocity for generic Kerr orbit  *)
(* Reference values from Zach Nasipak                                                               *)
(****************************************************************************************************)
VerificationTest[
    1-Through[Values[KerrGeoFourVelocity[0.9,10,0.1,0.5]][5]]/{1.2079169755393362`,0.01865126187702151`,-0.03265833538560572`,0.03280723490549257`}//Abs//Chop//Total
    ,
    0
    ,
    TestID -> "KerrGeoFourVelocity correctness test for the contravariant four-velocity for generic Kerr orbit"
]


(****************************************************************************************************)
(* KerrGeoFourVelocity correctness test for the covariant four-velocity for generic Kerr orbit  *)
(* Reference values from Zach Nasipak                                                               *)
(****************************************************************************************************)
VerificationTest[
    1-Through[Values[KerrGeoFourVelocity[0.2,8,0.7,0.4,"Covariant"->True]][20]]/{-0.9710424662201889`,-0.17158372820149176`,-3.373274660052031`,1.482176616261215`}//Abs//Chop//Total
    ,
    0
    ,
    TestID -> "KerrGeoFourVelocity correctness test for the covariant four-velocity for generic Kerr orbit"
]