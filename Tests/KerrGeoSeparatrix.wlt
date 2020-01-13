(****************************************************************)
(* KerrGeoSeparatrix documentation example 1                          *)
(****************************************************************)
VerificationTest[
    KerrGeoSeparatrix[0.9`20, 0.5`20, Cos[Pi/3]]
    ,
    4.34225968111274824428436590954469803095`19.
    ,
    TestID->"KerrGeoSeparatrix documentation example 1"
]

(****************************************************************)
(* KerrGeoSeparatrix Numerics 1                                 *)
(****************************************************************)
VerificationTest[
    KerrGeoSeparatrix[0.9, 0.5, Cos[Pi/3]]
    ,
    4.342259681112724`
    ,
    TestID -> "KerrGeoSeparatrix  Numerics 1"
    ,
    SameTest -> withinRoundoff 
]
