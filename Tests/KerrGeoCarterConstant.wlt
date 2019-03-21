(****************************************************************)
(* KerrGeoCarterConstant documentation example 1                *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[0.9`20, 10, 0.1`20, Cos[Pi/3]]
    ,
    9.73622324161324995959740006125691
    ,
    TestID->"KerrGeoCarterConstant documentation example 1"
]

(****************************************************************)
(* KerrGeoCarterConstant documentation example 2                *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[0, p, e, 1]
    ,
    0
    ,
    TestID->"KerrGeoCarterConstant documentation example 2"
]

(****************************************************************)
(* KerrGeoCarterConstant documentation example 3                *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[a, p, 0, 0]
    ,
    (p^2 (a^4 + 2 a^2 (-2 + p) p + p^4))/((a^2 + p^2) ((-3 + p) p^2 + 
   a^2 (1 + p)))
    ,
    TestID->"KerrGeoCarterConstant documentation example 3"
]
