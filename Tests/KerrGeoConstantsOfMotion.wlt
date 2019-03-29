(****************************************************************)
(* KerrGeoConstantsOfMotion with random numerical input         *)
(****************************************************************)
VerificationTest[
    Values[KerrGeoConstantsOfMotion[0.9`20, 10, 0.1`20, Cos[Pi/3]]]
    ,
    {
        0.95426997475311240549349843416842480591,
        1.79647169810973255845128705122101046531, 
        9.73622324161324995959740006125691526459
    }
    ,
    TestID->"KerrGeoConstantsOfMotion with random numerical input"
]

(****************************************************************)
(* KerrGeoConstantsOfMotion documentation example 2             *)
(****************************************************************)
VerificationTest[
    Values[KerrGeoConstantsOfMotion[0, p, e, 1]]
    ,
    {
        Sqrt[(-4 e^2 + (-2 + p)^2)/(p (-3 - e^2 + p))],
        p/Sqrt[-3 - e^2 + p], 
        0
    }
    ,
    TestID->"KerrGeoConstantsOfMotion documentation example 2"
]



