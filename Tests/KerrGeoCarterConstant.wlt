(****************************************************************)
(* KerrGeoCarterConstant formula for polar in Schw             *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[0, p, e, 1]
    ,
    0
    ,
    TestID -> "KerrGeoCarterConstant for polar orbit in Schwarzschild formula"
]

(****************************************************************)
(* KerrGeoCarterConstant formula for equatorial circular orbit  *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[a, p, 0, 0]
    ,
    (p^2 (a^4 + 2 a^2 (-2 + p) p + p^4))/((a^2 + p^2) ((-3 + p) p^2 + a^2 (1 + p)))
    ,
    TestID -> "KerrGeoCarterConstant for equatorial circular orbit formula"
]

(****************************************************************)
(* KerrGeoCarterConstant numerics                               *)
(****************************************************************)
VerificationTest[
    KerrGeoCarterConstant[0.3771034125400172, 0.20945249623845186, 0.49967208040572664, 0.014887423705454594]
    ,
    -1.9911772428885637` + 0.2263351638583924` I
    ,
    TestID -> "KerrGeoCarterConstant with random numerical input"
]