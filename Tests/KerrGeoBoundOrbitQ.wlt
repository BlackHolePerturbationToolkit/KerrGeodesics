(****************************************************************)
(* KerrGeoBoundOrbitQ False 1                                   *)
(****************************************************************)
VerificationTest[
    KerrGeodesics`SpecialOrbits`Private`KerrGeoBoundOrbitQ[0.3771034125400172, 0.20945249623845186, 0.49967208040572664, 0.014887423705454594]
    ,
    False
    ,
    TestID -> "KerrGeoBoundOrbitQ False 1"
]

(****************************************************************)
(* KerrGeoBoundOrbitQ False 2                                   *)
(****************************************************************)
VerificationTest[
    KerrGeodesics`SpecialOrbits`Private`KerrGeoBoundOrbitQ[0.604741601450882`, 3.5495012516854552`, 0.7295796239129471`, 0.874953488795457`]
    ,
    False
    ,
    TestID -> "KerrGeoBoundOrbitQ False 2"
]

(****************************************************************)
(* KerrGeoBoundOrbitQ True 1                                    *)
(****************************************************************)
VerificationTest[
    KerrGeodesics`SpecialOrbits`Private`KerrGeoBoundOrbitQ[0.3771034125400172, 10.20945249623845186, 0.49967208040572664, 0.014887423705454594]
    ,
    True
    ,
    TestID -> "KerrGeoBoundOrbitQ True 1"
]

(****************************************************************)
(* KerrGeoBoundOrbitQ True 2                                    *)
(****************************************************************)
VerificationTest[
    KerrGeodesics`SpecialOrbits`Private`KerrGeoBoundOrbitQ[0.974953488795457, 5.3771034125400172, 0.01967208040572664, 0.49967208040572664]
    ,
    True
    ,
    TestID -> "KerrGeoBoundOrbitQ True 2"
]