(* ::Package:: *)

(* ::Title:: *)
(*InitialConditions subpackage of KerrGeodesics*)


(* ::Chapter:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`InitialConditions`",
	{"KerrGeodesics`ConstantsOfMotion`",
	 "KerrGeodesics`OrbitalFrequencies`","KerrGeodesics`KerrGeoOrbit`"}];

KerrGeoInitOrbit::usage = "KerrGeoInitOrbit[a,{t0,r0,th0,ph0},u] returns a KerrGeoOrbitFunction[..] (which stores the orbital trajectory and parameters) corresponding to the initial positions and fourvelocity (which does not need to be normalized). 
One can also optionally provide {\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]} and {p,e,x} in other arguments.";

Begin["`Private`"];


(* ::Section:: *)
(*Kerr generic*)


(* ::Subsection::Closed:: *)
(*Auxiliary functions*)


gK[a_][t_,r_,\[Theta]_,\[Phi]_]:={{-1+(2 r)/(r^2+a^2 Cos[\[Theta]]^2),0,0,-((2 a r Sin[\[Theta]]^2)/(r^2+a^2 Cos[\[Theta]]^2))},{0,(r^2+a^2 Cos[\[Theta]]^2)/(a^2+(-2+r) r),0,0},{0,0,r^2+a^2 Cos[\[Theta]]^2,0},{-((2 a r Sin[\[Theta]]^2)/(r^2+a^2 Cos[\[Theta]]^2)),0,0,(Sin[\[Theta]]^2 ((a^2+r^2) (r^2+a^2 Cos[\[Theta]]^2)+2 a^2 r Sin[\[Theta]]^2))/(r^2+a^2 Cos[\[Theta]]^2)}};
\[ScriptCapitalK][a_][t_,r_,\[Theta]_,\[Phi]_]:={{(a^2 r Cos[\[Theta]]^2 (a^2+2 (-1+r) r+(a^2-2 r) Cos[2 \[Theta]]))/(r^2+a^2 Cos[\[Theta]]^2)^2,0,0,-((a r (a^4+2 r^4+a^2 r (-2+3 r)+a^2 (a^2+(-2+r) r) Cos[2 \[Theta]]) Sin[2 \[Theta]]^2)/(4 (r^2+a^2 Cos[\[Theta]]^2)^2))},{0,-((a^2 Cos[\[Theta]]^2 (r^2+a^2 Cos[\[Theta]]^2))/(a^2+(-2+r) r)),0,0},{0,0,r^4+a^2 r^2 Cos[\[Theta]]^2,0},{-((a r (a^4+2 r^4+a^2 r (-2+3 r)+a^2 (a^2+(-2+r) r) Cos[2 \[Theta]]) Sin[2 \[Theta]]^2)/(4 (r^2+a^2 Cos[\[Theta]]^2)^2)),0,0,(r Cos[\[Theta]]^2 (2 a^6-4 a^4 r+3 a^6 r+12 a^4 r^2+11 a^4 r^3+16 a^2 r^4+16 a^2 r^5+8 r^7+4 a^2 r (a^4+2 (-2+r) r^3+a^2 r (-2+3 r)) Cos[2 \[Theta]]+a^4 (-2+r) (a^2+(-2+r) r) Cos[4 \[Theta]]) Sin[\[Theta]]^2)/(8 (r^2+a^2 Cos[\[Theta]]^2)^2)}};
RRoots[a_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalQ]_]:=Evaluate@Solve[((r^2+a^2)\[ScriptCapitalE] - a \[ScriptCapitalL])^2 - (r^2-2 r +a^2)(r^2+(\[ScriptCapitalL] - a \[ScriptCapitalE])^2+\[ScriptCapitalQ])==0,r];
RDisc[a_,\[ScriptCapitalE]_,\[ScriptCapitalL]_,\[ScriptCapitalQ]_]:=-16 (-a^8 \[ScriptCapitalE]^4+a^10 \[ScriptCapitalE]^4+16 a^6 \[ScriptCapitalE]^6-16 a^8 \[ScriptCapitalE]^6-4 a^10 \[ScriptCapitalE]^6+62 a^8 \[ScriptCapitalE]^8+6 a^10 \[ScriptCapitalE]^8-72 a^8 \[ScriptCapitalE]^10-4 a^10 \[ScriptCapitalE]^10+27 a^8 \[ScriptCapitalE]^12+a^10 \[ScriptCapitalE]^12+4 a^7 \[ScriptCapitalE]^3 \[ScriptCapitalL]-4 a^9 \[ScriptCapitalE]^3 \[ScriptCapitalL]-96 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL]+100 a^7 \[ScriptCapitalE]^5 \[ScriptCapitalL]+16 a^9 \[ScriptCapitalE]^5 \[ScriptCapitalL]-428 a^7 \[ScriptCapitalE]^7 \[ScriptCapitalL]-24 a^9 \[ScriptCapitalE]^7 \[ScriptCapitalL]+540 a^7 \[ScriptCapitalE]^9 \[ScriptCapitalL]+16 a^9 \[ScriptCapitalE]^9 \[ScriptCapitalL]-216 a^7 \[ScriptCapitalE]^11 \[ScriptCapitalL]-4 a^9 \[ScriptCapitalE]^11 \[ScriptCapitalL]-6 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2+6 a^8 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2+240 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2-260 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2-21 a^8 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2+1274 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2+27 a^8 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2-1764 a^6 \[ScriptCapitalE]^8 \[ScriptCapitalL]^2-15 a^8 \[ScriptCapitalE]^8 \[ScriptCapitalL]^2+756 a^6 \[ScriptCapitalE]^10 \[ScriptCapitalL]^2+3 a^8 \[ScriptCapitalE]^10 \[ScriptCapitalL]^2+4 a^5 \[ScriptCapitalE] \[ScriptCapitalL]^3-4 a^7 \[ScriptCapitalE] \[ScriptCapitalL]^3-320 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3+360 a^5 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3+4 a^7 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3-2128 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL]^3+12 a^7 \[ScriptCapitalE]^5 \[ScriptCapitalL]^3+3276 a^5 \[ScriptCapitalE]^7 \[ScriptCapitalL]^3-20 a^7 \[ScriptCapitalE]^7 \[ScriptCapitalL]^3-1512 a^5 \[ScriptCapitalE]^9 \[ScriptCapitalL]^3+8 a^7 \[ScriptCapitalE]^9 \[ScriptCapitalL]^3-a^4 \[ScriptCapitalL]^4+a^6 \[ScriptCapitalL]^4+240 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4-280 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4+14 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4+2170 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4-45 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4-3780 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalL]^4+44 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalL]^4+1890 a^4 \[ScriptCapitalE]^8 \[ScriptCapitalL]^4-14 a^6 \[ScriptCapitalE]^8 \[ScriptCapitalL]^4-96 a \[ScriptCapitalE] \[ScriptCapitalL]^5+116 a^3 \[ScriptCapitalE] \[ScriptCapitalL]^5-12 a^5 \[ScriptCapitalE] \[ScriptCapitalL]^5-1372 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^5+24 a^5 \[ScriptCapitalE]^3 \[ScriptCapitalL]^5+2772 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL]^5-12 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL]^5-1512 a^3 \[ScriptCapitalE]^7 \[ScriptCapitalL]^5+16 \[ScriptCapitalL]^6-20 a^2 \[ScriptCapitalL]^6+3 a^4 \[ScriptCapitalL]^6+518 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6+9 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6-1260 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^6-26 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^6+756 a^2 \[ScriptCapitalE]^6 \[ScriptCapitalL]^6+14 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalL]^6-104 a \[ScriptCapitalE] \[ScriptCapitalL]^7-12 a^3 \[ScriptCapitalE] \[ScriptCapitalL]^7+324 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^7+20 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^7-216 a \[ScriptCapitalE]^5 \[ScriptCapitalL]^7-8 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL]^7+8 \[ScriptCapitalL]^8+3 a^2 \[ScriptCapitalL]^8-36 \[ScriptCapitalE]^2 \[ScriptCapitalL]^8+27 \[ScriptCapitalE]^4 \[ScriptCapitalL]^8-3 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^8-4 a \[ScriptCapitalE] \[ScriptCapitalL]^9+4 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^9+\[ScriptCapitalL]^10-\[ScriptCapitalE]^2 \[ScriptCapitalL]^10+a^8 \[ScriptCapitalQ]-a^10 \[ScriptCapitalQ]-20 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalQ]+19 a^8 \[ScriptCapitalE]^2 \[ScriptCapitalQ]+5 a^10 \[ScriptCapitalE]^2 \[ScriptCapitalQ]+48 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalQ]-28 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalQ]-98 a^8 \[ScriptCapitalE]^4 \[ScriptCapitalQ]-10 a^10 \[ScriptCapitalE]^4 \[ScriptCapitalQ]+192 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalQ]+170 a^8 \[ScriptCapitalE]^6 \[ScriptCapitalQ]+10 a^10 \[ScriptCapitalE]^6 \[ScriptCapitalQ]-252 a^6 \[ScriptCapitalE]^8 \[ScriptCapitalQ]-127 a^8 \[ScriptCapitalE]^8 \[ScriptCapitalQ]-5 a^10 \[ScriptCapitalE]^8 \[ScriptCapitalQ]+108 a^6 \[ScriptCapitalE]^10 \[ScriptCapitalQ]+35 a^8 \[ScriptCapitalE]^10 \[ScriptCapitalQ]+a^10 \[ScriptCapitalE]^10 \[ScriptCapitalQ]+40 a^5 \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]-44 a^7 \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]-192 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]+156 a^5 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]+268 a^7 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]-952 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]-540 a^7 \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]+1404 a^5 \[ScriptCapitalE]^7 \[ScriptCapitalL] \[ScriptCapitalQ]+452 a^7 \[ScriptCapitalE]^7 \[ScriptCapitalL] \[ScriptCapitalQ]-648 a^5 \[ScriptCapitalE]^9 \[ScriptCapitalL] \[ScriptCapitalQ]-136 a^7 \[ScriptCapitalE]^9 \[ScriptCapitalL] \[ScriptCapitalQ]-20 a^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+25 a^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-4 a^8 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+288 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-300 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-226 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+16 a^8 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+1920 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+541 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-24 a^8 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-3240 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-504 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+16 a^8 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+1620 a^4 \[ScriptCapitalE]^8 \[ScriptCapitalL]^2 \[ScriptCapitalQ]+164 a^6 \[ScriptCapitalE]^8 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-4 a^8 \[ScriptCapitalE]^8 \[ScriptCapitalL]^2 \[ScriptCapitalQ]-192 a \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]+244 a^3 \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]+40 a^5 \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]-2000 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3 \[ScriptCapitalQ]-80 a^5 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3 \[ScriptCapitalQ]+3960 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL]^3 \[ScriptCapitalQ]+40 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL]^3 \[ScriptCapitalQ]-2160 a^3 \[ScriptCapitalE]^7 \[ScriptCapitalL]^3 \[ScriptCapitalQ]+48 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-72 a^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+16 a^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-6 a^6 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+1120 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-156 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+18 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-2700 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+290 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-18 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+1620 a^2 \[ScriptCapitalE]^6 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-150 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalL]^4 \[ScriptCapitalQ]+6 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalL]^4 \[ScriptCapitalQ]-312 a \[ScriptCapitalE] \[ScriptCapitalL]^5 \[ScriptCapitalQ]+68 a^3 \[ScriptCapitalE] \[ScriptCapitalL]^5 \[ScriptCapitalQ]+972 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^5 \[ScriptCapitalQ]-188 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^5 \[ScriptCapitalQ]-648 a \[ScriptCapitalE]^5 \[ScriptCapitalL]^5 \[ScriptCapitalQ]+120 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL]^5 \[ScriptCapitalQ]+32 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-3 a^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-4 a^4 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-144 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]+48 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]+8 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]+108 \[ScriptCapitalE]^4 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-44 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-4 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^6 \[ScriptCapitalQ]-16 a \[ScriptCapitalE] \[ScriptCapitalL]^7 \[ScriptCapitalQ]+16 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^7 \[ScriptCapitalQ]+5 \[ScriptCapitalL]^8 \[ScriptCapitalQ]-a^2 \[ScriptCapitalL]^8 \[ScriptCapitalQ]-5 \[ScriptCapitalE]^2 \[ScriptCapitalL]^8 \[ScriptCapitalQ]+a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^8 \[ScriptCapitalQ]+8 a^4 \[ScriptCapitalQ]^2-12 a^6 \[ScriptCapitalQ]^2+4 a^8 \[ScriptCapitalQ]^2+48 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^2-44 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^2+24 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^2-16 a^8 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^2+206 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^2+22 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^2+24 a^8 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^2-324 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^2-68 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^2-16 a^8 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^2+162 a^4 \[ScriptCapitalE]^8 \[ScriptCapitalQ]^2+34 a^6 \[ScriptCapitalE]^8 \[ScriptCapitalQ]^2+4 a^8 \[ScriptCapitalE]^8 \[ScriptCapitalQ]^2-96 a \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^2+128 a^3 \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^2-44 a^5 \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^2-628 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]^2+88 a^5 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]^2+1188 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]^2-44 a^5 \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]^2-648 a^3 \[ScriptCapitalE]^7 \[ScriptCapitalL] \[ScriptCapitalQ]^2+48 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-84 a^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+35 a^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+4 a^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+686 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-255 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-12 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-1620 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+418 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+12 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2+972 a^2 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-198 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-4 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^2-312 a \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2+172 a^3 \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2+972 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2-436 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2-648 a \[ScriptCapitalE]^5 \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2+264 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL]^3 \[ScriptCapitalQ]^2+48 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-27 a^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-4 a^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-216 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2+144 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2+8 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2+162 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-114 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-4 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^2-24 a \[ScriptCapitalE] \[ScriptCapitalL]^5 \[ScriptCapitalQ]^2+24 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^5 \[ScriptCapitalQ]^2+10 \[ScriptCapitalL]^6 \[ScriptCapitalQ]^2-4 a^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]^2-10 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]^2+4 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^6 \[ScriptCapitalQ]^2+16 \[ScriptCapitalQ]^3-32 a^2 \[ScriptCapitalQ]^3+22 a^4 \[ScriptCapitalQ]^3-6 a^6 \[ScriptCapitalQ]^3+84 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^3-90 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^3+18 a^6 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^3-180 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^3+102 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^3-18 a^6 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^3+108 a^2 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^3-34 a^4 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^3+6 a^6 \[ScriptCapitalE]^6 \[ScriptCapitalQ]^3-104 a \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^3+92 a^3 \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^3+324 a \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]^3-228 a^3 \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]^3-216 a \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]^3+136 a^3 \[ScriptCapitalE]^5 \[ScriptCapitalL] \[ScriptCapitalQ]^3+32 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3-33 a^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3+4 a^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3-144 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3+144 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3-8 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3+108 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3-108 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3+4 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^3-16 a \[ScriptCapitalE] \[ScriptCapitalL]^3 \[ScriptCapitalQ]^3+16 a \[ScriptCapitalE]^3 \[ScriptCapitalL]^3 \[ScriptCapitalQ]^3+10 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^3-6 a^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^3-10 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^3+6 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^4 \[ScriptCapitalQ]^3+8 \[ScriptCapitalQ]^4-12 a^2 \[ScriptCapitalQ]^4+4 a^4 \[ScriptCapitalQ]^4-36 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^4+48 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^4-8 a^4 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^4+27 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^4-35 a^2 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^4+4 a^4 \[ScriptCapitalE]^4 \[ScriptCapitalQ]^4-4 a \[ScriptCapitalE] \[ScriptCapitalL] \[ScriptCapitalQ]^4+4 a \[ScriptCapitalE]^3 \[ScriptCapitalL] \[ScriptCapitalQ]^4+5 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^4-4 a^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^4-5 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^4+4 a^2 \[ScriptCapitalE]^2 \[ScriptCapitalL]^2 \[ScriptCapitalQ]^4+\[ScriptCapitalQ]^5-a^2 \[ScriptCapitalQ]^5-\[ScriptCapitalE]^2 \[ScriptCapitalQ]^5+a^2 \[ScriptCapitalE]^2 \[ScriptCapitalQ]^5);


(* ::Subsection:: *)
(*From initial data to p, e, x and phases: KerrGeoInit2Constants, KerrGeoInit2pex, KerrGeoInit2Phases*)


KerrGeoInit2Constants[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List]:=Module[{g=gK[a][t0,r0,\[Theta]0,\[Phi]0],\[ScriptCapitalK]=\[ScriptCapitalK][a][t0,r0,\[Theta]0,\[Phi]0],norm},
	norm=u . g . u;
	If[norm<0,
		Return[<|"E"->{-1,0,0,0} . g . (u/Sqrt[-norm]),"L"->{0,0,0,1} . g . (u/Sqrt[-norm]),"Q"->(u/Sqrt[-norm]) . \[ScriptCapitalK] . (u/Sqrt[-norm])|>],
		Message[KerrGeoInit2Constants::imnorm]; Return[$Failed]
	]
];
KerrGeoInit2Constants::imnorm = "Initial velocity is not timelike.";


KerrGeoInit2pex[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List,{\[ScriptCapitalE]_:Null,\[ScriptCapitalL]_:Null,\[ScriptCapitalQ]_:Null}]:=Module[{En,L,Q,rts,disc,r1,r2,zm},
	If[MemberQ[{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]},Null],{En,L,Q}=(Values@KerrGeoInit2Constants[a,{t0,r0,\[Theta]0,\[Phi]0},u]),{En,L,Q}={\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]}];
	disc=RDisc[a,En,L,Q];
	rts=Sort@Re@RRoots[a,En,L,Q][[All,1,2]];
	If[PossibleZeroQ[a],
		zm = Q/(L^2+Q),  
		zm = 1/(2 a^2 (1-En^2)) (Q +L^2+a^2 (1-En^2)-Sqrt[(Q+L^2+a^2 (1-En^2))^2-4 a^2 (1-En^2) Q])
	];
	Which[
		Q<0,
		Message[KerrGeoInit2pex::vortical],
		disc<0,
		Message[KerrGeoInit2pex::plunge],
		En>1,
		Message[KerrGeoInit2pex::unbound],
		True,
		r1=If[En<1,rts[[4]],rts[[1]]];
		r2=If[En<1,rts[[3]],rts[[4]]];
		Return[<|"p"->(2 r1 r2)/(r1+r2),"e"->(r1-r2)/(r1+r2),"x"->Sign[L] Sqrt[1-zm]|>]
	];
	$Failed
];
KerrGeoInit2pex[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List]:=KerrGeoInit2pex[a,{t0,r0,\[Theta]0,\[Phi]0},u,{}]; (*To make elegant notation for optional arguments*)
KerrGeoInit2pex::vortical = "Vortical (polar-cone) orbit detected";
KerrGeoInit2pex::plunge = "Plunge orbit detected";
KerrGeoInit2pex::unbound = "Unbound (escaping) orbit detected";


KerrGeoInit2Phases[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List,{\[ScriptCapitalE]_:Null,\[ScriptCapitalL]_:Null,\[ScriptCapitalQ]_:Null},{pp_:Null,ee_:Null,xx_:Null}] := Module[{zm,zp,mr,\[Lambda]r0,m\[Theta],\[Lambda]\[Theta]0,\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t,r1,r2,r3,r4,qr0,q\[Theta]0,\[Psi]r0,\[Psi]\[Theta]0,En,L,Q,p,e,x},
	If[MemberQ[{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]},Null],
		{En,L,Q}=(Values@KerrGeoInit2Constants[a,{t0,r0,\[Theta]0,\[Phi]0},u]),
		{En,L,Q}={\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]}
	];
	If[MemberQ[{pp,ee,xx},Null],
		{p,e,x} = (Values@KerrGeoInit2pex[a,{t0,r0,\[Theta]0,\[Phi]0},u,{En,L,Q}]),
		{p,e,x} = {pp,ee,xx}
	];
	{\[CapitalUpsilon]r,\[CapitalUpsilon]\[Theta],\[CapitalUpsilon]\[Phi],\[CapitalUpsilon]t} = Values[KerrGeodesics`OrbitalFrequencies`Private`KerrGeoMinoFrequencies[a,p,e,x]];
	{r1,r2,r3,r4} = KerrGeodesics`OrbitalFrequencies`Private`KerrGeoRadialRoots[a, p, e, x, En, Q];
	(*See Fujita & Hikida (0906.1420) section 4*)
	(*Preliminary definitions:*)
	\[Psi]r0 = ArcSin[Sqrt[(r1-r3)/(r1-r2) (r0-r2)/(r0-r3)]];
	mr = (r1-r2)/(r1-r3) (r3-r4)/(r2-r4); 
	(*Radial and polar phases*)
	\[Lambda]r0 = 1/Sqrt[1-En^2] 2/Sqrt[(r1-r3)(r2-r4)] EllipticF[\[Psi]r0,mr];
	If[PossibleZeroQ[a],
		zm = Q/(L^2+Q),  
		zm = 1/(2 a^2 (1-En^2)) (Q +L^2+a^2 (1-En^2)-Sqrt[(Q+L^2+a^2 (1-En^2))^2-4 a^2 (1-En^2) Q])
	];
	zp = 1/(2 a^2 (1-En^2)) (Q +L^2+a^2 (1-En^2)+Sqrt[(Q+L^2+a^2 (1-En^2))^2-4 a^2 (1-En^2) Q]);
	m\[Theta] = zm/zp;
	\[Psi]\[Theta]0 = ArcSin[ Cos[\[Theta]0]/Sqrt[zm]];
	\[Lambda]\[Theta]0 = Sign[L]/Sqrt[(Q +L^2+a^2 (1-En^2)+Sqrt[(Q+L^2+a^2 (1-En^2))^2-4 a^2 (1-En^2) Q])/2] EllipticF[\[Psi]\[Theta]0,m\[Theta]];
	If[u[[2]]>0,
		qr0 = \[Lambda]r0 \[CapitalUpsilon]r,
		qr0 = 2\[Pi] - \[Lambda]r0 \[CapitalUpsilon]r
	];
	Which[
		\[Theta]0==\[Pi]/2&&u[[3]]==0,q\[Theta]0=0,
		u[[3]]>0,q\[Theta]0 = \[Lambda]\[Theta]0 \[CapitalUpsilon]\[Theta],
		u[[3]]>=0, q\[Theta]0 = 2\[Pi] - \[Lambda]\[Theta]0 \[CapitalUpsilon]\[Theta];
	];
	(*Note that the current KerGeoOrbit implementation chooses integration constants so that at \[Lambda]=0 one has t0=qt0,\[Phi]=q\[Phi]0.*)
	{"\!\(\*SubscriptBox[SuperscriptBox[\(q\), \(t\)], \(0\)]\)"->t0,"\!\(\*SubscriptBox[SuperscriptBox[\(q\), \(r\)], \(0\)]\)"-> Re[qr0],"\!\(\*SubscriptBox[SuperscriptBox[\(q\), \(\[Theta]\)], \(0\)]\)"->Re[q\[Theta]0],"\!\(\*SubscriptBox[SuperscriptBox[\(q\), \(\[Phi]\)], \(0\)]\)"->\[Phi]0}
	(*Real parts used to get rid of imaginary noise near turning points.*)
];
KerrGeoInit2Phases[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List,{En_,L_,Q_}] :=KerrGeoInit2Phases[a,{t0,r0,\[Theta]0,\[Phi]0},u,{En,L,Q},{}] ;
KerrGeoInit2Phases[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List] :=KerrGeoInit2Phases[a,{t0,r0,\[Theta]0,\[Phi]0},u,{},{}] ;


(* ::Subsection::Closed:: *)
(*Old KerrGeoInitOrbit wrapper function*)


(*Old wrapper function, now the functionality is accessed directly from KerrGeoOrbit.*)
(*Clear[KerrGeoInitOrbit]
Options[KerrGeoInitOrbit] = {"Parametrization" -> "Mino", "Method" -> "FastSpec"};
SyntaxInformation[KerrGeoOrbit] = {"ArgumentsPattern"->{___,OptionsPattern[]}};
KerrGeoInitOrbit[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List,{\[ScriptCapitalE]_:Null,\[ScriptCapitalL]_:Null,\[ScriptCapitalQ]_:Null},{pp_:Null,ee_:Null,xx_:Null}] := Module[{initPhases,En,L,Q,p,e,x},
	(*If the optional arguments are supplemented, it is the responsibility of the user that they are consistent with the initial data!*)
	If[MemberQ[{\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]},Null],{En,L,Q}=(Values@KerrGeoInit2Constants[a,{t0,r0,\[Theta]0,\[Phi]0},u]),{En,L,Q}={\[ScriptCapitalE],\[ScriptCapitalL],\[ScriptCapitalQ]}];
	If[MemberQ[{pp,ee,xx},Null],{p,e,x} = (Values@KerrGeoInit2pex[a,{t0,r0,\[Theta]0,\[Phi]0},u,{En,L,Q}]),{p,e,x} = {pp,ee,xx}];
	initPhases = KerrGeoInit2Phases[a,{t0,r0,\[Theta]0,\[Phi]0},u,{En,L,Q},{p,e,x}];
	KerrGeoOrbit[a,p,e,x,initPhases,"Parametrization"->OptionValue["Parametrization"],"Method"->OptionValue["Method"]]
]
KerrGeoInitOrbit[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List,{En_,L_,Q_}] :=KerrGeoInitOrbit[a,{t0,r0,\[Theta]0,\[Phi]0},u,{En,L,Q},{}] ;
KerrGeoInitOrbit[a_,{t0_,r0_,\[Theta]0_,\[Phi]0_},u_List] :=KerrGeoInitOrbit[a,{t0,r0,\[Theta]0,\[Phi]0},u,{},{}]; (*To make elegant notation for optional arguments*)*)


(* ::Section::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
