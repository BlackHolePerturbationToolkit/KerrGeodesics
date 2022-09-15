(* ::Package:: *)

(* ::Title::Closed:: *)
(*NullGeodesics subpackage*)


(* ::Section:: *)
(*Define usage for public functions*)


BeginPackage["KerrGeodesics`NullGeodesics`"]

KerrNullEnergy::usage = "test fuction"

Begin["`Private`"];


(* ::Section:: *)
(*Constants of Motion*)


KerrNullEnergy[]:=Print["It works!"]


(* ::Title:: *)
(*Constants of Motion*)


(* ::Section:: *)
(*Spherical photon orbits constants of motion*)


(* ::Input::Initialization:: *)
KerrSphericalNullCoM[a_,rs_]:=Module[{M=1,\[Lambda]s,\[Eta]s},
\[Lambda]s=a+rs/a (rs-(2 (rs^2-2 M rs+a^2))/(rs-M))(*eqn 104*);
\[Eta]s=rs^3/a^2 ((4 M (rs^2-2 M rs+a^2))/(rs-M)^2-rs)(*eqn 108*);
 
{"\[Lambda]s"-> \[Lambda]s,"\[Eta]s"-> \[Eta]s}
] 


(* ::Input::Initialization:: *)
SNullCoM[a_,rs_,\[Theta]s_]:=Module[{M=1,\[Lambda]s,\[Eta]s},
\[Lambda]s=Sqrt[27/(1+Cot[\[Theta]s]^2)](*eqn 104*);
\[Eta]s=27-\[Lambda]s^2(*eqn 108*);
 
{"\[Lambda]s"-> \[Lambda]s,"\[Eta]s"-> \[Eta]s}
]


(* ::Input::Initialization:: *)
SchwarzSphericalNullCoM[a_,ts_,rs_,\[Theta]s_,\[Phi]s_,prs_,p\[Theta]s_,p\[Phi]s_]:=Module[{M=1,\[Lambda],\[Eta]},
\[Eta]=0;
\[Lambda]=(Sqrt[2] rs^(3/2))/Sqrt[-M+rs]; 
(*Print[rs^4-(rs^2-2M rs)(\[Eta]+\[Lambda]^2)];*)
{"\[Lambda]"-> \[Lambda],"\[Eta]"-> \[Eta]} 
]


(* ::Title::Closed:: *)
(*Radial Potential Roots*)


(* ::Input::Initialization:: *)
KerrNullGeoRadialRoots[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rt,NS,consts,A,B,C1,P,Q1,\[Omega]p,\[Omega]m,\[Xi]0,z,r1,r2,r3,r4,xp,xm,test1,test2,dis,t,\[Xi]1,\[Xi]2,x,y,R,\[CapitalTheta]},
(*consts=ConstantsOfMotion[a,ts,rs,\[Theta]s,\[Phi]s,prs,p\[Theta]s,p\[Phi]s];
{\[Lambda],\[Eta]}=Values[consts];*)
(*Print[Chop[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2),10^-10]];*)
NS=If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,NSolve[(rt^2+a^2-a \[Lambda])^2-(rt^2-2 M rt+a^2)(\[Eta]+(\[Lambda]-a)^2)==0,rt],"Error:negative Radial Potential"];
{r1,r2,r3,r4}={NS[[1,1,2]],NS[[2,1,2]],NS[[3,1,2]],NS[[4,1,2]]};
{"\!\(\*SubscriptBox[\(r\), \(1\)]\)"->r1,"\!\(\*SubscriptBox[\(r\), \(2\)]\)"->r2,"\!\(\*SubscriptBox[\(r\), \(3\)]\)"->r3,"\!\(\*SubscriptBox[\(r\), \(4\)]\)"->r4}]


(* ::Input::Initialization:: *)
zamo[a_]:=Module[{M=1,r},
r=M+2Sqrt[M^2-1/3 a^2]Cos[1/3 ArcCos[(M(M^2-a^2))/(M^2-1/3 a^2)^(3/2)]]
]

(* ::Title:: *)
(*SPHERICAL and NON-SPHERICAL PHOTON ORBITS*)


(* ::Chapter:: *)
(*Radial Motion*)


(* ::Section::Initialization:: *)
(*case One*)


(* ::Input::Initialization:: *)
case1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{ro,X1,k,Irs,F1,x1,roots,r1,r2,r3,r4,assoc,consts,\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x1=Sqrt[(rs-r2)/(rs-r1) (r3-r1)/(r3-r2)](*eqn B15*);
F1=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x1],k](*eqn B20*);
Irs=F1;
X1[\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 ( \[Tau]+\[Nu]r Irs)(*eqn B26*);
assoc=Association["Trajectory"-> ro,"RadialRoots"-> roots,"ConstantsofMotion"-> consts];
ro[\[Tau]_]:=(r2(r3-r1)-r1(r3-r2)JacobiSN[X1[\[Tau]],k]^2)/((r3-r1)-(r3-r2)JacobiSN[X1[\[Tau]],k]^2)(*eqn B27*);
case1Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]]


(* ::Input::Initialization:: *)
Format[case1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="case1Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
case1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
case1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Max\[Tau]1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,ro,X1,k,Irs,F1,x1,roots,r1,r2,r3,r4,assoc,consts,\[Nu]r,rp,\[Tau],tss,t3,tp,Is,x1s,\[Kappa]s,\[Tau]1,\[Tau]2},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
rp=M+Sqrt[M-a^2];
x1s=Sqrt[((r4-r2)(r3-rs))/((r3-r2)(r4-rs))]; 
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))]; 
\[Kappa]s=ArcSin[x1s];
Is=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[\[Kappa]s,k];  
\[Tau]=Floor[Re[1/(r1 r2-r2 r3-r1 r4+r3 r4) (Is r1 r2 \[Nu]r-Is r2 r3 \[Nu]r-Is r1 r4 \[Nu]r+Is r3 r4 \[Nu]r+2 Sqrt[(r1-r3) (r2-r4)] InverseJacobiSN[(Sqrt[r2-r4] Sqrt[r3-rp])/Sqrt[r2 r4-r3 r4-r2 rp+r3 rp],k])],1/1000]
]


(* ::Section::Initialization:: *)
(*case Two*)


(* ::Input::Initialization:: *)
case2[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{r1,r2,r3,r4,x2,k,F2,Irs,X2,r0,roots,assoc,consts,\[Nu]r,Kk,\[Tau],prd,\[Phi]},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
 x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
\[Phi]=ArcSin[x2];
F2=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*Sin[\[Phi]]CarlsonRF[Cos[\[Phi]]^2,1-k Sin[\[Phi]]^2,1]*)(*eqn B40*);
 Irs=F2; X2[\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 ( \[Tau]+ \[Nu]r Irs)(*eqn B45*);  
r0[\[Tau]_]:=(r4(r3-r1)-r3(r4-r1)JacobiSN[X2[\[Tau]],k]^2)/((r3-r1)-(r4-r1)JacobiSN[X2[\[Tau]],k]^2)(*eqn B46*); 
assoc=Association["Trajectory"->r0,"RadialRoots"-> roots,"ConstantsofMotion"-> consts]; 
case2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]]


(* ::Input::Initialization:: *)
\[Tau]r4[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],Irs,F2,\[Phi],k,x2,r1,r2,r3,r4,roots,\[Nu]r},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
 x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
\[Phi]=ArcSin[x2];
F2=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*Sin[\[Phi]]CarlsonRF[Cos[\[Phi]]^2,1-k Sin[\[Phi]]^2,1]*)(*eqn B40*);
 Irs=F2; 
\[Tau]=(-Irs r1 r2 \[Nu]r+Irs r2 r3 \[Nu]r+Irs r1 r4 \[Nu]r-Irs r3 r4 \[Nu]r)/(r1 r2-r2 r3-r1 r4+r3 r4)
] 


(* ::Input::Initialization:: *)
tauprac[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{r1,r2,r3,r4,x2,k,F2,Irs,X2,r0,roots,assoc,consts,\[Nu]r,Kk,\[Tau],prd,\[Phi],tau},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
 x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
\[Phi]=ArcSin[x2];
F2=2/Sqrt[(r3-r1)(r4-r2)](*EllipticF[ArcSin[x2],k]*)Sin[\[Phi]]CarlsonRF[Cos[\[Phi]]^2,1-k Sin[\[Phi]]^2,1](*eqn B40*);
 Irs=F2;
tau=(-Irs r1 r2 \[Nu]r+Irs r2 r3 \[Nu]r+Irs r1 r4 \[Nu]r-Irs r3 r4 \[Nu]r)/(r1 r2-r2 r3-r1 r4+r3 r4)(*eqn B46*)
]


(* ::Input::Initialization:: *)
Format[case2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="case2Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
case2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
case2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Max\[Tau]2p[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],\[Nu]r,roots,r1,r2,r3,r4,x2,k,F2,Irs},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
F2=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*eqn B40*);
 Irs=F2;
\[Tau]=Re[1/(r1 r2-r2 r3-r1 r4+r3 r4) (-Irs r1 r2 \[Nu]r+Irs r2 r3 \[Nu]r+Irs r1 r4 \[Nu]r-Irs r3 r4 \[Nu]r-2 Sqrt[(r1-r3) (r2-r4)] InverseJacobiSN[(Sqrt[r1-r3] Sqrt[-200+r4-rs])/Sqrt[-200 r1+r1 r3+200 r4-r3 r4-r1 rs+r4 rs],k])]]  


(* ::Input::Initialization:: *)
Max\[Tau]2nn[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,\[Tau],\[Nu]r,roots,r1,r2,r3,r4,x2,k,F2,Irs,\[Tau]p,\[Tau]4,\[Tau]p4},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rp=M+Sqrt[M-a^2];
x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
F2=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*eqn B40*);
 Irs=F2;
\[Tau]=Floor[Re[1/(r1 r2-r2 r3-r1 r4+r3 r4) (-Irs r1 r2 \[Nu]r+Irs r2 r3 \[Nu]r+Irs r1 r4 \[Nu]r-Irs r3 r4 \[Nu]r+2 Sqrt[(r1-r3) (r2-r4)] InverseJacobiSN[(Sqrt[r1-r3] Sqrt[r4-rp])/Sqrt[r1 r3-r3 r4-r1 rp+r4 rp],k])],1/1000]]  


(* ::Input::Initialization:: *)
Max\[Tau]2n[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],\[Nu]r,roots,r1,r2,r3,r4,x2,k,F2,Irs},
\[Nu]r=prs;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
F2=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*eqn B40*);
 Irs=F2;
\[Tau]=Re[1/(r1 r2-r2 r3-r1 r4+r3 r4) (-Irs r1 r2 \[Nu]r+Irs r2 r3 \[Nu]r+Irs r1 r4 \[Nu]r-Irs r3 r4 \[Nu]r+2 Sqrt[(r1-r3) (r2-r4)] InverseJacobiSN[(10 Sqrt[2] Sqrt[-r1+r3])/Sqrt[-200 r1+r1 r3+200 r4-r1 r4-r3 r4+r4^2],k])]]


(* ::Section::Initialization:: *)
(*Radial PNG*)


(* ::Input::Initialization:: *)
PrincipalNullGeodesics[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,ro,assoc,consts,\[Nu]r,\[CapitalTheta],r1,r2},
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
r1=I a Cos[\[Theta]s];
r2=-I a Cos[\[Theta]s];
\[Nu]r=prs;
(*Print[\[Eta]+(\[Lambda]-a)^2];*)
ro[\[Tau]_]:=If[a==0\[Or]\[Theta]s==\[Pi]/2,-((rs \[Nu]r)/(-\[Nu]r+rs \[Tau])),Re[(E^((r1 \[Tau])/\[Nu]r) r1 r2-E^((r2 \[Tau])/\[Nu]r) r1 r2+E^((r2 \[Tau])/\[Nu]r) r1 rs-E^((r1 \[Tau])/\[Nu]r) r2 rs)/(E^((r1 \[Tau])/\[Nu]r) r1-E^((r2 \[Tau])/\[Nu]r) r2-E^((r1 \[Tau])/\[Nu]r) rs+E^((r2 \[Tau])/\[Nu]r) rs)]]; 
assoc=Association["Trajectory"-> ro,"ConstantsofMotion"-> consts];
PrincipalNullGeodesicsFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[PrincipalNullGeodesicsFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="PrincipalNullGeodesicsFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
PrincipalNullGeodesicsFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
PrincipalNullGeodesicsFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Max\[Tau]PNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Tau],\[Nu]r,r1,r2,rp},
r1=I a Cos[\[Theta]s];
r2=-I a Cos[\[Theta]s];
\[Nu]r=prs;
rp=M+Sqrt[M-a^2];
\[Tau]=If[a==0\[Or]\[Theta]s==\[Pi]/2,If[\[Nu]r<0,Floor[\[Nu]r(-(1/rp)+1/rs),1/1000],(-1/200+1/rs)],If[\[Nu]r<0,Floor[Re[(\[Nu]r (Log[-(((r1-rp) (r2-rs))/((-r2+rp) (r1-rs)))]))/(r1-r2)],1/1000],Re[-((\[Nu]r (Log[((-200+r2) (r1-rs))/((-200+r1) (r2-rs))]))/(r1-r2))]]]   
]


(* ::Section::Initialization:: *)
(*case Three*)


(* ::Input::Initialization:: *)
case3[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{ro3,A,B,k3,X3\[Tau],Irs,F3s,x3s,assoc,roots,r1,r2,r3,r4,consts,\[Nu]r,D5},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
x3s=(1- (B(rs-r2))/(A (rs-r1)))/(1+(B(rs-r2))/(A(rs-r1)))(*eqn B58*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+\[Nu]r Irs)(*eqn B74*);
assoc=Association["Trajectory"-> ro3,"RadialRoots"-> roots,"ConstantsofMotion"-> consts];
ro3[\[Tau]_]:=((B r2-A r1)+(B r2 +A r1)JacobiCN[X3\[Tau][\[Tau]],k3])/((B-A)+(B+A)JacobiCN[X3\[Tau][\[Tau]],k3])(*eqn B75*);
case3Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[case3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="case3Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
case3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
case3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Max\[Tau]3n[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],roots,r1,r2,r3,r4,consts,\[Nu]r,A,B,k3,x3s,F3s,Irs,rp,M=1,tp,t2,tp2},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rp=M+Sqrt[M-a^2];
\[Nu]r=prs;
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
x3s=(1- (B(rs-r2))/(A (rs-r1)))/(1+(B(rs-r2))/(A(rs-r1)))(*eqn B58*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
tp=(-A B Irs \[Nu]r+Sqrt[A B] InverseJacobiCN[(A r1-B r2-A rp+B rp)/(A r1+B r2-A rp-B rp),k3])/(A B);
t2=-Irs \[Nu]r;
tp2=tp-t2;
\[Tau]=(*Floor[Re[-Irs \[Nu]r-tp2],1/100]*)(*Floor[Re[tp-2(tp-t2)],1/100]*)(*tp-2(tp-t2)*)Floor[Re[tp-2(tp-t2)],1/100]
]


(* ::Input::Initialization:: *)
Max\[Tau]3p[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],roots,r1,r2,r3,r4,consts,\[Nu]r,B,A,k3,x3s,F3s,Irs},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
x3s=(1- (B(rs-r2))/(A (rs-r1)))/(1+(B(rs-r2))/(A(rs-r1)))(*eqn B58*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
\[Tau]=Re[1/(A B) (-A B Irs \[Nu]r+Sqrt[A B] InverseJacobiCN[(-200 A+200 B+A r1-B r2-A rs+B rs)/(-200 A-200 B+A r1+B r2-A rs-B rs),k3])  ]   
]


(* ::Section::Initialization:: *)
(*case Four*)


(* ::Input::Initialization:: *)
case4[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{ro4,g0,X4\[Tau],Irs,F4s,k4,x4,a2,b1,b2,C,D,roots,r1,r2,r3,r4,assoc,consts,\[Nu]r},roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*); 
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b1=Re[(r3+r4)/2](*eqn B10*);
b2=Re[(r1+r2)/2](*eqn B11*);
x4=(rs+b1)/a2(*eqn B83*);
k4=( 4 C D)/(C+D)^2(*eqn B87*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2 -4 a2^2)](*eqn B88*);
F4s=2/(C+D) EllipticF[ArcTan[x4]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
X4\[Tau][\[Tau]_]:=(C+D)/2 (\[Nu]r \[Tau]+Irs)(*eqn B104*);
assoc=Association["Trajectory"-> ro4,"RadialRoots"-> roots,"ConstantsofMotion"-> consts];
ro4[\[Tau]_]:=-a2((g0-JacobiSC[X4\[Tau][\[Tau]],k4])/(1+g0 JacobiSC[X4\[Tau][\[Tau]],k4]))-b1(*eqn B109*);
case4Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[case4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="case4Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
case4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
case4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Max\[Tau]4p[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Tau],roots,r1,r2,r3,r4,consts,\[Nu]r,C,D,a2,b1,b2,x4,k4,g0,F4s,Irs,tss,tsp},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*); 
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b1=Re[(r3+r4)/2](*eqn B10*);
b2=Re[(r1+r2)/2](*eqn B11*);
x4=(rs-b2)/a2(*eqn B83*);   
k4=( 4 C D)/(C+D)^2(*eqn B87*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2 -4 a2^2)](*eqn B88*);
F4s=2/(C+D) EllipticF[ArcTan[x4]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
tss=(-C Irs-D Irs+2 InverseJacobiSC[(b1+a2 g0+rs)/(a2-b1 g0-g0 rs),k4])/((C+D) \[Nu]r);
tsp=(-C Irs-D Irs+2 InverseJacobiSC[(200+b1+a2 g0+rs)/(a2-200 g0-b1 g0-g0 rs),k4])/((C+D) \[Nu]r); 
\[Tau]=Re[tsp-tss ]  
]


(* ::Input::Initialization:: *)
Max\[Tau]4n1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,\[Tau],roots,r1,r2,r3,r4,consts,\[Nu]r,C,D,a2,b1,b2,x4,k4,g0,F4s,Irs,tp,tss,td},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rp=M+Sqrt[M-a^2];
\[Nu]r=prs;
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*); 
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b1=Re[(r3+r4)/2](*eqn B10*);
b2=Re[(r1+r2)/2](*eqn B11*);
x4=(rs+b1)/a2(*eqn B83*); 
k4=( 4 C D)/(C+D)^2(*eqn B87*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2 -4 a2^2)](*eqn B88*);
F4s=2/(C+D) EllipticF[ArcTan[x4]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
tp=(-C Irs-D Irs+2 InverseJacobiSC[(b1+a2 g0+rp)/(a2-b1 g0-g0 rp),k4])/((C+D) \[Nu]r)  ; 
tss=(-C Irs-D Irs+2 InverseJacobiSC[(b1+a2 g0+rs)/(a2-b1 g0-g0 rs),k4])/((C+D) \[Nu]r) ;
td=tp-tss;
\[Tau]=Re[td]
 ]


(* ::Input::Initialization:: *)
Max\[Tau]4n2[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,\[Tau],roots,r1,r2,r3,r4,consts,\[Nu]r,C,D,a2,b1,b2,x4,k4,g0,F4s,Irs,tp,tss,td},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rp=M+Sqrt[M-a^2];
\[Nu]r=prs;
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*); 
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b1=Re[(r3+r4)/2](*eqn B10*);
b2=Re[(r1+r2)/2](*eqn B11*);
x4=(rs+b1)/a2(*eqn B83*);
k4=( 4 C D)/(C+D)^2(*eqn B87*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2 -4 a2^2)](*eqn B88*);
F4s=2/(C+D) EllipticF[ArcTan[x4]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
tp=(-C Irs-D Irs+2 InverseJacobiSC[(b1+a2 g0+rp)/(a2-b1 g0-g0 rp),k4])/((C+D) \[Nu]r)  ; 
tss=(-C Irs-D Irs+2 InverseJacobiSC[(b1+a2 g0+rs)/(a2-b1 g0-g0 rs),k4])/((C+D) \[Nu]r) ;
td=tp-tss;
(*\[Tau]=Re[td]*)
\[Tau]=Re[tp ]  
 ]


(* ::Section::Initialization:: *)
(*Radial spherical motion*)


(* ::Input::Initialization:: *)
rSpherical[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,ro,assoc,consts,\[Nu]r,\[CapitalTheta],roots,r1,r2,r3,r4},
assoc=Association["Trajectory"-> ro,"ConstantsofMotion"-> consts];
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
ro[\[Tau]_]:=r4;
rSphericalFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[rSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="rSphericalFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
rSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
rSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Section::Initialization:: *)
(*Radial motion generalization code*)


(* ::Input::Initialization:: *)
RadialMotionGeneral[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Nu]r,rp,rm,r1,r2,r3,r4,roots,\[CapitalTheta],\[CapitalTheta]1,consts,R,constss,\[Lambda]s,\[Eta]s,k},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2));
If[(*r3\[Equal]r4==rs*)4 rs (a^2+rs^2-a \[Lambda])-(-2 M+2 rs) (\[Eta]+(-a+\[Lambda])^2)==0\[And](rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)==0,rSpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>r4,If[k<1\[And]k>0,Return[case2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["These values for the energy and angular momentum correspond to a critical orbit but the initial position does not. Please enter instead values for the energy and angular momentum which are very close -but not exactly equal- to the critical ones."]],If[rs<r3,Return[case1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[case3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[case4[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]]


(* ::Input::Initialization:: *)
RadialMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{consts,M=1},
If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,If[\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,PrincipalNullGeodesics[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s], RadialMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["Error:Negative Radial Potential"]]
]


(* ::Input::Initialization:: *)
RR[rs_]:=(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)


(* ::Input::Initialization:: *)
PolarSphericalPeriod[a_,rs_,\[Phi]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,G,\[CapitalDelta]\[Theta],up,um,consts},
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
G=4/Sqrt[-um a^2] EllipticK[up/um]]


(* ::Input::Initialization:: *)
Max\[Tau]General[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,rm,r1,r2,r3,r4,roots,\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
\[Nu]r=Sign[prs];
If[r3==r4,If[a==0,If[\[Theta]s==\[Pi]/2,2\[Pi],4/Sqrt[\[Eta]+\[Lambda]^2] EllipticK[0]],PolarSphericalPeriod[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,If[r4>rp,If[\[Nu]r>=0,Max\[Tau]2p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Max\[Tau]2n[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Max\[Tau]2nn[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Max\[Tau]1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],If[\[Nu]r>=0,Max\[Tau]3p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Max\[Tau]3n[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],If[\[Nu]r<0,If[a<0.8,Max\[Tau]4n1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Max\[Tau]4n2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Max\[Tau]4p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]
]  


(* ::Input::Initialization:: *)
Max\[Tau][a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{consts,M=1},
(*If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)\[GreaterEqual]0,If[\[Eta]==-a^2Cos[\[Theta]s]^4\[And]\[Lambda]\[Equal]a Sin[\[Theta]s]^2,Max\[Tau]PNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Max\[Tau]General[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],"Error:Negative Radial Potential"]*)If[\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,Max\[Tau]PNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Max\[Tau]General[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]
]


(* ::Chapter::Initialization:: *)
(*Angular Potential Roots*)


(* ::Input::Initialization:: *)
KerrNullGeoAngularRoots[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts},
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
\[Theta]1=ArcCos[Sqrt[up]](*eqn 20*);  
\[Theta]2=ArcCos[Sqrt[um]](*eqn 21*); 
\[Theta]3=ArcCos[-Sqrt[um]](*eqn 22*); 
\[Theta]4=ArcCos[-Sqrt[up]](*eqn 23*); 
{"\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\)"->\[Theta]1,"\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\)"->\[Theta]2,"\!\(\*SubscriptBox[\(\[Theta]\), \(3\)]\)"->\[Theta]3,"\!\(\*SubscriptBox[\(\[Theta]\), \(4\)]\)"->\[Theta]4}]


(* ::Input::Initialization:: *)
SNullGeoAngularRoots[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts},
up=\[Eta]/(\[Eta]+\[Lambda]^2); 
 um=\[Eta]/(\[Eta]+\[Lambda]^2); 
\[Theta]1=ArcCos[Sqrt[up]](*eqn 20*);  
\[Theta]2=ArcCos[Sqrt[um]](*eqn 21*); 
\[Theta]3=ArcCos[-Sqrt[um]](*eqn 22*); 
\[Theta]4=ArcCos[-Sqrt[up]](*eqn 23*); 
{"\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\)"->\[Theta]1,"\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\)"->\[Theta]2,"\!\(\*SubscriptBox[\(\[Theta]\), \(3\)]\)"->\[Theta]3,"\!\(\*SubscriptBox[\(\[Theta]\), \(4\)]\)"->\[Theta]4}]


(* ::Chapter:: *)
(*Polar Motion*)


(* ::Section::Initialization:: *)
(*Ordinary Motion (\[Eta] > 0)*)


(* ::Input::Initialization:: *)
ordinary[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta],u0},
roots=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[roots];
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2)(*eqn 19*);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
u0=-\[Eta]/a^2; 
\[Theta]p=ArcCos[-Sqrt[up]] (*eqn 28*);
\[Theta]m=ArcCos[Sqrt[up]];
Gcurl\[Theta]=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],-(1/Sqrt[- a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[u0]],u0],-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]](*eqn 29*);
assoc=Association["Trajectory"-> \[Theta]o,"AngularRoots"-> roots,"ConstantsofMotion"-> consts];
\[Theta]o[\[Tau]_]:=If[\[Theta]s==0 \[Or]\[Theta]s==\[Pi],(*-(-(\[Pi]/2)+JacobiAmplitude[(-Sqrt[\[Eta]] \[Tau]+\[Nu]\[Theta] EllipticF[\[Pi]/2-\[Theta]s,-(a^2/\[Eta])])/\[Nu]\[Theta],-(a^2/\[Eta])])*)Re[ArcCos[-\[Nu]\[Theta] Sqrt[u0]JacobiSN[Sqrt[-a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),u0]]],ArcCos[-\[Nu]\[Theta] Sqrt[up]JacobiSN[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um]]](*eqn 38*);
ordinaryFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[ordinaryFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinaryFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinaryFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinaryFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
KerrNullGeoAngularRootsx[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts},
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
\[Theta]1=ArcCos[Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2+Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2]](*eqn 20*);  
\[Theta]2=ArcCos[Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2-Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2]](*eqn 21*); 
\[Theta]3=ArcCos[-(Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2-Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2])](*eqn 22*); 
\[Theta]4=ArcCos[-(Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2+Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2])](*eqn 23*); 
{"\!\(\*SubscriptBox[\(\[Theta]\), \(1\)]\)"->\[Theta]1,"\!\(\*SubscriptBox[\(\[Theta]\), \(2\)]\)"->\[Theta]2,"\!\(\*SubscriptBox[\(\[Theta]\), \(3\)]\)"->\[Theta]3,"\!\(\*SubscriptBox[\(\[Theta]\), \(4\)]\)"->\[Theta]4}]


(* ::Input::Initialization:: *)
ordinaryx[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta],u0},
roots=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[roots];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2)(*eqn 19*);
up=Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2+Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2](*eqn 19*); 
 um=Sqrt[1-\[Eta]/a^2-\[Lambda]^2/a^2-Sqrt[4 a^2 \[Eta]+(-a^2+\[Eta]+\[Lambda]^2)^2]/a^2]/Sqrt[2](*eqn 19*); 
u0=\[Eta]/(\[Eta]+\[Lambda]^2); 
\[Theta]p=ArcCos[-Sqrt[up]] (*eqn 28*);
\[Theta]m=ArcCos[Sqrt[up]];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
Print[Sqrt[um^2-up^2]];
assoc=Association["Trajectory"-> \[Theta]o,"AngularRoots"-> roots,"ConstantsofMotion"-> consts];
\[Theta]o[\[Tau]_]:=(*ArcCos[um JacobiND[up a^2\[Tau],Sqrt[1-um^2/up^2]]]*)ArcCos[up JacobiCN[Sqrt[um^2-up^2] a^2 \[Tau],up^2/(um^2-up)]]; 
ordinaryxFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc] 
]     


(* ::Input::Initialization:: *)
Format[ordinaryxFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinaryxFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinaryxFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinaryxFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
ordinaryEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]o,assoc},
assoc=Association["Trajectory"-> \[Theta]o];
\[Theta]o[\[Tau]_]:=\[Theta]s;
ordinaryEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[ordinaryEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinaryEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinaryEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinaryEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Theta]PNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,ro,assoc,consts,\[Nu]r,\[CapitalTheta],\[Theta]o},
\[Nu]r=prs;
\[Theta]o[\[Tau]_]:=\[Theta]s;
assoc=Association["Trajectory"-> \[Theta]o,"ConstantsofMotion"-> consts];
\[Theta]PNGFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Theta]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Theta]PNGFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Theta]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Theta]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
KerrCelestialcordinates\[Alpha][a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{L,\[Alpha],\[Beta],consts,\[Theta],assoc},
\[Theta][\[Tau]_]:=PolarMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
\[Alpha][\[Tau]_]:=-\[Lambda] Csc[\[Theta][\[Tau]]];
assoc=Association["Trajectory"-> \[Alpha]];
KerrCelestialcordinates\[Alpha]Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
](*ovalles, results*)


(* ::Input::Initialization:: *)
Format[KerrCelestialcordinates\[Alpha]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="KerrCelestialcordinates\[Alpha]Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
KerrCelestialcordinates\[Alpha]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
KerrCelestialcordinates\[Alpha]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
KerrCelestialcordinates\[Beta][a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{L,\[Alpha],\[Beta],consts,\[Theta],assoc},
\[Theta][\[Tau]_]:=PolarMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
\[Beta][\[Tau]_]:=Sqrt[\[Eta]+a^2 Cos[\[Theta][\[Tau]]]^2-\[Lambda]^2 Cot[\[Theta][\[Tau]]]^2]; 
assoc=Association["Trajectory"-> \[Beta]];
KerrCelestialcordinates\[Beta]Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[KerrCelestialcordinates\[Beta]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="KerrCelestialcordinates\[Beta]Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
KerrCelestialcordinates\[Beta]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
KerrCelestialcordinates\[Beta]Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[PlusMinus]c9_:={c9,-c9}


(* ::Input::Initialization:: *)
newKDSPS[a_,M_,r_,\[Theta]_,\[CapitalLambda]_]:=Module[{ret,direct,\[Phi],A,p,q,as,b,c,d,B},
d=12 a^2 M;
c=-27 M^2;
b=18 M;
as=-3-4 a^2 \[CapitalLambda];
q=(2b^3-9 as b c+27as^2 d)/(27as^3);
p=(3as c-b^2)/(3as^2);
A=2Sqrt[-(p/3)]; 
\[Phi]=ArcCos[(3q)/(A p)];
B=-b/(3as); 
direct=B+A Cos[\[Phi]/3+(4\[Pi])/3];
ret=B+A Cos[\[Phi]/3];
{"direct"-> direct,"ret"-> ret}]


(* ::Section::Initialization:: *)
(*Schwarzchild case*)


(* ::Input::Initialization:: *)
ordinaryS[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta],u0,r1,r2,r3,r4,u,x},
assoc=Association["Trajectory"-> \[Theta]o];
\[Nu]\[Theta]=p\[Theta]s;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
u=Cos[\[Theta]s];
\[Theta]1=ArcCos[Sqrt[u0]];
\[Theta]2=ArcCos[-Sqrt[u0]];
(*x[\[Tau]_]:=\[Nu]\[Theta] Tan[Sqrt[\[Eta]/u0](\[Tau]+\[Nu]\[Theta] Sqrt[u0/\[Eta]]ArcTan[Cos[\[Theta]s]/Sqrt[u0]/Sqrt[1-Cos[\[Theta]s]/u0]])];
*)(*Print["a,b,c,d,e=",{Re[EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[u0]],0]]}];*)
\[Theta]o[\[Tau]_]:=If[\[Theta]s==\[Pi]\[Or]\[Theta]s==0,(*ArcCos[-\[Nu]\[Theta] JacobiSN[\[Tau]Sqrt[\[Eta]]-\[Nu]\[Theta] EllipticF[ArcSin[Cos[\[Theta]s]],0],0]]*)\[Nu]\[Theta] Sqrt[\[Eta]](\[Tau]+\[Nu]\[Theta] \[Theta]s/Sqrt[\[Eta]]),ArcCos[-\[Nu]\[Theta] Sqrt[u0]Sin[\[Tau] Sqrt[\[Eta]+\[Lambda]^2]-\[Nu]\[Theta] ArcSin[Cos[\[Theta]s]/Sqrt[u0]]]]](*ArcCos[-\[Nu]\[Theta] Sqrt[u0]JacobiSN[\[Tau]Sqrt[\[Eta]+\[Lambda]^2]-\[Nu]\[Theta] EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[u0]],0],0]]*)(*ArcCos[-\[Nu]\[Theta] Sqrt[u0]Sin[\[Tau]Sqrt[\[Eta]+\[Lambda]^2]-\[Nu]\[Theta] ArcSin[Cos[\[Theta]s]/Sqrt[u0]]]]*)(*eqn 38*);
ordinarySFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[ordinarySFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinarySFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinarySFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinarySFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
s[x_]:=-Sin[x]


(* ::Input::Initialization:: *)
ordinarySEq[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta],u0,r1,r2,r3,r4,u,x},
assoc=Association["Trajectory"-> \[Theta]o];
\[Nu]\[Theta]=p\[Theta]s;
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
\[Theta]1=ArcCos[Sqrt[u0]];
\[Theta]2=ArcCos[-Sqrt[u0]];
(*x[\[Tau]_]:=\[Nu]\[Theta] Tan[Sqrt[\[Eta]/u0](\[Tau]+\[Nu]\[Theta] Sqrt[u0/\[Eta]]ArcTan[Cos[\[Theta]s]/Sqrt[u0]/Sqrt[1-Cos[\[Theta]s]/u0]])];
*)(*Print["a,b,c,d,e=",{Re[EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[u0]],0]]}];*)
\[Theta]o[\[Tau]_]:=\[Theta]2(*eqn 38*);
ordinarySEqFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]  


(* ::Input::Initialization:: *)
Format[ordinarySEqFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinarySEqFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinarySEqFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinarySEqFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
ordinarySchwarzchild[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]1,\[Theta]2,u0},
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
\[Theta]1=ArcCos[Sqrt[u0]];
\[Theta]2=ArcCos[-Sqrt[u0]];
If[\[Theta]1==\[Theta]2,ordinarySEq[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],ordinaryS[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]
]


(* ::Input::Initialization:: *)
Format[ordinarySchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinarySchwarzchildFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinarySchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinarySchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Section::Initialization:: *)
(*Vortical Motion (\[Eta]<0)*)


(* ::Input::Initialization:: *)
vortical[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]v,Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,h,assoc,\[Theta]roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,roots,\[CapitalUpsilon],consts,\[Nu]\[Theta]},
roots=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[roots];
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
h=Sign[Cos[\[Theta]s]](*eqn 54*);
\[CapitalUpsilon]=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-(h/Sqrt[um a^2])EllipticF[\[CapitalUpsilon],1-up/um](*eqn 56*);
assoc=Association["Trajectory"-> \[Theta]v,"AngularRoots"-> roots];
\[Theta]v[\[Tau]_]:=ArcCos[Sqrt[um]h JacobiDN[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um]](*eqn 64*);
vorticalFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[vorticalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="vorticalFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
vorticalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
vorticalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Polar Non Spherical motion general code*)


(* ::Input::Initialization:: *)
PolarNonSphericalMotion[a_,ts_,rs_,\[Theta]s_,\[Phi]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,R,M=1},
R=(rs^2+a^2-a \[Lambda])^2-(rs^2-2M rs+a^2)(\[Eta]+(\[Lambda]-a)^2);
If[R>=0\[And]If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],(\[Eta]+a^2 Cos[\[Theta]s]^2)>=0,(\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2)>=0],If[a>0,Which[\[Eta]\[Element]PositiveReals,Return[ordinary[a,ts,rs,\[Theta]s,\[Phi]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinaryEquatorial[a,ts,rs,\[Theta]s,\[Phi]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]\[Element]NegativeReals,Return[vortical[a,ts,rs,\[Theta]s,\[Phi]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[ordinarySchwarzchild[a,ts,rs,\[Theta]s,\[Phi]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],"Error:Radial or Angular Potential is Negative"]
]


(* ::Section::Initialization:: *)
(*Polar Spherical motion*)


(* ::Input::Initialization:: *)
ordinarySpherical[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta]},
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2)(*eqn 19*);
up=If[a==0,\[Eta]/(\[Eta]+\[Lambda]^2),\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]](*eqn 19*); 
 um=If[a==0,\[Eta]/(\[Eta]+\[Lambda]^2),\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]](*eqn 19*); 
\[Theta]p=ArcCos[-Sqrt[up]] (*eqn 28*);
\[Theta]m=ArcCos[Sqrt[up]];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
assoc=Association["Trajectory"-> \[Theta]o,"ConstantsofMotion"-> consts];
\[Theta]o[\[Tau]_]:=If[\[Theta]s==0 \[Or]\[Theta]s==\[Pi],-(-(\[Pi]/2)-\[Nu]\[Theta] JacobiSN[Sqrt[\[Eta]](\[Tau]-\[Nu]\[Theta]/Sqrt[\[Eta]] EllipticF[\[Pi]/2-\[Theta]s,-(a^2/\[Eta])]),-(a^2/\[Eta])])-(-(\[Pi]/2)+JacobiAmplitude[(-Sqrt[\[Eta]] \[Tau]+\[Nu]\[Theta] EllipticF[\[Pi]/2-\[Theta]s,-(a^2/\[Eta])])/\[Nu]\[Theta],-(a^2/\[Eta])]),ArcCos[-\[Nu]\[Theta] Sqrt[up]JacobiSN[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um]]](*eqn 38*);
ordinarySphericalFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[ordinarySphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinarySphericalFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinarySphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinarySphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
ordinarySEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[CapitalDelta]\[Theta],up,um,Gcurl\[Theta],\[Theta]o,assoc,roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,consts,\[Nu]\[Theta]},
assoc=Association["Trajectory"-> \[Theta]o];
\[Theta]o[\[Tau]_]:=\[Theta]s;
ordinarySEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[ordinarySEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="ordinarySEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
ordinarySEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
ordinarySEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
vorticalSpherical[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]v,Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,h,assoc,\[Theta]roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Theta]p,\[Theta]m,roots,\[CapitalUpsilon],consts,\[Nu]\[Theta]},
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
h=Sign[Cos[\[Theta]s]](*eqn 54*);
\[CapitalUpsilon]=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-(h/Sqrt[um a^2])EllipticF[\[CapitalUpsilon],1-up/um](*eqn 56*);
assoc=Association["Trajectory"-> \[Theta]v,"ConstantsofMotion"-> consts];
\[Theta]v[\[Tau]_]:=ArcCos[Sqrt[um]h JacobiDN[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um]](*eqn 64*);
vorticalSphericalFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[vorticalSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="vorticalSphericalFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
vorticalSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
vorticalSphericalFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Polar Spherical motion general code*)


(* ::Input::Initialization:: *)
PolarSphericalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{\[Theta]roots,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,R,M=1},
R=(rs^2+a^2-a \[Lambda])^2-(rs^2-2M rs+a^2)(\[Eta]+(\[Lambda]-a)^2);
If[a==0,\[Theta]s,Which[\[Eta]>0,Return[ordinarySpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinarySEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]<0,Return[vorticalSpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]
]


(* ::Section::Initialization:: *)
(*Polar Motion general code*)


(* ::Input::Initialization:: *)
PolarMotionGeneralS[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[CapitalTheta],\[CapitalTheta]1,R,consts,roots,r1,r2,r3,r4,constss,\[Lambda]s,\[Eta]s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
If[(*r3\[Equal]r4*)4 rs (a^2+rs^2-a \[Lambda])-(-2 M+2 rs) (\[Eta]+(-a+\[Lambda])^2)==0\[And](rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)==0,If[a==0,ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Which[\[Eta]>0,Return[ordinarySpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinarySEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]<0,Return[vorticalSpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],If[a>0,Which[\[Eta]\[Element]PositiveReals,Return[ordinary[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinaryEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]\[Element]NegativeReals,Return[vortical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]
]


(* ::Input::Initialization:: *)
PolarMotionGeneralSS[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[CapitalTheta],\[CapitalTheta]1,R,consts,roots,r1,r2,r3,r4,constss,\[Lambda]s,\[Eta]s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
If[r3==r4,If[a==0,ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Which[\[Eta]>0,Return[ordinarySpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinarySEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]<0,Return[vorticalSpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],If[a>0,Which[\[Eta]\[Element]PositiveReals,Return[ordinary[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinaryEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]\[Element]NegativeReals,Return[vortical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]
]


(* ::Input::Initialization:: *)
PolarMotionGeneral[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[CapitalTheta],\[CapitalTheta]1,R,consts,roots,r1,r2,r3,r4,constss,\[Lambda]s,\[Eta]s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
If[(*r3\[Equal]r4*)4 rs (a^2+rs^2-a \[Lambda])-(-2 M+2 rs) (\[Eta]+(-a+\[Lambda])^2)==0\[And](rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)==0,If[a==0,ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Which[\[Eta]>0,Return[ordinary[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinarySEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]<0,Return[vorticalSpherical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],If[a>0,Which[\[Eta]\[Element]PositiveReals,Return[ordinary[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[ordinaryEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]\[Element]NegativeReals,Return[vortical[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]
]


(* ::Input::Initialization:: *)
PolarMotionS[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1},
If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,If[If[a==0,\[Eta]-\[Lambda]^2 Cot[\[Theta]s]^2>=0,If[\[Theta]s==0\[Or] \[Theta]s==\[Pi],\[Eta]+a^2 Cos[\[Theta]s]^2>=0,\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2>=0]],If[\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,\[Theta]PNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],PolarMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["Error:Negative Angular Potential"]],Print["Error:Negative Radial Potential"]]
](*Print produces alot of statements*)


(* ::Input::Initialization:: *)
PolarMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1},
If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,If[If[\[Theta]s==0\[Or] \[Theta]s==\[Pi],\[Eta]+a^2 Cos[\[Theta]s]^2>=0,\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2>=0],If[\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,\[Theta]PNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],PolarMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["Error:Negative Angular Potential"]],Print["Error:Negative Radial Potential"]]
]


(* ::Chapter:: *)
(*Azimuthal Motion*)


(* ::Section:: *)
(*Non spherical *)


(* ::Input:: *)
(*(*\[Phi]\[Eta]p1\[Rule] azimuthal motion for \[Eta] positive and case one, \[Phi]\[Eta]p2\[Rule] azimuthal motion for \[Eta] positive and case two, \[Phi]\[Eta]p3\[Rule] azimuthal motion for \[Eta] positive and case three, \[Phi]\[Eta]m3\[Rule] azimuthal motion for \[Eta] negative and case three, \[Phi]\[Eta]m4\[Rule] azimuthal motion for \[Eta] negative and case four*)*)


(* ::Subsection::Initialization:: *)
(*case one*)


(* ::Input::Initialization:: *)
\[Phi]\[Eta]p1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,Ips,Ipo,Imns,Imno,g,Intps,Intpo,Intms,Intmo,\[Alpha]2p,\[Alpha]2m,\[Alpha]12,x1s,us,uo,X1,k,\[Kappa]s,Irs,\[Psi]o,\[Psi]s,Imn,Ip,I\[Phi],\[Phi]o,assoc,\[Nu]r,roots,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,r1,r2,r3,r4,consts,\[Nu]\[Theta],rm,rp,G\[Phi],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi], um, up,\[CapitalDelta]\[Theta],\[CapitalPi]m,\[CapitalPi]p},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
x1s=Sqrt[((rs-r2)(r3-r1))/((rs-r1)(r3-r2))];
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))];
\[Kappa]s=ArcSin[x1s];
Irs=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[ArcSin[x1s],k];  
X1[\[Tau]_]:=Sqrt[(r4-r2)(r3-r1)]/2 (\[Tau]+\[Nu]r Irs);
\[Psi]s=ArcSin[x1s];
\[Psi]o[\[Tau]_]:= JacobiAmplitude[X1[\[Tau]],k];
\[CapitalPi]m[\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r2-r1)/((rm-r1)(rm-r2)) (EllipticPi[((rm-r1)(r3-r2))/((rm-r2)(r3-r1)),\[Psi]o[\[Tau]],k]-\[Nu]r EllipticPi[((rm-r1)(r3-r2))/((rm-r2)(r3-r1)),ArcSin[x1s],k]); 
\[CapitalPi]p[\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r2-r1)/((rp-r1)(rp-r2)) (EllipticPi[((rp-r1)(r3-r2))/((rp-r2)(r3-r1)),\[Psi]o[\[Tau]],k]-\[Nu]r EllipticPi[((rp-r1)(r3-r2))/((rp-r2)(r3-r1)),ArcSin[x1s],k]);
Imn[\[Tau]_]:=-\[Tau]/(rm-r1)-\[CapitalPi]m[\[Tau]](*eqn B30c*); 
Ip[\[Tau]_]:=-\[Tau]/(rp-r1)-\[CapitalPi]p[\[Tau]](*eqn B30c*); 
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots, "AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+\[Lambda] G\[Phi][\[Tau]]](*eqn 11*);
\[Phi]\[Eta]p1Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]\[Eta]p1Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*case two*)


(* ::Input::Initialization:: *)
\[Phi]\[Eta]p2[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*eqn B40*);
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+ \[Nu]r Irs)(*eqn B45*); 
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) (EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)), JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2],k])(*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) (EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2],k] )(*eqn B54*);
Imn[\[Tau]_]:=(-\[CapitalPi]m\[Tau][\[Tau]]-\[Tau]/(rm-r3))(*eqn B50*);   
Ip[\[Tau]_]:=(-\[CapitalPi]p\[Tau][\[Tau]]-\[Tau]/(rp-r3))(*eqn B50*);     
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*); 
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots, "AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[k>0\[And]k<1,If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+\[Lambda] G\[Phi][\[Tau]]],Print["k is not in the range (0,1)"]](*EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]*)(*eqn 11*);
\[Phi]\[Eta]p2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]\[Eta]p2Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*PN Geodesics*)


(* ::Input::Initialization:: *)
schwarzschild\[Phi]PNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,ro,assoc,consts,\[Nu]r,\[CapitalTheta],\[Phi]o},
\[Nu]r=prs;
\[Phi]o[\[Tau]_]:=0;
assoc=Association["Trajectory"-> \[Phi]o,"ConstantsofMotion"-> consts];
schwarzschild\[Phi]PNGFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[schwarzschild\[Phi]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="schwarzschild\[Phi]PNGFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
schwarzschild\[Phi]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
schwarzschild\[Phi]PNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Kerr\[Phi]PNGEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r,ro},
\[Nu]r=prs;
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
rp=M+Sqrt[M-a^2];
rm=M-Sqrt[M-a^2];
r1=I a Cos[\[Theta]s];
r2=-I a Cos[\[Theta]s];
\[Nu]r=Sign[prs];
ro[\[Tau]_]:=-((rs \[Nu]r)/(-\[Nu]r+rs \[Tau]));
Ims=1/(rs rm)-Log[rs]/rm^2+Log[rs-rm]/rm^2;   
Ips=1/(rs rp)-Log[rs]/rp^2+Log[rs-rp]/rp^2;      
Imn[\[Tau]_]:=\[Nu]r(1/(ro[\[Tau]] rm)-Log[ro[\[Tau]]]/rm^2+Log[ro[\[Tau]]-rm]/rm^2-Ims);    
Ip[\[Tau]_]:=\[Nu]r(1/(ro[\[Tau]] rp)-Log[ro[\[Tau]]]/rp^2+Log[ro[\[Tau]]-rp]/rp^2-Ips);   
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=I\[Phi][\[Tau]]+\[Lambda]  \[Tau](*eqn 11*);  
Kerr\[Phi]PNGEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[Kerr\[Phi]PNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="Kerr\[Phi]PNGEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
Kerr\[Phi]PNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
Kerr\[Phi]PNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Kerr\[Phi]PNGnonE[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r,ro},
\[Nu]r=prs;
\[Nu]\[Theta]=p\[Theta]s;
up=I Cos[\[Theta]s]^2; 
 um=I Cos[\[Theta]s]^2; 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
rp=M+Sqrt[M-a^2];
rm=M-Sqrt[M-a^2];
r1=I a Cos[\[Theta]s];
r2=-I a Cos[\[Theta]s];
\[Nu]r=Sign[prs];
ro[\[Tau]_]:=Re[(E^((r1 \[Tau])/\[Nu]r) r1 r2-E^((r2 \[Tau])/\[Nu]r) r1 r2+E^((r2 \[Tau])/\[Nu]r) r1 rs-E^((r1 \[Tau])/\[Nu]r) r2 rs)/(E^((r1 \[Tau])/\[Nu]r) r1-E^((r2 \[Tau])/\[Nu]r) r2-E^((r1 \[Tau])/\[Nu]r) rs+E^((r2 \[Tau])/\[Nu]r) rs)];

Ims=((r2-rm) Log[rs-r1]+(-r1+rm) Log[rs-r2]+(r1-r2) Log[rs-rm])/((r1-r2) (r1-rm) (r2-rm));   
Ips=((r2-rp) Log[rs-r1]+(-r1+rp) Log[rs-r2]+(r1-r2) Log[rs-rp])/((r1-r2) (r1-rp) (r2-rp));      
Imn[\[Tau]_]:=\[Nu]r(((r2-rm) Log[ro[\[Tau]]-r1]+(-r1+rm) Log[ro[\[Tau]]-r2]+(r1-r2) Log[ro[\[Tau]]-rm])/((r1-r2) (r1-rm) (r2-rm))-Ims);    
Ip[\[Tau]_]:=\[Nu]r(((r2-rp) Log[ro[\[Tau]]-r1]+(-r1+rp) Log[ro[\[Tau]]-r2]+(r1-r2) Log[ro[\[Tau]]-rp])/((r1-r2) (r1-rp) (r2-rp))-Ips);   
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[a==0,0,If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+(\[Lambda]  \[Tau])/Sin[\[Theta]s]^2]](*eqn 11*);  
Kerr\[Phi]PNGnonEFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[Kerr\[Phi]PNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="Kerr\[Phi]PNGnonEFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
Kerr\[Phi]PNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
Kerr\[Phi]PNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
AzimuthalPNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{},
If[a==0,schwarzschild\[Phi]PNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],If[\[Theta]s==\[Pi]/2,Kerr\[Phi]PNGEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Kerr\[Phi]PNGnonE[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]
]


(* ::Subsection::Initialization:: *)
(*case three*)


(* ::Input::Initialization:: *)
\[Phi]\[Eta]p3[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,R1p\[Tau],R1m\[Tau],F3\[Tau],Ips,Ims,\[CapitalPsi]\[Tau],f1p\[Tau],f1m\[Tau],F3s,x\[Tau],R1ps,R1ms,f1ps,f1ms,\[CapitalPsi]s,p1p,p1m,X3\[Tau],Irs,\[Alpha]p,\[Alpha]m,j,xs,x3s,k3,A,B,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,a1,b1,sk3,consts,\[Nu]r,\[Nu]\[Theta],kp,us,uo},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2)(*eqn 19*);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
a1=Sqrt[-((r4-r3)^2/4)](*eqn B10*);
b1=(r3+r4)/2(*eqn B10*);
x3s=(1- (B(rs-r2))/(A(rs-r1)))/(1+ (B(rs-r2))/(A (rs-r1)))(*eqn B58*);
j=k3^2;
kp=Sqrt[1-j];
\[Alpha]m=If[rm==r1\[Or]rm==r2,Print["Error:\[Alpha]m=1"],(B(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))](*eqn B66*);
\[Alpha]p=If[rp==r1\[Or]rp==r2,Print["Error:\[Alpha]p=1"],(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))](*eqn B66*); 
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
X3\[Tau][\[Tau]_]:=Sqrt[A B](  \[Tau]+ \[Nu]r Irs)(*eqn B74*);
F3\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs(*eqn B81 correction*); 
\[CapitalPsi]\[Tau][\[Tau]_]:= \[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3](*extracted from eqn B80*);
\[CapitalPsi]s=ArcCos[x3s](*extracted from eqn B70*);
us=EllipticF[\[CapitalPsi]s,j];
uo[\[Tau]_]:=EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j];
p1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
p1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
xs=(\[CapitalPsi]s-\[Pi]/2);
x\[Tau][\[Tau]_]:=(\[CapitalPsi]\[Tau][\[Tau]]-\[Pi]/2);
f1ms=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1ps=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1ms=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1ms+2HeavisideLambda[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1ps=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1ps+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
f1m\[Tau][\[Tau]_]:=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1p\[Tau][\[Tau]_]:=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1ms)(*eqn B70*);
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1ps)(*eqn B70*);
R1m\[Tau][\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1m\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1p\[Tau][\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1p\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1m\[Tau][\[Tau]])-\[Nu]r Ims (*eqn B80*);
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1p\[Tau][\[Tau]])-\[Nu]r Ips (*eqn B80*); 
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+\[Lambda] G\[Phi][\[Tau]]](*eqn 11*);
\[Phi]\[Eta]p3Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]\[Eta]p3Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Phi]\[Eta]m3[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,R1p\[Tau],R1m\[Tau],Ips,Ims,f1m\[Tau],f1p\[Tau],R1ms,R1ps,f1ms,f1ps,x\[Tau],xs,p1p,p1m,\[CapitalPsi]\[Tau],\[CapitalPsi]s,F3\[Tau],X3\[Tau],Irs,F3s,\[Alpha]p,\[Alpha]m,j,x3s,k3,A,B,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
\[Theta]p=ArcCos[Sqrt[um]];
\[Theta]m=ArcCos[Sqrt[up]];
h=Sign[Cos[\[Theta]s]](*eqn 54*);
\[CapitalUpsilon]=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-(h/Sqrt[um a^2])EllipticF[\[CapitalUpsilon],1-up/um];(*eqn 56*)
Gcurl\[Phi]=-(h/((1-um)Sqrt[um a^2]))EllipticPi[(up-um)/(1-um),\[CapitalUpsilon],1-up/um](*eqn 57*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
G\[Phi][\[Tau]_]:=Re[1/((1-um)Sqrt[um a^2]) EllipticPi[(up-um)/(1-um),\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurl\[Phi]](*eqn 67*); 
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
x3s=(1- (B(rs-r2))/(A(rs-r1)))/(1+ (B(rs-r2))/(A (rs-r1)))(*eqn B55*);
j=k3;
\[Alpha]m=(B(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))(*eqn B66*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+ \[Nu]r Irs)(*eqn B74*);
F3\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs(*eqn B81 correction*);
\[CapitalPsi]\[Tau][\[Tau]_]:= \[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3];
\[CapitalPsi]s=ArcCos[x3s];
p1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
p1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
xs=(\[CapitalPsi]s-Pi/2);
x\[Tau][\[Tau]_]:=(\[CapitalPsi]\[Tau][\[Tau]]-\[Pi]/2);
f1ms=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1ps=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1ms=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1ms+2HeavisideLambda[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1ps=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1ps+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
f1m\[Tau][\[Tau]_]:=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1p\[Tau][\[Tau]_]:=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1ms)(*eqn B70*);
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1ps)(*eqn B70*);
R1m\[Tau][\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1m\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1p\[Tau][\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1p\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1m\[Tau][\[Tau]])-\[Nu]r Ims (*eqn B80*);
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1p\[Tau][\[Tau]])-\[Nu]r Ips (*eqn B80*);
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+\[Lambda] G\[Phi][\[Tau]]](*eqn 11*);
\[Phi]\[Eta]m3Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]\[Eta]m3Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*case four*)


(* ::Input::Initialization:: *)
\[Phi]\[Eta]m4[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,S1\[Tau]p,S1\[Tau]m,\[CapitalPsi]\[Tau],j,Ips,Ims,f2\[Tau]p,f2\[Tau]m,S1sp,S1sm,\[CapitalPsi]s,F4\[Tau],f2sp,f2sm,X4\[Tau],Irs,x4p,x4m,P2p,P2m,x\[Tau],a2,b2,x4s,xs,F4s,gp,gm,rp,rm,C,D,g0,k4,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r,a1,b1,C2,D2,k},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
\[Theta]p=ArcCos[Sqrt[um]];
\[Theta]m=ArcCos[Sqrt[up]];
h=Sign[Cos[\[Theta]s]](*eqn 54*);
\[CapitalUpsilon]=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=Re[-(h/Sqrt[um a^2])EllipticF[\[CapitalUpsilon],1-up/um]](*eqn 56*);
Gcurl\[Phi]=Re[-(h/((1-um)Sqrt[um a^2]))EllipticPi[(up-um)/(1-um),\[CapitalUpsilon],1-up/um]](*eqn 57*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
G\[Phi][\[Tau]_]:=1/((1-um)Sqrt[um a^2]) EllipticPi[(up-um)/(1-um),\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 67*); 
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*);
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b2=Re[(r1+r2)/2](*eqn B11*);
a1=Re[Sqrt[-((r4-r3)^2/4)] ](*eqn B10*);
b1=Re[(r3+r4)/2 ](*eqn B10*);
C2=(a1-a2)^2+(b1-b2)^2;
D2=(a1+a2)^2+(b1-b2)^2;
k=D2/C2;
x4m=(rm-b2)/a2(*eqn B83*);
x4p=(rp-b2)/a2(*eqn B83*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2-4 a2^2)](*eqn B88*);
gp=(g0 x4p-1)/(g0+x4p)(*eqn B96*);
gm=(g0 x4m-1)/(g0+x4m)(*eqn B96*);
k4=Sqrt[(4 C D)/(C+D)^2](*eqn B87*);
j=k4(*eqn B92*);
P2m=Sqrt[(1+gm^2)/(1-j+gm^2)](*eqn B92*);
P2p=Sqrt[(1+gp^2)/(1-j+gp^2)](*eqn B92*); 
x4s=(rs-b2)/a2(*eqn B83*);
\[CapitalPsi]s=ArcTan[x4s]+ArcTan[g0](*eqn 100*);
f2sm=P2m/2 Log[Abs[(1-P2m)/(1+P2m) (1+P2m Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1-P2m Sqrt[1-j Sin[\[CapitalPsi]s]^2])]](*eqn B95*);
f2sp=P2p/2 Log[Abs[(1-P2p)/(1+P2p) (1+P2p Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1-P2p Sqrt[1-j Sin[\[CapitalPsi]s]^2])]](*eqn B95*);
xs=\[CapitalPsi]s-Pi/2(*eqn B92*);
S1sm=1/(1+gm^2) (EllipticF[\[CapitalPsi]s,j]+gm^2 EllipticPi[1+gm^2,\[CapitalPsi]s,j]-gm f2sm+2gm^2 HeavisideLambda[xs]EllipticPi[1+gm^2,j])(*eqn B92*);
S1sp=1/(1+gp^2) (EllipticF[\[CapitalPsi]s,j]+gp^2 EllipticPi[1+gp^2,\[CapitalPsi]s,j]-gp f2sp+2gp^2 HeavisideLambda[xs]EllipticPi[1+gp^2,j])(*eqn B92*);
F4s=2/(C+D) EllipticF[ArcTan[x4s]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
X4\[Tau][\[Tau]_]:=(C+D)/2 (\[Nu]r \[Tau]+Irs)(*eqn B104*);
\[CapitalPsi]\[Tau][\[Tau]_]:=JacobiAmplitude[X4\[Tau][\[Tau]],k4](*eqn B92*);
f2\[Tau]m[\[Tau]_]:=P2m/2 Log[Abs[(1-P2m)/(1+P2m) (1+P2m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1-P2m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])]](*eqn B95*);
f2\[Tau]p[\[Tau]_]:=P2p/2 Log[Abs[(1-P2p)/(1+P2p) (1+P2p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1-P2p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])]](*eqn B95*);
Ims=g0/(a2(1-g0 x4m)) (F4s-2/(C+D) ((1+g0^2)/(g0(g0+x4m)))S1sm)(*eqn B100*);
Ips=g0/(a2(1-g0 x4p)) (F4s-2/(C+D) ((1+g0^2)/(g0(g0+x4p)))S1sp)(*eqn B100*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-\[Pi]/2(*eqn B92*);
S1\[Tau]m[\[Tau]_]:=1/(1+gm^2) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+gm^2 EllipticPi[1+gm^2,\[CapitalPsi]\[Tau][\[Tau]],j]-gm f2\[Tau]m[\[Tau]]+2gm^2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[1+gm^2,j])(*eqn B92*);
S1\[Tau]p[\[Tau]_]:=1/(1+gp^2) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+gp^2 EllipticPi[1+gp^2,\[CapitalPsi]\[Tau][\[Tau]],j]-gp f2\[Tau]p[\[Tau]]+2gp^2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[1+gp^2,j])(*eqn B92*);
F4\[Tau][\[Tau]_]:=2 /(C+D) EllipticF[JacobiAmplitude[X4\[Tau][\[Tau]],k4],k4](*eqn B115*);
Imn[\[Tau]_]:=(\[Nu]r g0)/(a2(1-g0 x4m)) (F4\[Tau][\[Tau]]-2/(C+D) ((1+g0^2)/(g0(g0+x4m)))S1\[Tau]m[\[Tau]])-\[Nu]r Ims (*eqn B114*);
Ip[\[Tau]_]:=(\[Nu]r g0)/(a2(1-g0 x4p)) (F4\[Tau][\[Tau]]-2/(C+D) ((1+g0^2)/(g0(g0+x4p)))S1\[Tau]p[\[Tau]])-\[Nu]r Ips (*eqn B114*);
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o];
\[Phi]o[\[Tau]_]:=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],I\[Phi][\[Tau]],I\[Phi][\[Tau]]+\[Lambda] G\[Phi][\[Tau]]](*eqn 11*);
\[Phi]\[Eta]m4Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]\[Eta]m4Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*\[Phi] equatorial *)


(* ::Subsubsection::Initialization:: *)
(*\[Phi]1*)


(* ::Input::Initialization:: *)
\[Phi]1Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,\[CapitalPi]p\[Tau],F1\[Tau],\[CapitalPi]m\[Tau],Ims,X1\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F1s,Irs,x1s,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r,aa,b,c,d,g,\[Alpha]2p,\[Alpha]2m,\[Alpha]12,\[Kappa]s,X1,\[Psi]s,\[Psi]o,us,uo,Imno,Imns,Ipo,Ips},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
aa=r4;
b=r3;
c=r2;
d=r1;
g=2/Sqrt[(r4-r2)(r3-r1)];
\[Alpha]2p=((rp-r4)(r3-r2))/((rp-r3)(r4-r2));
\[Alpha]2m=((rm-r4)(r3-r2))/((rm-r3)(r4-r2));
\[Alpha]12=(r3-r2)/(r4-r2);
x1s=Sqrt[((r4-r2)(r3-rs))/((r3-r2)(r4-rs))];
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))];
\[Kappa]s=ArcSin[x1s];
Irs=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[\[Kappa]s,k];  
X1[\[Tau]_]:=Sqrt[(r4-r2)(r3-r1)]/2 (\[Tau]+\[Nu]r Irs);
\[Psi]s=ArcSin[x1s];
\[Psi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X1[\[Tau]],k];
us=EllipticF[\[Psi]s,k];
uo[\[Tau]_]:=EllipticF[\[Psi]o[\[Tau]],k];
Imns=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]s,k]+(rm-b)EllipticF[\[Psi]s,k]);
Imno[\[Tau]_]:=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]o[\[Tau]],k]+(rm-b)EllipticF[\[Psi]o[\[Tau]],k]);
Ips=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]s,k]+(rp-b)EllipticF[\[Psi]s,k]);
Ipo[\[Tau]_]:=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]o[\[Tau]],k]+(rp-b)EllipticF[\[Psi]o[\[Tau]],k]);
Imn[\[Tau]_]:=-\[Nu]r(Imno[\[Tau]]-Imns)(*eqn B30c*);
Ip[\[Tau]_]:=-\[Nu]r(Ipo[\[Tau]]-Ips)(*eqn B30c*); 
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots, "AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=I\[Phi][\[Tau]]+\[Lambda] \[Tau](*eqn 11*);
\[Phi]1EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]1EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization:: *)
(*\[Phi]2*)


(* ::Input::Initialization:: *)
\[Phi]2Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2],k](*eqn B40*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2],k](*eqn B43*);
Irs=F2s; 
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Nu]r \[Tau]+ Irs)(*eqn B45*); 
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*); 
F2\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B51*);
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*); 
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*);   
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);   
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);    
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots, "AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=I\[Phi][\[Tau]]+\[Lambda] \[Tau](*eqn 11*); 
\[Phi]2EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]2EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization:: *)
(*\[Phi]3*)


(* ::Input::Initialization:: *)
\[Phi]3Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,I\[Phi],Ip,Imn,R1p\[Tau],R1m\[Tau],Ips,Ims,f1m\[Tau],f1p\[Tau],R1ms,R1ps,f1ms,f1ps,x\[Tau],xs,p1p,p1m,\[CapitalPsi]\[Tau],\[CapitalPsi]s,F3\[Tau],X3\[Tau],Irs,F3s,\[Alpha]p,\[Alpha]m,j,x3s,k3,A,B,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=p\[Theta]s;
\[Nu]r=prs;
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
x3s=(A(rs-r1)-B (rs-r2))/(A (rs-r1)+B (rs-r2))(*eqn B55*);
j=k3^2;
\[Alpha]m=(B(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))(*eqn B66*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+\[Nu]r Irs)(*eqn B74*);
F3\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B81*);
\[CapitalPsi]\[Tau][\[Tau]_]:=\[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3];
\[CapitalPsi]s=ArcCos[x3s];
p1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
p1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
xs=(\[CapitalPsi]s-Pi/2);
x\[Tau][\[Tau]_]:=(\[CapitalPsi]\[Tau][\[Tau]]-Pi/2);
f1ms=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1ps=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(p1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1ms=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1ms+2HeavisideLambda[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1ps=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1ps+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
f1m\[Tau][\[Tau]_]:=p1m/2 Log[Abs[(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1p\[Tau][\[Tau]_]:=p1p/2 Log[Abs[(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(p1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1ms)(*eqn B70*);
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1ps)(*eqn B70*);
R1m\[Tau][\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1m\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1p\[Tau][\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1p\[Tau][\[Tau]]+2HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1m\[Tau][\[Tau]])-\[Nu]r Ims (*eqn B80*);
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1p\[Tau][\[Tau]])-\[Nu]r Ips (*eqn B80*);
I\[Phi][\[Tau]_]:=(2 M a)/(rp-rm) ((rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])(*eqn B2*);
assoc=Association["Trajectory"-> \[Phi]o,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=I\[Phi][\[Tau]]+\[Lambda] \[Tau] (*eqn 11*);
\[Phi]3EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]3EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Schwarzchild case*)


(* ::Input::Initialization:: *)
\[Phi]Schwarzchild1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,consts,\[Nu]\[Theta],\[Nu]r,u0,\[Theta]o,\[Phi]\[Tau],n,\[Phi]s,u\[Tau],us},
\[Nu]\[Theta]=p\[Theta]s;
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
n=u0;
\[Phi]\[Tau][\[Tau]_]:=(*JacobiAmplitude[\[Tau]Sqrt[\[Eta]+\[Lambda]^2]-\[Nu]\[Theta] ArcSin[Cos[\[Theta]s]/Sqrt[u0]],0]*)\[Tau] Sqrt[\[Eta]+\[Lambda]^2]-\[Nu]\[Theta] ArcSin[Cos[\[Theta]s]/Sqrt[u0]];
\[Phi]s=ArcSin[Cos[\[Theta]s]/Sqrt[u0]];
(*\[Theta]o[\[Tau]_]:=ordinarySchwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];*)
us=Cos[\[Theta]s];
G\[Phi][\[Tau]_]:=(*Sqrt[u0/\[Eta]]EllipticPi[u0,\[Phi]\[Tau][\[Tau]],0]+\[Nu]\[Theta]Sqrt[u0/\[Eta]]EllipticPi[u0,\[Phi]s,0]*)Sqrt[u0/\[Eta]](ArcTanh[Sqrt[-1+u0] Tan[\[Phi]\[Tau][\[Tau]]]]/Sqrt[-1+u0]+(\[Pi] Round[Re[\[Phi]\[Tau][\[Tau]]]/\[Pi]])/Sqrt[1-u0])+\[Nu]\[Theta] Sqrt[u0/\[Eta]](ArcTanh[Sqrt[-1+u0] Tan[\[Phi]s]]/Sqrt[-1+u0]+(\[Pi] Round[Re[\[Phi]s]/\[Pi]])/Sqrt[1-u0])(*Eqn 26b*);    
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],0,(*If[\[Nu]\[Theta]==0,\[Lambda]/Sin[\[Theta]s]^2\[Tau],\[Lambda] G\[Phi][\[Tau]]]*)\[Lambda] G\[Phi][\[Tau]]](*eqn 11*);
\[Phi]Schwarzchild1Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
X[rs_]:=rs^4-(rs^2-2 rs)(\[Eta]+\[Lambda]^2)


(* ::Input::Initialization:: *)
Format[\[Phi]Schwarzchild1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]Schwarzchild1Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]Schwarzchild1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]Schwarzchild1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Phi]Schwarzchild2[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,consts,\[Nu]\[Theta],\[Nu]r,u0,\[Theta]o,\[Phi]\[Tau],n,\[Phi]s,u\[Tau],us},
\[Nu]\[Theta]=p\[Theta]s;    
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=\[Lambda] \[Tau]/Sin[\[Theta]s]^2(*eqn 11*);
\[Phi]Schwarzchild2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]  


(* ::Input::Initialization:: *)
Format[\[Phi]Schwarzchild2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]Schwarzchild2Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]Schwarzchild2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]Schwarzchild2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Phi]Schwarzchild[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{u0,\[Theta]1,\[Theta]2},
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
\[Theta]1=ArcCos[Sqrt[u0]];
\[Theta]2=ArcCos[-Sqrt[u0]];
If[\[Theta]1==\[Theta]2,\[Phi]Schwarzchild2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Phi]Schwarzchild1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]


(* ::Subsection::Initialization:: *)
(*Azimuthal Motion non spherical motion general code*)


(* ::Input::Initialization:: *)
AzimuthalNonSphericalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,rm,roots,r1,r2,r3,r4,consts,\[Nu]r,R},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
R=(rs^2+a^2-a \[Lambda])^2-(rs^2-2M rs+a^2)(\[Eta]+(\[Lambda]-a)^2);
If[R>=0\[And]If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],(\[Eta]+a^2 Cos[\[Theta]s]^2)>=0,(\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2)>=0],If[a>0,Which[\[Eta]\[Element]PositiveReals,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>r4,Return[\[Phi]\[Eta]p2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<r3,If[\[Nu]r<0,Return[\[Phi]\[Eta]p1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]](*,Return[\[Phi]\[Eta]p145[a,ts,rs,\[Theta]s,\[Phi]s,prs,p\[Theta]s,p\[Phi]s]]*)]]],Return[\[Phi]\[Eta]p3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]==0,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>r4,Return[\[Phi]2Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<r3,If[\[Nu]r<0,Return[\[Phi]1Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Return[\[Phi]145Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]],Return[\[Phi]3Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]\[Element]NegativeReals,If[r2 \[Element] Reals,
If[r4 \[Element]Complexes,Return[\[Phi]\[Eta]m3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[\[Phi]\[Eta]m4[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],"Error:Radial or Angular Potential is Negative"]
]


(* ::Input::Initialization:: *)



(* ::Section::Initialization:: *)
(*Azimuthal Spherical motion*)


(* ::Input::Initialization:: *)
\[Phi]Spherical\[Eta]p[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,consts,\[Nu]\[Theta],\[Nu]r,zamo},
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2);
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
zamo=M+2 Sqrt[M^2-a^2/3] Cos[1/3 ArcCos[(M (M^2-a^2))/(M^2-a^2/3)^(3/2)]];
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[(*\[Lambda]\[Equal]0*)rs==zamo\[Or]\[Theta]s==0\[Or]\[Theta]s==\[Pi],(a/(rs^2-2 M rs+a^2) (rs^2+a^2-a \[Lambda])-a)\[Tau],((a/(rs^2-2 M rs+a^2) (rs^2+a^2-a \[Lambda])-a)\[Tau]+\[Lambda] G\[Phi][\[Tau]])](*eqn 11*);
\[Phi]Spherical\[Eta]pFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]Spherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]Spherical\[Eta]pFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]Spherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]Spherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Phi]Spherical\[Eta]m[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,consts,\[Nu]\[Theta],\[Nu]r,zamo},
\[Nu]\[Theta]=p\[Theta]s;
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
\[Theta]p=ArcCos[Sqrt[um]];
\[Theta]m=ArcCos[Sqrt[up]];
h=Sign[Cos[\[Theta]s]](*eqn 54*);
\[CapitalUpsilon]=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-(h/Sqrt[um a^2])EllipticF[\[CapitalUpsilon],1-up/um];(*eqn 56*)
Gcurl\[Phi]=-(h/((1-um)Sqrt[um a^2]))EllipticPi[(up-um)/(1-um),\[CapitalUpsilon],1-up/um](*eqn 57*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
G\[Phi][\[Tau]_]:=1/((1-um)Sqrt[um a^2]) EllipticPi[(up-um)/(1-um),\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 67*); 
zamo=M+2 Sqrt[M^2-a^2/3] Cos[1/3 ArcCos[(M (M^2-a^2))/(M^2-a^2/3)^(3/2)]];
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[rs==zamo\[Or]\[Theta]s==0\[Or]\[Theta]s==\[Pi],(a/(rs^2-2 M rs+a^2) (rs^2+a^2-a \[Lambda])-a)\[Tau],((a/(rs^2-2 M rs+a^2) (rs^2+a^2-a \[Lambda])-a)\[Tau]+\[Lambda] G\[Phi][\[Tau]])](*eqn 11*);
\[Phi]Spherical\[Eta]mFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]Spherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]Spherical\[Eta]mFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]Spherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]Spherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
\[Phi]SEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,consts,\[Nu]\[Theta],\[Nu]r},
\[Nu]\[Theta]=p\[Theta]s;
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=(a/(rs^2-2 M rs+a^2) (rs^2+a^2-a \[Lambda])-a)\[Tau]+\[Lambda] \[Tau](*eqn 11*);
\[Phi]SEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[\[Phi]SEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]SEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]SEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]SEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
Clear[\[Phi]SSph]


(* ::Input::Initialization:: *)
\[Phi]SSph[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,\[Phi]o,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[CapitalUpsilon]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,\[CapitalUpsilon],\[Theta]p,\[Theta]m,h,assoc,consts,\[Nu]\[Theta],\[Nu]r,u0,\[Theta]o},
\[Nu]\[Theta]=p\[Theta]s;
u0=\[Eta]/(\[Eta]+\[Lambda]^2);
assoc=Association["Trajectory"->\[Phi]o,"ConstantsofMotion"-> consts];
\[Phi]o[\[Tau]_]:=If[\[Lambda]==0,0,\[Lambda]/(*Sin[\[Theta]s]^2*)1 \[Tau]](*eqn 11*);
\[Phi]SSphFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]   


(* ::Input::Initialization:: *)
Format[\[Phi]SSphFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="\[Phi]SSphFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
\[Phi]SSphFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
\[Phi]SSphFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Azimuthal Spherical motion general code*)


(* ::Input::Initialization:: *)
AzimuthalSphericalMotion[a_,rs_,\[Theta]s_,prs_,p\[Theta]s_,p\[Phi]s_]:=Module[{M=1,rp,rm,roots,r1,r2,r3,r4,consts,\[Lambda],\[Eta],R},
R=(rs^2+a^2-a \[Lambda])^2-(rs^2-2M rs+a^2)(\[Eta]+(\[Lambda]-a)^2);
If[R==0\[And]If[\[Theta]s==0\[Or]\[Theta]s==\[Pi],(\[Eta]+a^2 Cos[\[Theta]s]^2)>=0,(\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2)>=0],Which[\[Eta]\[Element]PositiveReals,\[Phi]Spherical\[Eta]p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]==0,\[Phi]SEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]\[Element]NegativeReals,\[Phi]Spherical\[Eta]m[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],"Error:Radial Potential \[NotEqual] 0 or Angular Potential is negative"]
]


(* ::Section::Initialization:: *)
(*Azimuthal Motion General Code*)


(* ::Input::Initialization:: *)
AzimuthalMotionGeneralS[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,roots,r1,r2,r3,r4,consts,\[Nu]r,constss,\[Lambda]s,\[Eta]s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
If[(*r3\[Equal]r4==rs*)4 rs (a^2+rs^2-a \[Lambda])-(-2 M+2 rs) (\[Eta]+(-a+\[Lambda])^2)==0\[And](rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)==0,If[a==0,\[Phi]SSph[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Which[\[Eta]\[Element]PositiveReals,\[Phi]Spherical\[Eta]p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]==0,\[Phi]SEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]\[Element]NegativeReals,\[Phi]Spherical\[Eta]m[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],If[a>0,Which[\[Eta]\[Element]PositiveReals,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[\[Phi]\[Eta]p2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[\[Phi]\[Eta]p1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]\[Eta]p3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]==0,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[\[Phi]2Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[\[Phi]1Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]3Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]\[Element]NegativeReals,If[r2 \[Element] Reals,
If[r4 \[Element]Complexes,Return[\[Phi]\[Eta]m3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[\[Phi]\[Eta]m4[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]]


(* ::Input::Initialization:: *)
AzimuthalMotionGeneral[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,roots,r1,r2,r3,r4,consts,\[Nu]r,constss,\[Lambda]s,\[Eta]s,k},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=prs;
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2));
If[r3==r4==rs(*4 rs (a^2+rs^2-a \[Lambda])-(-2 M+2 rs) (\[Eta]+(-a+\[Lambda])^2)==0\[And](rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)==0*),If[a==0,\[Phi]Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],Which[\[Eta]\[Element]PositiveReals,\[Phi]Spherical\[Eta]p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]==0,\[Phi]SEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],\[Eta]\[Element]NegativeReals,\[Phi]Spherical\[Eta]m[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],If[a>0,Which[\[Eta]\[Element]PositiveReals,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>r4,If[k<1\[And]k>0,Return[\[Phi]\[Eta]p2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],"These values for the energy and angular momentum correspond to a critical orbit but the initial position does not. Please enter instead values for the energy and angular momentum which are very close -but not exactly equal- to the critical ones."],If[rs<=r3,Return[\[Phi]\[Eta]p1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]\[Eta]p3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]==0,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[\[Phi]2Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[\[Phi]1Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]3Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]\[Element]NegativeReals,If[r2 \[Element] Reals,
If[r4 \[Element]Complexes,Return[\[Phi]\[Eta]m3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[\[Phi]\[Eta]m4[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[\[Phi]Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]]


(* ::Input::Initialization:: *)
AzimuthalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{consts,M=1},
If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,If[If[\[Theta]s==0\[Or] \[Theta]s==\[Pi],\[Eta]+a^2 Cos[\[Theta]s]^2>=0,\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2>=0],If[(*\[Eta]+(\[Lambda]-a)^2\[Equal]0*)\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,AzimuthalPNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],AzimuthalMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["Error:Negative Angular Potential"]],Print["Error:Negative Radial Potential"]]
]


(* ::Input::Initialization:: *)
AzimuthalMotions[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{consts,M=1},
If[(*\[Eta]+(\[Lambda]-a)^2\[Equal]0*)\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,AzimuthalPNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],AzimuthalMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]
]


(* ::Chapter:: *)
(*Temporal Motion*)


(* ::Subsection::Initialization:: *)
(*case one*)


(* ::Input::Initialization:: *)
t\[Eta]p1[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,Ips,Ipo,Imns,Imno,g,Intps,Intpo,Intms,Intmo,\[Alpha]2p,\[Alpha]2m,\[Alpha]12,x1s,us,uo,X1,k,\[Kappa]s,Irs,\[Psi]o,\[Psi]s,Imn,Ip,It,to1,assoc,\[Nu]r,roots,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,r1,r2,r3,r4,consts,\[Nu]\[Theta],rm,rp,G\[Phi],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi], um, up,\[CapitalDelta]\[Theta],Gt,Ep1,Ep, Gcurlt,I0,aa,b,c,d,\[Alpha]1,\[Alpha],I1o,I1s,I1,I2,I2o,I2s,V1o,V1s,V2o,V2s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])(*eqn 32*);
Gcurlt=(2 up)/Sqrt[-um a^2] Ep1 (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-Gcurlt ](*eqn 48*);
aa=r4;
b=r3;
c=r2;
d=r1;
\[Alpha]1=Sqrt[(aa(b-c))/(b(aa-c))];
\[Alpha]=Sqrt[(b-c)/(aa-c)];
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
g=2/Sqrt[(r4-r2)(r3-r1)];
\[Alpha]2p=((rp-r4)(r3-r2))/((rp-r3)(r4-r2));
\[Alpha]2m=((rm-r4)(r3-r2))/((rm-r3)(r4-r2));
\[Alpha]12=(r3-r2)/(r4-r2);
x1s=Sqrt[((r4-r2)(r3-rs))/((r3-r2)(r4-rs))];
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))];
\[Kappa]s=ArcSin[x1s];
Irs=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[\[Kappa]s,k];  
X1[\[Tau]_]:=Sqrt[(r4-r2)(r3-r1)]/2 (\[Tau]+\[Nu]r Irs);
\[Psi]s=ArcSin[x1s];
\[Psi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X1[\[Tau]],k];
us=EllipticF[\[Psi]s,k];
uo[\[Tau]_]:=EllipticF[\[Psi]o[\[Tau]],k];
Imns=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]s,k]+(rm-b)EllipticF[\[Psi]s,k]);
Imno[\[Tau]_]:=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]o[\[Tau]],k]+(rm-b)EllipticF[\[Psi]o[\[Tau]],k]);
Ips=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]s,k]+(rp-b)EllipticF[\[Psi]s,k]);
Ipo[\[Tau]_]:=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]o[\[Tau]],k]+(rp-b)EllipticF[\[Psi]o[\[Tau]],k]);
Imn[\[Tau]_]:=-\[Nu]r(Imno[\[Tau]]-Imns)(*eqn B30c*);
Ip[\[Tau]_]:=-\[Nu]r(Ipo[\[Tau]]-Ips)(*eqn B30c*);
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[us,k]+(k^2-\[Alpha]^2)us+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]-(\[Alpha]^4 JacobiSN[us,k]JacobiCN[us,k]JacobiDN[us,k])/(1-\[Alpha]^2 JacobiSN[us,k]^2)); 
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[uo[\[Tau]],k]+(k^2-\[Alpha]^2)uo[\[Tau]]+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]-(\[Alpha]^4 JacobiSN[uo[\[Tau]],k]JacobiCN[uo[\[Tau]],k]JacobiDN[uo[\[Tau]],k])/(1-\[Alpha]^2 JacobiSN[uo[\[Tau]],k]^2));
V1s=EllipticPi[\[Alpha]^2,\[Psi]s,k];
V1o[\[Tau]_]:=EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k];
I2s=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 us+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2o[\[Tau]_]:=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 uo[\[Tau]]+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s);
I1s=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]+\[Alpha]1^2 us);
I1o[\[Tau]_]:=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]+\[Alpha]1^2 uo[\[Tau]]);
I1[\[Tau]_]:=\[Nu]r (I1o[\[Tau]]-I1s);
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to1,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to1[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]p1Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t\[Eta]p1Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t\[Eta]p1Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*case two*)


(* ::Input::Initialization:: *)
t\[Eta]p2x[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],E2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,I2o,V1s,V1o,V2s,V2o,\[CapitalPsi]o,\[CapitalPsi]s,\[Alpha],\[Alpha]1,aa,b,c,d,y,uo,us,g,u1,u2,x2ss,x2o,ro,V14,V24,\[CapitalPsi]4,x24,I24,I2s4,I24o},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=Re[1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])](*eqn 32*);
Gcurlt=Re[(2 up)/Sqrt[-um a^2] Ep1] (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-\[Nu]\[Theta] Gcurlt] (*eqn 48*);
k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]EllipticE[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B52*);
F2\[Tau][\[Tau]_]:=(*\[Nu]r \[Tau]+Irs*)2/Sqrt[(r3-r1)(r4-r2)] EllipticF[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k]; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
aa=r4;
b=r3;
c=r2;
d=r1;
y=rs;
\[Alpha]1=Sqrt[((aa-d)b)/((b-d)aa)];
\[Alpha]=Sqrt[(aa-d)/(b-d)];
g=1/Sqrt[(aa-c)(b-d)];
x2ss=Sqrt[((b-d)(y-aa))/((aa-d)(y-b))];
x24=Sqrt[(r4-r4)/(r4-r3) (r3-r1)/(r4-r1)];  
ro[\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
x2o[\[Tau]_]:=Sqrt[((b-d)(ro[\[Tau]]-aa))/((aa-d)(ro[\[Tau]]-b))];
\[CapitalPsi]s=ArcSin[x2s];
\[CapitalPsi]4=ArcSin[x24];
\[CapitalPsi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k];
V1o[\[Tau]_]:=EllipticPi[\[CapitalPsi]o[\[Tau]],\[Alpha]^2,k];
V1s=EllipticPi[\[CapitalPsi]s,\[Alpha]^2,k];
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]o[\[Tau]],k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1o[\[Tau]]-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiCN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiDN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]^2));    
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]s,k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]s,k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1s-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]s,k],k] JacobiCN[EllipticF[\[CapitalPsi]s,k],k] JacobiDN[EllipticF[\[CapitalPsi]s,k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]s,k],k]^2));    
I2o[\[Tau]_]:=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]o[\[Tau]],k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]); 
I2s=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]s,k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s); 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]p2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
t\[Eta]p2original[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],E2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,H2\[Tau],R,ro\[Tau],rop\[Tau]},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=Re[1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])](*eqn 32*);
Gcurlt=Re[(2 up)/Sqrt[-um a^2] Ep1] (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-\[Nu]\[Theta] Gcurlt] (*eqn 48*);
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]EllipticE[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B52*);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
ro\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
rop\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]'[\[Tau]];
H2\[Tau][\[Tau]_]:=rop\[Tau][\[Tau]]/(ro\[Tau][\[Tau]]-r3)-\[Nu]r Sqrt[Rs]/(rs-r3); 
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r3-r2)](EllipticE[JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticE[ArcSin[x2s],k]);
I2[\[Tau]_]:=H2\[Tau][\[Tau]]-((r1-r4)+(r2-r3))/2-E2\[Tau][\[Tau]]; 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]p2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
t\[Eta]p2[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],E2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,H2\[Tau],R,ro\[Tau],rop\[Tau]},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=Re[1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])](*eqn 32*);
Gcurlt=Re[(2 up)/Sqrt[-um a^2] Ep1] (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-\[Nu]\[Theta] Gcurlt] (*eqn 48*);
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]EllipticE[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B52*);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
ro\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
rop\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]'[\[Tau]];
H2\[Tau][\[Tau]_]:=rop\[Tau][\[Tau]]/(ro\[Tau][\[Tau]]-r3)-\[Nu]r Sqrt[Rs]/(rs-r3); 
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r3-r2)](EllipticE[JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticE[ArcSin[x2s],k]);
I2[\[Tau]_]:=H2\[Tau][\[Tau]]-((r1-r4)+(r2-r3))/2-E2\[Tau][\[Tau]]; 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=-((a^4+2 a^2 rp^2+rp^4-a^3 \[Lambda]-a rp^2 \[Lambda])/(rm-rp))Ip[\[Tau]]-(a^4+2 a^2 rm^2+rm^4-a^3 \[Lambda]-a rm^2 \[Lambda])/(-rm+rp) Imn[\[Tau]]+(a^2+rm^2+rm rp+rp^2)I0[\[Tau]]+(rm+rp)I1[\[Tau]]+I2[\[Tau]];
(*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=(*It[\[Tau]]+a^2Gt[\[Tau]]*)It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]p2Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t\[Eta]p2Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t\[Eta]p2Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*case three*)


(* ::Input::Initialization:: *)
t\[Eta]p3[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_]:=Module[{M=1,to3,It,Ip,Imn,I0,I1,I2,R1\[Tau]p,R1\[Tau]m,x\[Tau],\[CapitalPsi]\[Tau],Ips,Ims,\[CapitalPi]1\[Tau],I1s,\[CapitalPi]2\[Tau],I2s,f1\[Tau]p,f1\[Tau]m,P1p,P1m,j,X3\[Tau],R1sp,R1sm,\[CapitalPsi]s,xs,R10\[Tau],\[CapitalPi]1s,R20\[Tau],f1sp,f1sm,F3\[Tau],Irs,\[CapitalPi]2s,f10\[Tau],P10,R10s,R20s,f10s,F3s,\[Alpha]p,\[Alpha]m,\[Alpha]0,x3s,A,B,k3,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])(*eqn 32*);
Gcurlt=Re[(2 up)/Sqrt[-um a^2] Ep1] (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-\[Nu]\[Theta] Gcurlt ](*eqn 48*);
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
x3s=(A(rs-r1)-B(rs-r2))/(A(rs-r1)+B(rs-r2))(*eqn B55*);
\[Alpha]0=(B+A)/(B-A)(*eqn B58*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*);
\[Alpha]m=(B(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))(*eqn B66*); 
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
j=k3(*eqn B80 and B62*);
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+\[Nu]r Irs)(*eqn B74*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-Pi/2(*eqn B62*);
\[CapitalPsi]\[Tau][\[Tau]_]:=\[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3](*eqn B80 and B62*);
P10=Sqrt[(\[Alpha]0^2-1)/(j+(1-j)\[Alpha]0^2)](*eqn B65*);
f10\[Tau][\[Tau]_]:=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
R10\[Tau][\[Tau]_]:=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0 f10\[Tau][\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20\[Tau][\[Tau]_]:=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]\[Tau][\[Tau]],j]-(\[Alpha]0 Sin[\[CapitalPsi]\[Tau][\[Tau]]]Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]\[Tau][\[Tau]]])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10\[Tau][\[Tau]](*eqn B64*);
\[CapitalPsi]s=ArcCos[x3s](*eqn B70*);
xs=\[CapitalPsi]s-Pi/2(*eqn B62*);
f10s=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R10s=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]s,j]-\[Alpha]0 f10s+2 HeavisideLambda[xs]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20s=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]s,j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]s,j]-(\[Alpha]0 Sin[\[CapitalPsi]s]Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]s])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10s(*eqn B64*);
\[CapitalPi]2s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20s(*eqn B72*);
\[CapitalPi]1s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10s(*eqn B72*);
P1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
P1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
f1sm=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1sp=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1sm=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1sm+2HeavisideLambda[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1sp=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1sp+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
I2s=((B r2+A r1)/(B+A))^2 F3s+2((B r2+A r1)/(B+A))\[CapitalPi]1s+Sqrt[A B]\[CapitalPi]2s(*eqn B69*);
\[CapitalPi]2\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20\[Tau][\[Tau]](*eqn B82*);
I1s=((B r2+A r1)/(B+A))F3s+\[CapitalPi]1s(*eqn B68*);
\[CapitalPi]1\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10\[Tau][\[Tau]](*eqn B82*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1sm)(*eqn B70*); 
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1sp)(*eqn B70*);
f1\[Tau]m[\[Tau]_]:=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1\[Tau]p[\[Tau]_]:=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
F3\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B81*);
R1\[Tau]m[\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1\[Tau]m[\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1\[Tau]p[\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1\[Tau]p[\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
I2[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))^2 F3\[Tau][\[Tau]]+2((B r2+A r1)/(B+A))\[CapitalPi]1\[Tau][\[Tau]]+Sqrt[A B]\[CapitalPi]2\[Tau][\[Tau]]-I2s)(*eqn B79*);
I1[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))F3\[Tau][\[Tau]]+\[CapitalPi]1\[Tau][\[Tau]]-I1s)(*eqn B78*);
I0[\[Tau]_]:=\[Tau];
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1\[Tau]m[\[Tau]])-\[Nu]r Ims(*eqn B80*); 
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1\[Tau]p[\[Tau]])-\[Nu]r Ips(*eqn B80*);
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to3];
to3[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]p3Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t\[Eta]p3Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t\[Eta]p3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
t\[Eta]m3[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to3,It,Ip,Imn,I0,I1,I2,R1\[Tau]p,R1\[Tau]m,x\[Tau],\[CapitalPsi]\[Tau],Ips,Ims,\[CapitalPi]1\[Tau],I1s,\[CapitalPi]2\[Tau],I2s,f1\[Tau]p,f1\[Tau]m,P1p,P1m,j,X3\[Tau],R1sp,R1sm,\[CapitalPsi]s,xs,R10\[Tau],\[CapitalPi]1s,R20\[Tau],f1sp,f1sm,F3\[Tau],Irs,\[CapitalPi]2s,f10\[Tau],P10,R10s,R20s,f10s,F3s,\[Alpha]p,\[Alpha]m,\[Alpha]0,x3s,A,B,k3,Gt,\[CapitalUpsilon]\[Tau],Gcurlt,Gcurl\[Theta],\[CapitalUpsilon]s,\[Theta]p,\[Theta]m,h,\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
\[Theta]p=ArcCos[Sqrt[um]](*eqn 55*);
\[Theta]m=ArcCos[Sqrt[up]](*eqn 55*);
h=Sign[Cos[\[Theta]s]];
\[CapitalUpsilon]s=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=Re[-h/Sqrt[um a^2] EllipticF[\[CapitalUpsilon]s,1-up/um]](*eqn 56*);
Gcurlt=Re[-h Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]s,1-up/um]](*eqn 58*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
Gt[\[Tau]_]:=Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurlt (*eqn 68*);
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
x3s=(A(rs-r1)-B(rs-r2))/(A(rs-r1)+B(rs-r2))(*eqn B55*);
\[Alpha]0=(B+A)/(B-A)(*eqn B58*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*); 
\[Alpha]m=(A(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))(*eqn B66*);
k3=Sqrt[((A+B)^2-(r2-r1)^2)/(4 A B)](*eqn B59*);
j=k3(*eqn B80 and B62*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-Pi/2(*eqn B62*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
X3\[Tau][\[Tau]_]:=Sqrt[A B](\[Tau]+\[Nu]r  Irs)(*eqn B74*);
\[CapitalPsi]\[Tau][\[Tau]_]:=\[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3](*eqn B80 and B62*);
P10=Sqrt[(\[Alpha]0^2-1)/(j+(1-j)\[Alpha]0^2)](*eqn B65*);
f10\[Tau][\[Tau]_]:=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
R10\[Tau][\[Tau]_]:=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0 f10\[Tau][\[Tau]]+2 HeavisidePi[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20\[Tau][\[Tau]_]:=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]\[Tau][\[Tau]],j]-(\[Alpha]0 Sin[\[CapitalPsi]\[Tau][\[Tau]]]Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]\[Tau][\[Tau]]])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10\[Tau][\[Tau]](*eqn B64*);
\[CapitalPsi]s=ArcCos[x3s](*eqn B70*);
xs=\[CapitalPsi]s-Pi/2(*eqn B62*);
f10s=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R10s=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]s,j]-\[Alpha]0 f10s+2 HeavisidePi[xs]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20s=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]s,j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]s,j]-(\[Alpha]0 Sin[\[CapitalPsi]s]Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]s])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10s(*eqn B64*);
\[CapitalPi]2s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20s(*eqn B72*);
\[CapitalPi]1s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10s(*eqn B72*);
P1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
P1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
f1sm=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1sp=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1sm=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1sm+2HeavisidePi[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1sp=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1sp+2HeavisidePi[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
I2s=((B r2+A r1)/(B+A))^2 F3s+2((B r2+A r1)/(B+A))\[CapitalPi]1s+Sqrt[A B]\[CapitalPi]2s(*eqn B69*);
\[CapitalPi]2\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20\[Tau][\[Tau]](*eqn B82*);
I1s=((B r2+A r1)/(B+A))F3s+\[CapitalPi]1s(*eqn B68*);
\[CapitalPi]1\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10\[Tau][\[Tau]](*eqn B82*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r1)-A(rm-r1)) R1sm)(*eqn B70*);
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r1)-A(rp-r1)) R1sp)(*eqn B70*);
f1\[Tau]m[\[Tau]_]:=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1\[Tau]p[\[Tau]_]:=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
F3\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B81*);
R1\[Tau]m[\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1\[Tau]m[\[Tau]]+2 HeavisideTheta[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1\[Tau]p[\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1\[Tau]p[\[Tau]]+2 HeavisideTheta[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
I2[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))^2 F3\[Tau][\[Tau]]+2((B r2+A r1)/(B+A))\[CapitalPi]1\[Tau][\[Tau]]+Sqrt[A B]\[CapitalPi]2\[Tau][\[Tau]]-I2s)(*eqn B79*);
I1[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))F3\[Tau][\[Tau]]+\[CapitalPi]1\[Tau][\[Tau]]-I1s)(*eqn B78*);
I0[\[Tau]_]:=\[Tau];
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1\[Tau]m[\[Tau]])-\[Nu]r Ims(*eqn B80*);
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1\[Tau]p[\[Tau]])-\[Nu]r Ips(*eqn B80*);
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to3,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to3[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]m3Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t\[Eta]m3Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t\[Eta]m3Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*case four*)


(* ::Input::Initialization:: *)
t\[Eta]m4[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,]:=Module[{M=1,to3,It,Ip,Imn,I0,I1,I2,S1\[Tau]p,S1\[Tau]m,Ips,Ims,\[CapitalPi]1\[Tau],\[CapitalPi]2\[Tau],I2s,I1s,F4\[Tau],\[CapitalPsi]\[Tau],f2\[Tau]p,f2\[Tau]m,x\[Tau],S1sp,S1sm,xs,\[CapitalPsi]s,F4s,S10\[Tau],S20\[Tau],\[CapitalPi]1s,\[CapitalPi]2s,X4\[Tau],Irs,S20s,S10s,f0\[Tau],p2p,p2m,p20,x4p,x4m,gp,gm,g0,f2sm,f2sp,f0s,\[CapitalPsi]1\[Tau],x1\[Tau],j,D2,k,C2,C,D,a1,b1,x4s,a2,b2,k4,Gt,\[CapitalUpsilon]\[Tau],Gcurlt,Gcurl\[Theta],\[CapitalUpsilon]s,\[Theta]p,\[Theta]m,h,\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
\[Theta]p=ArcCos[Sqrt[um]](*eqn 55*);
\[Theta]m=ArcCos[Sqrt[up]](*eqn 55*);
h=Sign[Cos[\[Theta]s]];
\[CapitalUpsilon]s=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-h/Sqrt[um a^2] EllipticF[\[CapitalUpsilon]s,1-up/um](*eqn 56*);
Gcurlt=-h Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]s,1-up/um](*eqn 58*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
Gt[\[Tau]_]:=Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurlt (*eqn 68*);
b2=Re[(r1+r2)/2](*eqn B11*);
a2=Re[Sqrt[-((r2-r1)^2/4)]](*eqn B11*);
b1=Re[(r3+r4)/2](*eqn B10*);
a1=Re[Sqrt[-((r4-r3)^2/4)]](*eqn B10*);
x4s=Re[(rs-b2)/a2](*eqn B83*);
g0=Sqrt[(4 a2^2-(C-D)^2)/((C+D)^2-4 a2^2)](*eqn B88*);
\[CapitalPsi]s=ArcTan[x4s]+ArcTan[g0](*from eqn B100,B102...*);
xs=\[CapitalPsi]s-\[Pi]/2;
D=Re[Sqrt[(r3-r2)(r4-r1)]](*eqn B85*);
C=Re[Sqrt[(r3-r1)(r4-r2)]](*eqn B85*);
k4=(4 C D)/(C+D)^2 (*eqn B87*);
F4s=2/(C+D) EllipticF[ArcTan[x4s]+ArcTan[g0],k4](*eqn B101*);
Irs=F4s;
X4\[Tau][\[Tau]_]:=(C+D)/2 ( \[Nu]r \[Tau]+ Irs)(*eqn B104*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-\[Pi]/2;
x4m=(rm-b2)/a2(*eqn B83*);
x4p=(rp-b2)/a2(*eqn B83*);
gm=(g0 x4m-1)/(g0+x4m)(*eqn B96*);
gp=(g0 x4p-1)/(g0+x4p)(*eqn B96*);
\[CapitalPsi]\[Tau][\[Tau]_]:= JacobiAmplitude[X4\[Tau][\[Tau]],k4];
j=k4;
p20=(*Sqrt[(1+g0^2)/(1-j+g0^2)]*)Sqrt[1/((1+g0^2)(1-j+g0^2))](*eqn B95*); 
p2m=Sqrt[(1+gm^2)/(1-j+gm^2)](*eqn B95*); 
p2p=Sqrt[(1+gp^2)/(1-j+gp^2)](*eqn B95*);
f0s=p20/2 Log[Abs[(1-p20)/(1+p20) (1+p20 Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1-p20 Sqrt[1-j Sin[\[CapitalPsi]s]^2])]] (*eqn B95*); 
f2sm=p2m/2 Log[Abs[(1-p2m)/(1+p2m) (1+p2m Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1-p2m Sqrt[1-j Sin[\[CapitalPsi]s]^2])]] (*eqn B95*);
f2sp=p2p/2 Log[Abs[(1-p2p)/(1+p2p) (1+p2p Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1-p2p Sqrt[1-j Sin[\[CapitalPsi]s]^2])]] (*eqn B95*); 
S10s=1/(1+g0^2) (EllipticF[\[CapitalPsi]s,j]+g0^2 EllipticPi[1+g0^2,\[CapitalPsi]s,j]-gp f0s+2 g0^2 HeavisideLambda[xs]EllipticPi[1+g0^2,j]) (*eqn B92*);
S20s=-(1/((1+g0^2)(1-j+g0^2)))((1-j)EllipticF[\[CapitalPsi]s,j]+g0^2 EllipticE[\[CapitalPsi]s,j]+(g0^2 Sqrt[1-j Sin[\[CapitalPsi]s]^2](g0-Tan[\[CapitalPsi]s]))/(1+g0 Tan [\[CapitalPsi]s])-g0^3)+(1/(1+g0^2)+(1-j)/(1-j+g0^2))S10s (*eqn B94*);
\[CapitalPi]2s=2/(C+D) (a2/g0 (1+g0^2))^2 S20s (*eqn B116*);
\[CapitalPi]1s=2/(C+D) (a2/g0 (1+g0^2))S10s (*eqn B102*);
f0\[Tau][\[Tau]_]:=p20/2 Log[Abs[(1-p20)/(1+p20) (1+p20 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1-p20 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])]] (*eqn B95*); 
S10\[Tau][\[Tau]_]:=1/(1+g0^2) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+g0^2 EllipticPi[1+g0^2,\[CapitalPsi]\[Tau][\[Tau]],j]-gp f0\[Tau][\[Tau]]+2 g0^2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[1+g0^2,j]) (*eqn B92*);
S20\[Tau][\[Tau]_]:=-(1/((1+g0^2)(1-j+g0^2)))((1-j)EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+g0^2 EllipticE[\[CapitalPsi]\[Tau][\[Tau]],j]+(g0^2 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2](g0-Tan[\[CapitalPsi]\[Tau][\[Tau]]]))/(1+g0 Tan [\[CapitalPsi]\[Tau][\[Tau]]])-g0^3)+(1/(1+g0^2)+(1-j)/(1-j+g0^2))S10\[Tau][\[Tau]] (*eqn B94*);
S1sm=1/(1+gm^2) (EllipticF[\[CapitalPsi]s,j]+gm^2 EllipticPi[1+gm^2,\[CapitalPsi]s,j]-gm f2sm+2gm^2 HeavisideLambda[xs]EllipticPi[1+gm^2,j])(*eqn B92*);
S1sp=1/(1+gp^2) (EllipticF[\[CapitalPsi]s,j]+gp^2 EllipticPi[1+gp^2,\[CapitalPsi]s,j]-gp f2sp+2gp^2 HeavisideLambda[xs]EllipticPi[1+gp^2,j])(*eqn B92*);
\[CapitalPsi]\[Tau][\[Tau]_]:=JacobiAmplitude[X4\[Tau][\[Tau]],k4];
f2\[Tau]m[\[Tau]_]:=p2m/2 Log[Abs[(1-p2m)/(1+p2m) (1+p2m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1-p2m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])]] (*eqn B95*);
f2\[Tau]p[\[Tau]_]:=p2p/2 Log[Abs[(1-p2p)/(1+p2p) (1+p2p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1-p2p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])]] (*eqn B95*); 
F4\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B115*);
I2s=(a2/g0-b1)^2 F4s-2(a2/g0-b1)\[CapitalPi]1s+\[CapitalPi]2s (*eqn B99*);
I1s=(a2/g0-b1)F4s-\[CapitalPi]1s (*eqn B98*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/(C+D) (a2/g0 (1+g0^2))^2 S20\[Tau][\[Tau]] (*eqn B116*);
\[CapitalPi]1\[Tau][\[Tau]_]:=2/(C+D) (a2/g0 (1+g0^2))S10\[Tau][\[Tau]] (*eqn B116*);
Ims=g0/(a2(1-g0 x4m)) (F4s-2/(C+D) ((1+g0^2)/(g0(g0+x4m)))S1sm)(*eqn B100*);
Ips=g0/(a2(1-g0 x4p)) (F4s-2/(C+D) ((1+g0^2)/(g0(g0+x4p)))S1sp)(*eqn B100*);
D2=Re[(a1+a2)^2+(b1-b2)^2](*eqn B85*);
C2=Re[(a1-a2)^2+(b1-b2)^2](*eqn B85*);
k=D2/C2 (*eqn B86*);
\[CapitalPsi]1\[Tau][\[Tau]_]:=JacobiAmplitude[X4\[Tau][\[Tau]],k4];
x1\[Tau][\[Tau]_]:=\[CapitalPsi]1\[Tau][\[Tau]]-\[Pi]/2;
S1\[Tau]m[\[Tau]_]:=1/(1+gm^2) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+gm^2 EllipticPi[1+gm^2,\[CapitalPsi]\[Tau][\[Tau]],j]-gm f2\[Tau]m[\[Tau]]+2gm^2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[1+gm^2,j])(*eqn B92*);
S1\[Tau]p[\[Tau]_]:=1/(1+gp^2) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]+gp^2 EllipticPi[1+gp^2,\[CapitalPsi]\[Tau][\[Tau]],j]-gp f2\[Tau]p[\[Tau]]+2gp^2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[1+gp^2,j])(*eqn B92*);
I2[\[Tau]_]:=\[Nu]r((a2/g0-b1)^2 F4\[Tau][\[Tau]]-2(a2/g0-b1)\[CapitalPi]1\[Tau][\[Tau]]+\[CapitalPi]2\[Tau][\[Tau]]-I2s) (*eqn B112*);
I1[\[Tau]_]:=\[Nu]r((a2/g0-b1)F4\[Tau][\[Tau]]-\[CapitalPi]1\[Tau][\[Tau]]-I1s) (*eqn B112*);
I0[\[Tau]_]:=\[Tau];
Imn[\[Tau]_]:=(\[Nu]r g0)/(a2(1-g0 x4m)) (F4\[Tau][\[Tau]]-2/(C+D) ((1+g0^2)/(g0(g0+x4m)))S1\[Tau]m[\[Tau]])-\[Nu]r Ims (*eqn B114*);
Ip[\[Tau]_]:=(\[Nu]r g0)/(a2(1-g0 x4p)) (F4\[Tau][\[Tau]]-2/(C+D) ((1+g0^2)/(g0(g0+x4p)))S1\[Tau]p[\[Tau]])-\[Nu]r Ips (*eqn B114*);
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to3,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to3[\[Tau]_]:=It[\[Tau]]+a^2 Gt[\[Tau]](*eqn 12*);
t\[Eta]m4Function[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t\[Eta]m4Function["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t\[Eta]m4Function[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*t Principal Null Geodesics*)


(* ::Input::Initialization:: *)
schwarzschildtPNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to,It,Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,consts,\[Nu]\[Theta],\[Nu]r,ro,I0,I1,I2},
\[Nu]r=Sign[prs];
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
rp=M+Sqrt[M-a^2];
rm=M-Sqrt[M-a^2];
\[Nu]r=Sign[prs];
ro[\[Tau]_]:=-((rs \[Nu]r)/(-\[Nu]r+rs \[Tau]));

Ims=-(1/(2 rs^2));   
Ips=1/(rs rp)-Log[rs]/rp^2+Log[rs-rp]/rp^2;      
Imn[\[Tau]_]:=\[Nu]r(-(1/(2 ro[\[Tau]]^2))-Ims);     
Ip[\[Tau]_]:=\[Nu]r(1/(ro[\[Tau]] rp)-Log[ro[\[Tau]]]/rp^2+Log[ro[\[Tau]]-rp]/rp^2-Ips);   
I2[\[Tau]_]:=\[Nu]r(ro[\[Tau]]-rs);
I1[\[Tau]_]:=\[Nu]r(Log[ro[\[Tau]]]-Log[rs]);
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] ; 
assoc=Association["Trajectory"-> to];
to[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 11*);   
schwarzschildtPNGFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[schwarzschildtPNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="schwarzschildtPNGFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
schwarzschildtPNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
schwarzschildtPNGFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
KerrtPNGEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to,It,Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,consts,\[Nu]\[Theta],\[Nu]r,ro,I0,I1,I2},
\[Nu]r=Sign[prs];
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
rp=M+Sqrt[M-a^2];
rm=M-Sqrt[M-a^2];
\[Nu]r=Sign[prs];
ro[\[Tau]_]:=-((rs \[Nu]r)/(-\[Nu]r+rs \[Tau]));

Ims=1/(rs rm)-Log[rs]/rm^2+Log[rs-rm]/rm^2;   
Ips=1/(rs rp)-Log[rs]/rp^2+Log[rs-rp]/rp^2;      
Imn[\[Tau]_]:=\[Nu]r(1/(ro[\[Tau]] rm)-Log[ro[\[Tau]]]/rm^2+Log[ro[\[Tau]]-rm]/rm^2-Ims);    
Ip[\[Tau]_]:=\[Nu]r(1/(ro[\[Tau]] rp)-Log[ro[\[Tau]]]/rp^2+Log[ro[\[Tau]]-rp]/rp^2-Ips);   
I2[\[Tau]_]:=\[Nu]r(ro[\[Tau]]-rs);
I1[\[Tau]_]:=\[Nu]r(Log[ro[\[Tau]]]-Log[rs]);
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] ; 
assoc=Association["Trajectory"-> to];
to[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 11*);   
KerrtPNGEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[KerrtPNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="KerrtPNGEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
KerrtPNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
KerrtPNGEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
KerrtPNGnonE[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to,It,Ip,Imn,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],F2\[Tau],Ips,Ims,X2\[Tau],\[CapitalPi]ps,\[CapitalPi]ms,F2s,Irs,x2,k,rp,rm,r1,r2,r3,r4,G\[Phi],up,um,\[CapitalDelta]\[Theta],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi],roots,assoc,\[Nu]\[Theta],\[Nu]r,ro,I0,I1s,I1,I2s,I2},
\[Nu]r=Sign[prs];
\[Nu]\[Theta]=Sign[p\[Theta]s];
up=I Cos[\[Theta]s]^2; 
 um=I Cos[\[Theta]s]^2; 
Gcurl\[Phi]=-(1/Sqrt[-um a^2])EllipticPi[up,ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 30*);
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
G\[Phi][\[Tau]_]:=1/Sqrt[-um a^2] EllipticPi[up,\[Psi]\[Tau][\[Tau]],up/um]-\[Nu]\[Theta] Gcurl\[Phi](*eqn 47*);
(*\[Lambda]=a Sin[\[Theta]s]^2;
\[Eta]=-a^2Cos[\[Theta]s]^4;*)
rp=M+Sqrt[M-a^2];
rm=M-Sqrt[M-a^2];
r1=I a Cos[\[Theta]s];
r2=-I a Cos[\[Theta]s];
\[Nu]r=Sign[prs];
ro[\[Tau]_]:=Re[(E^((r1 \[Tau])/\[Nu]r) r1 r2-E^((r2 \[Tau])/\[Nu]r) r1 r2+E^((r2 \[Tau])/\[Nu]r) r1 rs-E^((r1 \[Tau])/\[Nu]r) r2 rs)/(E^((r1 \[Tau])/\[Nu]r) r1-E^((r2 \[Tau])/\[Nu]r) r2-E^((r1 \[Tau])/\[Nu]r) rs+E^((r2 \[Tau])/\[Nu]r) rs)];

Ims=((r2-rm) Log[rs-r1]+(-r1+rm) Log[rs-r2]+(r1-r2) Log[rs-rm])/((r1-r2) (r1-rm) (r2-rm));   
Ips=((r2-rp) Log[rs-r1]+(-r1+rp) Log[rs-r2]+(r1-r2) Log[rs-rp])/((r1-r2) (r1-rp) (r2-rp));      
Imn[\[Tau]_]:=\[Nu]r(((r2-rm) Log[ro[\[Tau]]-r1]+(-r1+rm) Log[ro[\[Tau]]-r2]+(r1-r2) Log[ro[\[Tau]]-rm])/((r1-r2) (r1-rm) (r2-rm))-Ims);    
Ip[\[Tau]_]:=\[Nu]r(((r2-rp) Log[ro[\[Tau]]-r1]+(-r1+rp) Log[ro[\[Tau]]-r2]+(r1-r2) Log[ro[\[Tau]]-rp])/((r1-r2) (r1-rp) (r2-rp))-Ips); 
I2s=rs+(r1^2 Log[rs-r1])/(r1-r2)-(r2^2 Log[rs-r2])/(r1-r2);
I2[\[Tau]_]:=\[Nu]r((ro[\[Tau]]+(r1^2 Log[ro[\[Tau]]-r1])/(r1-r2)-(r2^2 Log[ro[\[Tau]]-r2])/(r1-r2))-I2s); 
I1s=(r1 Log[rs-r1]-r2 Log[rs-r2])/(r1-r2); 
I1[\[Tau]_]:=\[Nu]r ((r1 Log[ro[\[Tau]]-r1]-r2 Log[ro[\[Tau]]-r2])/(r1-r2)-I1s);  
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] ;
assoc=Association["Trajectory"-> to];
to[\[Tau]_]:=It[\[Tau]]+a^2 Sin[\[Theta]s]^2 \[Tau](*eqn 11*);  
KerrtPNGnonEFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[KerrtPNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="KerrtPNGnonEFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
KerrtPNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
KerrtPNGnonEFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
TemporalPNG[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{},
If[a==0,schwarzschildtPNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],If[\[Theta]s==\[Pi]/2,KerrtPNGEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],KerrtPNGnonE[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]
]


(* ::Subsection::Initialization:: *)
(*t equatorial*)


(* ::Subsubsection::Initialization:: *)
(*t1*)


(* ::Input::Initialization:: *)
t1Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_]:=Module[{M=1,Ips,Ipo,Imns,Imno,g,Intps,Intpo,Intms,Intmo,\[Alpha]2p,\[Alpha]2m,\[Alpha]12,x1s,us,uo,X1,k,\[Kappa]s,Irs,\[Psi]o,\[Psi]s,Imn,Ip,It,to1,assoc,\[Nu]r,roots,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,r1,r2,r3,r4,consts,\[Nu]\[Theta],rm,rp,G\[Phi],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi], um, up,\[CapitalDelta]\[Theta],Gt,Ep1,Ep, Gcurlt,I0,aa,b,c,d,\[Alpha]1,\[Alpha],I1o,I1s,I1,I2,I2o,I2s,V1o,V1s,V2o,V2s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
aa=r4;
b=r3;
c=r2;
d=r1;
\[Alpha]1=Sqrt[(aa(b-c))/(b(aa-c))];
\[Alpha]=Sqrt[(b-c)/(aa-c)];
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
g=2/Sqrt[(r4-r2)(r3-r1)];
\[Alpha]2p=((rp-r4)(r3-r2))/((rp-r3)(r4-r2));
\[Alpha]2m=((rm-r4)(r3-r2))/((rm-r3)(r4-r2));
\[Alpha]12=(r3-r2)/(r4-r2);
x1s=Sqrt[((r4-r2)(r3-rs))/((r3-r2)(r4-rs))];
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))];
\[Kappa]s=ArcSin[x1s];
Irs=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[\[Kappa]s,k];  
X1[\[Tau]_]:=Sqrt[(r4-r2)(r3-r1)]/2 (\[Tau]+\[Nu]r Irs);
\[Psi]s=ArcSin[x1s];
\[Psi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X1[\[Tau]],k];
us=EllipticF[\[Psi]s,k];
uo[\[Tau]_]:=EllipticF[\[Psi]o[\[Tau]],k];
Imns=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]s,k]+(rm-b)EllipticF[\[Psi]s,k]);
Imno[\[Tau]_]:=2/((rm-aa)(rm-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rm-a))/((aa-c)(rm-b)),\[Psi]o[\[Tau]],k]+(rm-b)EllipticF[\[Psi]o[\[Tau]],k]);
Ips=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]s,k]+(rp-b)EllipticF[\[Psi]s,k]);
Ipo[\[Tau]_]:=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]o[\[Tau]],k]+(rp-b)EllipticF[\[Psi]o[\[Tau]],k]);
Imn[\[Tau]_]:=-\[Nu]r(Imno[\[Tau]]-Imns)(*eqn B30c*);
Ip[\[Tau]_]:=-\[Nu]r(Ipo[\[Tau]]-Ips)(*eqn B30c*);
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[us,k]+(k^2-\[Alpha]^2)us+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]-(\[Alpha]^4 JacobiSN[us,k]JacobiCN[us,k]JacobiDN[us,k])/(1-\[Alpha]^2 JacobiSN[us,k]^2)); 
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[uo[\[Tau]],k]+(k^2-\[Alpha]^2)uo[\[Tau]]+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]-(\[Alpha]^4 JacobiSN[uo[\[Tau]],k]JacobiCN[uo[\[Tau]],k]JacobiDN[uo[\[Tau]],k])/(1-\[Alpha]^2 JacobiSN[uo[\[Tau]],k]^2));
V1s=EllipticPi[\[Alpha]^2,\[Psi]s,k];
V1o[\[Tau]_]:=EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k];
I2s=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 us+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2o[\[Tau]_]:=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 uo[\[Tau]]+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s);
I1s=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]+\[Alpha]1^2 us);
I1o[\[Tau]_]:=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]+\[Alpha]1^2 uo[\[Tau]]);
I1[\[Tau]_]:=\[Nu]r (I1o[\[Tau]]-I1s);
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to1,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to1[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 12*);
t1EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_,assoc_]]:="t1EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t1EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization:: *)
(*t2*)


(* ::Input::Initialization:: *)
t2Equatorialx[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],E2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,I2o,V1s,V1o,V2s,V2o,\[CapitalPsi]o,\[CapitalPsi]s,\[Alpha],\[Alpha]1,aa,b,c,d,y,uo,us,g,u1,u2,x2ss,x2o,ro,V14,V24,\[CapitalPsi]4,x24,I24,I2s4,I24o},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]EllipticE[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B52*);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
aa=r4;
b=r3;
c=r2;
d=r1;
y=rs;
\[Alpha]1=Sqrt[((aa-d)b)/((b-d)aa)];
\[Alpha]=Sqrt[(aa-d)/(b-d)];
g=1/Sqrt[(aa-c)(b-d)];
x2o[\[Tau]_]:=Sqrt[((b-d)(ro[\[Tau]]-aa))/((aa-d)(ro[\[Tau]]-b))];
\[CapitalPsi]s=ArcSin[x2s];
\[CapitalPsi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k];
V1o[\[Tau]_]:=EllipticPi[\[CapitalPsi]o[\[Tau]],\[Alpha]^2,k];
V1s=EllipticPi[\[CapitalPsi]s,\[Alpha]^2,k];
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]o[\[Tau]],k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1o[\[Tau]]-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiCN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiDN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]^2));    
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]s,k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]s,k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1s-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]s,k],k] JacobiCN[EllipticF[\[CapitalPsi]s,k],k] JacobiDN[EllipticF[\[CapitalPsi]s,k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]s,k],k]^2));    
I2o[\[Tau]_]:=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]o[\[Tau]],k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]); 
I2s=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]s,k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s); 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 12*);
t2EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
t2Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,H2\[Tau],R,ro\[Tau],rop\[Tau],E2\[Tau]},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r3-r2)](EllipticE[JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticE[ArcSin[x2s],k]);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
ro\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
rop\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]'[\[Tau]];
H2\[Tau][\[Tau]_]:=rop\[Tau][\[Tau]]/(ro\[Tau][\[Tau]]-r3)-\[Nu]r Sqrt[Rs]/(rs-r3); 
I2[\[Tau]_]:=H2\[Tau][\[Tau]]-((r1-r4)+(r2-r3))/2-E2\[Tau]; 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 12*);
t2EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t2EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t2EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization:: *)
(*t3*)


(* ::Input::Initialization:: *)
t3Equatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to3,It,Ip,Imn,I0,I1,I2,R1\[Tau]p,R1\[Tau]m,x\[Tau],\[CapitalPsi]\[Tau],Ips,Ims,\[CapitalPi]1\[Tau],I1s,\[CapitalPi]2\[Tau],I2s,f1\[Tau]p,f1\[Tau]m,P1p,P1m,j,X3\[Tau],R1sp,R1sm,\[CapitalPsi]s,xs,R10\[Tau],\[CapitalPi]1s,R20\[Tau],f1sp,f1sm,F3\[Tau],Irs,\[CapitalPi]2s,f10\[Tau],P10,R10s,R20s,f10s,F3s,\[Alpha]p,\[Alpha]m,\[Alpha]0,x3s,A,B,k3,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
x3s=(A(rs-r1)-B(rs-r2))/(A(rs-r1)+B(rs-r2))(*eqn B55*);
\[Alpha]0=(B+A)/(B-A)(*eqn B58*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*);
\[Alpha]m=(B(rm-r2)+A(rm-r1))/(B(rm-r2)-A(rm-r1))(*eqn B66*); 
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
j=k3(*eqn B80 and B62*);
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+\[Nu]r Irs)(*eqn B74*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-Pi/2(*eqn B62*);
\[CapitalPsi]\[Tau][\[Tau]_]:=\[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3](*eqn B80 and B62*);
P10=Sqrt[(\[Alpha]0^2-1)/(j+(1-j)\[Alpha]0^2)](*eqn B65*);
f10\[Tau][\[Tau]_]:=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
R10\[Tau][\[Tau]_]:=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0 f10\[Tau][\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20\[Tau][\[Tau]_]:=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]\[Tau][\[Tau]],j]-(\[Alpha]0 Sin[\[CapitalPsi]\[Tau][\[Tau]]]Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]\[Tau][\[Tau]]])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10\[Tau][\[Tau]](*eqn B64*);
\[CapitalPsi]s=ArcCos[x3s](*eqn B70*);
xs=\[CapitalPsi]s-Pi/2(*eqn B62*);
f10s=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R10s=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]s,j]-\[Alpha]0 f10s+2 HeavisideLambda[xs]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20s=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]s,j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]s,j]-(\[Alpha]0 Sin[\[CapitalPsi]s]Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]s])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10s(*eqn B64*);
\[CapitalPi]2s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20s(*eqn B72*);
\[CapitalPi]1s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10s(*eqn B72*);
P1m=Sqrt[(\[Alpha]m^2-1)/(j+(1-j)\[Alpha]m^2)](*eqn B65*);
P1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
f1sm=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1m Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
f1sp=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1sm=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]s,j]-\[Alpha]m f1sm+2HeavisideLambda[xs]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1sp=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1sp+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
I2s=((B r2+A r1)/(B+A))^2 F3s+2((B r2+A r1)/(B+A))\[CapitalPi]1s+Sqrt[A B]\[CapitalPi]2s(*eqn B69*);
\[CapitalPi]2\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20\[Tau][\[Tau]](*eqn B82*);
I1s=((B r2+A r1)/(B+A))F3s+\[CapitalPi]1s(*eqn B68*);
\[CapitalPi]1\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10\[Tau][\[Tau]](*eqn B82*);
Ims=-(1/(B(rm-r2)+A(rm-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1sm)(*eqn B70*); 
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1sp)(*eqn B70*);
f1\[Tau]m[\[Tau]_]:=P1m/2 Log[Abs[(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1m Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
f1\[Tau]p[\[Tau]_]:=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
F3\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B81*);
R1\[Tau]m[\[Tau]_]:=1/(1-\[Alpha]m^2) (EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]m f1\[Tau]m[\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]m^2/(\[Alpha]m^2-1),j])(*eqn B62*);
R1\[Tau]p[\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1\[Tau]p[\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
I2[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))^2 F3\[Tau][\[Tau]]+2((B r2+A r1)/(B+A))\[CapitalPi]1\[Tau][\[Tau]]+Sqrt[A B]\[CapitalPi]2\[Tau][\[Tau]]-I2s)(*eqn B79*);
I1[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))F3\[Tau][\[Tau]]+\[CapitalPi]1\[Tau][\[Tau]]-I1s)(*eqn B78*);
I0[\[Tau]_]:=\[Tau];
Imn[\[Tau]_]:=-(\[Nu]r/(B(rm-r2)+A(rm-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rm-r2)-A(rm-r1)) R1\[Tau]m[\[Tau]])-\[Nu]r Ims(*eqn B80*); 
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1\[Tau]p[\[Tau]])-\[Nu]r Ips(*eqn B80*);
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]]-rm(rm-(a \[Lambda])/(2 M))Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to3];
to3[\[Tau]_]:=It[\[Tau]]+a^2 \[Tau](*eqn 12*);
t3EquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t3EquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t3EquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*t Schwarzchild*)


(* ::Subsubsection::Initialization::Closed:: *)
(*t1*)


(* ::Input::Initialization:: *)
t1Schwarzchild[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,p\[Phi]s_]:=Module[{M=1,Ips,Ipo,Imns,Imno,g,Intps,Intpo,Intms,Intmo,\[Alpha]2p,\[Alpha]2m,\[Alpha]12,x1s,us,uo,X1,k,\[Kappa]s,Irs,\[Psi]o,\[Psi]s,Imn,Ip,It,to1,assoc,\[Nu]r,roots,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,r1,r2,r3,r4,consts,\[Nu]\[Theta],rm,rp,G\[Phi],\[Psi]\[Tau],Gcurl\[Theta],Gcurl\[Phi], um, up,\[CapitalDelta]\[Theta],Gt,Ep1,Ep, Gcurlt,I0,aa,b,c,d,\[Alpha]1,\[Alpha],I1o,I1s,I1,I2,I2o,I2s,V1o,V1s,V2o,V2s},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
aa=r4;
b=r3;
c=r2;
d=r1;
\[Alpha]1=Sqrt[(aa(b-c))/(b(aa-c))];
\[Alpha]=Sqrt[(b-c)/(aa-c)];
rm=M-Sqrt[M^2-a^2];
rp=M+Sqrt[M^2-a^2];
g=2/Sqrt[(r4-r2)(r3-r1)];
\[Alpha]12=(r3-r2)/(r4-r2);
x1s=Sqrt[((r4-r2)(r3-rs))/((r3-r2)(r4-rs))];
k=Sqrt[((r3-r2)(r4-r1))/((r4-r2)(r3-r1))];
\[Kappa]s=ArcSin[x1s];
Irs=2/Sqrt[(r4-r2)(r3-r1)] EllipticF[\[Kappa]s,k];  
X1[\[Tau]_]:=Sqrt[(r4-r2)(r3-r1)]/2 (\[Tau]+\[Nu]r Irs);
\[Psi]s=ArcSin[x1s];
\[Psi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X1[\[Tau]],k];
us=EllipticF[\[Psi]s,k];
uo[\[Tau]_]:=EllipticF[\[Psi]o[\[Tau]],k];
Ips=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]s,k]+(rp-b)EllipticF[\[Psi]s,k]);
Ipo[\[Tau]_]:=2/((rp-aa)(rp-b)Sqrt[(aa-c)(b-d)]) ((b-aa)EllipticPi[((b-c)(rp-a))/((aa-c)(rp-b)),\[Psi]o[\[Tau]],k]+(rp-b)EllipticF[\[Psi]o[\[Tau]],k]);
Ip[\[Tau]_]:=-\[Nu]r(Ipo[\[Tau]]-Ips)(*eqn B30c*);
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[us,k]+(k^2-\[Alpha]^2)us+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]-(\[Alpha]^4 JacobiSN[us,k]JacobiCN[us,k]JacobiDN[us,k])/(1-\[Alpha]^2 JacobiSN[us,k]^2)); 
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[uo[\[Tau]],k]+(k^2-\[Alpha]^2)uo[\[Tau]]+(2\[Alpha]^2 k^2+2\[Alpha]^2-\[Alpha]^4-3k^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]-(\[Alpha]^4 JacobiSN[uo[\[Tau]],k]JacobiCN[uo[\[Tau]],k]JacobiDN[uo[\[Tau]],k])/(1-\[Alpha]^2 JacobiSN[uo[\[Tau]],k]^2));
V1s=EllipticPi[\[Alpha]^2,\[Psi]s,k];
V1o[\[Tau]_]:=EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k];
I2s=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 us+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2o[\[Tau]_]:=(b^2 g)/\[Alpha]^4 (\[Alpha]1^4 uo[\[Tau]]+2\[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s);
I1s=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]s,k]+\[Alpha]1^2 us);
I1o[\[Tau]_]:=(b g)/\[Alpha]^2 ((\[Alpha]^2-\[Alpha]1^2)EllipticPi[\[Alpha]^2,\[Psi]o[\[Tau]],k]+\[Alpha]1^2 uo[\[Tau]]);
I1[\[Tau]_]:=\[Nu]r (I1o[\[Tau]]-I1s);
I0[\[Tau]_]:=\[Tau];
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp)Ip[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to1,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to1[\[Tau]_]:=It[\[Tau]](*eqn 12*);
t1SchwarzchildFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t1SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t1SchwarzchildFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t1SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t1SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization::Closed:: *)
(*t2*)


(* ::Input::Initialization:: *)
t2Shwarzchildx[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],E2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,I2o,V1s,V1o,V2s,V2o,\[CapitalPsi]o,\[CapitalPsi]s,\[Alpha],\[Alpha]1,aa,b,c,d,y,uo,us,g,u1,u2},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]EllipticE[\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B52*);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
aa=r4;
b=r3;
c=r2;
d=r1;
y=rs;
\[Alpha]1=Sqrt[((aa-d)b)/((b-d)aa)];
\[Alpha]=Sqrt[(aa-d)/(b-d)];
g=1/Sqrt[(aa-c)(b-d)];
\[CapitalPsi]s=ArcSin[x2s];
\[CapitalPsi]o[\[Tau]_]:=\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k];
V1o[\[Tau]_]:=EllipticPi[\[CapitalPsi]o[\[Tau]],\[Alpha]^2,k];
V1s=EllipticPi[\[CapitalPsi]s,\[Alpha]^2,k];
V2o[\[Tau]_]:=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]o[\[Tau]],k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1o[\[Tau]]-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiCN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k] JacobiDN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]o[\[Tau]],k],k]^2));    
V2s=1/(2(\[Alpha]^2-1)(k^2-\[Alpha]^2)) (\[Alpha]^2 EllipticE[EllipticF[\[CapitalPsi]s,k],k]+(k^2-\[Alpha]^2)EllipticF[\[CapitalPsi]s,k]+(2 \[Alpha]^2 k^2+2 \[Alpha]^2-\[Alpha]^4-3 k^2)V1s-(\[Alpha]^4 JacobiSN[EllipticF[\[CapitalPsi]s,k],k] JacobiCN[EllipticF[\[CapitalPsi]s,k],k] JacobiDN[EllipticF[\[CapitalPsi]s,k],k])/(1-\[Alpha]^2 JacobiSN[EllipticF[\[CapitalPsi]s,k],k]^2));    
I2o[\[Tau]_]:=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]o[\[Tau]],k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1o[\[Tau]]+(\[Alpha]^2-\[Alpha]1^2)^2 V2o[\[Tau]]); 
I2s=aa^2 g 1/\[Alpha]^4 (\[Alpha]1^4 EllipticF[\[CapitalPsi]s,k]+2 \[Alpha]1^2 (\[Alpha]^2-\[Alpha]1^2)V1s+(\[Alpha]^2-\[Alpha]1^2)^2 V2s);
I2[\[Tau]_]:=\[Nu]r (I2o[\[Tau]]-I2s); 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp)Ip[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]](*eqn 12*);
t2ShwarzchildxFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
t2Shwarzchild[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to2,It,Ip,Imn,I0,I1,I2,\[CapitalPi]p\[Tau],\[CapitalPi]m\[Tau],Ips,Ims,\[CapitalPi]2\[Tau],I1s,\[CapitalPi]ps,\[CapitalPi]ms,F2s,E2s,F2\[Tau],\[CapitalPi]2s,Ir,I2s,X2\[Tau],Irs,x2s,k,Rs,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,consts,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,\[Nu]\[Theta],\[Nu]r,H2\[Tau],R,ro\[Tau],rop\[Tau],E2\[Tau]},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
rootsA=KerrNullGeoAngularRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4}=Values[rootsA];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
 k=((r3-r2)(r4-r1))/((r3-r1)(r4-r2))(*eqn B13*);
x2s=Sqrt[(rs-r4)/(rs-r3) (r3-r1)/(r4-r1)](*eqn B35*);
F2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticF[ArcSin[x2s],k](*eqn B40*); 
Irs=F2s;
X2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r4-r2)]/2 (\[Tau]+\[Nu]r Irs)(*eqn B45*);  
E2\[Tau][\[Tau]_]:=Sqrt[(r3-r1)(r3-r2)](EllipticE[JacobiAmplitude[X2\[Tau][\[Tau]],k],k]-\[Nu]r EllipticE[ArcSin[x2s],k]);
F2\[Tau][\[Tau]_]:=\[Nu]r \[Tau]+Irs; 
E2s=Sqrt[(r3-r1)(r4-r2)]EllipticE[ArcSin[x2s],k](*eqn B41*);
\[CapitalPi]2s=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),ArcSin[x2s],k](*eqn B42*);
\[CapitalPi]ms=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
\[CapitalPi]ps=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),ArcSin[x2s],k](*eqn B43*);
I1s=r3 F2s+(r4-r3)\[CapitalPi]2s (*eqn B37*);
Rs=(rs-r1)(rs-r2)(rs-r3)(rs-r4)(*eqn B8*);
\[CapitalPi]2\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] EllipticPi[(r4-r1)/(r3-r1),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B53*);
Ims=-\[CapitalPi]ms-F2s/(rm-r3)(*eqn B39*); 
Ips=-\[CapitalPi]ps-F2s/(rp-r3)(*eqn B39*);    
\[CapitalPi]m\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rm-r3)(rm-r4)) EllipticPi[((rm-r3)(r4-r1))/((rm-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
\[CapitalPi]p\[Tau][\[Tau]_]:=2/Sqrt[(r3-r1)(r4-r2)] (r4-r3)/((rp-r3)(rp-r4)) EllipticPi[((rp-r3)(r4-r1))/((rp-r4)(r3-r1)),\[Nu]r JacobiAmplitude[X2\[Tau][\[Tau]],k],k](*eqn B54*);
ro\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s][\[Tau]];
rop\[Tau][\[Tau]_]:=RadialMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]'[\[Tau]];
H2\[Tau][\[Tau]_]:=rop\[Tau][\[Tau]]/(ro\[Tau][\[Tau]]-r3)-\[Nu]r Sqrt[Rs]/(rs-r3); 
I2[\[Tau]_]:=H2\[Tau][\[Tau]]-((r1-r4)+(r2-r3))/2-E2\[Tau]; 
I1[\[Tau]_]:=\[Nu]r(r3 F2\[Tau][\[Tau]]+(r4-r3)\[CapitalPi]2\[Tau][\[Tau]]-I1s)(*eqn B48*);
Ir[\[Tau]_]:=\[Tau];
I0[\[Tau]_]:=Ir[\[Tau]](*eqn B1*);
Imn[\[Tau]_]:=\[Nu]r(-\[CapitalPi]m\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rm-r3)-Ims)(*eqn B50*); 
Ip[\[Tau]_]:=\[Nu]r(-\[CapitalPi]p\[Tau][\[Tau]]-F2\[Tau][\[Tau]]/(rp-r3)-Ips)(*eqn B50*);  
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp)Ip[\[Tau]]-rm(rm)Imn[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to2,"RadialRoots"-> roots,"AngularRoots"-> rootsA,"ConstantsofMotion"-> consts];
to2[\[Tau]_]:=It[\[Tau]](*eqn 12*);
t2ShwarzchildFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t2ShwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t2ShwarzchildFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t2ShwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t2ShwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsubsection::Initialization:: *)
(*t3*)


(* ::Input::Initialization:: *)
t3Schwarzchild[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,to3,It,Ip,Imn,I0,I1,I2,R1\[Tau]p,R1\[Tau]m,x\[Tau],\[CapitalPsi]\[Tau],Ips,Ims,\[CapitalPi]1\[Tau],I1s,\[CapitalPi]2\[Tau],I2s,f1\[Tau]p,f1\[Tau]m,P1p,P1m,j,X3\[Tau],R1sp,R1sm,\[CapitalPsi]s,xs,R10\[Tau],\[CapitalPi]1s,R20\[Tau],f1sp,f1sm,F3\[Tau],Irs,\[CapitalPi]2s,f10\[Tau],P10,R10s,R20s,f10s,F3s,\[Alpha]p,\[Alpha]m,\[Alpha]0,x3s,A,B,k3,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[Nu]r=Sign[prs]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
A=Re[Sqrt[(r3-r2)(r4-r2)]](*eqn B57*);
B=Re[Sqrt[(r3-r1)(r4-r1)]](*eqn B57*);
x3s=(A(rs-r1)-B(rs-r2))/(A(rs-r1)+B(rs-r2))(*eqn B55*);
\[Alpha]0=(B+A)/(B-A)(*eqn B58*);
\[Alpha]p=(B(rp-r2)+A(rp-r1))/(B(rp-r2)-A(rp-r1))(*eqn B66*);
k3=((A+B)^2-(r2-r1)^2)/(4 A B)(*eqn B59*);
j=k3(*eqn B80 and B62*);
X3\[Tau][\[Tau]_]:=Sqrt[A B]( \[Tau]+\[Nu]r Irs)(*eqn B74*);
x\[Tau][\[Tau]_]:=\[CapitalPsi]\[Tau][\[Tau]]-Pi/2(*eqn B62*);
\[CapitalPsi]\[Tau][\[Tau]_]:=\[Nu]r JacobiAmplitude[X3\[Tau][\[Tau]],k3](*eqn B80 and B62*);
P10=Sqrt[(\[Alpha]0^2-1)/(j+(1-j)\[Alpha]0^2)](*eqn B65*);
f10\[Tau][\[Tau]_]:=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P10 Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
R10\[Tau][\[Tau]_]:=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0 f10\[Tau][\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20\[Tau][\[Tau]_]:=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]\[Tau][\[Tau]],j]-(\[Alpha]0 Sin[\[CapitalPsi]\[Tau][\[Tau]]]Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]\[Tau][\[Tau]]])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10\[Tau][\[Tau]](*eqn B64*);
\[CapitalPsi]s=ArcCos[x3s](*eqn B70*);
xs=\[CapitalPsi]s-Pi/2(*eqn B62*);
f10s=P10/2 Log[Abs[(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P10 Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R10s=1/(1-\[Alpha]0^2) (EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),\[CapitalPsi]s,j]-\[Alpha]0 f10s+2 HeavisideLambda[xs]EllipticPi[\[Alpha]0^2/(\[Alpha]0^2-1),j])(*eqn B62*);
R20s=1/(\[Alpha]0^2-1) (EllipticF[\[CapitalPsi]s,j]-\[Alpha]0^2/(j+(1-j)\[Alpha]0^2) (EllipticE[\[CapitalPsi]s,j]-(\[Alpha]0 Sin[\[CapitalPsi]s]Sqrt[1-j Sin[\[CapitalPsi]s]^2])/(1+\[Alpha]0 Cos[\[CapitalPsi]s])))+1/(j+(1-j)\[Alpha]0^2) (2 j-\[Alpha]0^2/(\[Alpha]0^2-1))R10s(*eqn B64*);
\[CapitalPi]2s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20s(*eqn B72*);
\[CapitalPi]1s=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10s(*eqn B72*);
P1p=Sqrt[(\[Alpha]p^2-1)/(j+(1-j)\[Alpha]p^2)](*eqn B65*);
f1sp=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]+Sin[\[CapitalPsi]s])/(P1p Sqrt[1-j Sin[\[CapitalPsi]s]^2]-Sin[\[CapitalPsi]s])]](*eqn B65*);
R1sp=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]s,j]-\[Alpha]p f1sp+2HeavisideLambda[xs]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
F3s=1/Sqrt[A B] EllipticF[ArcCos[x3s],k3](*eqn B71*);
Irs=F3s;
I2s=((B r2+A r1)/(B+A))^2 F3s+2((B r2+A r1)/(B+A))\[CapitalPi]1s+Sqrt[A B]\[CapitalPi]2s(*eqn B69*);
\[CapitalPi]2\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^2 R20\[Tau][\[Tau]](*eqn B82*);
I1s=((B r2+A r1)/(B+A))F3s+\[CapitalPi]1s(*eqn B68*);
\[CapitalPi]1\[Tau][\[Tau]_]:=((2 (r2-r1)Sqrt[A B])/(B^2-A^2))^1 R10\[Tau][\[Tau]](*eqn B82*);
Ips=-(1/(B(rp-r2)+A(rp-r1)))((B+A)F3s+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1sp)(*eqn B70*);
f1\[Tau]p[\[Tau]_]:=P1p/2 Log[Abs[(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]+Sin[\[CapitalPsi]\[Tau][\[Tau]]])/(P1p Sqrt[1-j Sin[\[CapitalPsi]\[Tau][\[Tau]]]^2]-Sin[\[CapitalPsi]\[Tau][\[Tau]]])]](*eqn B65*);
F3\[Tau][\[Tau]_]:=(\[Nu]r \[Tau]+ Irs)(*eqn B81*);
R1\[Tau]p[\[Tau]_]:=1/(1-\[Alpha]p^2) (EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),\[CapitalPsi]\[Tau][\[Tau]],j]-\[Alpha]p f1\[Tau]p[\[Tau]]+2 HeavisideLambda[x\[Tau][\[Tau]]]EllipticPi[\[Alpha]p^2/(\[Alpha]p^2-1),j])(*eqn B62*);
I2[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))^2 F3\[Tau][\[Tau]]+2((B r2+A r1)/(B+A))\[CapitalPi]1\[Tau][\[Tau]]+Sqrt[A B]\[CapitalPi]2\[Tau][\[Tau]]-I2s)(*eqn B79*);
I1[\[Tau]_]:=\[Nu]r(((B r2+A r1)/(B+A))F3\[Tau][\[Tau]]+\[CapitalPi]1\[Tau][\[Tau]]-I1s)(*eqn B78*);
I0[\[Tau]_]:=\[Tau];
Ip[\[Tau]_]:=-(\[Nu]r/(B(rp-r2)+A(rp-r1)))((B+A)F3\[Tau][\[Tau]]+(2(r2-r1)Sqrt[A B])/(B(rp-r2)-A(rp-r1)) R1\[Tau]p[\[Tau]])-\[Nu]r Ips(*eqn B80*);
It[\[Tau]_]:=(2 M)^2/(rp-rm) (rp(rp-(a \[Lambda])/(2 M))Ip[\[Tau]])+(2 M)^2 I0[\[Tau]]+(2M)I1[\[Tau]]+I2[\[Tau]] (*eqn B3*);
assoc=Association["Trajectory"-> to3];
to3[\[Tau]_]:=It[\[Tau]](*eqn 12*);
t3SchwarzchildFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[t3SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="t3SchwarzchildFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
t3SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
t3SchwarzchildFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Temporal motion non spherical general code*)


(* ::Input::Initialization:: *)
TemporalNonSphericalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,rm,roots,r1,r2,r3,r4,consts,\[Nu]r},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
\[Nu]r=Sign[prs];
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
If[a>0,Which[\[Eta]\[Element]PositiveReals,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[t\[Eta]p2[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[t\[Eta]p1[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[t\[Eta]p3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]==0,If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[t2Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[t1Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[t3Equatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],\[Eta]\[Element]NegativeReals,If[r2 \[Element] Reals,
If[r4 \[Element]Complexes,Return[t\[Eta]m3[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]],Return[t\[Eta]m4[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],If[r2 \[Element] Reals,
If[r4 \[Element] Reals,
If[rs>=r4,Return[t2Shwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],If[rs<=r3,Return[t1Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]],Return[t3Schwarzchild[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]]]
]


(* ::Section:: *)
(*Temporal spherical Motion*)


(* ::Input::Initialization:: *)
tSpherical\[Eta]p[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,t,Gt,Ep,Gcurlt,Ep1,\[Psi]\[Tau],Gcurl\[Theta],\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
Gcurl\[Theta]=-(1/Sqrt[-um a^2])EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um](*eqn 29*);
\[Psi]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[-um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),up/um](*eqn 46*);
Ep1=1/(2 up/um) (EllipticE[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um]-EllipticF[ArcSin[Cos[\[Theta]s]/Sqrt[up]],up/um])(*eqn 32*);
Gcurlt=Re[(2 up)/Sqrt[-um a^2] Ep1] (*eqn 31*);
Ep[\[Tau]_]:=(EllipticE[\[Psi]\[Tau][\[Tau]],up/um]-EllipticF[\[Psi]\[Tau][\[Tau]],up/um])/(2 up/um)(*eqn 32*);
Gt[\[Tau]_]:=Re[-((2 up)/Sqrt[-um a^2])Ep[\[Tau]]-\[Nu]\[Theta] Gcurlt ](*eqn 48*);
assoc=Association["Trajectory"-> t];
t[\[Tau]_]:=((rs^2+a^2)/(rs^2-2M rs+a^2) (rs^2+a^2-a \[Lambda])+a \[Lambda])\[Tau]+a^2 Gt[\[Tau]](*eqn 12*);
tSpherical\[Eta]pFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[tSpherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="tSpherical\[Eta]pFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
tSpherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
tSpherical\[Eta]pFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
tSpherical\[Eta]m[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,t,Gt,\[CapitalUpsilon]\[Tau],Gcurlt,Gcurl\[Theta],\[CapitalUpsilon]s,\[Theta]p,\[Theta]m,h,\[CapitalDelta]\[Theta],up,um,rp,rm,roots,r1,r2,r3,r4,assoc,rootsA,\[Theta]1,\[Theta]2,\[Theta]3,\[Theta]4,consts,\[Nu]\[Theta],\[Nu]r},
\[Nu]\[Theta]=Sign[p\[Theta]s];
\[CapitalDelta]\[Theta]=1/2 (1-(\[Eta]+\[Lambda]^2)/a^2); 
up=\[CapitalDelta]\[Theta]+Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2](*eqn 19*); 
 um=\[CapitalDelta]\[Theta]-Sqrt[\[CapitalDelta]\[Theta]^2+\[Eta]/a^2]; 
rp=M+Sqrt[M^2-a^2];
rm=M-Sqrt[M^2-a^2];
\[Theta]p=ArcCos[Sqrt[um]](*eqn 55*);
\[Theta]m=ArcCos[Sqrt[up]](*eqn 55*);
h=Sign[Cos[\[Theta]s]];
\[CapitalUpsilon]s=ArcSin[Sqrt[(Cos[\[Theta]s]^2-um)/(up-um)]](*eqn 59*);
Gcurl\[Theta]=-h/Sqrt[um a^2] EllipticF[\[CapitalUpsilon]s,1-up/um](*eqn 56*);
Gcurlt=-h Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]s,1-up/um](*eqn 58*);
\[CapitalUpsilon]\[Tau][\[Tau]_]:=JacobiAmplitude[Sqrt[um a^2](\[Tau]+\[Nu]\[Theta] Gcurl\[Theta]),1-up/um](*eqn 66*);
Gt[\[Tau]_]:=Sqrt[um/a^2]EllipticE[\[CapitalUpsilon]\[Tau][\[Tau]],1-up/um]-\[Nu]\[Theta] Gcurlt (*eqn 68*);
assoc=Association["Trajectory"-> t,"ConstantsofMotion"-> consts];
t[\[Tau]_]:=((rs^2+a^2)/(rs^2-2M rs+a^2) (rs^2+a^2-a \[Lambda])+a \[Lambda])\[Tau]+a^2 Gt[\[Tau]](*eqn 12*);
tSpherical\[Eta]mFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[tSpherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="tSpherical\[Eta]mFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
tSpherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
tSpherical\[Eta]mFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Input::Initialization:: *)
tSEquatorial[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,t,assoc,consts,\[Nu]\[Theta],\[Nu]r},
assoc=Association["Trajectory"-> t,"ConstantsofMotion"-> consts];
t[\[Tau]_]:=((rs^2+a^2)/(rs^2-2M rs+a^2) (rs^2+a^2-a \[Lambda])+a \[Lambda])\[Tau]+a^2 \[Tau](*eqn 12*);
tSEquatorialFunction[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s,assoc]
]


(* ::Input::Initialization:: *)
Format[tSEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_]]:="tSEquatorialFunction["<>ToString[a]<>","<>ToString[rs]<>","<>ToString[\[Theta]s]<>","<>ToString[\[Lambda]]<>","<>ToString[\[Eta]]<>","<>ToString[prs]<>","<>ToString[p\[Theta]s]<>",<<>>]";
tSEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][\[Tau]_/;StringQ[\[Tau]]==False]:= assoc["Trajectory"][\[Tau]]
tSEquatorialFunction[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_,assoc_][y_?StringQ] := assoc[y]


(* ::Subsection::Initialization:: *)
(*Temporal Spherical motion general code*)


(* ::Input::Initialization:: *)
TemporalSphericalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{M=1,rp,rm,roots,r1,r2,r3,r4,consts},
Which[\[Eta]>0,Return[tSpherical\[Eta]p[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]==0,Return[tSEquatorial[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],\[Eta]<0,Return[tSpherical\[Eta]m[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]
]


(* ::Section:: *)
(* Temporal motion general code*)


(* ::Input::Initialization:: *)
TemporalMotionGeneral[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=
Module[{roots, r1,r2,r3,r4},
roots=KerrNullGeoRadialRoots[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s];
{r1,r2,r3,r4}=Values[roots];
If[r3==r4,TemporalSphericalMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],TemporalNonSphericalMotion[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]]]


(* ::Input::Initialization:: *)
TemporalMotion[a_,rs_,\[Theta]s_,\[Lambda]_,\[Eta]_,prs_,p\[Theta]s_]:=Module[{consts,M=1},
If[(rs^2+a^2-a \[Lambda])^2-(rs^2-2 M rs+a^2)(\[Eta]+(\[Lambda]-a)^2)>=0,If[If[\[Theta]s==0\[Or] \[Theta]s==\[Pi],\[Eta]+a^2 Cos[\[Theta]s]^2>=0,\[Eta]+a^2 Cos[\[Theta]s]^2-\[Lambda]^2 Cot[\[Theta]s]^2>=0],If[\[Eta]==-a^2 Cos[\[Theta]s]^4\[And]\[Lambda]==a Sin[\[Theta]s]^2,TemporalPNG[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s],TemporalMotionGeneral[a,rs,\[Theta]s,\[Lambda],\[Eta],prs,p\[Theta]s]],Print["Error:Negative Angular Potential"]],Print["Error:Negative Radial Potential"]]
]


(* ::Title::Closed:: *)
(*Close the package*)


End[];

EndPackage[];
