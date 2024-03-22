function [aLong,bLong,aLD,bLD] = Lin(Ref,Phys,D)
% function [aLong,bLong,aLD,bLD] = Lin(Ref,Phys,D)
%
% Lin Takes reference flight conditions,
% aircraft physical characteristics,
% nondimensional derivatives; and returns
% dimensional a and b matrices.
%
% Inputs:
% Ref is a 6-vector of reference flight conditions,
% Ref(1)=density (slugs/ft^3)
% Ref(2)=TAS (ft/s)
% Ref(3)=Mach
% Ref(4)=CL (trim)
% Ref(5)=CD (trim)
% Ref(6)=Gamma, flight path angle (radians)
% Phys is an 9-vector of physical characteristics,
% Phys(1) = Weight (pounds)
% Phys(2) = Ixx (slug-ft^2)
% Phys(3) = Iyy (slug-ft^2)
% Phys(4) = Izz (slug-ft^2)
% Phys(5) = Ixz (slug-ft^2)
% Phys(6) = Area (ft^2)
% Phys(7) = Span (ft)
% Phys(8) = Chord (ft)
% Phys(9) = Thrust angle (radians)
% D is a 29-vector of nondimensional derivatives
% D(1)=CLAlpha
% D(2)=CDAlpha
% D(3)=CmAlpha
% D(4)=CLAlphaDot
% D(5)=CmAlphaDot
% D(6)=CLq
% D(7)=Cmq
% D(8)=CLM
% D(9)=CDM
% D(10)=CmM
% D(11)=CLDeltaM
% D(12)=CDDeltaM
% D(13)=CMDeltaM
% D(14)=CTV
% D(15)=CTDeltaT
% D(16)=CyBeta
% D(17)=ClBeta
% D(18)=CnBeta
% D(19)=Clp
% D(20)=Cnp
% D(21)=Cyp
% D(22)=Clr
% D(23)=Cnr
% D(24)=Cyr
% D(25)=ClDeltaL
% D(26)=CnDeltaL
% D(27)=ClDeltaN
% D(28)=CnDeltaN
% D(29)=CyDeltaN
%
% Outputs
% aLong and aLD are 4x4 longitudinal
% and lateral/directional system matrices
% bLong and bLD are 4x2 longitudinal
% and lateral/directional control matrices
gee=32.174;
mass=Phys(1)/gee;
q=0.5*Ref(1)*Ref(2)^2;
qS=q*Phys(6);
qSb=qS*Phys(7);
qSc=qS*Phys(8);
V=Ref(2);
mV=mass*V;
qSoV=qS/V;
qSoM=qS/mass;
CW=Phys(1)/qS;
CosEps=cos(Phys(9));
SinEps=sin(Phys(9));
CosGam=cos(Ref(6));
SinGam=sin(Ref(6));
CT=CW*sin(Ref(6)-Phys(9))+Ref(5)*CosGam+Ref(4)*SinGam;
aLong(4,1)=0;
aLong(4,2)=0;
aLong(4,3)=1;
aLong(4,4)=0;
Xu=qSoV*(2*CW*SinGam-2*CT*CosEps-Ref(3)*D(9));
Tu=qSoV*(2*CT+D(14));
Tu=0;
aLong(1,1)=Xu+Tu*CosEps;
Xw=qSoV*(Ref(4)-D(2));
aLong(1,2)=Xw;
aLong(1,3)=0;
aLong(1,4)=-Phys(1)*CosGam;
aLong(1,:)=aLong(1,:)/mass;
Zu=-qSoV*(2*Ref(4)+Ref(3)*D(8));
aLong(2,1)=Zu+Tu*SinEps;
Zw=-qSoV*(Ref(5)+D(1));
aLong(2,2)=Zw;
Zq=-qSc*D(6)/(2*V);
aLong(2,3)=Zq+mass*V;
aLong(2,4)=-Phys(1)*SinGam;
m_ZwDot=mass+qSc*D(4)/(2*V^2);
aLong(2,:)=aLong(2,:)/m_ZwDot;
MwDot=qSc*Phys(8)*D(5)/(2*V^2);
Mu=Ref(3)*qSc*D(10)/V;
aLong(3,1)=Mu+MwDot*aLong(2,1);
Mw=qSc*D(3)/V;
aLong(3,2)=Mw+MwDot*aLong(2,2);
Mq=qSc*Phys(8)*D(7)/(2*V);
aLong(3,3)=Mq+MwDot*aLong(2,3);
aLong(3,4)=MwDot*aLong(2,4);
aLong(3,:)=aLong(3,:)/Phys(3);
bLong(1,1)=qSoM*D(15)*CosEps;
bLong(1,2)=-qSoM*D(12);
bLong(2,1)=qS*D(15)*SinEps/m_ZwDot;
bLong(2,2)=-qSoM*D(11)/m_ZwDot;
bLong(3,1)=0;
MDeltaM=qSc*D(13);
bLong(3,2)=(MDeltaM+MwDot*bLong(2,2))/Phys(3);
bLong(4,1)=0;
bLong(4,2)=0;
qSboV=qSoV*Phys(7);
Yv=qSoV*D(16);
Yp=qSboV*D(21)/2;
Yr=qSboV*D(24)/2;
aLD(1,1)=Yv;
aLD(1,2)=Yp;
aLD(1,3)=Yr-mass*Ref(2);
aLD(1,4)=Phys(1)*CosGam;
aLD(1,:)=aLD(1,:)/mass;
Lv=qSboV*D(17);
Lp=qSboV*Phys(7)*D(19)/2;
Lr=qSboV*Phys(7)*D(22)/2;
Nv=qSboV*D(18);
Np=qSboV*Phys(7)*D(20)/2;
Nr=qSboV*Phys(7)*D(23)/2;
aLD(2,1)=Phys(4)*Lv+Phys(5)*Nv;
aLD(2,2)=Phys(4)*Lp+Phys(5)*Np;
aLD(2,3)=Phys(4)*Lr+Phys(5)*Nr;
aLD(2,4)=0;
aLD(3,1)=Phys(5)*Lv+Phys(2)*Nv;
aLD(3,2)=Phys(5)*Lp+Phys(2)*Np;
aLD(3,3)=Phys(5)*Lr+Phys(2)*Nr;
aLD(3,4)=0;
aLD(2:3,:)=aLD(2:3,:)/(Phys(2)*Phys(4)-Phys(5)^2);
aLD(4,1)=0;
aLD(4,2)=1;
aLD(4,3)=SinGam/CosGam;
aLD(4,4)=0;
bLD(1,1)=0;
bLD(1,2)=qSoM*D(29);
LdL=qSb*D(25);
NdL=qSb*D(26);
LdN=qSb*D(27);
NdN=qSb*D(28);
bLD(2,1)=Phys(4)*LdL+Phys(5)*NdL;
bLD(2,2)=Phys(4)*LdN+Phys(5)*NdN;
bLD(3,1)=Phys(5)*LdL+Phys(2)*NdL;
bLD(3,2)=Phys(5)*LdN+Phys(2)*NdN;
bLD(2:3,1:2)=bLD(2:3,1:2)/(Phys(2)*Phys(4)-Phys(5)^2);
bLD(4,1)=0;
bLD(4,2)=0;
