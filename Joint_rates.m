clc
clear
close all;
 
%% Input
k=[0;0;1];
 
 
%% ~--------------Dimensions of the current prototype of the McGill SMG
 
phi=0;
l0=0.6;
l1=0.06;
l2=0.6329;
l3=0.3;
l44=0.0621;
l5=0.0254;
l4=l44+l5; % (5.6)thesis
h=0.03;
L=0.01987;
pdat=[0.2;0.1;0.2];
phidat=0.3;
%% ~~--------------------given coordinates--------------
x=0.300;
y=0;
zp=0.8422;
%% 
 
%% Limb I II inverse
tetI1=atand(y/x);
tetI4=tetI1+phi;
tetII1=atand((y-l5*sind(phi))/(x+l5*cos(phi)-l0));
 
N=-2*zp*l2;
rI=sqrt(x^2+y^2);
rII=sqrt((x+l5*cosd(phi)-l0)^2+(y-l5*sind(phi))^2);
 
QI=rI^2+zp^2+l1^2+l2^2-l3^2+2*l4*l1-2*rI*l1-2*l4*rI+l4^2;
QII=rII^2+zp^2+l1^2+l2^2-l3^2+2*l44*l1-2*rII*l1-2*l44*rII+l44^2;
MI=2*l1*l2+2*l4*l2-2*rI*l2;
MII=2*l1*l2+2*l44*l2-2*rII*l2;
tetI2=2*atand(-(N-sqrt(N^2-QI^2+MI^2))/(QI-MI));
tetII2=2*atand(-(N-sqrt(N^2-QII^2+MII^2))/(QII-MII));
tetI3=asind((zp-l2*sind(tetI2))/l3);
tetII3=asind((zp-l2*sind(tetII2))/l3);
%% Method for joint rates
f1=[cos(tetI2+pi/2);sin(tetI2+pi/2);0];
f2=[cos(tetII2+pi/2);sin(tetII2+pi/2);0];
aI1=l1*[cos(tetI1);sin(tetI1);0];
aI2=l2*[cos(tetI2)*cos(tetI1);cos(tetI2)*sin(tetI1);sin(tetI2)];
aI3=l3*[cos(tetI3)*cos(tetI1);cos(tetI3)*sin(tetI1);sin(tetI3)];
aI4=l4*[cos(tetI1);sin(tetI1);0];
 
aII1=l1*[cos(tetII1);sin(tetII1);0];
aII2=l2*[cos(tetII2)*cos(tetII1);cos(tetII2)*sin(tetII1);sin(tetII2)];
aII3=l3*[cos(tetII3)*cos(tetII1);cos(tetII3)*sin(tetII1);sin(tetII3)];
aII4=l44*[cos(tetII1);sin(tetII1);0];
aII5=l5*[-cos(phi);sin(phi);0];
PI=aI1+aI2+aI3+aI4;
PII=aII1+aII2+aII3+aII4+l0*[1;0;0];
RII4=aII1+aII2+aII3+aII4;
RI14=aI1+aI2+aI3+aI4;
RII14=aII1+aII2+aII3+aII4;
RI23=aI2+aI3;
RII23=aII2+aII3;
v1=cross(cross(f1,aI3),cross(k,RI14));
v2=cross(cross(f2,aII3),cross(k,RII14));
DeltaI=cross(cross(f1,aI3),cross(f1,aI2));
DeltaII=cross(cross(f2,aII3),cross(f2,aII2));
B=[ v1,DeltaI,zeros(3,1),zeros(3,1);zeros(3,1),zeros(3,1),v2,DeltaII];
sayI= cross(f1,aI3);
sayII= cross(f2,aII3);
SAYI= [ 0,-sayI(3),sayI(2);sayI(3),0,-sayI(1);-sayI(2),sayI(1),0];
SAYII= [ 0,-sayII(3),sayII(2);sayII(3),0,-sayII(1);-sayII(2),sayII(1),0];

A=[ SAYI ,zeros(3,1) ; SAYII, SAYII*cross(k,aII5)];
p=[pdat;phidat];
thetadat= pinv(B)*A*p;
tet_dot_11=thetadat(1);
tet_dot_12=thetadat(2);
tet_dot_21=thetadat(3);
tet_dot_22=thetadat(4);
%% Results
fprintf(' /theta_{dot}_I1 is %d \n',double(thetadat(1)))
fprintf(' /theta_{dot}_I2 is %d \n',double(thetadat(2)))
fprintf(' /theta_{dot}_II1 is %d \n',double(thetadat(3)))
fprintf(' /theta_{dot}_II2 is %d \n',double(thetadat(4)))