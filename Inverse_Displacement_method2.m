clc
clear
close all;
 
%% ~~~Input~~~
k=[0;0;1];
 
%tetI1=0;
%tetII1=pi;
%tetI2=0.53*pi;
%tetII2=0.53*pi;
%tetI3=pi/4;
%tetII3=pi/4;
%% ~~~~~dimensions of the current prototype of the McGill SMG~~~~ 
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
%% ~~~given coordinates~~
x=0.300;
y=0;
zp=0.8422;
%% 

%% ~~~Limb I II~~~
tetI1=atand(y/x);
tetI4=tetI1+phi;
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
tet_II_1=atand((y-l5*sind(phi))/(x+l5*cos(phi)-l0));
%% ~~~Results~~
fprintf(' theta_I_1 is %d (deg) \n',tetI1)
fprintf(' theta_I_2 is %d (deg) \n',tetI2)
fprintf(' theta_I_3 is %d (deg) \n',tetI3)
fprintf(' theta_I_4 is %d (deg) \n',tetI4)
fprintf(' theta_II_1 is %d (deg) \n',tet_II_1)
fprintf(' theta_II_2 is %d (deg) \n',tetII2)
fprintf(' theta_II_3 is %d (deg) \n',tetII3)
