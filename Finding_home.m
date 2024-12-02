clc
clear
close all;
 
k=[0;0;1];

%% Given angles
tetI1=0;
tetII1=pi;
tetI2=0.53*pi;
tetII2=0.53*pi;
tetI3=pi/4;
tetII3=pi/4;
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
 
%% Geometry
 
%-------------------Eq.5----------------------
%phi=tetII4-tetII4+180;
%phi=tetI4-tetI1;   % Eq.5
 
%-------------------Eq.1----------------------
aI1=l1*[cos(tetI1);sin(tetI1);0];
aI2=l2*[cos(tetI2)*cos(tetI1);cos(tetI2)*sin(tetI1);sin(tetI2)];
aI3=l3*[cos(tetI3)*cos(tetI1);cos(tetI3)*sin(tetI1);sin(tetI3)];
aI4=l4*[cos(tetI1);sin(tetI1);0];
 
aII1=l1*[cos(tetII1);sin(tetII1);0];
aII2=l2*[cos(tetII2)*cos(tetII1);cos(tetII2)*sin(tetII1);sin(tetII2)];
aII3=l3*[cos(tetII3)*cos(tetII1);cos(tetII3)*sin(tetII1);sin(tetII3)];
aII4=l44*[cos(tetII1);sin(tetII1);0];
aII5=l5*[-cos(phi);sin(phi);0];
 
 
%-------------------Eq.1----------------------
PI=aI1+aI2+aI3+aI4;
PII=aII1+aII2+aII3+aII4+l0*[1;0;0];
 
 
%-------------------Eq.4----------------------
xII=(l1+l2*cos(tetII2)+l3*cos(tetII3)+l44)*cos(tetII1)-l5*cos(phi)+l0;
yII=(l1+l2*cos(tetII2)+l3*cos(tetII3)+l44)*sin(tetII1)+l5*sin(phi);
zpII=l2*sin(tetII2)+l3*sin(tetII3);
 
%-------------------Eq.3----------------------
xI=(l1+l2*cos(tetI2)+l3*cos(tetI3)+l4)*cos(tetI1);
yI=(l1+l2*cos(tetI2)+l3*cos(tetI3)+l4)*sin(tetI1);
zpI=l2*sin(tetI2)+l3*sin(tetI3);

%% Results
fprintf(' The x of home posture is %d. \n',xI);
fprintf(' The y of home posture is %d. \n',yI);
fprintf(' The z of home posture is %d. \n',zpI);

