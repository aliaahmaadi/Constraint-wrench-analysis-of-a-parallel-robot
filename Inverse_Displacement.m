clc
clear
close all; 
%% Input
k=[0;0;1];
%tetI1=0;
%tetII1=pi;
%tetI2=0.53*pi;
%tetII2=0.53*pi;
%tetI3=pi/4;
%tetII3=pi/4;
%--------------dimensions of the current prototype of the McGill SMG
phi=0;
l0=0.6;
l1=0.06;
l2=0.6329;
l3=0.3;
l44=0.0621;
l5=0.0254;
l4=l44+l5; 
h=0.03;
L=0.01987;
%--------------------given coordinates--------------
x=0.300;
y=0;
zp=0.8422;
%% 
syms tetI1 tetII1 tetI2 tetII2 tetI3 tetII3
%% Limb I
eq(1)=x==(l1+l2*cosd(tetI2)+l3*cosd(tetI3)+l4)*cosd(tetI1);
eq(2)=y==(l1+l2*cosd(tetI2)+l3*cosd(tetI3)+l4)*sind(tetI1);
eq(3)=zp==l2*sind(tetI2)+l3*sind(tetI3);
s=(solve(eq));
double(s.tetI1(1));
double(s.tetI2(1));
double(s.tetI3(1));

fprintf(' theta_I_1 is %d (deg) \n',double(s.tetI1(1)))
fprintf(' theta_I_2 is %d (deg) \n',double(s.tetI2(1)))
fprintf(' theta_I_3 is %d (deg) \n',double(s.tetI3(1)))
%% Limb II
eqII(1)=x==(l1+l2*cosd(tetII2)+l3*cosd(tetII3)+l44)*cosd(tetII1)-l5*cosd(phi)+l0;
eqII(2)=y==(l1+l2*cosd(tetII2)+l3*cosd(tetII3)+l44)*sind(tetII1)+l5*sind(phi);
eqII(3)=zp==l2*sind(tetII2)+l3*sind(tetII3);
sII=(solve(eqII));
double(sII.tetII1(1));
double(sII.tetII2(1));
double(sII.tetII3(1));
fprintf(' theta_II_1 is %d (deg) \n',double(sII.tetII1(1)))
fprintf(' theta_II_2 is %d (deg) \n',double(sII.tetII2(1)))
fprintf(' theta_II_3 is %d (deg) \n',double(sII.tetII3(1)))