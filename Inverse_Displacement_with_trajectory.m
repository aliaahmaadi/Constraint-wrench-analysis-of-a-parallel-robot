clc
clear

%% given values
k=[0;0;1];
pdat=[0.1;0.1;0.1];
phidat=0.1;
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
aII6=l44/3;
aII7=l5/2;
aI5=l4/3;
a=0.12;
x=0.300;
y0=0;
zp=0.8422;
z0=zp;

iter=2000;
t=linspace(0,5.8,iter);
time=t;

%% Method
for i=1:iter
    if t(i)<0.2
        z=10*t(i)+z0;
        y=0+y0;
    elseif t(i)>=0.2 && t(i)<0.4
        z=-12.5*(t(i)-0.2)^2+5*(t(i)-0.2)+2+z0;
        y=12.5*(t(i)-0.2)^2+y0;
    elseif t(i)>=0.4 && t(i)<2.2
        z=2.5+z0;
        y=5*(t(i)-0.4)+0.5+y0;
    elseif t(i)>=2.2 && t(i)<2.4
        z=-12.5*(t(i)-2.2)^2+2.5+z0;
        y=-12.5*(t(i)-2.2)^2+5*(t(i)-2.2)+9.5+y0;
    elseif t(i)>=2.4 && t(i)<2.8
        z=-5*(t(i)-2.4)+2+z0;
        y=10+y0;
    elseif t(i)>=2.8 && t(i)<3.2
        y=10+y0;
        z=5*(t(i)-2.8)+z0;
    elseif t(i)>=3.2 && t(i)<3.4
        z=-12.5*(t(i)-3.2)^2+5*(t(i)-3.2)+2+z0;
        y=-12.5*(t(i)-3.2)^2+10+y0;
    elseif t(i)>=3.4 && t(i)<5.2
        z=2.5+z0;
        y=-5*(t(i)-3.4)+9.5+y0;
    elseif t(i)>=5.2 && t(i)<5.4
        z=-12.5*(t(i)-5.2)^2+2.5+z0;
        y= 12.5*(t(i)-5.2)^2-5*(t(i)-5.2)+0.5+y0;
    elseif t(i)>=5.4 && t(i)<=5.6
        y=y0;
        z=-10*(t(i)-5.4)+2+z0;
    end

%phi=-0.0969*(t(i))^3 + 0.7252*(t(i))^2 - 0.5881*t(i) - 1.5478;
%phi=phi*180/pi;
%phidat=3*(-0.0969)*(t(i))^2+2*(t(i))*0.7252-0.5881;    
    
tet_I_1(i)=atand(y/x);
tet_I_4(i)=tet_I_1(i)+phi;
N=-2*zp*l2;
rI=sqrt(x^2+y^2);
rII=sqrt((x+l5*cosd(phi)-l0)^2+(y-l5*sind(phi))^2);
 
QI=rI^2+zp^2+l1^2+l2^2-l3^2+2*l4*l1-2*rI*l1-2*l4*rI+l4^2;
QII=rII^2+zp^2+l1^2+l2^2-l3^2+2*l44*l1-2*rII*l1-2*l44*rII+l44^2;
MI=2*l1*l2+2*l4*l2-2*rI*l2;
MII=2*l1*l2+2*l44*l2-2*rII*l2;
tet_I_2(i)=2*atand(-(N-sqrt(N^2-QI^2+MI^2))/(QI-MI));
tet_II_2(i)=2*atand(-(N-sqrt(N^2-QII^2+MII^2))/(QII-MII));
tet_I_3(i)=asind((zp-l2*sind(tet_I_2(i)))/l3);
tet_II_3(i)=asind((zp-l2*sind(tet_II_2(i)))/l3);
tet_II_1(i)=atand((y-l5*sind(phi))/(x+l5*cos(phi)-l0));
tet_II_4(i)=phi+tet_II_1(i)-180;
%%  ~~~~ joint rates ~~~~
f1=[cosd(tet_I_2(i)+90);sind(tet_I_2(i)+90);0];
f2=[cos(tet_II_2(i)+90);sind(tet_II_2(i)+90);0];
aI1=l1*[cosd(tet_I_1(i));sind(tet_I_1(i));0];
aI2=l2*[cosd(tet_I_2(i))*cosd(tet_I_1(i));cosd(tet_I_2(i))*sind(tet_I_1(i));sind(tet_I_2(i))];
aI3=l3*[cosd(tet_I_3(i))*cosd(tet_I_1(i));cosd(tet_I_3(i))*sind(tet_I_1(i));sind(tet_I_3(i))];
aI4=l4*[cosd(tet_I_1(i));sind(tet_I_1(i));0];

aII1=l1*[cosd(tet_II_1(i));sind(tet_II_1(i));0];
aII2=l2*[cosd(tet_II_2(i))*cosd(tet_II_1(i));cosd(tet_II_2(i))*sind(tet_II_1(i));sind(tet_II_2(i))];
aII3=l3*[cosd(tet_II_3(i))*cosd(tet_II_1(i));cosd(tet_II_3(i))*sind(tet_II_1(i));sind(tet_II_3(i))];
aII4=l44*[cosd(tet_II_1(i));sind(tet_II_1(i));0];
aII5=l5*[-cosd(phi);sind(phi);0];
PI=aI1+aI2+aI3+aI4;
PII=aII1+aII2+aII3+aII4+l0*[1;0;0];
RI14=aI1+aI2+aI3+aI4;
RII4=aII1+aII2+aII3+aII4;
RI23=aI2+aI3;
RII23=aII2+aII3;
v1=cross(cross(f1,aI3),cross(k,RI14));
v2=cross(cross(f2,aII3),cross(k,RII4));
p=[pdat;phidat];
DeltaI=cross(cross(f1,aI3),cross(f1,aI2));
DeltaII=cross(cross(f2,aII3),cross(f2,aII2));
B=[ v1,DeltaI,zeros(3,1),zeros(3,1);zeros(3,1),zeros(3,1),v2,DeltaII];
sayI= cross(f1,aI3);
sayII= cross(f2,aII3);
SAYI= [ 0,-sayI(3),sayI(2);sayI(3),0,-sayI(1);-sayI(2),sayI(1),0];
SAYII= [ 0,-sayII(3),sayII(2);sayII(3),0,-sayII(1);-sayII(2),sayII(1),0];

A=[ SAYI ,zeros(3,1) ; SAYII, SAYII*cross(k,aII5)];
tetadat= pinv(B)*A*p;
tet_I1_dot(i)=tetadat(1);
tet_I2_dot(i)=tetadat(2);
tet_II1_dot(i)=tetadat(3);
tet_II2_dot(i)=tetadat(4);

psi1=cross(f1,aI3);
psi2=cross(f1,aII3);
C=[zeros(3,1) -psi1 -cross(k,aII5) psi2 ;-1 0 1 0];
g=RII4+aII5;
D=[cross(k,RI14) cross(f1,aI2) -cross(k,g) -cross(f2,aII5) ;1 0 1 0];
DL=(pinv(D));
thetadat_P=DL*C*tetadat;
tet_I3_dot(i)=double(thetadat_P(1));
tet_I4_dot(i)=double(thetadat_P(2));
tet_II3_dot(i)=double(thetadat_P(3));
tet_II4_dot(i)=double(thetadat_P(4));
  
yy(i)=y;
zpp(i)=z;
fprintf(' Calculated for z= %d & y= %d. \n',z,y)    

%%
if i>1
    tdot11_test(i)=0.01*((tet_I_1(i)-tet_I_1(i-1))/(t(i)-t(i-1)));
    tdot12_test(i)=0.01*((tet_I_2(i)-tet_I_2(i-1))/(t(i)-t(i-1)));
    tdot13_test(i)=0.01*((tet_I_3(i)-tet_I_3(i-1))/(t(i)-t(i-1)));
    tdot14_test(i)=0.01*((tet_I_4(i)-tet_I_4(i-1))/(t(i)-t(i-1)));
    tdot21_test(i)=0.01*((tet_II_1(i)-tet_II_1(i-1))/(t(i)-t(i-1)));
    tdot22_test(i)=0.01*((tet_II_2(i)-tet_II_2(i-1))/(t(i)-t(i-1)));
    tdot23_test(i)=0.01*((tet_II_3(i)-tet_II_3(i-1))/(t(i)-t(i-1)));
    tdot24_test(i)=0.01*((tet_II_4(i)-tet_II_4(i-1))/(t(i)-t(i-1)));
end

if i>1
    tet_I1_2dot(i)=0.01*(( tdot11_test(i)- tdot11_test(i-1))/(t(i)-t(i-1)));
    tet_I2_2dot(i)=0.01*(( tdot12_test(i)- tdot12_test(i-1))/(t(i)-t(i-1)));
    tet_I3_2dot(i)=0.01*(( tdot13_test(i)- tdot13_test(i-1))/(t(i)-t(i-1)));
    tet_I4_2dot(i)=0.01*(( tdot14_test(i)- tdot14_test(i-1))/(t(i)-t(i-1)));
    tet_II1_2dot(i)=0.01*(( tdot21_test(i)- tdot21_test(i-1))/(t(i)-t(i-1)));
    tet_II2_2dot(i)=0.01*(( tdot22_test(i)- tdot22_test(i-1))/(t(i)-t(i-1)));
    tet_II3_2dot(i)=0.01*(( tdot23_test(i)- tdot23_test(i-1))/(t(i)-t(i-1)));
    tet_II4_2dot(i)=0.01*(( tdot24_test(i)- tdot24_test(i-1))/(t(i)-t(i-1)));
end 

end
%%

%% Plots for displacement
figure 
plot(time,real(tet_II_1), 'LineWidth' , 2)
title('\theta_I_I_1')
ylabel('\theta_I_I_1 (deg)')
xlabel('time (s)')
title('Plot for \theta_I_I_1 ')
axis([0,5.8 -180 180])
grid on

figure
plot(time,real(tet_II_2), 'LineWidth' , 2)
title('\theta_I_I_2')
ylabel('\theta_I_2 (deg)')
xlabel('time (s)')
title('Plot for \theta_I_I_2 ')
axis([0,5.8 -180 180])
grid on

figure
plot(time,real(tet_II_3), 'LineWidth' , 2)
title('\theta_I_I_3')
ylabel('\theta_I_I_3 (deg)')
title('Plot for \theta_I_I_3 ')
xlabel('time (s)')
axis([0,5.8 -180 180])
grid on

figure
plot(time,real(tet_I_1), 'LineWidth' , 2)
ylabel('\theta_I_1 (deg)')
title('Plot for \theta_I_1 ')
xlabel('time (s)')
title('\theta_I_1')
axis([0,5.8 -180 180])
grid on

figure
plot(time,real(tet_I_2), 'LineWidth' , 2)
ylabel('\theta_I_2 (deg)')
title('Plot for \theta_I_2 ')
xlabel('time (s)')
title('\theta_I_2')
axis([0,5.8 -180 180])
grid on

figure
plot(time,real(tet_I_3), 'LineWidth' , 2)
ylabel('\theta_I_I_3 (deg)')
title('Plot for \theta_I_3 ')
xlabel('time (s)')
title('\theta_I_3')
axis([0,5.8 -180 180])
grid on
%% Plots FOR J_Rates..
figure
plot(time,(tet_II1_dot), 'LineWidth' , 2)
title('$\dot{\theta_{II1}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_I1_dot), 'LineWidth' , 2)
title('$\dot{\theta_{I1}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_II2_dot), 'LineWidth' , 2)
title('$\dot{\theta_{II2}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_II3_dot), 'LineWidth' , 2)
title('$\dot{\theta_{II3}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_I2_dot), 'LineWidth' , 2)
title('$\dot{\theta_{I2}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_I3_dot), 'LineWidth' , 2)
title('$\dot{\theta_{I3}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_I4_dot), 'LineWidth' , 2)
title('$\dot{\theta_{I4}}$', 'Interpreter','latex')
grid on

figure
plot(time,(tet_II4_dot), 'LineWidth' , 2)
title('$\dot{\theta_{II4}}$', 'Interpreter','latex')
grid on
close all
%% ALL
plot(time,real(tet_II_1), 'LineWidth' , 2)
hold on
plot(time,real(tet_II_2), 'LineWidth' , 2)
hold on
plot(time,real(tet_II_3), 'LineWidth' , 2)
hold on
plot(time,real(tet_I_1), 'LineWidth' , 2)
hold on
plot(time,real(tet_I_2), 'LineWidth' , 2)
hold on
plot(time,real(tet_I_3), 'LineWidth' , 2)
ylabel('\theta (deg)')
xlabel('time (s)')
title('\theta')
axis([0,5.8 -180 180])
title('Plot for all \theta ')
grid on
legend('\theta_{II1}','\theta_{II2}','\theta_{II3}','\theta_{I1}','\theta_{I2}','\theta_{I3}','Location', 'Best')


%% 2DOTS
figure
plot(time,(tet_II1_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{II1}}$', 'Interpreter','latex')
title('The dot  $\dot{\theta_{II1}}$ with finite diff ')
xlim([0 5.6])
grid on
figure
plot(time,(tet_I1_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{I1}}$', 'Interpreter','latex')
title('The dot  $\dot{\theta_{I1}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tet_I2_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{I2}}$', 'Interpreter','latex')
title('The dot  $\dot{\theta_{I2}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tet_II2_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{II2}}$', 'Interpreter','latex')
title('The  dot $\dot{\theta_{II2}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tet_II3_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{II3}}$', 'Interpreter','latex')
title('The dot  $\dot{\theta_{II3}}$ with finite diff ')
xlim([0 5.6])
grid on


figure
plot(time,(tet_I3_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{I3}}$', 'Interpreter','latex')
title('The dot $\dot{\theta_{I3}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tet_I4_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{I4}}$', 'Interpreter','latex')
title('The dot $\dot{\theta_{I4}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tet_II4_2dot), 'LineWidth' , 2)
title('$\dot{\theta_{II4}}$', 'Interpreter','latex')
title('The dot  $\dot{\theta_{I4}}$ with finite diff ')

xlim([0 5.6])
grid on
%% DOTS
figure
plot(time,(tdot21_test), 'LineWidth' , 2)
title('$\dot{\theta_{II1}}$', 'Interpreter','latex')
title('The $\dot{\theta_{II1}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot11_test), 'LineWidth' , 2)
title('$\dot{\theta_{I1}}$', 'Interpreter','latex')
title('The $\dot{\theta_{I1}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot22_test), 'LineWidth' , 2)
title('$\dot{\theta_{II2}}$', 'Interpreter','latex')
title('The $\dot{\theta_{II2}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot23_test), 'LineWidth' , 2)
title('$\dot{\theta_{II3}}$', 'Interpreter','latex')
title('The $\dot{\theta_{II3}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot12_test), 'LineWidth' , 2)
title('$\dot{\theta_{I2}}$', 'Interpreter','latex')
title('The $\dot{\theta_{I2}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot13_test), 'LineWidth' , 2)
title('$\dot{\theta_{I3}}$', 'Interpreter','latex')
title('The $\dot{\theta_{I3}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot14_test), 'LineWidth' , 2)
title('$\dot{\theta_{I4}}$', 'Interpreter','latex')
title('The $\dot{\theta_{I4}}$ with finite diff ')
xlim([0 5.6])
grid on

figure
plot(time,(tdot24_test), 'LineWidth' , 2)
title('$\dot{\theta_{II4}}$', 'Interpreter','latex')
title('The $\dot{\theta_{II4}}$ with finite diff ')

xlim([0 5.6])
grid on

%% Plot for trajectory

figure 
plot(yy,zpp,'r', 'LineWidth' , 2)
title('The Test Trajectory')
axis([-0.5,10.5 0 4])
ylabel('z (cm)')
xlabel('y (cm)')
legend('The Test Trajectory line' )
