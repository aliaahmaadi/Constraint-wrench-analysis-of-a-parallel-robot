clc
clear
close all
%% 
k=[0;0;1];
pdat=[0;0;0];
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

iter=10000;
t=linspace(0,5.8,iter);
time=t;
fprintf(' Number of iteration is %d.',iter)
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

phi=-0.0969*(t(i))^3 + 0.7252*(t(i))^2 - 0.5881*t(i) - 1.5478;
phi=phi*180/pi;
phidat=3*(-0.0969)*(t(i))^2+2*(t(i))*0.7252-0.5881;   
yy(i)=y;
zpp(i)=z;
if i>1
pdat=[0;(yy(i)-yy(i-1))/(t(i)-t(i-1));(zpp(i)-zpp(i-1))/(t(i)-t(i-1))];
ypp(i)=(yy(i)-yy(i-1))/(t(i)-t(i-1));
zppp(i)=(zpp(i)-zpp(i-1))/(t(i)-t(i-1));
end

tet_I1(i)=atand(y/x);
tet_I4(i)=tet_I1(i)+phi;
N=-2*zp*l2;
rI=sqrt(x^2+y^2);
rII=sqrt((x+l5*cosd(phi)-l0)^2+(y-l5*sind(phi))^2);
 
QI=rI^2+zp^2+l1^2+l2^2-l3^2+2*l4*l1-2*rI*l1-2*l4*rI+l4^2;
QII=rII^2+zp^2+l1^2+l2^2-l3^2+2*l44*l1-2*rII*l1-2*l44*rII+l44^2;
MI=2*l1*l2+2*l4*l2-2*rI*l2;
MII=2*l1*l2+2*l44*l2-2*rII*l2;
tet_I2(i)=2*atand(-(N-sqrt(N^2-QI^2+MI^2))/(QI-MI));
tet_II2(i)=2*atand(-(N-sqrt(N^2-QII^2+MII^2))/(QII-MII));
tet_I3(i)=asind((zp-l2*sind(tet_I2(i)))/l3);
tet_II3(i)=asind((zp-l2*sind(tet_II2(i)))/l3);
tet_II1(i)=atand((y-l5*sind(phi))/(x+l5*cos(phi)-l0));
tet_II4(i)=phi+tet_II1(i)-180;

%%  ~~~~ joint rates ~~~~
f1=[cosd(tet_I2(i)+90);sind(tet_I2(i)+90);0];
f2=[cos(tet_II2(i)+90);sind(tet_II2(i)+90);0];
aI1=l1*[cosd(tet_I1(i));sind(tet_I1(i));0];
aI2=l2*[cosd(tet_I2(i))*cosd(tet_I1(i));cosd(tet_I2(i))*sind(tet_I1(i));sind(tet_I2(i))];
aI3=l3*[cosd(tet_I3(i))*cosd(tet_I1(i));cosd(tet_I3(i))*sind(tet_I1(i));sind(tet_I3(i))];
aI4=l4*[cosd(tet_I1(i));sind(tet_I1(i));0];

aII1=l1*[cosd(tet_II1(i));sind(tet_II1(i));0];
aII2=l2*[cosd(tet_II2(i))*cosd(tet_II1(i));cosd(tet_II2(i))*sind(tet_II1(i));sind(tet_II2(i))];
aII3=l3*[cosd(tet_II3(i))*cosd(tet_II1(i));cosd(tet_II3(i))*sind(tet_II1(i));sind(tet_II3(i))];
aII4=l44*[cosd(tet_II1(i));sind(tet_II1(i));0];
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

%% saving theta 
tI1(i)= tet_I1(i);tI2(i)= tet_I2(i);tI3(i)= tet_I3(i);tI4(i)= tet_I4(i);
tII1(i)= tet_II1(i);tII2(i)= tet_II2(i);tII3(i)= tet_II3(i);tII4(i)= tet_II4(i);

if i>1
    tet_I1_dot(i)=((tet_I1(i)-tet_I1(i-1))/(t(i)-t(i-1)));
    tet_I2_dot(i)=((tet_I2(i)-tet_I2(i-1))/(t(i)-t(i-1)));
    tet_I3_dot(i)=((tet_II2(i)-tet_II2(i-1))/(t(i)-t(i-1)));
    tet_I4_dot(i)=((tet_I4(i)-tet_I4(i-1))/(t(i)-t(i-1)));
    tet_II1_dot(i)=((tet_II1(i)-tet_II1(i-1))/(t(i)-t(i-1)));
    tet_II2_dot(i)=((tet_II2(i)-tet_II2(i-1))/(t(i)-t(i-1)));
    tet_II3_dot(i)=((tet_II3(i)-tet_II3(i-1))/(t(i)-t(i-1)));
    tet_II4_dot(i)=((tet_II4(i)-tet_II4(i-1))/(t(i)-t(i-1)));
end

%% saving theta dot 

tI1d(i)= tet_I1_dot(i);tI2d(i)= tet_I2_dot(i);tI3d(i)= tet_I3_dot(i);tI4d(i)= tet_I4_dot(i);
tII1d(i)= tet_II1_dot(i);tII2d(i)= tet_II2_dot(i);tII3d(i)= tet_II3_dot(i);tII4d(i)= tet_II4_dot(i);
fprintf(' Calculated for z= %d & y= %d. \n',z,y)

%% 
tet_I1_2dot(1)=0;tet_I2_2dot(1)=0;tet_I3_2dot(1)=0;tet_I4_2dot(1)=0;
tet_II1_2dot(1)=0;tet_II2_2dot(1)=0;tet_II3_2dot(1)=0;tet_II4_2dot(1)=0;
if i>1
    tet_I1_2dot(i)=(( tI1d(i)- tI1d(i-1))/(t(i)-t(i-1)));
    tet_I2_2dot(i)=(( tI2d(i)- tI2d(i-1))/(t(i)-t(i-1)));
    tet_I3_2dot(i)=(( tI3d(i)- tI3d(i-1))/(t(i)-t(i-1)));
    tet_I4_2dot(i)=(( tI4d(i)- tI4d(i-1))/(t(i)-t(i-1)));
    tet_II1_2dot(i)=(( tII1d(i)- tII1d(i-1))/(t(i)-t(i-1)));
    tet_II2_2dot(i)=(( tII2d(i)- tII2d(i-1))/(t(i)-t(i-1)));
    tet_II3_2dot(i)=(( tII3d(i)- tII3d(i-1))/(t(i)-t(i-1)));
    tet_II4_2dot(i)=(( tII4d(i)- tII4d(i-1))/(t(i)-t(i-1)));
end 
qd=[tet_I1_dot(i);tet_I2_dot(i);tet_I3_dot(i);tet_I4_dot(i);tet_II1_dot(i);tet_II2_dot(i);tet_II3_dot(i);tet_II4_dot(i)];
qdd=[tet_I1_2dot(i);tet_I2_2dot(i);tet_I3_2dot(i);tet_I4_2dot(i);tet_II1_2dot(i);tet_II2_2dot(i);tet_II3_2dot(i);tet_II4_2dot(i)];
%% ~~~~~~~~DERIVATION OF THE TWIST SHAPING MATRICES~~~~~~~~~~



%--------------------------for T1
T1=[zeros(2,1),zeros(2,7);1,zeros(7,1).';zeros(3,1),zeros(3,7)];
T1_dot=zeros(6,8);

%--------------------------for T2
t41=-l1*sind(tet_I1(i))-l2*sind(tet_I1(i))*sind(tet_I2(i))+0.47*a*sind(tet_I1(i));
t42=l2*cosd(tet_I1(i))*cosd(tet_I2(i));
t51=l1*cosd(tet_I1(i))+l2*cosd(tet_I1(i))*sind(tet_I2(i))-0.47*a*cosd(tet_I1(i));
t52=l2*sind(tet_I1(i))*cosd(tet_I2(i));
t62=-l2*sind(tet_I2(i));

t41_dot=-l1*tet_I1_dot(i)*cosd(tet_I1(i))-l2*tet_I1_dot(i)*cosd(tet_I1(i))*sind(tet_I2(i))-l2*tet_I2_dot(i)*sind(tet_I1(i))*cosd(tet_I2(i))+0.47*a*tet_I1_dot(i)*cosd(tet_I1(i));
t42_dot=-l2*tet_I1_dot(i)*sind(tet_I1(i))*cosd(tet_I2(i))-l2*tet_I2_dot(i)*cosd(tet_I1(i))*sind(tet_I2(i));
t51_dot=-l1*tet_I1_dot(i)*sind(tet_I1(i))-l2*tet_I1_dot(i)*sind(tet_I1(i))*sind(tet_I2(i))+l2*tet_I2_dot(i)*cosd(tet_I1(i))*cosd(tet_I2(i))+0.47*a*tet_I1_dot(i)*sind(tet_I1(i));
t52_dot=l2*tet_I1_dot(i)*cosd(tet_I1(i))*cosd(tet_I2(i))-l2*tet_I2_dot(i)*sind(tet_I1(i))*sind(tet_I2(i));
t62_dot=-l2*tet_I2_dot(i)*cosd(tet_I2(i));


T2=[0,0,zeros(6,1).';
    0,0,zeros(6,1).';
    1,0,zeros(6,1).';
    t41,t42,zeros(6,1).';
    t51,t52,zeros(6,1).';
    0,t62,zeros(6,1).'];

T2_dot=[0,0,zeros(6,1).';
    0,0,zeros(6,1).';
    0,0,zeros(6,1).';
    t41_dot,t42_dot,zeros(6,1).';
    t51_dot,t52_dot,zeros(6,1).';
    0,t62_dot,zeros(6,1).'];


%--------------------------for T3
t41=-l1*sind(tet_I1(i))-l2*sind(tet_I1(i))*sind(tet_I2(i))-l3*sind(tet_I1(i))*cosd(tet_I3(i))-aI5*sind(tet_I1(i));
t42=l2*cosd(tet_I1(i))*cosd(tet_I2(i));
t43=-l3*cosd(tet_I1(i))*sind(tet_I3(i));
t51=l1*cosd(tet_I1(i))+l2*cosd(tet_I1(i))*sind(tet_I2(i))+l3*cosd(tet_I1(i))*cosd(tet_I3(i))+aI5*cosd(tet_I1(i));
t52=l2*sind(tet_I1(i))*cosd(tet_I2(i));
t53=-l3*sind(tet_I1(i))*sind(tet_I3(i));
t62=-l2*sind(tet_I2(i));
t63=-l3*cosd(tet_I3(i));

t41_dot=-l1*tet_I1_dot(i)*cosd(tet_I1(i))-l2*tet_I1_dot(i)*cosd(tet_I1(i))*sind(tet_I2(i))-l2*tet_I2_dot(i)*sind(tet_I1(i))*cosd(tet_I2(i))-...
    l3*tet_I1_dot(i)*cosd(tet_I1(i))*cosd(tet_I3(i))+l3*tet_I3_dot(i)*sind(tet_I1(i))*sind(tet_I3(i))-aI5*tet_I1_dot(i)*cosd(tet_I1(i));
t42_dot=-l2*tet_I1_dot(i)*sind(tet_I1(i))*cosd(tet_I2(i))-l2*tet_I2_dot(i)*cosd(tet_I1(i))*sind(tet_I2(i));

t51_dot=-l1*tet_I1_dot(i)*sind(tet_I1(i))-l2*tet_I1_dot(i)*sind(tet_I1(i))*sind(tet_I2(i))+l2*tet_I2_dot(i)*cosd(tet_I1(i))*cosd(tet_I2(i))-...
    l3*tet_I1_dot(i)*sind(tet_I1(i))*cosd(tet_I3(i))-l3*tet_I3_dot(i)*cosd(tet_I1(i))*sind(tet_I3(i))-aI5*tet_I1_dot(i)*sind(tet_I1(i));
t52_dot=l2*tet_I1_dot(i)*cosd(tet_I1(i))*cosd(tet_I2(i))-l2*tet_I2_dot(i)*sind(tet_I1(i))*sind(tet_I2(i));
t53_dot=-l3*tet_I1_dot(i)*cosd(tet_I1(i))*sind(tet_I3(i))-l3*tet_I3_dot(i)*sind(tet_I1(i))*cosd(tet_I3(i));
t62_dot=-l2*tet_I2_dot(i)*cosd(tet_I2(i));
t63_dot=l3*tet_I3_dot(i)*sind(tet_I3(i));
t43_dot=l3*tet_I1_dot(i)*sind(tet_I1(i))*sind(tet_I3(i))-l3*tet_I3_dot(i)*cosd(tet_I3(i))*cosd(tet_I1(i));


T3=[0,0,0,zeros(5,1).';
    0,0,0,zeros(5,1).';
    1,0,0,zeros(5,1).';
    t41,t42,t43,zeros(5,1).';
    t51,t52,t53,zeros(5,1).';
    0,t62,t63,zeros(5,1).'];

T3_dot=[0,0,0,zeros(5,1).';
    0,0,0,zeros(5,1).';
    0,0,0,zeros(5,1).';
    t41_dot,t42_dot,t43_dot,zeros(5,1).';
    t51_dot,t52_dot,t53_dot,zeros(5,1).';
    0,t62_dot,t63_dot,zeros(5,1).'];
tet_I3_dot=tet_I3_dot(i);
tet_I4_dot=tet_I4_dot(i);
tet_II3_dot=tet_II3_dot(i);
tet_II4_dot=tet_II4_dot(i);
tet_I1_dot=tet_I1_dot(i);
tet_I2_dot=tet_I2_dot(i);
tet_II1_dot=tet_II1_dot(i);
tet_II2_dot=tet_II2_dot(i);

tet_I1=tet_I1(i);
tet_I4=tet_I4(i);
tet_I2=tet_I2(i);
tet_II2=tet_II2(i);
tet_I3=tet_I3(i);
tet_II3=tet_II3(i);
tet_II1=tet_II1(i);
tet_II4=tet_II4(i);


%--------------------------for T4
t41=-l1*sind(tet_I1)-l2*sind(tet_I1)*sind(tet_I2)-l3*sind(tet_I1)*cosd(tet_I3)-l4*cosd(tet_I1)+aII7*sind(tet_II1+tet_II4);
t42=l2*cosd(tet_I1)*cosd(tet_I2);
t43=-l3*cosd(tet_I1)*sind(tet_I3);
t44=aII7*sind(tet_II1+tet_II4);
t51=l1*cosd(tet_I1)+l2*cosd(tet_I1)*sind(tet_I2)+l3*cosd(tet_I1)*cosd(tet_I3)+l4*cosd(tet_I1)-aII7*cosd(tet_II1+tet_II4);
t52=l2*sind(tet_I1)*cosd(tet_I2);
t53=-l3*sind(tet_I1)*sind(tet_I3);
t54=-aII7*cosd(tet_II1+tet_II4);
t62=-l2*sind(tet_I2);
t63=-l3*cos(tet_I3);

t41_dot=-l1*tet_I1_dot*cosd(tet_I1)-l2*tet_I1_dot*cosd(tet_I1)*sind(tet_I2)-l2*tet_I2_dot*sind(tet_I1)*cosd(tet_I2)-...
    l3*tet_I1_dot*cosd(tet_I1)*cosd(tet_I3)+l3*tet_I3_dot*sind(tet_I1)*sind(tet_I3)+l4*tet_I1_dot*sind(tet_I1)+aII7*(tet_I1_dot+tet_I4_dot)*cosd(tet_II1+tet_II4);
t42_dot=-l2*tet_I1_dot*sind(tet_I1)*cosd(tet_I2)-l2*tet_I2_dot*cosd(tet_I1)*sind(tet_I2);
t43_dot=l3*tet_I1_dot*sind(tet_I1)*sind(tet_I3)-l3*tet_I3_dot*cosd(tet_I3)*cosd(tet_I1);
t44_dot=aII7*(tet_I1_dot+tet_I4_dot)*cosd(tet_II1+tet_II4);

t51_dot=-l1*tet_I1_dot*sind(tet_I1)+l2*tet_I1_dot*cosd(tet_I1)*cosd(tet_I2)-l2*tet_I2_dot*sind(tet_I1)*sind(tet_I2)-...
    l3*tet_I1_dot*sind(tet_I1)*cosd(tet_I3)-l3*tet_I3_dot*cosd(tet_I1)*sind(tet_I3)-l4*tet_I1_dot*sind(tet_I1)+aII7*(tet_I1_dot+tet_I4_dot)*sind(tet_II1+tet_II4);

t52_dot=l2*tet_I1_dot*cosd(tet_I1)*cosd(tet_I2)-l2*tet_I2_dot*sind(tet_I1)*sind(tet_I2);
t53_dot=-l3*tet_I1_dot*cosd(tet_I1)*sind(tet_I3)-l3*tet_I3_dot*sind(tet_I1)*cosd(tet_I3);
t54_dot=aII7*(tet_I1_dot+tet_I4_dot)*sind(tet_II1+tet_II4);
t62_dot=-l2*tet_I2_dot*cosd(tet_I2);
t63_dot=l3*tet_I3_dot*sind(tet_I3);



T4=[0,0,0,0,zeros(4,1).';
    0,0,0,0,zeros(4,1).';
    1,0,0,0,zeros(4,1).';
    t41,t42,t43,t44,zeros(4,1).';
    t51,t52,t53,t54,zeros(4,1).';
    0,t62,t63,0,zeros(4,1).'];

T4_dot=[0,0,0,0,zeros(4,1).';
    0,0,0,0,zeros(4,1).';
    0,0,0,0,zeros(4,1).';
    t41_dot,t42_dot,t43_dot,t44_dot,zeros(4,1).';
    t51_dot,t52_dot,t53_dot,t54_dot,zeros(4,1).';
    0,t62_dot,t63_dot,0,zeros(4,1).'];


%--------------------------for T5
t41=-l1*sind(tet_II1)-l2*sind(tet_II1)*sind(tet_II2)-l3*sind(tet_II1)*cosd(tet_II3)-aII6*sind(tet_II1);
t42=l2*cosd(tet_II1)*cosd(tet_II2);
t43=-l3*cosd(tet_II1)*sind(tet_II3);
t51=l1*cosd(tet_II1)+l2*cosd(tet_II1)*sind(tet_II2)+l3*cosd(tet_II1)*cosd(tet_II3)+aII6*cosd(tet_II1);
t52=l2*sind(tet_II1)*cosd(tet_II2);
t53=-l3*sind(tet_II1)*sind(tet_II3);
t62=-l2*sind(tet_II2);
t63=-l3*cosd(tet_II3);

t41_dot=-l1*tet_II1_dot*cosd(tet_II1)-l2*tet_II1_dot*cosd(tet_II1)*sind(tet_II2)-l2*tet_II2_dot*sind(tet_II1)*cosd(tet_II2)-...
    l3*tet_II1_dot*cosd(tet_II1)*cosd(tet_II3)+l3*tet_II3_dot*sind(tet_II1)*sind(tet_II3)-aII6*tet_II1_dot*cosd(tet_II1);
t42_dot=-l2*tet_II1_dot*sind(tet_II1)*cosd(tet_II2)-l2*tet_II2_dot*cosd(tet_II1)*sind(tet_II2);
t43_dot=l3*tet_II1_dot*sind(tet_II1)*sind(tet_II3)-l3*tet_II3_dot*cosd(tet_II3)*cosd(tet_II1);

t51_dot=-l1*tet_II1_dot*sind(tet_II1)-l2*tet_II1_dot*sind(tet_II1)*sind(tet_II2)+l2*tet_II2_dot*cosd(tet_II1)*cosd(tet_II2)-...
    l3*tet_II1_dot*sind(tet_II1)*cosd(tet_II3)-l3*tet_II3_dot*cosd(tet_II1)*sind(tet_II3)+aII6*tet_II1_dot*sind(tet_II1);
t52_dot=l2*tet_II1_dot*cosd(tet_II1)*cosd(tet_II2)-l2*tet_II2_dot*sind(tet_II1)*sind(tet_II2);
t53_dot=-l3*tet_II1_dot*cosd(tet_II1)*sind(tet_II3)-l3*tet_II3_dot*sind(tet_II1)*cosd(tet_II3);
t62_dot=-l2*tet_II2_dot*cosd(tet_II2);
t63_dot=l3*tet_II3_dot*sind(tet_II3);

T5=[zeros(5,1).',0,0,0;
    zeros(5,1).',0,0,0;
    zeros(5,1).',1,0,0;
    zeros(5,1).',t41,t42,t43;
    zeros(5,1).',t51,t52,t53;
    zeros(5,1).',0,t62,t63];

T5_dot=[zeros(5,1).',0,0,0;
    zeros(5,1).',0,0,0;
    zeros(5,1).',0,0,0;
    zeros(5,1).',t41_dot,t42_dot,t43_dot;
    zeros(5,1).',t51_dot,t52_dot,t53_dot;
    zeros(5,1).',0,t62_dot,t63_dot];

%--------------------------for T2
t41=-l1*sind(tet_II1)-l2*sind(tet_II1)*sind(tet_II2)+0.47*a*sind(tet_II1);
t42=l2*cosd(tet_II1)*cosd(tet_II2);
t51=l1*cosd(tet_II1)+l2*cosd(tet_II1)*sind(tet_II2)-0.47*a*cosd(tet_II1);
t52=l2*sind(tet_II1)*cosd(tet_II2);
t62=-l2*sind(tet_II2);

t41_dot=-l1*tet_II1_dot*cosd(tet_II1)-l2*tet_II1_dot*cosd(tet_II1)*sind(tet_II2)-l2*tet_II2_dot*sind(tet_II1)*cosd(tet_II2)+0.47*a*tet_II1_dot*cosd(tet_II1);
t42_dot=-l2*tet_II1_dot*sind(tet_II1)*cosd(tet_I2)-l2*tet_II2_dot*cosd(tet_II1)*sind(tet_II2);
t51_dot=-l1*tet_II1_dot*sind(tet_II1)-l2*tet_II1_dot*sind(tet_II1)*sind(tet_II2)+l2*tet_II2_dot*cosd(tet_II1)*cosd(tet_II2)+0.47*a*tet_II1_dot*sind(tet_II1);
t52_dot=l2*tet_II1_dot*cosd(tet_II1)*cosd(tet_II2)-l2*tet_II2_dot*sind(tet_II1)*sind(tet_II2);
t62_dot=-l2*tet_II2_dot*cosd(tet_II2);

T6=[zeros(6,1).',0,0;
    zeros(6,1).',0,0;
    zeros(6,1).',1,0;
    zeros(6,1).',t41,t42;
    zeros(6,1).',t51,t52;
    zeros(6,1).',0,t62];

T6_dot=[zeros(6,1).',0,0;
    zeros(6,1).',0,0;
    zeros(6,1).',0,0;
    zeros(6,1).',t41_dot,t42_dot;
    zeros(6,1).',t51_dot,t52_dot;
    zeros(6,1).',0,t62_dot];

%--------------------------for T7
T7=[zeros(2,7),zeros(1,2).';zeros(7,1).',0;zeros(2,7),zeros(1,2).';zeros(1,7),1];
T7_dot=zeros(6,8);


%~~~~~~~~~~~~~~~~~~~twist-shaping matrices Ti of all rigid bodies,
T=[T1;T2;T3;T4;T5;T6;T7];
T_dot=[T1_dot;T2_dot;T3_dot;T4_dot;T5_dot;T6_dot;T7_dot];


%% ~~~~The Angular Velocity Dyads~~~~~
omeg_1=[0,0,tet_I1].';
omeg_2=omeg_1;omeg_3=omeg_1;

omeg_5=[0,0,tet_II1].';
omeg_6=omeg_5;omeg_7=omeg_5;

omeg_4=[0,0,tet_I1-tet_II1].';

W_1=[CPM(omeg_1),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_2=[CPM(omeg_2),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_3=[CPM(omeg_3),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_4=[CPM(omeg_4),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_5=[CPM(omeg_5),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_6=[CPM(omeg_6),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W_7=[CPM(omeg_7),zeros(3,3);
     zeros(3,3),zeros(3,3)];
W = blkdiag(W_1,W_2,W_3,W_4,W_5,W_6,W_7);
%% ~~~~~~~  Screw matrix
R_I_1=[cosd(tet_I1),-sind(tet_I1),0;
       sind(tet_I1),cosd(tet_I1),0;
       0,0,1];
R_II_1=[cosd(tet_II1),-sind(tet_II1),0;
       sind(tet_II1),cosd(tet_II1),0;
       0,0,1];

R_I_2=[cosd(tet_I2),0,sind(tet_I2);
       0,1,0;
       -sind(tet_I2),0,cosd(tet_I2)];
R_II_2=[cosd(tet_II2),0,sind(tet_II2);
       0,1,0;
       -sind(tet_II2),0,cosd(tet_II2)];   
R_I_3=[cosd(tet_I3) 0 sind(tet_I3) ;
            0 1 0 ;
            -sind(tet_I3) 0 cosd(tet_I3)];  %Y
R_II_3=[cosd(tet_II3) 0 sind(tet_II3) ;
            0 1 0 ;
            -sind(tet_II3) 0 cosd(tet_II3)];

e_I_1=-k;
Ro_I_1=[0;0;0.2];

e_I_2=R_I_1*[0;1;0];
u_I_P=R_I_1*R_I_2*[l2;0;0];
gama_I_1=R_I_1*[0;0;0.2];
Ro_I_2=R_I_1*(-[0.47*a;0;0]+[l1;0;0]);

e_I_3=R_I_1*[0;1;0];
u_I_D=R_I_1*R_I_3*[l3;0;0];
gama_I_2=R_I_1*[0.47*a;0;0];
Ro_I_3=R_I_1*([aI5;0;0]+[l1;0;0]);

Ro_M=R_I_1*[aII7;0;tan(45)*aII7];

e_II_4=-k;
gama_I_3=R_I_1*[aII7;0;tan(45)*aII7];
Ro_II_3=R_II_1*([aII6;0;0]+[l1;0;0]);

e_II_3=R_II_1*R_II_3*[0;1;0];
u_II_D=(R_II_1*R_II_3*[l3;0;0])-([l0;0;0]);
gama_II_3=R_II_1*([aII6;0;0]+[l1;0;0]);
Ro_II_2=R_II_1*(-[0.47*a;0;0]+[l1;0;0]);

e_II_2=R_II_1*[0;1;0];
u_II_P=(R_II_1*R_II_2*[l2;0;0])-([l0;0;0]);
gama_II_2=R_II_1*([0.47*a;0;0]+[l1;0;0]);
Ro_II_1=R_II_1*[0.47*a;0;0];


e_II_1=R_II_1*(-k);
gama_II_1=-[0;0;0.2];


%% Part B Screw matrix


s_I_r=[e_I_1;zeros(3,1)];
s_II_r=[e_II_1;zeros(3,1)];

s_I_piP=[zeros(3,1); cross(e_I_2,u_I_P)];
s_II_piP=[zeros(3,1);cross(e_II_2,u_II_P)];

s_I_pid=[zeros(3,1); cross(e_I_3,u_I_D)];
s_II_pid=[zeros(3,1); cross(e_II_3,u_II_D)];

S_I_s=[e_I_3,e_I_2,e_I_1;
       zeros(3,1),zeros(3,1),zeros(3,1)]; 

s_II_rm=[e_II_1;zeros(3,1)];

S = blkdiag(s_I_r,s_I_piP,s_I_pid,S_I_s,s_II_rm,s_II_pid,s_II_piP,s_II_r);

%% ~~~~~  actuator wrench-shaping matrix
a_I=sqrt( (2*a/3)^2 + (a/3)^2 )*[1;0;1];
b_I=sqrt( (a/3)^2 + (a/3)^2 )*[-1;0;1];
a_II=sqrt( (2*a/3)^2 + (a/3)^2 )*[-1;0;1];
b_II=sqrt( (a/3)^2 + (a/3)^2 )*[1;0;1];
A_I=CPM(a_I);A_II=CPM(a_II);   
B_I=CPM(b_I);B_II=CPM(b_II);     
   
n_I=1/(2*l2)*[R_I_1*(A_I+B_I)*R_I_2;2*R_I_1*R_I_2]*[e_I_2];
n_II=1/(2*l2)*[R_II_1*(A_II+B_II)*R_II_2;2*R_II_1*R_II_2]*[e_II_2];

   
   
TA=[s_I_r,zeros(6,1),zeros(6,1),zeros(6,1);
    zeros(6,1),n_I,zeros(6,1),zeros(6,1);
    zeros(18,1),zeros(18,1),zeros(18,1),zeros(18,1);
    zeros(6,1),zeros(6,1),zeros(6,1),n_II;
    zeros(6,1),zeros(6,1),s_II_r,zeros(6,1)];

%% ~~~ Derivation of the constraint matrix ~~~



E_I_1=CPM(e_I_1);
%E_I_4=CPM(e_I_4);
E_II_1=CPM(e_II_1);
E_II_4=CPM(e_II_4);


R_I_1=CPM(Ro_I_1);
R_I_2=CPM(Ro_I_2);
R_I_3=CPM(Ro_I_3);
R_M=CPM(Ro_M);
R_II_3=CPM(Ro_II_3);
R_II_2=CPM(Ro_II_2);
R_II_1=CPM(Ro_II_1);

L_I_P=CPM(cross(e_I_2,u_I_P));
L_I_D=CPM(cross(e_I_3,u_I_D));
L_I_3=CPM(cross(e_I_3,u_I_D));
L_II_D=CPM(cross(e_II_3,u_II_D));
L_II_P=CPM(cross(e_II_2,u_II_P));

D_I_1=CPM(gama_I_1+u_I_P);
D_I_2=CPM(gama_I_2+u_I_P);
D_I_3=CPM(gama_I_3+u_I_D);
D_II_3=CPM(gama_II_3+u_II_D);
D_II_2=CPM(gama_II_2+u_II_P);
D_II_1=CPM(gama_II_1);
D_M=CPM(gama_I_3);




K11=[E_I_1,zeros(3,3);R_I_1,eye(3)];
K21=[-eye(3),zeros(3,3);L_I_P*D_I_1,-L_I_P];
K22=[eye(3),zeros(3,3);L_I_P*R_I_2,L_I_P];
K32=[-eye(3),zeros(3,3);L_I_D*D_I_2,-L_I_D];
K33=[eye(3),zeros(3,3);L_I_3*R_I_3,L_I_3];
K43=[zeros(3,3),zeros(3,3);D_I_3,-eye(3)];
K44=[zeros(3,3),zeros(3,3);R_M,eye(3)];
K54=[-E_II_4,zeros(3,3);D_M,-eye(3)];
K55=[E_II_4,zeros(3,3);R_II_3,eye(3)];
K65=[-eye(3),zeros(3,3);L_II_D*D_II_3,-L_II_D]; %%% 
K66=[eye(3),zeros(3,3);L_II_D*R_II_2,L_II_D]; %%% 
K76=[-eye(3),zeros(3,3);L_II_P*D_II_2,-L_II_P]; %%% 
K77=[eye(3),zeros(3,3);L_II_P*R_II_1,L_II_P]; %%% 
K87=[-E_II_1,zeros(3,3);D_II_1,-eye(3)]; %%% 

z6=zeros(6,6);
K=[K11,z6,z6,z6,z6,z6,z6;   
   K21,K22,z6,z6,z6,z6,z6;
    z6,K32,K33,z6,z6,z6,z6;
     z6,z6,K43,K44,z6,z6,z6;
      z6,z6,z6,K54,K55,z6,z6;
      z6,z6,z6,z6,K65,K66,z6;
       z6,z6,z6,z6,z6,K76,K77;
        z6,z6,z6,z6,z6,z6,K87];

%% %% Ws
m11=1.887; m12=0.5;  
m21=1.887; m22=0.5; 
m31=1.887; m32=0.5;  
mp = 0.09;
cc=zeros(3,3);
g = 9.81;
Ws = zeros(42,1);

Ws(6) = -m11*g;
Ws(12) = -m12*g;
Ws(18) = -m21*g;
Ws(24) = -m22*g;
Ws(30) = -m31*g;
Ws(36) = -m32*g;
Ws(42) = -mp*g;
%% the Inertia Dyads
M=0.1*eye(42);
eta = (M*T*qdd + (W*M*T+M*T_dot)*qd-Ws);

%% Compact system of 6n + d linear equations,

KTAS=[K' TA ; S' zeros(10,4)];
ETA=[eta;zeros(10,1)];
A = pinv(KTAS)*ETA;



Mx_I1(i)=A(1);
My_I1(i)=A(2);
Mz_I1(i)=A(3);
fx_I1(i)=A(4);
fy_I1(i)=A(5);
fz_I1(i)=A(6);

Mx_I2(i)=A(7);
My_I2(i)=A(8);
Mz_I2(i)=A(9);
fx_I2(i)=A(10);
fy_I2(i)=A(11);
fz_I2(i)=A(12);

Mx_I3(i)=A(13);
My_I3(i)=A(14);
Mz_I3(i)=A(15);
fx_I3(i)=A(16);
fy_I3(i)=A(17);
fz_I3(i)=A(18);

Mx_I4(i)=A(19);
My_I4(i)=A(20);
Mz_I4(i)=A(21);
fx_I4(i)=A(22);
fy_I4(i)=A(23);
fz_I4(i)=A(24);

Mx_II1(i)=A(25);
My_II1(i)=A(26);
Mz_II1(i)=A(27);
fx_II1(i)=A(28);
fy_II1(i)=A(29);
fz_II1(i)=A(30);

Mx_II2(i)=A(31);
My_II2(i)=A(32);
Mz_II2(i)=A(33);
fx_II2(i)=A(34);
fy_II2(i)=A(35);
fz_II2(i)=A(36);

Mx_II3(i)=A(37);
My_II3(i)=A(38);
Mz_II3(i)=A(39);
fx_II3(i)=A(40);
fy_II3(i)=A(41);
fz_II3(i)=A(42);

Mx_II4(i)=A(43);
My_II4(i)=A(44);
Mz_II4(i)=A(45);
fx_II4(i)=A(46);
fy_II4(i)=A(47);
fz_II4(i)=A(48);


ta1(i)=A(49);
ta2(i)=A(50);
ta3(i)=A(51);
ta4(i)=A(52);
end
%% 

figure
plot(time,real(fx_I1), 'LineWidth' , 2)
hold on
plot(time,real(fy_I1), 'LineWidth' , 2)
hold on
plot(time,real(fz_I1), 'LineWidth' , 2)
xlim ([0.5 5.1])
legend('f_x','f_y','f_z')
ylabel('Force')
xlabel('Time (s)')
title('Forces at the first revolute joints of the limb I',' along the test trajectory')
grid on

figure
plot(time,real(Mx_I1), 'LineWidth' , 2)
hold on
plot(time,real(My_I1), 'LineWidth' , 2)
hold on
plot(time,real(Mz_I1), 'LineWidth' , 2)
xlim ([0.5 5.1])
legend('M_x','M_y','M_z')
ylabel('Moment')
xlabel('Time (s)')
title('Moments at first revolute joints of the limb I',' along the test trajectory')

grid on
%%
figure
plot(time,real(ta1), 'LineWidth' , 2)
hold on
plot(time,real(ta2),'--', 'LineWidth' , 2)
legend('\tau_I_1','\tau_I_2')
xlim([0.45 5.15])
ylabel('Moment')
xlabel('Time (s)')
%ylim([-500,500])
title('The time-histories of the actuator torques')
grid on
%% 
figure
plot(time,real(ta3), 'LineWidth' , 2)
hold on
plot(time,real(ta4),'--', 'LineWidth' , 2)
legend('\tau_I_I_1','\tau_I_I_2')
xlim([0.45 5.15])
ylabel('Moment')
xlabel('Time (s)')
%ylim([-500,500])
title('The time-histories of the actuator torques')
grid on






