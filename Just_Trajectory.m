clc
clear

%% ~~~~Given values~~~~
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

t=linspace(0,5.8,2000);
time=t;

%% ~~~Trajectory Function~~~
tic
for i=1:2000
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

yy(i)=y;
zpp(i)=z;
fprintf(' Calculated for z= %d & y= %d. \n',z,y)    
toc
end

%% ~~~Plot for trajectory~~~

figure 
plot(yy,zpp,'r', 'LineWidth' , 2)
title('The Test Trajectory')
axis([-0.5,10.5 0 4])
ylabel('z (cm)')
xlabel('y (cm)')
legend('The Test Trajectory line' )
grid on