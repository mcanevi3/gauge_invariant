clear;clc;

alpha=pi/6;
beta=pi/12;
Dw=[1;2];
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(Dw);

t=0:0.01:1;
x=sin(10*t);
w=2*ones(size(t));

y=[sin(x+alpha);sin(x+beta)]+Dw*w;
z=[1,0]*Pival*y;

A=sqrt((Pival(1,1)*cos(alpha)+Pival(1,2)*cos(beta))^2+(Pival(1,1)*sin(alpha)+Pival(1,2)*sin(beta))^2);
theta=atan2(Pival(1,1)*sin(alpha)+Pival(1,2)*sin(beta),Pival(1,1)*cos(alpha)+Pival(1,2)*cos(beta));

x1=asin(z/A)-theta;
x2=pi-asin(z/A)-theta;

figure(1);clf;hold on;grid on;
plot(t,x,'k','LineWidth',3,'DisplayName','x');
plot(t,x1,'r--','LineWidth',1.5,'DisplayName','x_1');
plot(t,x2,'b--','LineWidth',1.5,'DisplayName','x_2');

figure(2);clf;hold on;grid on;
plot(t,y(1,:),'k','LineWidth',3,'DisplayName','x');
plot(t,sin(x+alpha),'m--','LineWidth',3,'DisplayName','x');
plot(t,sin(x1+alpha),'r','LineWidth',1.5,'DisplayName','x_1');
plot(t,sin(x2+alpha),'b','LineWidth',1.5,'DisplayName','x_1');

