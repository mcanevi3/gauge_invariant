clear;clc;

syms x real;

fx=-x^2;
Bw=1;
hx=[sin(x+pi/2);sin(x)];
Dw=[1;2];

%% linearization 
Ax=diff(fx,x);
Ax 

Cx=diff(hx,x);
Cx

%% projection 
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(Cx);

Pival=proj(Dw);
Pival*hx




% funPi11=matlabFunction(Pival(1,1));
% funPi12=matlabFunction(Pival(1,2));
% funPi22=matlabFunction(Pival(2,2));
% 
% x0=4;
% t=0:0.01:1;
% xval=x0*exp(-4*t)+2;
% xhatval=(x0/2)*exp(-4*t)+2;
% 
% figure(1);clf;
% subplot(3,1,1);cla;hold on;grid on;
% plot(t,funPi11(xval),'k','LineWidth',2);
% plot(t,funPi11(xhatval),'b','LineWidth',2);
% 
% subplot(3,1,2);cla;hold on;grid on;
% plot(t,funPi12(xval),'k','LineWidth',2);
% plot(t,funPi12(xhatval),'b','LineWidth',2);
% 
% subplot(3,1,3);cla;hold on;grid on;
% plot(t,funPi22(xval),'k','LineWidth',2);
% plot(t,funPi22(xhatval),'b','LineWidth',2);
% 


