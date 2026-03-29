clear;clc;

Ct=2.2*1e-3;
Lt=1.8*1e-3;
Rt=0.2;
A=[0,1/Ct;-1/Lt,-Rt/Lt];
Bw=[-1/Ct;0];

C=[1,0;2,1;0,1;0,1];
Cperp=null(C');
Dw=Cperp*[1;1]+C*[1;1];
Dn=rand(4,1);

n=size(A,1);
p=size(C,1);

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';

M=[Dw Dn];

t=0:0.001:5;
x1=sin(2*t);
x2=cos(2*t);
w=t.^2-4*t;
n=randn(size(t));
n=n/var(n);
n=n-mean(n);

phi=C*[x1;x2]+Dw*w+Dn*n;
% proj(M)*(Dw*w+Dn*n)
% proj(M)*(C*[x1;x2])
figure(1);clf;hold on;grid on;legend("show");
xlabel("time(sec)");

xhat=pinv(proj(M)*C)*(proj(M)*phi);
plot(t,phi(1,:),'m','LineWidth',1,'DisplayName','y_1');
plot(t,phi(2,:),'c','LineWidth',1,'DisplayName','y_2');
plot(t,phi(3,:),'g','LineWidth',1,'DisplayName','y_3');
plot(t,phi(4,:),'k','LineWidth',1,'DisplayName','y_4');
plot(t,xhat(1,:),'r','LineWidth',2,'DisplayName','xhat_1');
plot(t,xhat(2,:),'b','LineWidth',2,'DisplayName','xhat_2');


exportgraphics(gcf,"../img/ex2_results.pdf",'ContentType',"vector");



M=[C Dn];
what=pinv(proj(M)*Dw)*(proj(M)*phi);
figure(2);clf;hold on;grid on;legend("show");
xlabel("time(sec)");

plot(t,w,'k','LineWidth',3,'DisplayName','w');
plot(t,what,'r','LineWidth',2,'DisplayName','what');
exportgraphics(gcf,"../img/ex2_results2.pdf",'ContentType',"vector");
