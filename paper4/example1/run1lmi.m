clear;clc;
%%%"Location","best"

%% System model
Ct=2.2*1e-3;
Lt=1.8*1e-3;
Rt=0.2;
A=[0,1/Ct;-1/Lt,-Rt/Lt];
Bw=[-1/Ct;0];
C=[1,0;1,0;0,1];
Cperp=null(C');
Dw=Cperp+C*[1;1];

n=size(A,1);
p=size(C,1);
r=size(Bw,2);
%% Projection
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);

Pival_Dw=proj(Dw);

%% LMI
alpha = 100;   % desired stability margin
% (P*A-Y*C)+(P*A-Y*C)'+2*alpha*P<=0;

Cz=[1,0];
Cz=eye(2);
nz=size(Cz,1);

cvx_begin sdp quiet;
variable P(n,n) symmetric
variable Y(n,p) 
variable Z(nz,nz)
variable gmma
variable gmma_h2
minimize(gmma+gmma_h2)
subject to
P >= 1e-6*eye(n);
[(P*A-Y*C)+(P*A-Y*C)', P*Bw-Y*Dw,Cz';
(P*Bw-Y*Dw)',-gmma*eye(r),zeros(r,nz);
Cz, zeros(r,nz)',-gmma*eye(nz)
]<= 0;

[ (P*A-Y*C)+(P*A-Y*C)',  (P*Bw-Y*Dw);
  (P*Bw-Y*Dw)',         -eye(r) ] <= 0;

% 2. Trace constraint for H2
[ P,   Cz';
    Cz, Z ] >= 0;
    
trace(Z) <= gmma_h2;

cvx_end

% Recover observer gain
L = P\Y;
disp("L");
disp(L);
disp("eig(A-LC)");
disp(eig(A-L*C));
disp("Bw-L*Dw");
disp(Bw-L*Dw);

%% simulate
sim("model1");
%% results
results=load('results.mat');
results=results.data;

t=results(1,:);
x=results(2:3,:);
xluen=results(4:5,:);
xgi=results(6:7,:);
w=results(8,:);
wgi=results(9,:);
wluen=results(10,:);

results=load('outputs.mat');
results=results.data;
y=results(2:4,:);
yluen=results(5:7,:);
ygi=results(8:10,:);

color_luen='r--';
color_gi='c--';
figure(1);clf;
subplot(2,1,1);cla;hold on;grid on;xlabel("time(s)");ylabel("x_1(t)");title("State 1");legend("show");ax1=gca;
subplot(2,1,2);cla;hold on;grid on;xlabel("time(s)");ylabel("x_2(t)");title("State 2");legend("show");ax2=gca;
plot(ax1,t,x(1,:),'k','LineWidth',3,'DisplayName','state');
plot(ax1,t,xluen(1,:),color_luen,'LineWidth',2,'DisplayName','luenberger');
plot(ax1,t,xgi(1,:),color_gi,'LineWidth',2,'DisplayName','projection');

plot(ax2,t,x(2,:),'k','LineWidth',3,'DisplayName','state');
plot(ax2,t,xluen(2,:),color_luen,'LineWidth',2,'DisplayName','luenberger');
plot(ax2,t,xgi(2,:),color_gi,'LineWidth',2,'DisplayName','projection');

exportgraphics(gcf,"../img/ex1_states1.pdf",'ContentType',"vector");

figure(2);clf;hold on;grid on;
xlabel("time(s)");ylabel("w(t)");title("Disturbance");legend("show");
plot(t,w,'k','LineWidth',3,'DisplayName','disturbance');
plot(t,wluen,color_luen,'LineWidth',2,'DisplayName','luenberger');
plot(t,wgi,color_gi,'LineWidth',2,'DisplayName','projection');

exportgraphics(gcf,"../img/ex1_disturbance.pdf",'ContentType',"vector");

figure(3);clf;
subplot(3,1,1);cla;hold on;grid on;xlabel("time(s)");ylabel("y_1(t)");title("Output 1");legend("show");ax1=gca;
subplot(3,1,2);cla;hold on;grid on;xlabel("time(s)");ylabel("y_2(t)");title("Output 2");legend("show");ax2=gca;
subplot(3,1,3);cla;hold on;grid on;xlabel("time(s)");ylabel("y_3(t)");title("Output 3");legend("show");ax3=gca;
plot(ax1,t,y(1,:),'k','LineWidth',3,'DisplayName','y');
plot(ax1,t,yluen(1,:),color_luen,'LineWidth',2,'DisplayName','y luenberger');
plot(ax1,t,ygi(1,:),color_gi,'LineWidth',2,'DisplayName','y gi');

plot(ax2,t,y(2,:),'k','LineWidth',3,'DisplayName','y');
plot(ax2,t,yluen(2,:),color_luen,'LineWidth',2,'DisplayName','y luenberger');
plot(ax2,t,ygi(2,:),color_gi,'LineWidth',2,'DisplayName','y gi');

plot(ax3,t,y(3,:),'k','LineWidth',3,'DisplayName','y');
plot(ax3,t,yluen(3,:),color_luen,'LineWidth',2,'DisplayName','y luenberger');
plot(ax3,t,ygi(3,:),color_gi,'LineWidth',2,'DisplayName','y gi');

exportgraphics(gcf,"../img/ex1_outputs.pdf",'ContentType',"vector");
