clear;clc;


%% simulate
sim("model2");
%% results
results=load('res_extended.mat');
results=results.data;

t=results(1,:);
x=results(2,:);
xluen=results(3,:);
xgi=results(4,:);
w=results(5,:);
wluen=results(6,:);
wgi=results(7,:);

y=results(8:9,:);
yluen=results(10:11,:);
ygi=results(12:13,:);

color_luen='r';
color_gi='b';
figure(1);clf;hold on;grid on;xlabel("time(s)");ylabel("x(t)");title("State");legend("show","Location","best");ax1=gca;
plot(ax1,t,x,'k','LineWidth',3,'DisplayName','system');
plot(ax1,t,xluen,color_luen,'LineWidth',2,'DisplayName','luenberger');
plot(ax1,t,xgi,color_gi,'LineWidth',2,'DisplayName','gi');

exportgraphics(gcf,"../img/ex2_states1.pdf",'ContentType',"vector");

figure(2);clf;hold on;grid on;
xlabel("time(s)");ylabel("w(t)");title("Disturbance");legend("show","Location","best");
plot(t,w,'k','LineWidth',3,'DisplayName','w');
plot(t,wluen,color_luen,'LineWidth',2,'DisplayName','w luenberger');
plot(t,wgi,color_gi,'LineWidth',2,'DisplayName','w gi');

exportgraphics(gcf,"../img/ex2_disturbance.pdf",'ContentType',"vector");

figure(3);clf;
subplot(2,1,1);cla;hold on;grid on;xlabel("time(s)");ylabel("y_1(t)");title("Output 1");legend("show");ax1=gca;
subplot(2,1,2);cla;hold on;grid on;xlabel("time(s)");ylabel("y_2(t)");title("Output 2");legend("show");ax2=gca;
plot(ax1,t,y(1,:),'k','LineWidth',3,'DisplayName','y');
plot(ax1,t,yluen(1,:),color_luen,'LineWidth',2,'DisplayName','y luenberger');
plot(ax1,t,ygi(1,:),color_gi,'LineWidth',2,'DisplayName','y gi');

plot(ax2,t,y(2,:),'k','LineWidth',3,'DisplayName','y');
plot(ax2,t,yluen(2,:),color_luen,'LineWidth',2,'DisplayName','y luenberger');
plot(ax2,t,ygi(2,:),color_gi,'LineWidth',2,'DisplayName','y gi');


exportgraphics(gcf,"../img/ex2_outputs.pdf",'ContentType',"vector");
