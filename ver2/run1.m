clear;clc;

A_=[1,2;-3,-4];
B1_=zeros(2,1);
B2_=zeros(2,1);
C1_=[1,0];
D1_=zeros(1,1);
D2_=ones(1,1);
C2_=[1,0];

A=[A_,zeros(2,1);zeros(1,2),zeros(1,1)];
B1=[B1_;zeros(1,1)];
B2=[B2_;zeros(1,1)];
C1=[C1_,eye(1)];
D1=D1_;
D2=D2_;
C2=[C2_,zeros(1,1)];

n=size(A,1);
m=size(C1,1);
cvx_clear;
cvx_begin sdp quiet;
variable P(n,n) symmetric;
variable W(n,m);
variable gm;
minimize(gm);
[
    A'*P+P*A+W*C1+C1'*W' , P*B2+W*D2 , C2';
    (P*B2+W*D2)' , -gm*eye(m) , zeros(m,m);
    C2 , zeros(m,m) , -gm*eye(m)
]<=0;
A'*P+P*A+W*C1+C1'*W'+2*2*P<=0;
P>=0;
gm>=0;
cvx_end;

if strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved')

    L=inv(P)*W;
    if norm(L,inf)>1e5
        disp("L too large");
%        return;
    end
    disp("eig(A+L*C1):");
    disp(eig(A+L*C1));

    x0=[1;2;3;];
    simout=sim("model_pi.slx");
    
    data=simout.simout;
    t=simout.tout;
    x=data(:,1:3);
    y=data(:,4);
    xhat=data(:,5:7);
    yhat=data(:,8);
    
    figure(1);clf;
    subplot(2,2,1);cla;hold on;grid on;title("x_1");ax1=gca;
    subplot(2,2,2);cla;hold on;grid on;title("x_2");ax2=gca;
    subplot(2,2,3);cla;hold on;grid on;title("b");ax3=gca;
    subplot(2,2,4);cla;hold on;grid on;title("y");ax4=gca;
    plot(ax1,t,x(:,1),'k','LineWidth',2);
    plot(ax1,t,xhat(:,1),'b','LineWidth',2);
    
    plot(ax2,t,x(:,2),'k','LineWidth',2);
    plot(ax2,t,xhat(:,2),'b','LineWidth',2);
    
    plot(ax3,t,x(:,3),'k','LineWidth',2);
    plot(ax3,t,xhat(:,3),'b','LineWidth',2);
    
    plot(ax4,t,y,'k','LineWidth',2);
    plot(ax4,t,yhat,'b','LineWidth',2);
else
    disp("cvx status:"+string(cvx_status));

end

