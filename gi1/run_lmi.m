clear;clc;

% n states, m controls,p outputs, k noise inputs
n=1;m=1;p=3; k=3;
A=[-2,zeros(n,p);zeros(p,n),zeros(p,p)]; % (n+p)x(n+p) 
Bu=[ones(n,m);zeros(p,m)]; % (n+p)xm
Bw=[zeros(n,k);zeros(p,k)]; %(n+p)xk
C=[ones(p,n),eye(p)]; %px(n+p)
Du=zeros(p,m); %pxm
Dw=zeros(p,k); %pxk
Dw=eye(k);

Phio=[C;C*A;C*A*A];

Ce=eye(n+p);

cvx_clear;
cvx_begin sdp quiet;
variable P(n+p,n+p) symmetric;
variable Y(n+p,p);
variable W(p,p) symmetric;
minimize(trace(W));
A'*P+P*A-Y*C-(Y*C)'+Ce'*Ce<=0*-1e-3*eye(n+p);
[W,Y';Y,P]>=0;
A'*P+P*A-Y*C-(Y*C)'+2*1*P<=0;
cvx_end;

if strcmp(cvx_status,"Solved") || strcmp(cvx_status,"Inaccurate/Solved")
    L=inv(P)*Y;
    
    Gz=ss(A-L*C,-L,Ce,zeros(n+p,p));
    Gz=tf(Gz);

    mu=sqrt(trace(W));
    disp("eig(A-LC):");
    disp(eig(A-L*C));
    
    disp("mu:"+string(mu));
    disp("||Gz||_2:"+string(norm(Gz,2)));


    B1=Bu;
    B2=Bw;
    C1=C;
    D1=Du;
    D2=Dw;
    C2=Ce;
    L=-L;

    x0=[4;5;-24;13];
    simout=sim("model_pi.slx");

    data=simout.simout;
    t=simout.tout;
    x=data(:,1:n+p);
    y=data(:,n+p+1:n+p+p);
    xhat=data(:,n+p+p+1:n+p+p+n+p);
    yhat=data(:,n+p+p+n+p+1:end);

    figure(1);clf;
    subplot(1,3,1);cla;hold on;grid on;legend("show");
    plot(t,x(:,1:n),'k','LineWidth',2,'DisplayName','x');
    plot(t,xhat(:,1:n),'LineWidth',2,'DisplayName','xhat');

    subplot(1,3,2);cla;hold on;grid on;legend("show");
    plot(t,x(:,n+1:n+p),'k','LineWidth',2,'DisplayName','b');
    plot(t,xhat(:,n+1:n+p),'LineWidth',2,'DisplayName','bhat');

    subplot(1,3,3);cla;hold on;grid on;legend("show");
    plot(t,y,'k','LineWidth',2,'DisplayName','y');
    plot(t,yhat,'LineWidth',2,'DisplayName','yhat');

else
    disp("Status:"+string(cvx_status));
    return;
end