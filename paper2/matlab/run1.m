clear;clc;

%% System model
A=[1,2;-3,-4];
Bw=[1;2];
C=[1,0;2,0;1,2];

Cperp=null(C');
Dw=Cperp+C*[1;1];

n=size(A,1);
p=size(C,1);
%% Projection
proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(C);

%% Disturbance model
Aw=0;
Cw=1;
nw=size(Aw,1);
pw=size(Cw,1);

%% Augmented model
Aa=[A,Bw*Cw;zeros(nw,n),Aw];
Ca=[C,Dw*Cw];

La=place(Aa',Ca',[-10,-10,-20])'



