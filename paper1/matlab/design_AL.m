clear;clc;

A=[1,2;-3,-4];
Bw=[1;1];
C=[2,1;4,2];
Dw=[1;1];
Cz=[1,0;0,1];

% Augmented Luenberger design
A_aug=[A,Bw;zeros(size(Bw')),zeros(size(Bw,2),size(Bw,2))];
C_aug=[C,Dw];

L_aug=place(A_aug',C_aug',[-3,-4,-5])'
eig(A_aug-L_aug*C_aug)
