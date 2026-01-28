clear;clc;

C=ones(3,1);

eig(C*C')

Pim=eye(3)-inv(3)*C*C'
