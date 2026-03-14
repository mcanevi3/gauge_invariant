clear;clc;

A=[1,2;-3,-4];
Bw=[2;1];
C=[1,0;1,0;0,1];
Dw=[1;2;3];

L=sym('L',[2,3]);

temp=Bw-L*Dw;
sol=solve(temp,L(:,1));
Lval=subs(L,sol);

syms s;
temp=det(s*eye(2)-(A-Lval*C));
prob=coeffs(temp,s,'all')==conv([1 4],[1 5]);
sol=solve(prob,L(:,2));

Lval2=subs(Lval,sol);
funL=matlabFunction(Lval2);

Lval3=funL(1,1);

disp("L:");
disp(Lval3);
disp("LC:");
disp(Lval3*C);
disp("||LC||:"+string(norm(Lval3*C,inf)));
% disp("Eig:");
% disp(eig(A-Lval3*C));
% disp("Bw-L*Dw");
% disp(Bw-Lval3*Dw);

Lval3=funL(10,12);

disp("L:");
disp(Lval3);
disp("LC:");
disp(Lval3*C);
disp("||LC||:"+string(norm(Lval3*C,inf)));
% disp("Eig:");
% disp(eig(A-Lval3*C));
% disp("Bw-L*Dw");
% disp(Bw-Lval3*Dw);
