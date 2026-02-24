clear;clc;

syms x real;
%x=1;
Cx=[cos(x);sin(x)];

proj=@(C)eye(size(C,1))-C*pinv(C'*C)*C';
Pival=proj(Cx);
Pival=simplify(Pival)

expr=Pival*Cx;
%simplify(expr)

simplify(pinv(Pival*[1;2]))

