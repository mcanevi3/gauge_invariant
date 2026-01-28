clear;clc;
%% system parameters
A=[1,2;-3,-4];
B=[1;0];
C=[1,0];

%% \dot{e}=(A-L*C)-L*b pd(s)=s^2+6*s+9 -> A*A+6*A+9I
Phio=[C;C*A];
coef=conv([1,3],[1,3]); % (s+3)(s+3)
L=[0,1]*inv(Phio')*(coef(1)*A'*A'+coef(2)*A'+coef(3)*eye(2));
L=L';
disp("L:");
disp(L');
disp("eig(A-LC)");disp(eig(A-L*C));

%% augmented system
Aa=[A,zeros(2,1);zeros(1,2),zeros(1,1)];
Ca=[C,0*eye(1)];
Phioa=[Ca;Ca*Aa;Ca*Aa*Aa];
if(rank(Phioa)==size(Phioa,1))
    disp("Observable");
    coef=conv([1,3],conv([1,3],[1,3])); % (s+3)(s+3)(s+3)
    La=[0,0,1]*inv(Phioa')*(coef(1)*Aa'*Aa'*Aa'+coef(2)*Aa'*Aa'+coef(3)*Aa'+coef(4)*eye(3));
    La=La';
    disp("La:");
    disp(La');
    disp("eig(Aa-LaCa)");disp(eig(Aa-La*Ca));
else
    disp("Not observable det(Phioa)="+string(det(Phioa)));
end
disp("**********");
if true
    res=Aa-La*Ca;
    fprintf("Aa:");disp(res);
    fprintf("eig:");disp(eig(res));
    res=Aa-[L;0]*Ca;
    fprintf("Aa:");disp(res);
    fprintf("eig:");disp(eig(res));

end