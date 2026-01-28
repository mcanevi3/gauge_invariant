clear;clc;

yi=@(x,b,n)x+b+n;

x0=10;
b0=[0.3,-0.1,0.0]';
mu=0.5;

% n=[[0.02,-0.01,0.01]',[-0.01,0.02,-0.01]',[0.01,-0.02,0.01]',[-0.02,0.01,0.01]'];
% y=[[10.32,9.89,10.01]',[10.29,9.92,9.99]',[10.31,9.88,10.01]',[10.28,9.91,10.00]'];

nsample=4;
n=randn(3,nsample);n=n-mean(n);n=n./std(n);n=n*0.01;
y=x0+b0+n;

x=zeros(1,nsample+1);x(:,1)=x0;
bhat=zeros(3,nsample+1);
r=zeros(3,nsample+1);
rmean=zeros(1,nsample+1);
for k=1:nsample
    r(:,k)=y(:,k)-x(k)-bhat(:,k);
    rmean(k)=mean(r(:,k));
    x(k+1)=x(k)+mu*rmean(k);
    bhat(:,k+1)=bhat(:,k)+mu*r(:,k);
end
disp("b:");
disp(b0');
disp("bhat:");
disp(bhat(:,end)');
