
%% export
t=out.tout;
result=out.result;

figure(1);clf;hold on;grid on;
xlabel("Zaman(s)");
ylabel("\theta");
title("Sonuclar");
plot(t,result,'k','LineWidth',2);
print("result1.pdf","-dpdf","-vector");
