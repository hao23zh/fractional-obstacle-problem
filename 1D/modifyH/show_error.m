% uniform
subplot(1,2,1)
plot(x1,y1,'o-r',x1,y1(1)+(x1-x1(1))*0.3,'--r')
hold on
plot(x2,y2,'*-k',x2,y2(1)+(x2-x2(1))*0.6,'--k')
plot(x3,y3,'s-b',x3,y3(1)+(x3-x3(1))*0.9,'--b')

xlabel("log(N^{-1})")
ylabel("log(||u_h-u_I||_{L_\infty})")
title("uniform")
legend('s=0.3','0.3','s=0.6','0.6','s=0.9','0.9')
hold off

% graded
subplot(1,2,2)
plot(x4,y4,'o-r',x4,y4(1)+(x4-x4(1))*0.5,'--r')
hold on
% hold on
% p = polyfit(x4,y4,1);disp(p(1))
% x14 = linspace(-8.5,-3.5);
% z4 = polyval(p,x14);
% plot(x14,z4,':r')

plot(x5,y5,'*-k',x5,y5(1)+(x5-x5(1))*1.4,'--k')
% p = polyfit(x5,y5,1); disp(p(1))
% x14 = linspace(-8.5,-3.5);
% z5 = polyval(p,x14);
% plot(x14,z5,':k')

plot(x6,y6,'s-b',x6,y6(1)+(x6-x6(1))*1,'--b')
% p = polyfit(x6,y6,1); disp(p(1))
% x14 = linspace(-8.5,-3.5);
% z6 = polyval(p,x14);
% plot(x14,z6,':b')

xlabel("log(N^{-1})")
ylabel("log(||u_h-u_I||_{L_\infty})")
title("graded")
%legend('s=0.3','0.3','0.28','s=0.6','0.6','0.43','s=0.9','0.55','0.35')
legend('s=0.3','0.5','s=0.6','1.4','s=0.9','1')
hold off
