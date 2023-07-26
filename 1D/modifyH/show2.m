% s = 0.3, 0.6, 0.9
hold on
plot(node03,psih03,'k','linewidth',2);
plot(node03,uh3,Color=[0 0.4470 0.7410], Linewidth=1)
plot(node03,uh32,':',Color=[0 0.4470 0.7410], Linewidth=1.5)
plot(node06,uh6,Color=[0.4660 0.6740 0.1880], Linewidth=1)
plot(node06,uh62,':',Color=[0.4660 0.6740 0.1880], Linewidth=1.5)
plot(node09,uh9,Color=[0.6350 0.0780 0.1840], Linewidth=1)
plot(node09,uh92,':',Color=[0.6350 0.0780 0.1840], Linewidth=1.5)
legend('\psi','s=0.3 improve','s=0.3','s=0.6 improve','s=0.6','s=0.9 improve','s=0.9')
xlabel('x'); %ylabel('');
hold off