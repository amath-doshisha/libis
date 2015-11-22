%function TEST_plot_double
clear all
close all

x=0:0.1:10;
y=sin(x);
z=sin(x)+0.5*sin(2*x)+0.25*sin(4*x);

S = get(0, 'ScreenSize');
figure('Position',S);
plot(x,y,'k-o',x,z,'k-x','LineWidth',3,'MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor',[1 1 1]);
grid on;
grid minor;
xlabel('x','FontName','Times','FontSize',24,'FontAngle','italic','FontWeight','bold','Color','k');
ylabel('y','FontName','Times','FontSize',24,'FontAngle','italic','FontWeight','bold','Color','k');
title('Waves','FontName','Times','FontSize',24,'FontAngle','normal','FontWeight','bold','Color','k');
