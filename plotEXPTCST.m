close all;

figure('position',[0 0 1000 600]);

p1 = plot (dv2,'-r','linewidth',3); m1 = "l_c = 0.4, k_c = 2";
%  m1 = "Case 1";
hold on;
p3 = plot (dv5,'-g','linewidth',3); m3 = "l_c = 0.9, k_c = 2";
% m3 = "Case 2";
p5 = plot (dv8,'-.b','linewidth',3); m5 = "l_c = 1, k_c = 2";
% m5 = "Case 3";
p7 = plot (dv11,'-.c','linewidth',3); m7 = "l_c = 2, k_c = 2";
% m7 = "Case 4";


% p1 = plot (dv2,'-.c','linewidth',3); m1 = "Case 4";
% hold on;
% p3 = plot (dv5,'-.b','linewidth',3); m3 = "Case 3";
% p5 = plot (dv8,'-g','linewidth',3); m5 = "Case 2";
% p7 = plot (dv11,'-r','linewidth',3); m7 = "Case 1";


legend([p1,p3,p5,p7],[m1,m3,m5,m7]);
% legend([p7,p5,p3,p1],[m7,m5,m3,m1]);

grid on;
% set(gca,'fontsize',25,'FontWeight', 'bold') 
set(gca,'fontsize',25) 
% 
xlabel('Time', 'Fontsize', 25,'FontAngle','italic', 'FontWeight', 'bold');
ylabel('Closed Loop ||Ec||_F', 'Fontsize', 25, 'FontAngle','italic', 'FontWeight', 'bold');
xlim([0,iter]);

figure('position',[0 0 1000 600]);

p2 = plot (dv3,'-r','linewidth',3); m2 = "\it{l}_o = 0.04, k_o = 2";
% m2 = "Case 1";
hold on;
p4 = plot (dv6,'-g','linewidth',3); m4 = "\it{l}_o = 0.07, k_o = 2";
% m4 = "Case 2";
p6 = plot (dv9,'-.b','linewidth',3); m6 = "\it {l}_o = 0.1, k_o = 2";
% m6 = "Case 3";
p8 = plot (dv12,'-.c','linewidth',3); m8 = "\it {l}_o = 0.2, k_o = 2";
% m8 = "Case 4";

% p2 = plot (dv3,'-.c','linewidth',3); m2 = "Case 4";
% hold on;
% p4 = plot (dv6,'-.b','linewidth',3); m4 = "Case 3";
% p6 = plot (dv9,'-g','linewidth',3); m6 = "Case 2";
% p8 = plot (dv12,'-r','linewidth',3); m8 = "Cass 1";


legend([p2,p4,p6,p8],[m2,m4,m6,m8]);
% legend([p8,p6,p4,p2],[m8,m6,m4,m1]);


grid on;
% set(gca,'fontsize',25,'FontWeight', 'bold') 
set(gca,'fontsize',25) 
% 
xlabel('Time', 'Fontsize', 25,'FontAngle','italic', 'FontWeight', 'bold');
ylabel('||Eo||_F', 'Fontsize', 25, 'FontAngle','italic', 'FontWeight', 'bold');
xlim([0,iter]);

figure('position',[0 0 1400 600]);
d1 = plot3(px1,py1,pz1,'-r','linewidth',3); k1= "l_o = 0.03, k = 0";

hold on
d2 = plot3(px2,py2,pz2,'-g','linewidth',3); k2 = "l_o = 0.004, k = 0.1";
d3 = plot3(px3,py3,pz3,'-.b','linewidth',3); k3 = "l_o = 0.007, k = 0.1";
d4 = plot3(px4,py4,pz4,'-.c','linewidth',3); k4 = "l_o = 0.01, k = 0.1";

d5 = plot3([T(1,4),T(1,4)+1*T(1,1)],[T(2,4),T(2,4)+1*T(1,2)],[T(3,4),T(3,4)+1*T(1,3)],'k','linewidth',7); k5 = "T";
     
    
plot3([T(1,4),T(1,4)+1*T(2,1)],[T(2,4),T(2,4)+1*T(2,2)],[T(3,4),T(3,4)+1*T(2,3)],'k','linewidth',7);
plot3([T(1,4),T(1,4)+1*T(3,1)],[T(2,4),T(2,4)+1*T(3,2)],[T(3,4),T(3,4)+1*T(3,3)],'k','linewidth',7);
hold on
for j = 1:length(y)
      plot3(y(1,j),y(2,j),y(3,j),'.','Color','r','MarkerSize',40);
 end;  
% for j = 1:landmarks
%       plot3(y(1,j),y(2,j),y(3,j),'filled','Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1);
% end; 


legend([d1,d2,d3,d4,d5],[k1,k2,k3,k4,k5]);

grid on;
% set(gca,'fontsize',25,'FontWeight', 'bold') 
set(gca,'fontsize',30) 
