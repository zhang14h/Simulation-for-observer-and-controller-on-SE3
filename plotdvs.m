% observabilitytestc;
% visualwithoutobsdamping;
% visualwithobs;
% visualwithobsdgamma;
% visualwithobslargergamma;

figure;
% p1 = plot(dv1,'-r'); m1 = "OB1";
hold on;
% p2 = plot(dv2,'-.g'); m2 = "OB2";
hold on;
 p3 = plot(dv3,'--b'); m3 = "\eta = 4,K_o = 1.5*I";
hold on;
 p4 = plot(dv4,'-k'); m4 = "\eta = 6,K_o = 1.5*I";
hold on;
p5 = plot(dv5,'-.m'); m5 = "\eta = 8,K_o = 1.5*I";
% legend([p1,p5],[m1,m5]); 
% legend([p1,p2,p4],[m1,m2,m4]);
legend([p3,p4,p5],[m3,m4,m5]);
% legend([p1,p2,p3,p4,p5],[m1,m2,m3,m4,m5]);
% title('w = 0.00, sr = 0 and no detection error');
xlabel('Time');
ylabel('Estimation error ||Eo||_F');
figure;
plot (num,'-r','linewidth',1);
figure
plot (robs,'-b','linewidth',1);
figure
for i = 1:iter
    px(i) = SE_3t(1,4,i);
    py(i) = SE_3t(2,4,i);
    pz(i) = SE_3t(3,4,i);
    hpx(i) = SE_3s(1,4,i);
    hpy(i) = SE_3s(2,4,i);
    hpz(i) = SE_3s(3,4,i);
end;    
 plot3(px,py,pz);
%  hold on;
%  plot3(hpx,hpy,hpz);
%  xlim([-3,3]);
%  ylim([-3,3]);
%  zlim([-3,3])
 figure
 for j = 1:landmarks
            
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
                    hold on;
%                    
end;
hold on;
% plot3(px,py,pz);
%  xlim([-3,3]);
%  ylim([-3,3]);
%  zlim([-3,3])