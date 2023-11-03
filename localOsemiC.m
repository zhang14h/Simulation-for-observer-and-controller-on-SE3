% figure()
for i = 1:length(localo_semic)
    SE3_s(:,:,i) = [rotation(localo_semic(i,4),localo_semic(i,5),localo_semic(i,6)) [localo_semic(i,1); localo_semic(i,2); localo_semic(i,3)]; [0 0 0] 1];
    SE3_t(:,:,i) = [rotation(localo_semic(i,10),localo_semic(i,11),localo_semic(i,12)) [localo_semic(i,7); localo_semic(i,8); localo_semic(i,9)]; [0 0 0] 1];
    
    dv4(i) = norm(SE3_s(:,:,i)*inv(T)-eye(4),'fro');  %T - hatX
    dv5(i) = norm(SE3_t(:,:,i)*inv(T)-eye(4),'fro');  %T - X
    dv6(i) = norm(SE3_s(:,:,i)*inv(SE3_t(:,:,i))-eye(4),'fro');
    
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
   
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];

%     plot3([TPs(1,i),TPs(1,i)+1*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+1*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+1*TRs(1,3,i)],'--r');
%     hold on 
%     xlim([-4,4]);
%     ylim([-4,4]);
%     zlim([-1,1]);
%     plot3([TPs(1,i),TPs(1,i)+1*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+1*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+1*TRs(2,3,i)],'--g');
%     plot3([TPs(1,i),TPs(1,i)+1*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+1*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+1*TRs(3,3,i)],'--b');
%     plot3([TPt(1,i),TPt(1,i)+1*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+1*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+1*TRt(1,3,i)],'color',[10 10 10]/255);
%     plot3([TPt(1,i),TPt(1,i)+1*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+1*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+1*TRt(2,3,i)],'color',[100 100 100]/255);
%     plot3([TPt(1,i),TPt(1,i)+1*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+1*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+1*TRt(3,3,i)],'color',[200 200 200]/255);
%     plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5);
%     hold off  
%     pause(0.01)
end;  

L2_5 = sqrt(sum(dv5*dv5.')/length(localo_semic));
Linf_5 = max(dv5);

L2_6 = sqrt(sum(dv6*dv6.')/length(localo_semic));
Linf_6 = max(dv6);


figure();
for i = 1:length(localo_semic)
    px2(i) = SE3_t(1,4,i);
    py2(i) = SE3_t(2,4,i);
    pz2(i) = SE3_t(3,4,i)+0.6;
    hpx2(i) = SE3_s(1,4,i);
    hpy2(i) = SE3_s(2,4,i);
    hpz2(i) = SE3_s(3,4,i)+0.6;
end;    

plot3(px2,py2,pz2,'b');  
xlim([-5,5]);
ylim([-3,3]);
zlim([-1.5,1.5]);
hold on
for j = 1:length(y)
      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',7);
 end; 
plot3(hpx2,hpy2,hpz2,'-.k');


figure()  
  p3 = plot (dv5,'-r','linewidth',1); m3 = "Ec with Kc = 0";
%   title('attitude and position control performance with Kc = 0');
%  
%   txt = {'dV = ' ,k2};
%   text(iter/2,10,txt);
  grid
  hold on
  p4 = plot (dv6,'-.g','linewidth',1); m4 = "Eo with Ko = 0";
  title('Closed-loop Performance ');
  xlabel('Time');
  ylabel('Error');
  legend([p3,p4],[m3,m4]);