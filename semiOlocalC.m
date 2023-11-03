% figure()
for i = 1:length(semio_localc)
    SE3_s(:,:,i) = [rotation(semio_localc(i,4),semio_localc(i,5),semio_localc(i,6)) [semio_localc(i,1); semio_localc(i,2); semio_localc(i,3)]; [0 0 0] 1];
    SE3_t(:,:,i) = [rotation(semio_localc(i,10),semio_localc(i,11),semio_localc(i,12)) [semio_localc(i,7); semio_localc(i,8); semio_localc(i,9)]; [0 0 0] 1];
    
    dv7(i) = norm(SE3_s(:,:,i)*inv(T)-eye(4),'fro');  %T - hatX
    dv8(i) = norm(SE3_t(:,:,i)*inv(T)-eye(4),'fro');  %T - X
    dv9(i) = norm(SE3_s(:,:,i)*inv(SE3_t(:,:,i))-eye(4),'fro');
    
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

L2_8 = sqrt(sum(dv8*dv8.')/length(semio_localc));
Linf_8 = max(dv8);

L2_9 = sqrt(sum(dv9*dv9.')/length(semio_localc));
Linf_9 = max(dv9);

figure();
for i = 1:length(semio_localc)
    px3(i) = SE3_t(1,4,i);
    py3(i) = SE3_t(2,4,i);
    pz3(i) = SE3_t(3,4,i)+0.6;
    hpx3(i) = SE3_s(1,4,i);
    hpy3(i) = SE3_s(2,4,i);
    hpz3(i) = SE3_s(3,4,i)+0.6;
end;    

plot3(px3,py3,pz3,'b');  
xlim([-5,5]);
ylim([-3,3]);
zlim([-1.5,1.5]);
hold on
for j = 1:length(y)
      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',7);
 end; 
plot3(hpx3,hpy3,hpz3,'-.k');


figure()  
  p5 = plot (dv8,'-r','linewidth',1); m5 = "Ec with Kc = 0";
%   title('attitude and position control performance with Kc = 0');
%  
%   txt = {'dV = ' ,k2};
%   text(iter/2,10,txt);
  grid
  hold on
  p6 = plot (dv9,'-.g','linewidth',1); m6 = "Eo with Ko = 0";
  title('Closed-loop Performance ');
  xlabel('Time');
  ylabel('Error');
  legend([p5,p6],[m5,m6]);