Ut(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
Ui(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
SE3_s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];
SE3_t(:,:,1) = [R_t(:,:,1) P_t;
                [0,0,0]            1     ];

P =  [3 1 1 1; 2 3 1 1; 3 1 2 1; 1 1 1 1];
P=  20*eye(4)+ 20*[0 1 0 0 ; 0 0 1 0; 0 0 0 0; 0 0 0 0];
% figure
for i = 1:iter;
     TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)]; 
   
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];
    
    
    error4(:,:,i) =  SE3_s(:,:,i)*inv(T);
    error6(:,:,i) = SE3_t(:,:,i)*inv(T);
    error7(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    
%     dv4(i) = norm(T-SE3_s(:,:,i),'fro');  %T - hatX
%     dv5(i) = norm(T-SE3_t(:,:,i),'fro');  %T - X
%     dv6(i) = norm(SE3_s(:,:,i)-SE3_t(:,:,i),'fro');
    
    dv4(i) = norm(SE3_s(:,:,i)*inv(T)-eye(4),'fro');  %T - hatX
    dv5(i) = norm(SE3_t(:,:,i)*inv(T)-eye(4),'fro');  %T - X
    dv6(i) = norm(SE3_s(:,:,i)*inv(SE3_t(:,:,i))-eye(4),'fro');
    
    omega1(:,:,i)=zeros(4,4);
    omega2(:,:,i)=zeros(4,4);
    omega3(:,:,i)=zeros(4,4);
    
    omega(:,:,i) =   4*projection(1/4*((error4(:,:,i) - error4(:,:,i).')));
%     omega4(:,:,i) =   0.5*inv(SE3_s(:,:,i))*projection(1/4*((error4(:,:,i) - eye(4))-(error4(:,:,i) - eye(4)).'))*SE3_s(:,:,i);
     omega4(:,:,i) =   inv(T)*projection(1/4*(P*error4(:,:,i)-error4(:,:,i).'*P.'))*T;
    
      
%     SE3_t(:,:,70) =[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 0];
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
%     Ui(:,:,i) =  -inv(SE3_s(:,:,i))*(omega(:,:,i))*SE3_s(:,:,i);
    Ui(:,:,i) =  -0.1*inv(T)*(omega(:,:,i))*(T)-0.1*omega4(:,:,i);
    Ut(:,:,i) =  Ui(:,:,i)+[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ];
%     U(:,:,i) =  -(omega(:,:,i))-[delta_hat_t(:,:,i) v_delta(:,1,i);[0,0,0]    0 ]-omega4(:,:,i);
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(Ut(:,:,i));
    error5(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
    for p = 1:lengthy
       cframept(:,p) =  inv(SE3_t(:,:,i)) *  [y(1,p) y(2,p) y(3,p) 1].';
       omega2(:,:,i) = omega2(:,:,i)+0.003*(SE3_s(:,:,i)*(cframept(:,p)+deltay(:,p,i))-[y(1,p) y(2,p) y(3,p) 1].')*transpose([y(1,p) y(2,p) y(3,p) 1].');
       
    end;
    omega_t(:,:,i) =   projection(1/4*((omega2(:,:,i) - omega2(:,:,i).')));
%     se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_t(:,:,i)*SE3_s(:,:,i)-projection (1.1*[1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1]*1/4*((omega3(:,:,i) - omega3(:,:,i).')));
    se3_s = Ui(:,:,i)-inv(SE3_s(:,:,i))*omega_t(:,:,i)*SE3_s(:,:,i); 
  
    SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s);
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(TRs(:,:,i));
    roll_s(:,i) =  [roll_s1(i),roll_s2(i),roll_s3(i)];
    
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(1,3,i)],'--r');
%     hold on 
%     
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(2,3,i)],'--g');
%     plot3([TPs(1,i),TPs(1,i)+5*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+5*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+5*TRs(3,3,i)],'--b');
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(1,3,i)],'color',[10 10 10]/255);
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(2,3,i)],'color',[100 100 100]/255);
%     plot3([TPt(1,i),TPt(1,i)+5*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+5*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+5*TRt(3,3,i)],'color',[200 200 200]/255);
%     plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5);
%     hold off  
%     pause(0.01)

end;


L2_5 = sqrt(sum(dv5*dv5.')/iter);
Linf_5 = max(dv5);

L2_6 = sqrt(sum(dv6*dv6.')/iter);
Linf_6 = max(dv6);



figure();
for i = 1:iter
    px2(i) = SE3_t(1,4,i);
    py2(i) = SE3_t(2,4,i);
    pz2(i) = SE3_t(3,4,i);
    hpx2(i) = SE3_s(1,4,i);
    hpy2(i) = SE3_s(2,4,i);
    hpz2(i) = SE3_s(3,4,i);
end;    
 plot3(px2,py2,pz2,'b');  
 hold on
for j = 1:length(y)
      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',7);
 end;     
 plot3(hpx2,hpy2,hpz2,'-.k');

% figure()  
%   p1 = plot (dv2,'-r','linewidth',1); m1 = " observer in (8)";
%   hold on;
%   p2 = plot (dv5,'--g','linewidth',1); m2 = "observer with k = 1.5";
%   legend([p1,p2],[m1,m2]);
%   title('attitude and position control performance');
%   xlabel('Time');
%   ylabel('Tracking error');
% %   txt = {'dV = ' ,k2};
% %   text(iter/2,10,txt);
%   grid
% 
% figure()  
%   p3 = plot (dv3,'-r','linewidth',1); m3 = "observer in (8)";
%   hold on;
%   p4 = plot (dv6,'--g','linewidth',1); m4 = "observer with k = 1.5";
%   legend([p3,p4],[m3,m4]);
%   title('Performance of the observer in (8) and the proposed observers');
%   xlabel('Time');
%   ylabel('Estimation error');
% %   txt = {' = ' ,k3};
% %   text(iter/2,3,txt);
%   grid
  