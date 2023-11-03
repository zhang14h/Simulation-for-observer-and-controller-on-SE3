newobserver;
for i = 1:iter;
   
   
%     dV(i) = real(visiondistance(R_t(:,:,i),R_s(:,:,i),P_t,P_s));
%     R_t(:,:,i+1) = exp_t(:,:,i)*R_t(:,:,i);
%     
%     error(:,:,i) = R_s(:,:,i).'*(R_t(:,:,i)+0.02*randi(1,3));
%     
%     proj= 0.5*1/2*((error(:,:,i) - error(:,:,i).'));
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0]+proj;
%     if A == 0
%       exp_s = 1;
%     else   
% %        a_s = sqrt(A(2,1)^2+A(3,1)^2+A(3,2)^2);
% %        A = A/a_s;
% %        exp_s = [1 0 0;0 1 0;0 0 1] + A*sin(a_s)+ A^2*(1-cos(a_s));
%      exp_s = expm(A);
%     end;
%     
%     R_s(:,:,i+1) = exp_s*R_s(:,:,i);
%     
    

% end of so3 estimation;    
    
    se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
                  [0,0,0]        0   ];
    
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(se3_t(:,:,i));
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)];
    SE3_tm(:,:,i) = SE3_t(:,:,i)*SE3_mo(:,:,i);
    
    TRtm(:,:,i) = [SE3_tm(1,1,i) SE3_tm(1,2,i) SE3_tm(1,3,i);
                 SE3_tm(2,1,i) SE3_tm(2,2,i) SE3_tm(2,3,i);
                 SE3_tm(3,1,i) SE3_tm(3,2,i) SE3_tm(3,3,i)];
             
    TPtm(:,i) = [SE3_tm(1,4,i) SE3_tm(2,4,i) SE3_tm(3,4,i)];
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
             
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];  
    dv2(i) = real(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
    error4(:,:,i) =  SE3_s(:,:,i)*inv(SE3_tm(:,:,i));
     omega1(:,:,i)=zeros(4,4);
    omega2(:,:,i)=zeros(4,4);
    for p = 1:14
       omega1(:,:,i) = omega1(:,:,i)+0.01*(error4(:,:,i)*y(:,p)-y(:,p))*transpose(y(:,p));
       
    end;
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
    omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
    se3_s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-inv(SE3_s(:,:,i))*omega(:,:,i)*SE3_s(:,:,i)); 
    SE3_s(:,:,i+1) = SE3_s(:,:,i)*expm(se3_s(:,:,i));
    [roll_s1(i),roll_s2(i),roll_s3(i)] = derotation(TRs(:,:,i));
    roll_s(:,i) =  [roll_s1(i),roll_s2(i),roll_s3(i)];
%     
    [roll_t1(i),roll_t2(i),roll_t3(i)] = derotation(TRt(:,:,i));
    roll_t(:,i) = [roll_t1(i),roll_t2(i),roll_t3(i)];
    
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(1,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(1,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(1,3,i)],'r');
%     hold on 
%     xlim([-20 20]);
%     ylim([-20,20]);
%     zlim([-20,20]);
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(2,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(2,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(2,3,i)],'g');
%     plot3([TPs(1,i),TPs(1,i)+10*TRs(3,1,i)],[TPs(2,i),TPs(2,i)+10*TRs(3,2,i)],[TPs(3,i),TPs(3,i)+10*TRs(3,3,i)],'b');
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(1,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(1,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(1,3,i)],'color',[10 10 10]/255);
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(2,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(2,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(2,3,i)],'color',[100 100 100]/255);
%     plot3([TPt(1,i),TPs(1,i)+10*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(3,3,i)],'color',[200 200 200]/255);
%     hold off  
%     pause(0.05)
%     
end 
k2 = 0;
for i=1:iter
k2 = k2+(dv2(i));
end;
p1 = plot (dv1,'-r','linewidth',1); m1 = "observer in (8)";

  hold on;
  p2 = plot (dv2,'--g','linewidth',1); m2 = "observer with k = 1.5";
  
  hold on;
%   p3 = plot (dv3,'.b','linewidth',1); m3 = "observer with k = 3";
  
   legend ([p1,p2],[m1,m2]); 
% legend ([p1,p3],[m1,m3]);
%  title('Performance of the proposed observerscwith k = 1.5 and k =3');
 title('Performance of the observer in (8) and the proposed observers');
 grid
  xlabel('Time')
  ylabel('Estimation error')