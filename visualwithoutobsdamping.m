% close (2);
% figure('position',[0 0 1000 800]);
Q = 0.01 * [1 4 3 1; 4 3 1 1; 3 1 2 1; 1 1 1 3];
Q=  Q+1.5*eye(4);
for i = 1:iter
        se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
                  [0,0,0]        0   ];
                     
        SE_3t(:,:,i+1) = SE_3t(:,:,i)*expm(se3_t(:,:,i));
        TRt(:,:,i) = [SE_3t(1,1,i) SE_3t(1,2,i) SE_3t(1,3,i);
                 SE_3t(2,1,i) SE_3t(2,2,i) SE_3t(2,3,i);
                 SE_3t(3,1,i) SE_3t(3,2,i) SE_3t(3,3,i)];
             
        TPt(:,i) = [SE_3t(1,4,i) SE_3t(2,4,i) SE_3t(3,4,i)];
        
        TRs(:,:,i) = [SE_3s(1,1,i) SE_3s(1,2,i) SE_3s(1,3,i);
                 SE_3s(2,1,i) SE_3s(2,2,i) SE_3s(2,3,i);
                 SE_3s(3,1,i) SE_3s(3,2,i) SE_3s(3,3,i)];
             
        TPs(:,i) = [SE_3s(1,4,i) SE_3s(2,4,i) SE_3s(3,4,i)];  
%         
        error4(:,:,i) =  SE_3s(:,:,i)*inv(SE_3t(:,:,i));
        
%         vd3(i) = abs(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
        dv2(i) = norm(SE_3t(:,:,i)-SE_3s(:,:,i),'fro');
        omega1(:,:,i)=zeros(4,4);
        omega2(:,:,i)=zeros(4,4);
        
%         sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,Rad,1);
        counter = 0;
        p = 0;
%         obs3(:,i) = 0;
        if i >= 400 && i<= 700 
          for j = 1:landmarks
           cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
           cframeps(:,j) =  inv(SE_3s(:,:,i)) *  y(:,j);
           
           
            cframept(:,10) = [4;1.5;-1;1];
            cframept(:,150) = [3.5;-2;1;1];
            cframept(:,100) = [2;-2;1.5;1];
            cframept(:,91) = [3;-1.5;-1;1];
             cframept(:,11) = [2.5;-2;-1;1];
            cframept(:,67) = [2.8;-2;-1.8;1];
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
                
%                 if cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                    omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*transpose(y(:,j));             
                    omega2(:,:,i) = omega2(:,:,i)+alpha*0.3*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*transpose(y(:,j))/landmarks;
%                     omega2(:,:,i) = omega2(:,:,i)+alpha*exp(0.8-(diss-dist))*(cframept(:,j))*transpose(y(:,j));
                    counter = counter+1;
%
                yprime(:,counter,i) = [y(:,j)];
                robs(i) = rank(yprime(:,:,i));
             else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%                    
%                     
             end;
          end; 
        else
          for j = 1:landmarks
            cframeps(:,j) =  inv(SE_3s(:,:,i)) *  y(:,j);
            cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
            
            
%             focalmprime(:,j) = inv(SE_3s(:,:,i)) * focallm;
            
             if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
                
%                 if cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                    omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*transpose(y(:,j));             
                    omega2(:,:,i) = omega2(:,:,i)+alpha*0.3*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*transpose(y(:,j))/landmarks;
%                     omega2(:,:,i) = omega2(:,:,i)+alpha*exp(0.8-(diss-dist))*(cframept(:,j))*transpose(y(:,j));
                    counter = counter+1;
%
                yprime(:,counter,i) = [y(:,j)];
                robs(i) = rank(yprime(:,:,i));
             else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%                    
%                     
             end;
            
%           obs3(:,i) = obs3(:,i) + dis(j);
          end;
        end;
       num(i) = counter;
        
        omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
        
%         se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
%                   [0,0,0]        0   ]-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i))-inv(SE_3s(:,:,i))*projection (1*[1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1]*1/4*((omega2(:,:,i) - omega2(:,:,i).')))*SE_3s(:,:,i); 
              
        se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i))-inv(SE_3s(:,:,i))*projection (1/4*((Q*omega2(:,:,i) - omega2(:,:,i).'*Q.')))*SE_3s(:,:,i); 
              
        SE_3s(:,:,i+1) = SE_3s(:,:,i)*expm(se_3s(:,:,i));
        
%       sensorshape(SE_3s(:,:,i),focall,focalh,hov,vov,Rad,2);
%       hold off;
%       pause(0.001)
%       
        
end;
% figure('position',[0 0 1000 800]);
% sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,Rad,1);
% for j = 1:landmarks
%             cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
%             deltay(:,j,i) = snr*[0.3*cos(i*j/2) 0.2*sin(i*j/5) 0.1*sin(i*j/4) 0].';
%             
%             
%             if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%          
%                       plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
%                
%             else 
% %                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%                   
%             end;
% end;            
% hold on;
% %        if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
% %                       plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
% %             else 
% %                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
% %        end;
% plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,1,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,1,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,1,i)],'-','color',[10 10 10]/255);
% plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,2,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,2,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,2,i)],'-','color',[100 100 100]/255);
% plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,3,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,3,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,3,i)],'-','color',[200 200 200]/255);

L2_2 = sqrt(sum(dv2*dv2.'));
Linf_2 = max(dv2);
display (L2_2);
display (Linf_2);
figure;
plot (dv2,':k','linewidth',1);
% figure;
% plot (robs,'-r','linewidth',1);