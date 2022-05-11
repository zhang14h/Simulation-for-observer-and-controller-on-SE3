% visualwithoutobs;
% close (4);
figure();
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
        error4(:,:,i) =  SE_3t(:,:,i)*inv(SE_3s(:,:,i));
        dv4(i) =  abs(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
%         dv4(i) = norm(error4(:,:,i)-eye(4),'fro');
        omega1(:,:,i)=zeros(4,4);
        omega2(:,:,i)=zeros(4,4);
        
        sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,1);
        counter = 0;
        p = 0;
        obs4(i)=0;
        if i >= 400 && i<= 500 
          for j = 1:landmarks
           cframeps(:,j) =  inv(SE_3s(:,:,i)) *  y(:,j);
           cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
%            cframept(:,17) = [4;1.5;-1;1];
%             cframept(:,50) = [3.5;-2;1;1];
%             cframept(:,21) = [2;-2;1.5;1];
%             cframept(:,33) = [3;-1.5;-1;1];
%             cframept(:,42) = [2.5;-2;-1;1];
%             cframept(:,19) = [2.8;-2;-1.8;1];
             if cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                
                if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                    counter = counter+1;
                    p(counter) = j;
                    dis(j) = sqrt((cframept(1,j)-(focall+focalh)/2)^2+ cframept(2,j)^2+ cframept(3,j)^2+0.01);
                else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','b','MarkerSize',5); %avoid points missing display
                    dis(j) = 0;
                end;
             else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','b','MarkerSize',5);
                    dis(j) = 0;
                    
             end;
            
%            obs4(:,i) = obs4(:,i) + dis(j);
           end;
          
        else
          for j = 1:landmarks
            cframeps(:,j) =  inv(SE_3s(:,:,i)) *  y(:,j);
            cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
            
             if cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                
                if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                    counter = counter+1;
                    p(counter) = j;
                    dis(j) = sqrt((cframept(1,j)-(focall+focalh)/2)^2+ cframept(2,j)^2+ cframept(3,j)^2+0.01);
                else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','b','MarkerSize',5); %avoid points missing display
                    dis(j) = 0;
                end;
             else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','b','MarkerSize',5);
                    dis(j) = 0;
%                     
             end;
            
%           obs4(:,i) = obs4(:,i) + dis(j);
          end;
        end;
        obs4(:,i) = counter;
        for k = 1:counter;  
                    omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3t(:,:,i)*cframept(:,p(k))-error4(:,:,i)*y(:,p(k)))*transpose(y(:,p(k)));
                    omega2(:,:,i) = omega2(:,:,i)+0.003*(SE_3t(:,:,i)*cframept(:,p(k))-error4(:,:,i)*y(:,p(k)))*transpose(y(:,p(k)))/norm(y(:,p(k))*y(:,p(k)).',2);
        end; 
        omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
%         if counter > 3;
% %    
%             se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
%                   [0,0,0]        0   ]-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i));
%               
%         else
%             se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
%                   [0,0,0]        0   ]);
%         end;
        se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i))-projection (3*[1 0 0 0;0 1 0 0;0 0 1 0; 0 0 0 1]*1/4*((omega2(:,:,i).' - omega2(:,:,i))));
        SE_3s(:,:,i+1) = SE_3s(:,:,i)*expm(se_3s(:,:,i));
        
      sensorshape(SE_3s(:,:,i),focall,focalh,hov,vov,2);
      hold off;
      pause(0.001)
      
        
end;

L2_2 = std(dv2);
Linf_2 = max(dv2);
display (L2_2);
display (Linf_2);


figure();
% plot (dv1,'-r','linewidth',1);
% hold on;
plot (dv4,':k','linewidth',1);
hold on;
plot (obs4,'-r','linewidth',1);
