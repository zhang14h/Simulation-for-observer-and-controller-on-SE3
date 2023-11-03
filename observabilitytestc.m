


          
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
       
        error4(:,:,i) =  SE_3s(:,:,i)*inv(SE_3t(:,:,i));
        
%         vd1(i) = abs(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
        dv1(i) = norm(SE_3t(:,:,i)-SE_3s(:,:,i),'fro');
        
        omega1(:,:,i)=zeros(4,4);
        omega2(:,:,i)=zeros(4,4);
        
%         sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,Rad,1);
%         obs1(:,i) = 0;
        counter = 0;
        if i >= 400 && i<= 700 
           
           for j = 1:landmarks
              cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
              deltay(:,j,i) = snr*[0.3*cos(i*j/2) 0.2*sin(i*j/5) 0.1*sin(i*j/4) 0].';
              
%             

            cframept(:,10) = [4;1.5;-1;1];
            cframept(:,150) = [3.5;-2;1;1];
            cframept(:,100) = [2;-2;1.5;1];
            cframept(:,91) = [3;-1.5;-1;1];
             cframept(:,11) = [2.5;-2;-1;1];
            cframept(:,67) = [2.8;-2;-1.8;1];
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%                     
%             if cframept(1,j)^2 + cframept(2,j)^2 + cframept(3,j)^2 <= Rad^2; 
%                       plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                      omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*y(:,j).';    
%                      
                      counter = counter+1;
                      yprime(:,counter,i) = [y(:,j)];
                      
                      robs(i) = rank(yprime(:,:,i));
            else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
                     
%                     
            end;
%           
           end;
        else
           for j = 1:landmarks
            cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
            deltay(:,j,i) = snr*[0.3*cos(i*j/2) 0.2*sin(i*j/5) 0.1*sin(i*j/4) 0].';
            
            
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%                     
%             if cframept(1,j)^2 + cframept(2,j)^2 + cframept(3,j)^2 <= Rad^2; 
%                       plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                      omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3s(:,:,i)*(cframept(:,j)+deltay(:,j,i))-y(:,j))*y(:,j).';    
%                       
                      counter = counter+1;
                      yprime(:,counter,i) = [y(:,j)];
                      
                      robs(i) = rank(yprime(:,:,i));
            else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
                     
%                     
            end;
%            
           end;
           
        end;
        num(i) = counter;
        omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
 
        se_3s(:,:,i) =  ([w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i)); 
        SE_3s(:,:,i+1) = SE_3s(:,:,i)*expm(se_3s(:,:,i));
        

%         
%       plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,1,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,1,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,1,i)],'-','color',[10 10 10]/255);
% 
%       plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,2,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,2,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,2,i)],'-','color',[100 100 100]/255);
%       plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,3,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,3,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,3,i)],'-','color',[200 200 200]/255);
%       hold off;
%       pause(0.001)

        
end;

figure('position',[0 0 1000 800]);
sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,Rad,1);
hold on;
for j = 1:landmarks
            cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
            deltay(:,j,i) = snr*[0.3*cos(i*j/2) 0.2*sin(i*j/5) 0.1*sin(i*j/4) 0].';
            
            
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
         
                      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
               
            else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
                  
            end;
end;            
%        if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%                       plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
%             else 
%                     plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%        end;
plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,1,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,1,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,1,i)],'-','color',[10 10 10]/255);
plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,2,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,2,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,2,i)],'-','color',[100 100 100]/255);
plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,3,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,3,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,3,i)],'-','color',[200 200 200]/255);


L2_1 = sqrt(sum(dv1*dv1.'));
Linf_1 = max(dv1);
display (L2_1);
display (Linf_1);
figure;
plot (dv1,':k','linewidth',1);
% figure;
% plot (robs,'-r','linewidth',1);
