clc
clear
close all

iter = 1000;
alpha = 0.004;

landmarks = 5^3;

for i = 1:iter
 
  v_delta(:,1,i) = 0.0*[1*cos(i/200),1*sin(i/200+pi/6),0.5*cos(i/200)*sin(i/50+pi/7)] ;
  
  delta1_t(i) = 0.0*pi/200;
  delta2_t(i) = 0*pi/200;
  delta3_t(i) = 0.0*pi/200;
  
  
  if delta1_t(i) == 0 & delta2_t(i) == 0 & delta3_t(i) ==0;
     
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].';
     delta_hat_t(:,:,i) = hatoperate ([delta1_t(i) delta2_t(i) delta3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = 1;
  else
   
   
     delta_t(:,1,i) = [delta1_t(i) delta2_t(i) delta3_t(i)].'; %body-frame angular velocisy
     
     delta_hat_t(:,:,i) = [0 -delta_t(3,1,i) delta_t(2,1,i);delta_t(3,1,i) 0 -delta_t(1,1,i);-delta_t(2,1,i) delta_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_deltat(:,:,i) = expm(delta_hat_t(:,:,i));
  end;
  
end;



hov=3*pi/4;  %horizontal fieled of view in radians
vov=2*pi/3;  %vertical field of view in radians
focall = 1; 
focalh = 8;  %parameters for camera
Rad = 10;


for i = 1: landmarks
    y(:,i) = [randi([-20 20]); randi([-20 20]); randi([-20 20]); 1 ];
    
end;
j = -20;
k = -20;
h = -20;

% for i = 1: landmarks
%     
%        if h > 20
%           h = -20;
%           k = k +40/(nthroot(landmarks,3)-1);
%        end;   
%        if k >  20;
%           k = -20;
%           j = j + 40/(nthroot(landmarks,3)-1);
%        end;
%     y(:,i) = [j; k; h; 1 ];
%     h = h +40/(nthroot(landmarks,3)-1);
% end;
%sensor infomration
% P_t = [10;10;10];
 P_t = [-1;-2;-2];

fei_t = 0;
theta_t = 0;   %sensor angle in radians
miu_t = 0;

R_t(:,:,1) = rotation(fei_t,theta_t,miu_t);
SE_3t(:,:,1) = [R_t(:,:,1) P_t; [0 0 0] 1];

%observer information
fei_s =0;
theta_s =0;
miu_s = 0;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3

R_s = rotation(fei_s,theta_s,miu_s);

P_s = [0;0;0];

   
SE_3s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];
figure('position',[0 0 1000 800]);

Ut(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
Ui(:,:,1) = [ 0 0 0 0; 0 0 0 0 ; 0 0 0 0; 0 0 0 0 ];
R = rotation( 0.4, 0.7 , 0.1 );
P = [5 4 1 ];
T = [R P.'; [ 0 0 0] 1];

for i = 1:iter
        
                     
        
        TRt(:,:,i) = [SE_3t(1,1,i) SE_3t(1,2,i) SE_3t(1,3,i);
                 SE_3t(2,1,i) SE_3t(2,2,i) SE_3t(2,3,i);
                 SE_3t(3,1,i) SE_3t(3,2,i) SE_3t(3,3,i)];
             
        TPt(:,i) = [SE_3t(1,4,i) SE_3t(2,4,i) SE_3t(3,4,i)];
        
        TRs(:,:,i) = [SE_3s(1,1,i) SE_3s(1,2,i) SE_3s(1,3,i);
                 SE_3s(2,1,i) SE_3s(2,2,i) SE_3s(2,3,i);
                 SE_3s(3,1,i) SE_3s(3,2,i) SE_3s(3,3,i)];
             
        TPs(:,i) = [SE_3s(1,4,i) SE_3s(2,4,i) SE_3s(3,4,i)];  
       
        error4(:,:,i) =  SE_3t(:,:,i)*inv(SE_3s(:,:,i));
        error5(:,:,i) = SE_3s(:,:,i)*inv(T);
        error7(:,:,i) =  SE_3t(:,:,i)*inv(SE_3s(:,:,i));
        
        
        dv1(i) = abs(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
%         dv1(i) = norm(error4(:,:,i)-eye(4),'fro');
        dv2(i) = real(visiondistance(R,TRt(:,:,i),P.',TPt(:,i)));
        Ui(:,:,i) =  -inv(T)*1*projection(1/4*((error5(:,:,i) - error5(:,:,i).')))*T;
        Ut(:,:,i) = Ui(:,:,i)+ [delta_hat_t(:,:,i) v_delta(:,1,i);
                  [0,0,0]        0   ];
        SE_3t(:,:,i+1) = SE_3t(:,:,i)*expm(Ut(:,:,i));
        omega1(:,:,i)=zeros(4,4);
        omega2(:,:,i)=zeros(4,4);
        
        sensorshape(SE_3t(:,:,i),focall,focalh,hov,vov,Rad,1);
        obs1(:,i) = 0;
        counter = 0;
        if i >= 400 && i<= 500 
            
           for j = 1:landmarks
             cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
%             cframept(:,17) = [4;1.5;-1;1];
%             cframept(:,50) = [3.5;-2;1;1];
%             cframept(:,21) = [2;-2;1.5;1];
%             cframept(:,33) = [3;-1.5;-1;1];
%             cframept(:,42) = [2.5;-2;-1;1];
%             cframept(:,17) = [2.8;-2;-1.8;1];
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%                     
%             if cframept(1,j)^2 + cframept(2,j)^2 + cframept(3,j)^2 <= Rad^2; 
                      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                      omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3t(:,:,i)*cframept(:,j)-error4(:,:,i)*y(:,j))*transpose(y(:,j));    
%                      
                      counter = counter+1;
                      yprime(:,counter,i) = [y(:,j)];
                      robs(i) = rank(yprime(:,:,i));
            else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%                     
            end;
%           
           end;
        else
           for j = 1:landmarks
            cframept(:,j) =  inv(SE_3t(:,:,i)) *  y(:,j);
            
            if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
%                     
%             if cframept(1,j)^2 + cframept(2,j)^2 + cframept(3,j)^2 <= Rad^2; 
                      plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5);
                      omega1(:,:,i) = omega1(:,:,i)+alpha*(SE_3t(:,:,i)*cframept(:,j)-error4(:,:,i)*y(:,j))*transpose(y(:,j));    
%                       
                      counter = counter+1;
                      yprime(:,counter,i) = [y(:,j)];
                      robs(i) = rank(yprime(:,:,i));
            else 
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','k','MarkerSize',5);
%                     
            end;
%            
           end;
           
        end;
        obs1(:,i) = counter;
        omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
 
        
        se_3s(:,:,i) = Ui(:,:,i)-inv(SE_3s(:,:,i))*omega(:,:,i)*SE_3s(:,:,i); 
        SE_3s(:,:,i+1) = SE_3s(:,:,i)*expm(se_3s(:,:,i));
        

%         
      plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,1,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,1,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,1,i)],'-','color',[10 10 10]/255);

      plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,2,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,2,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,2,i)],'-','color',[100 100 100]/255);
      plot3([SE_3s(1,4,i),SE_3s(1,4,i)+1*SE_3s(1,3,i)],[SE_3s(2,4,i),SE_3s(2,4,i)+1*SE_3s(2,3,i)],[SE_3s(3,4,i),SE_3s(3,4,i)+1*SE_3s(3,3,i)],'-','color',[200 200 200]/255);
      hold off;
      pause(0.001)

        
end;
