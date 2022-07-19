clc
clear
close all
a0 =1;
iter = 2000;
for i = 1:iter
    if rem(i,100) <=40;
       q(i) = 0.05*pi/5250;
    else
       q(i) = 0;
    end;
  w1_t(i) =0;
  w2_t(i) = 2*pi/500;
  w3_t(i) = 2*pi/500;
  v_t(:,1,i) = 0.3*[1*pi/2000*cos(2*i*pi/2000) 1*pi/iter*sin(2*i*pi/2000)  1*pi/iter*sin(2*i*pi/2000)].';
%   v_t(:,1,i) = 0.3*[1*pi*5/2000 1*pi*5/2000  1*pi*5/2000].';
  v_delta(:,1,i) = 11*0.01*[1*cos(i/2000),1*sin(i/3000+pi/6),0.5*cos(i/1000)*sin(i/7000+pi/7)] ;

  
  delta1_t(i) = 5*0.01*q(i);
  delta2_t(i) = 5*0.01*pi/6000;
  delta3_t(i) = 5*0.03*pi/8000;
  
  
  if w1_t(i) == 0 & w2_t(i) == 0 & w3_t(i) ==0;
     
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].';
     w_hat_t(:,:,i) = hatoperate ([w1_t(i) w2_t(i) w3_t(i)].'); %lefs-righs invarians coefficiens,so3;
     exp_t(:,:,i) = 1;
  else
   
   
     w_t(:,1,i) = [w1_t(i) w2_t(i) w3_t(i)].'; %body-frame angular velocisy
     
     w_hat_t(:,:,i) = [0 -w_t(3,1,i) w_t(2,1,i);w_t(3,1,i) 0 -w_t(1,1,i);-w_t(2,1,i) w_t(1,1,i) 0]; %lefs-righs invarians coefficiens,so3;
     exp_t(:,:,i) = expm(w_hat_t(:,:,i));
  end;
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

fei_t= 0;
theta_t= 0.3;
miu_t= 0.2;
% fei_t = 0;
% theta_t = pi;
% miu_t = pi;
rpy = [fei_t;
       theta_t;
       miu_t];

R_t(:,:,1) = rotation(fei_t,theta_t,miu_t); %SO3 


px_t = 0.1;
py_t = 0.2;
pz_t = 0.5;
P_t= [px_t;
      py_t;
      pz_t];
   

         
SE3_t(:,:,1) = [R_t(:,:,1) P_t;
                [0,0,0]            1     ];


% SE3_t(:,:,1) = [0.1313  -0.0230    0.9911 2; 0.3377   0.9410   -0.0230 2; -0.9320    0.3377    0.1313 2; 0 0 0 1];  

fei_s = 0;
theta_s = 0;
miu_s = 0;

R_s(:,:,1) = rotation(fei_s,theta_s,miu_s); %SO3
% R_s(:,:,1) = [-R_t(1,1,1) -R_t(1,2,1) -R_t(1,3,1);
%               -R_t(2,1,1) -R_t(2,2,1) -R_t(2,3,1)
%               -R_t(3,1,1) -R_t(3,2,1) -R_t(3,3,1)];

px_s =0;
py_s =0;
pz_s =0;
P_s= [px_s;
       py_s;
       pz_s];

v_s(:,1) = [0,0,0].';
   
SE3_s(:,:,1) = [R_s(:,:,1) P_s;
               [0,0,0]      1     ];



y(:,1) = [2;2;2;1]; y(:,2) = [4; 1; 1; 1]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;1]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;1]; 
y(:,7) = [2; 1; 2; 1]; y(:,8) = [3;4;5;1] ; y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

% y(:,1) = [2;2;2;1]; y(:,2) = [4;4;4;1]; y(:,3) = [1; 2;3;1]; y(:,4) = [3;2;1;1]; y(:,5) = [-1;3;1;1]; y(:,6) = [2;2;3;1]; 
% y(:,7) = [2; 1; 2; 1]; y(:,8) = [3;4;5;1] ; y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
% y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [1.9;3.8;5.7;1]; 
% y(:,7) = [2.2;4.4;6.6;1]; y(:,8) = [3;5;4;1] ;  y(:,9) = [ 5;3;1;1] ; y(:,10) = [7;8;9;1]; y(:,11) = [2;5;7;1]; y(:,12) = [3;6;8;1];
% y(:,13) = [1;5;9;1]; y(:,14) = [3;7;8;1];

%  y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [1.9;3.8;5.7;1]; y(:,7) = [2.2;4.4;6.6;1]; 


% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1]; y(:,5) = [1.7;3.4;5.1;1]; y(:,6) = [-1;-2;-3;1];  y(:,7) = [-5;-10;-15;1]; 

% y(:,1) = [1;2;3;1]; y(:,2) = [1.2;2.4;3.6;1]; y(:,3) = [1.5;3.5;5;1]; y(:,4) = [3;7;8;1]; y(:,5) = [3;5;4;1];  y(:,6) = [3.6;6;4.8;1] ;y(:,7) = [7.2;13;10.4;1];

%  y(:,1) = [1;2;3;1]; y(:,2) = [2;7;1;1]; y(:,3) = [1.5;3.5;5;1]; y(:,4) = [3;7;8;1]; y(:,5) = [3;5;4;1];  y(:,6) = [5;4;1;1] ;y(:,7) = [7.2;13;10.4;1];

% y(:,1) = [1;2;3;1]; y(:,2) = [2;7;1;1]; y(:,3) = [5;4;1;1]; y(:,4) = [3;5;4;1];
% y(:,1) = [0;0;1;1]; y(:,2) = [0;0;2;1]; y(:,3) = [0;0;3;1]; y(:,4) = [0;0;4;1];
% y(:,1) = [1;2;3;1];  y(:,2) = [5;4;1;1]; y(:,3) = [1.3;2.6;3.9;1]; y(:,4) = [2.6;5.2;7.8;1];
% y(:,1) = [0;0;1;1]; y(:,2) = [0;0;2;1]; y(:,3) = [0;0;3;1]; y(:,4) = [0;0;4;1]; y(:,5) = [0;0;-1;1]; y(:,6) = [0;0;-2;1]; y(:,7) = [0;0;-3;1]; 

%  y(:,1) = [0;1;0;1]; y(:,2) = [0;2;0;1]; y(:,3) = [0;-1;0;1]; y(:,4) = [0;3;0;1];
% y(:,5) = [0;4;1;1]; y(:,6) = [0;4;2;1]; y(:,7) = [0;4;3;1]; 
 
% y(:,1) = [1;1;1;1]; y(:,2) = [1;1;1;1]; y(:,3) = [3;2;1;1]; y(:,4)= [2;1;3;1]; y(:,5)=[1;3;2;1] ;
% y(:,1) = [6;3;2;1]; y(:,2) = [1;2;1;1];  y(:,3)= [5;8;1;1]; y(:,4)= [2;1;1;1];
% y(:,1) = [1;0;3;1]; y(:,2) = [1;0;2;1];  y(:,3)= [2;0;5;1]; y(:,4)= [2;0;1;1];
% y(:,1) = [1;0;4;1]; y(:,2) = [2;0;3;1];  y(:,3)= [5;0;1;1]; y(:,4)= [-5;0;6;1];
%  y(:,4)= [2;1;0;1];
%  y(:,4) = [7;-1;2;2];
%   y(:,1) = [1;2;4;1]; y(:,2) = [2;4;8;1];   y(:,3)= [3;6;12;1];
% y(:,1) = [5;4;7;1];  y(:,2) = [1;2;3;1]; 
% y(:,3) = [2;0;3;1];
% y(:,1) = [3;0;0;1] ; y(:,2) = [0;2;0;1]; y(:,3) = [0;0;3;1];


% end of so3 estimation;   
[rownum,lengthy]=size(y); 
for i = 1:iter;    
    se3_t(:,:,i) =  [w_hat_t(:,:,i)+delta_hat_t(:,:,i) v_t(:,1,i)+v_delta(:,1,i);
                  [0,0,0]        0   ];
    
    SE3_t(:,:,i+1) = SE3_t(:,:,i)*expm(se3_t(:,:,i));
    TRt(:,:,i) = [SE3_t(1,1,i) SE3_t(1,2,i) SE3_t(1,3,i);
                 SE3_t(2,1,i) SE3_t(2,2,i) SE3_t(2,3,i);
                 SE3_t(3,1,i) SE3_t(3,2,i) SE3_t(3,3,i)];
             
    TPt(:,i) = [SE3_t(1,4,i) SE3_t(2,4,i) SE3_t(3,4,i)];
    
    TRs(:,:,i) = [SE3_s(1,1,i) SE3_s(1,2,i) SE3_s(1,3,i);
                 SE3_s(2,1,i) SE3_s(2,2,i) SE3_s(2,3,i);
                 SE3_s(3,1,i) SE3_s(3,2,i) SE3_s(3,3,i)];
             
    TPs(:,i) = [SE3_s(1,4,i) SE3_s(2,4,i) SE3_s(3,4,i)];  
%   
    
    error5(:,:,i) =  SE3_t(:,:,i)*inv(SE3_s(:,:,i));
%     dv1(i) = real(visiondistance(TRt(:,:,i),TRs(:,:,i),TPt(:,i),TPs(:,i)));
    dv1(i) = norm(error5(:,:,i)-eye(4),'fro');
    omega1(:,:,i)=zeros(4,4);
    deltay(:,:,i) = zeros(4,14);
    for p = 1:14
       deltay(:,p,i) = [0*cos(i*p/2) 0*sin(i*p/5) 0*sin(i*p/4) 0].';   
       omega1(:,:,i) = omega1(:,:,i)+0.01*(-error5(:,:,i)*y(:,p)-SE3_s(:,:,i)*deltay(:,p,i)+y(:,p))*transpose(y(:,p));
       
    end;
%     A = w_t(:,1,i)+0.04*rand(1);
%     A = [0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
    omega(:,:,i) =   projection(1/4*((omega1(:,:,i) - omega1(:,:,i).')));
    
    
%     v_s(:,i) = v_t(:,1,i)+0.1*rand(1)+ 0.3.*(TPt(:,i)+0.03*rand(1)-TPs(:,i));
    se3_s(:,:,i) =  [w_hat_t(:,:,i) v_t(:,1,i);
                  [0,0,0]        0   ]-1*inv(SE3_s(:,:,i))*omega(:,:,i)*SE3_s(:,:,i); 
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
%     plot3([TPt(1,i),TPt(1,i)+10*TRt(3,1,i)],[TPt(2,i),TPt(2,i)+10*TRt(3,2,i)],[TPt(3,i),TPt(3,i)+10*TRt(3,3,i)],'color',[200 200 200]/255);
%     plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5)
%     hold off  
%     pause(0.01)
    
end 
k1 = 0;
for i=1:iter
k1 = k1+(dv1(i));
end;
display (k1);
figure()
 plot3([TPs(1,iter),TPs(1,iter)+10*TRs(1,1,iter)],[TPs(2,iter),TPs(2,iter)+10*TRs(1,2,iter)],[TPs(3,iter),TPs(3,iter)+10*TRs(1,3,iter)],'r');
    hold on 
    xlim([-20 20]);
    ylim([-20,20]);
    zlim([-20,20]);
    plot3([TPs(1,iter),TPs(1,iter)+10*TRs(2,1,iter)],[TPs(2,iter),TPs(2,iter)+10*TRs(2,2,iter)],[TPs(3,iter),TPs(3,iter)+10*TRs(2,3,iter)],'g');
    plot3([TPs(1,iter),TPs(1,iter)+10*TRs(3,1,iter)],[TPs(2,iter),TPs(2,iter)+10*TRs(3,2,iter)],[TPs(3,iter),TPs(3,iter)+10*TRs(3,3,iter)],'b');
    plot3([TPt(1,iter),TPt(1,iter)+10*TRt(1,1,iter)],[TPt(2,iter),TPt(2,iter)+10*TRt(1,2,iter)],[TPt(3,iter),TPt(3,iter)+10*TRt(1,3,iter)],'color',[10 10 10]/255);
    plot3([TPt(1,iter),TPt(1,iter)+10*TRt(2,1,iter)],[TPt(2,iter),TPt(2,iter)+10*TRt(2,2,iter)],[TPt(3,iter),TPt(3,iter)+10*TRt(2,3,iter)],'color',[100 100 100]/255);
    plot3([TPt(1,iter),TPt(1,iter)+10*TRt(3,1,iter)],[TPt(2,iter),TPt(2,iter)+10*TRt(3,2,iter)],[TPt(3,iter),TPt(3,iter)+10*TRt(3,3,iter)],'color',[200 200 200]/255);
    plot3(y(1,:),y(2,:),y(3,:),'o','Color','k','MarkerSize',5)
    hold off  
    pause(0.05)
 
figure()
 
  plot (dv1);
  title('observer and proposed controlled output');
  grid