hold off;
hold on;
sensorshape(SE_3s(:,:,1),focall,focalh,hov,vov,2);
sensorshape(SE_3t(:,:,1),focall,focalh,hov,vov,1);
counter = 0;
p = 0;
for j = 1:landmarks
           cframeps(:,j) =  inv(SE_3s(:,:,1)) *  y(:,j);
           cframept(:,j) =  inv(SE_3t(:,:,1)) *  y(:,j);
%            cframept(:,17) = [4;1.5;-1;1];
%             cframept(:,50) = [3.5;-2;1;1];
%             cframept(:,21) = [2;-2;1.5;1];
%             cframept(:,33) = [3;-1.5;-1;1];
%             cframept(:,42) = [2.5;-2;-1;1];
%             cframept(:,19) = [2.8;-2;-1.8;1];
             if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2) & cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                    plot3(y(1,j),y(2,j),y(3,j),'o','Color','g','MarkerSize',5);
                    counter = counter+1;
                    p(counter) = j;
                    dis(j) = sqrt((cframept(1,j)-(focall+focalh)/2)^2+ cframept(2,j)^2+ cframept(3,j)^2+0.01);
             else   if cframept(1,j) >= focall & cframept(1,j) <= focalh & abs(cframept(2,j)) <= cframept(1,j)*tan(hov/2) & abs(cframept(3,j)) <= cframept(1,j)*tan(vov/2);
                            plot3(y(1,j),y(2,j),y(3,j),'o','Color','r','MarkerSize',5); %avoid points missing display
%                     dis(j) = 0;
                     
                    else if cframeps(1,j) >= focall & cframeps(1,j) <= focalh & abs(cframeps(2,j)) <= cframeps(1,j)*tan(hov/2) & abs(cframeps(3,j)) <= cframeps(1,j)*tan(vov/2);
                            plot3(y(1,j),y(2,j),y(3,j),'o','Color','b','MarkerSize',5); %avoid points missing display
                         
                         else 
                             plot3(y(1,j),y(2,j),y(3,j),'o','Color','w','MarkerSize',5);
                             dis(j) = 0;
                    
                         end;
                    end;
            
             end;
%            obs3(:,i) = obs3(:,i) + dis(j);    
 end;
 
plot3(y(1,13),y(2,13),y(3,13),'o','Color','m','MarkerSize',5);
plot3(y(1,46),y(2,46),y(3,46),'o','Color','k','MarkerSize',5);