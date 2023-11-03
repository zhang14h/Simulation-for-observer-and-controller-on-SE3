close all;
clc;
clear;

landmarks = 20;

figure()
for i = 1: landmarks
    y(:,i) = [randi([-10 10]); randi([-10 10])];
        
end;

[miny,pos] = mink(y(2,:),1);
hull(:,1)= y(:,pos);
i = 1;
plot(y(1,:),y(2,:),'o','Color','b','MarkerSize',5); 
while 1
    if (hull(1,i) == y(1,pos) && hull(2,i) == y(2,pos) && i >2 ); 
        break;
    end;
    for j = 1:landmarks
        for k = 1:landmarks
          p(j,k) = det([hull(1,i) hull(2,i) 1; y(1,j) y(2,j) 1; y(1,k) y(2,k) 1]);
        end;
    end;
    ct = 1;
    for j = 1:landmarks
          
          if min(p(j,:)) >= 0 && max(p(j,:)) > 0 && not(hull(1,i) == y(1,j) && hull(2,i) == y(2,j));
             hull(:,i+1) = y(:,j);
          end    
        
    end;
    
    i = i+1;
    
 end

hold on;
plot(hull(1,:),hull(2,:),'--','Color','r','MarkerSize',3);

