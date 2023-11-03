close all;
clc;
clear;

landmarks = 20;


for i = 1: landmarks
    y(:,i) = [randi([-8 8]); randi([-8 8])];
        
end;
[miny,pos] = mink(y(2,:),1);



[ysorted(1,:),I] = sort(y(1,:));
ysorted(2,:) = y(2,I);
for i = 1:landmarks-1
    if ysorted(1,i) == ysorted(1,i+1)
       if ysorted(2,i) > ysorted(2,i+1) 
          tep = ysorted(2,i);
          ysorted(2,i) = ysorted(2,i+1);
          ysorted(2,i+1) = tep;
       end;
    end;
end;
hullu(:,1) = ysorted(:,1) ;
hullu(:,2) = ysorted(:,2) ;
% hull(:,3) = ysorted(:,3) ;
if det([ysorted(1,3) ysorted(2,3) 1; hullu(1,2) hullu(2,2) 1; hullu(1,1) hullu(2,1) 1]) < 0 
    hullu (:,2) = ysorted(:,3);
    i = 4;
else  
    i = 3;
end;

j = 3;



 while i <= landmarks
      temp1 = hullu(:,j-1);
      temp2 = hullu(:,j-2); 
     
      
      while det([ysorted(1,i) ysorted(2,i) 1; temp1(1,:) temp1(2,:) 1; temp2(1,:) temp2(2,:) 1]) <=0
              
              j = j-1;
              hullu(:,j) = [];
              temp1 = hullu(:,j-1);
              if j - 2 == 0
                 temp2 = [hullu(1,1) ; hullu(2,1)-0.0001];
              else   
                 temp2 = hullu(:,j-2);
              end;   
                 
      end;
      hullu(:,j) = ysorted(:,i);     
      j = j+1;
      i = i+1;
end; 
for i = 1:landmarks-1
    if ysorted(1,i) == ysorted(1,i+1)
       if ysorted(2,i) < ysorted(2,i+1) 
          tep = ysorted(2,i);
          ysorted(2,i) = ysorted(2,i+1);
          ysorted(2,i+1) = tep;
       end;
    end;
end;
hulll(:,1) = ysorted(:,1) ;
hulll(:,2) = ysorted(:,2) ;
if det([ysorted(1,3) ysorted(2,3) 1; hulll(1,2) hulll(2,2) 1; hulll(1,1) hulll(2,1) 1]) > 0 
    hulll (:,2) = ysorted(:,3);
    i = 4;
else  
    i = 3;
end;

j = 3;
while i <= landmarks
      temp1 = hulll(:,j-1);
      temp2 = hulll(:,j-2); 
     
      
      while det([ysorted(1,i) ysorted(2,i) 1; temp1(1,:) temp1(2,:) 1; temp2(1,:) temp2(2,:) 1]) >=0
              
              j = j-1;
              hulll(:,j) = [];
              temp1 = hulll(:,j-1);
              if j - 2 == 0
                 temp2 = [hulll(1,1) ; hulll(2,1)+0.0001];
              else   
                 temp2 = hulll(:,j-2);
              end;   
                 
      end;
      hulll(:,j) = ysorted(:,i);     
      j = j+1;
      i = i+1;
end; 
figure();

plot(y(1,:),y(2,:),'o','Color','b','MarkerSize',5); 
hold on;
plot(hullu(1,:),hullu(2,:),'--','Color','r','MarkerSize',3);
hold on;
plot(hulll(1,:),hulll(2,:),'--','Color','r','MarkerSize',3);

[k,av] = convhull(y.');
figure;
plot(y(1,:),y(2,:),'o','Color','b','MarkerSize',5);
hold on
plot(y(1,k),y(2,k))
