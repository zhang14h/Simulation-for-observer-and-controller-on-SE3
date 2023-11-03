close all;
clc;
clear;

landmarks = 20;


for i = 1: landmarks
    y(:,i) = [randi([-8 8]); randi([-8 8])];
        
end;
[miny,pos] = mink(y(2,:),1);
hull(:,1)= y(:,pos);


plot(y(1,:),y(2,:),'o','Color','b','MarkerSize',5); 
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
if ysorted(:,1) == y(:,pos)
  hull(:,2) = ysorted(:,2);
  hull(:,3) = ysorted(:,3);
  i = 4;
else 
   if ysorted(:,2) == y(:,pos)
     hull(:,2) = ysorted(:,1);
     hull(:,3) = ysorted(:,3);
     i = 4;
   else
     hull(:,2) = ysorted(:,1);
     hull(:,3) = ysorted(:,2);
     i = 3;
   end;
end;


  


j = 4;

cout = 1;

 while i <= landmarks
      temp1 = hull(:,j-1);
      temp2 = hull(:,j-2); 
     
      if not(ysorted(1,i) == y(1,pos) &&  ysorted(2,i) == y(2,pos))
         while det([ysorted(1,i) ysorted(2,i) 1; temp1(1,:) temp1(2,:) 1; temp2(1,:) temp2(2,:) 1]) <0
              j = j-1;
              temp1 = hull(:,j-1);
              temp2 = hull(:,j-2);
         end;
         hull(:,j) = ysorted(:,i);
            
      else 
         i = i+1;
         hull(:,j) = ysorted(:,i);
         
      end;
      j = j+1;
      i = i+1;
end; 


figure();

plot(y(1,:),y(2,:),'o','Color','b','MarkerSize',5); 
hold on;
plot(hull(1,:),hull(2,:),'--','Color','r','MarkerSize',3);


