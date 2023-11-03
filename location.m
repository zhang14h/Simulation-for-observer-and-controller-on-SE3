close all;
clc;
clear;

landmarks = 20;

X = [randi([-8 8]); randi([-8 8])];
for i = 1: landmarks
    y(:,i) = [randi([-8 8]); randi([-8 8])];
        
end;


[k,av] = convhull(y.');
hull = [y(1,k);y(2,k)];
figure;
plot(X(1),X(2),'o','Color','k','MarkerSize',5);
hold on
plot(hull(1,:),hull(2,:))

if det([hull(1,1) hull(2,1) 1;
       hull(1,2) hull(2,2) 1;
       X(1) X(2) 1]) < 0 | det([hull(1,1) hull(2,1) 1;
       hull(1,length(k)-1) hull(2,length(k)-1) 1;
       X(1) X(2) 1]) > 0
    disp ('X is outside')
end

min = 2;
max = length(k)-1;

while max - min > 1
      mid = ceil((max+min)/2);
      if det([hull(1,1) hull(2,1) 1;
         hull(1,mid) hull(2,mid) 1;
         X(1) X(2) 1]) < 0
         max = mid;
      else
         min = mid; 
      end;   
end;

if det([hull(1,min) hull(2,min) 1;
       hull(1,max) hull(2,max) 1;
       X(1) X(2) 1]) > 0 
    disp ('X is inside')
else
    disp ('X is outside')
    
end