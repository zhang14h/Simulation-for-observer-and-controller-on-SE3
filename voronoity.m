clear;
clc;
close all;

%width, height and number of points
W = 100;
L = 200;
p =10;
fdistance = 100;
%% method 1
%generating points
tic;
for i = 1: p
    l(:,i) = [randi([0 L]); randi([0 W])];       
end;

% delaunay triangulation
tri = delaunay(l(1,:).',l(2,:).');
tr = triangulation(tri,l(1,:).',l(2,:).');

%center of the delaunay triangles
c = tr.circumcenter();

%determine the running time for special case(point less than 3)
if length(tri) <= 3
   iter = rank(tri);
else 
   iter = length(tri);
end;  

%determine the triagnles that share edges
for i = 1:iter
    for j = 1:iter
          v = double(ismember(tri(i,:),tri(j,:)));
          if sum(v) == 2
             e(j,i) = 1;
          else  
             e(j,i) = 0; 
          end;   
    end; 
end;

%determine the inner edges of the voronoi partition by connecting the
%neighbour center of the circles
counter = 1;
for i = 1:iter
    for j = i:iter  
          if e(j,i) == 1
             edge(counter,:) = [c(i,1) c(i,2) c(j,1) c(j,2)];
             counter = counter + 1;
          end;
    end; 
end;


%determine the outer edges
%the convex hull
h = convhull(l(1,:),l(2,:));

%find the edge of the convex hull and its corresponding center of circle
for i = 1:length(h)-1 
    for j = 1:iter
        v = ismember([h(i) h(i+1)],tri(j,:));
        v = double(v);
        if sum(v) == 2
           oedge(i,:) = [h(i) h(i+1) j]; 
        end;    
    end;    
end;  

%establish coneection between the midpoint of the edge and the center of the circle 
outedge = [(l(1,oedge(:,1))+l(1,oedge(:,2))).'/2 (l(2,oedge(:,1))+l(2,oedge(:,2))).'/2 c(oedge(:,3),:)];
[rp,M] = max([l(1,oedge(:,1)).' l(1,oedge(:,2)).'],[],2);
if length(oedge) <= 3
   iter1 = rank(oedge);
else 
   iter1 = length(oedge);
end;  

%find the right point of the edge
for i = 1:iter1
    sidep(i,:) = [l(1,oedge(i,M(i))) l(2,oedge(i,M(i)))];
end;    
%determine whether it is right turn or left turn and draw the line
for i = 1:iter1
    r(i) = det([c(oedge(i,3),1)  c(oedge(i,3),2) 1
               l(1,oedge(i,2)) l(2,oedge(i,2)) 1;
               l(1,oedge(i,1))  l(2,oedge(i,1)) 1]);
    if r(i) >= 0
       r(i) = 1; %1 is outside
       outedge(i,:) = [outedge(i,3) outedge(i,4) outedge(i,3)+fdistance*(outedge(i,3)-outedge(i,1)) outedge(i,4)+fdistance*(outedge(i,4)-outedge(i,2))];
    else
       r(i) = 0; %0 is inside
       outedge(i,:) = [outedge(i,3) outedge(i,4) outedge(i,1)+fdistance*(outedge(i,1)-outedge(i,3)) outedge(i,2)+fdistance*(outedge(i,2)-outedge(i,4))];
    end; 
    hold on;
end;   
toc;
% plot

plot(l(1,:),l(2,:),'o','Color','b','MarkerSize',2);
hold on
plot(c(:,1),c(:,2),'o','Color','r','MarkerSize',2)
hold on
if p==3
     plot([outedge(:,1) outedge(:,3)].',[outedge(:,2) outedge(:,4)].','color','k');
else     
     plot([edge(:,1) edge(:,3)].',[edge(:,2) edge(:,4)].','color','k');
     plot([outedge(:,1) outedge(:,3)].',[outedge(:,2) outedge(:,4)].','color','k');
end;     
hold on
plot([outedge(:,1) outedge(:,3)].',[outedge(:,2) outedge(:,4)].','color','k');
xlim([0 L]);
ylim([0 W]);

%% method 2
% matrix creation
tic;
dis = ones(W,L,p);
rows = ones(W,L).*(1:L);
cols = ones(W,L).*(1:W).';
vorons = zeros(W,L);
max = sqrt(L^2+W^2)*ones(W,L); 
%distance calculation
for i = 1:p
   dis(:,:,i) = sqrt((rows-l(1,i)*ones(W,L)).^2+(cols-l(2,i)*ones(W,L)).^2);
end;   
%voronoi partition
for i = 1:p
    temp = dis(:,:,i);
    csind = find(temp<=max);
    max(csind) = temp(csind);
    vorons(csind) = i;
end;  
toc;
figure
mesh(1:L,1:W,vorons);