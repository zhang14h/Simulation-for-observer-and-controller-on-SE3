clc;
clear
counter = 1;

for i = 0:pi/10:2*pi-0.01
    for j = 0:pi/10:2*pi-0.01
        for k = 0:pi/10:2*pi-0.01
            dis(:,:,counter) = rotation(i,j,k);
            
            if dis(1,1,counter) == 1 && dis(2,2,counter) == 1 && dis(3,3,counter) == 1
               display (counter);
               display (i),display (j),display (k);
%                display (j);
%                display (k);
            end;   
            if dis(1,1,counter) == 1 && dis(2,2,counter) == -1 && dis(3,3,counter) == -1
               display (counter);
               display (i),display (j),display (k);
%                display (j);
%                display (k);
            end;   
            if dis(1,1,counter) == -1 && dis(2,2,counter) == 1 && dis(3,3,counter) == -1
               display (counter);
               display (i),display (j),display (k);
%                display (j);
%                display (k);
            end;   
            if dis(1,1,counter) == -1 && dis(2,2,counter) == -1 && dis(3,3,counter) == 1
               display (counter);
               display (i),display (j),display (k);
%                display (j);
%                display (k);
            end;  
         counter = counter+1; 
        end;
    end;
end    