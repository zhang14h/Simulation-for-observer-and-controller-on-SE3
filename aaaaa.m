clc;
clear;
i = 1;
for theta = 0:pi/4:pi
    for fei = 0:pi/4:pi
        for miu = 0:pi/4:pi
            
          S(:,:,i) = rotation (theta,fei,miu);
          deter(i) = det(S(:,:,i));
          i = i+1;
            
        end;
    end;
end;