function [fei,theta,miu] = derotation(A)
      
      miu = atan (A(2,1)/A(1,1));
      theta = atan(-A(3,2)/sqrt(A(3,2)^2+A(3,3)^2));
      fei = atan(A(3,2)/A(3,3));
 end