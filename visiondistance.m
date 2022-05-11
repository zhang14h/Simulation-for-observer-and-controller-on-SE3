function dV = visiondistance(a,b,c,d)
         
         dF = sqrt(6-trace(2*transpose(a)*b));
         oF = 1/(1-dF/sqrt(8));
         eF = (sqrt(transpose(c-d)*(c-d)));
         if eF == 0;
            eF = 0.01;
         end;
         if oF ==0;
            of = 0.01;
         end;
         dV = oF * eF;
end
