clear;
counter = 0;
% for i = 0:7
%     for j = 0:7
%         for k = 0:7
%             for l = 0:7
%                 for m = 0:7
%                     for n = 0:7
%                         counter = counter +1; 
%                         y = [i l;j m;k n;1 1]*[i l;j m;k n;1 1].';
%                         ei(:,1) = eig(y); 
%                         if (ei(1,1) == ei(2,1)  && rank(y) == 2 )|| (ei(1,1) == ei(3,1) && rank(y) == 2 )|| (ei(1,1) == ei(4,1)&& rank(y) == 2) || (ei(2,1) == ei(3,1)&& rank(y) == 2) || (ei(2,1) == ei(4,1)&& rank(y) == 2) ||(ei(3,1) == ei(4,1)&& rank(y) == 2)
%                           
%                            display ([i l;j m;k n;1 1]);
%                            display( y(:,:,1));
%                            display( ei(:,1));
%                         end;   
%                     end;
%                 end;
%             end;
%         end;
%     end;
% end;    
    

% for i = 0:7
%     for j = 0:7
%         for k = 0:7
%             for l = 0:7
%                 for m = 0:7
%                     for n = 0:7
%                         for o = 0:7
%                             for p = 0:7
%                                 for q = 0:7
%                                    counter = counter +1; 
%                                    y = [i l o;j m p;k n 1;1 1 1]*[i l o;j m p;k n q;1 1 1].';
%                                    ei(:,1) = eig(y); 
%                                    if (ei(1,1) == ei(2,1)  && rank(y) >= 3 )|| (ei(1,1) == ei(3,1) && rank(y) >= 3 )|| (ei(1,1) == ei(4,1)&& rank(y) >= 3) || (ei(2,1) == ei(3,1)&& rank(y) >= 3) || (ei(2,1) == ei(4,1)&& rank(y) >= 3) ||(ei(3,1) == ei(4,1)&& rank(y) >= 2)
%                                         display ([i l o;j m p;k n q;1 1 1]);
%                                         display( y(:,:,1));
%                                         display( ei(:,1));
%                                    end;
%                                 end;   
%                             end;
%                         end;           
%                     end;
%                 end;
%             end;
%         end;
%     end;
% end;    

% for i = 1:4
%     for j = 1:4
%         for k = 1:4
%             for l = 1:4
%                 for m =1:4
%                     for n = 1:4
%                         for o = 1:4
%                             for p = 1:4
%                                 for q = 1:4
%                                     for r =1:4
%                                         for s = 1:4
%                                             for t = 1:4
%                                                  counter = counter +1; 
%                                                  a = [i l o r;j m p s;k n 1 t;1 1 1 1];
%                                                  y = a*a.';
%                                                  ei(:,1) = eig(y); 
%                                                  if (ei(1,1) == ei(2,1)  && rank(y) >= 3 )|| (ei(1,1) == ei(3,1) && rank(y) >= 3 )|| (ei(1,1) == ei(4,1)&& rank(y) >= 3) || (ei(2,1) == ei(3,1)&& rank(y) >= 3) || (ei(2,1) == ei(4,1)&& rank(y) >= 3) ||(ei(3,1) == ei(4,1)&& rank(y) >= 3)
%                                                         display(a);
%                                                         display(y);
%                                                         display(ei(:,1)); 
%                                                         display(rank(y));
%                                                  end;
%                                             end;
%                                         end;
%                                     end;    
%                                 end;   
%                             end;
%                         end;           
%                     end;
%                 end;
%             end;
%         end;
%     end;
% end; 

for i = 1:4
    for j = 1:4
        for k = 1:4
            for l = 1:4
                for m =1:4 
                    for n = 1:4
                        for o = 1:4
                            for p = 1:4
                                for q = 1:4
                                    for r =1:4
                                        for s = 1:4
                                            for t = 1:4
                                                for u = 1:4
                                                    for v = 1:4
                                                        for w = 1:4
                                                               counter = counter +1; 
                                                               a = [i l o r u;j m p s v;k n 1 t w;1 1 1 1 1];
                                                               y = a*a.';
                                                               ei(:,1) = eig(y); 
                                                               if (ei(1,1) == ei(2,1)  && rank(y) >= 4 )|| (ei(1,1) == ei(3,1) && rank(y) >= 4 )|| (ei(1,1) == ei(4,1)&& rank(y) >= 4) || (ei(2,1) == ei(3,1)&& rank(y) >= 4) || (ei(2,1) == ei(4,1)&& rank(y) >= 4) ||(ei(3,1) == ei(4,1)&& rank(y) >= 4)
                                                                   display(a);
                                                                   display(y);
                                                                   display(ei(:,1)); 
                                                                   display(rank(y));
                                                               end;
                                                        end;
                                                    end;
                                                end;    
                                            end;
                                        end;
                                    end;    
                                end;   
                            end;
                        end;           
                    end;
                end;
            end;
        end;
    end;
end; 
    