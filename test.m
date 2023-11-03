Q =  0.01 * [1 2 3 1; 4 3 1 1; 2 1 2 1; 1 2 1 3];
Q=  Q+1*eye(4);
Q = Q*Q;

s = 1/2*projection(SE_3t(:,:,12)-SE_3t(:,:,12).')
m = 1/2*projection(100*Q*Q.'*SE_3t(:,:,12)-100*SE_3t(:,:,12).'*Q*Q.')
t = 1/2*projection(10*Q*SE_3t(:,:,12)-10*SE_3t(:,:,12)*Q.')

trace(m.'*SE_3t(:,:,12))
trace(m.'*s)
trace(t.'*t)