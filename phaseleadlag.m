num = 1600*[1 2*0.19 0.19^2+11.15^2];
den = [0.0101*0.0103 0 0.0101*1.37+0.0103*1.37 0 0];
g = tf(num,den);
bode (g)