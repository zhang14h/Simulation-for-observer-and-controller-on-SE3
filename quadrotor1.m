theta = 0;
fei = 0;
miu = 0;



R11 = cos(theta)*cos(miu);
R12 = sin(fei)*sin(theta)*cos(miu) - cos(fei)*sin(miu);
R13 = cos(miu)*sin(theta)*cos(fei) + sin(miu)*sin(fei);
R21 = cos(theta)*sin(miu);                                                                                    
R22 = sin(miu)*sin(theta)*sin(fei) + cos(fei)*cos(miu);
R23 = cos(fei)*sin(theta)*sin(miu) - sin(fei)*cos(miu);
R31 = -sin(theta);
R32 = cos(theta)*sin(fei);
R33 = cos(theta)*cos(fei);



Jx = 0.033;
Jy = 0.033;
Jz = 0.061;
g = 9.8;
m = 0.25;
kt =3.1*10^(-7); %thrust coefficient
kd = 1.12*10^(-7);
l=0.25;
