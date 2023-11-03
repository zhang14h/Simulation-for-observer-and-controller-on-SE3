theta = 1.2;
fei = 0.7;
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

SO3= [R11 R12 R13;
      R21 R22 R23;
      R31 R32 R33];

