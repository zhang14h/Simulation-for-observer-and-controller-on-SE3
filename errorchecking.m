for i = 1:iter
    e1 = 1/2*projection((eye(4)+P)*error6(:,:,i)-error6(:,:,i).'*(eye(4)+P).');
    ne1(i) = norm(e1,'fro');
    e2 = 1/2*projection((eye(4)+P)*((error7(:,:,i))-eye(4))*error6(:,:,i)-error6(:,:,i).'*((error7(:,:,i))-eye(4)).'*(eye(4)+P).');
    ne2(i) = norm(e2,'fro');
    et1(i)=ne1(i)-ne2(i);
end;    

figure()
plot(et1);