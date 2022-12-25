function dQ = nonLinear_model_Observer3(~,y,F,L_Matrix3)
m1 = 100; m2 = 100; M = 1000; L1 = 20; L2 = 10; g = 9.81;
X = y(1);
dX = y(2);
t1 = y(3);
dt1 = y(4);
t2 = y(5);
dt2 = y(6);
dQ=zeros(6,1);
Y3 = [X; 0; t2];
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
sum = L_Matrix3*(Y3-C3*y);
dQ(1) = dX + sum(1);
dQ(2) = (F-((m1*sin(t1)*cos(t1))+(m2*sin(t2)*cos(t2)))*g - (L1*m1*(dQ(3)^2)*sin(t1)) - (L2*m2*(dQ(5)^2)*sin(t2)))/(m1+m2+M-(m1*(cos(t1)^2))-(m2*(cos(t2)^2)))+sum(2);
dQ(3) = dt1+sum(3);
dQ(4) = ((cos(t1)*dQ(2)-g*sin(t1))/L1) + sum(4);
dQ(5) = dt2 + sum(5);
dQ(6) = (cos(t2)*dQ(2)-g*sin(t2))/L2 + sum(6);
end