clc; close; clear;
syms q1 q2 q3 q4 q5 d1 d3 d5 L2 q1p q2p q3p

A_0_0 = eye(4,4);
A_0_1 = NMQ_LDGR(0, q1, d1, 0, -90);
A_1_2 = NMQ_LDGR(0, q2, 0, 0, 90);
A_2_3 = NMQ_LDGR(0, 0, d3+q3, 0, 0);

A_0_2 = A_0_1*A_1_2;
A_0_3 = A_0_1*A_1_2*A_2_3; %Esto es la matriz de Denavit 

Z_0_0 = transpose(A_0_0(1:3,3));
Z_0_1 = transpose(A_0_1(1:3,3));
Z_0_2 = transpose(A_0_2(1:3,3));

P_0_3 = (A_0_3(1:3,4)) - (A_0_0(1:3,4));  %%P_0_3 = transpose(A_0_3(1:3,4)) -transpose(A_0_0(1:3,4));
P_1_3 = (A_0_3(1:3,4)) - (A_0_1(1:3,4));


JG1(1:3,1) = transpose(cross(Z_0_0, P_0_3));
JG1(4:6,1) = transpose(Z_0_0);

JG2(1:3,1) = transpose(simplify(cross(Z_0_1, P_1_3)));
JG2(4:6,1) = transpose(Z_0_1);

JG3(1:3,1) = transpose(Z_0_2);
JG3(4:6,1) = zeros(3,1);

JG = [JG1,JG2,JG3];
Theta = [q1p;q2p;q3p];
Pv =  JG*Theta;

Vx = Pv(1,1);
Wx =Pv(4,1);

Vy = Pv(2,1);
Wy =Pv(5,1);

Vz = Pv(3,1);
Wz =Pv(6,1);
disp('Perfil de velocidades')
disp(Vx)
disp(Vy)
disp(Vz)
disp(Wx)
disp(Wy)
disp(Wz)
function [A] = NMQ_LDGR(T,q,d,a,alpha)

    A = [cosd(T)*cos(q) - sind(T)*sin(q), (-sind(T)*cos(q) - sin(q)*cosd(T))*(cosd(alpha)), (sind(T)*cos(q) + sin(q)*cosd(T))*(sind(alpha)), a*(cosd(T)*cos(q) - sind(T)*sin(q)); 
         sind(T)*cos(q) + sin(q)*cosd(T), (cosd(T)*cos(q) - sind(T)*sin(q))*cosd(alpha), (-cosd(T)*cos(q) + sind(T)*sin(q))*sind(alpha), a*(sind(T)*cos(q) + sin(q)*cosd(T)); 
                           0,                                 sind(alpha),                                cosd(alpha),                                   d; 
                           0,                                     0,                                          0,                                       1];    
end