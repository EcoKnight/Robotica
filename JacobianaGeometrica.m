clc; close; clear;
% Definimos las variables simbólicas
syms q1 q2 q3 q1p q2p q3p d1 d3 

A_0_1 = NMQ_LDGR(0, q1, d1, 0, -90);
A_1_2 = NMQ_LDGR(0, q2, 0, 0, 90);
A_2_3 = NMQ_LDGR(0, 0, d3+q3, 0, 0);
A_0_3 = A_0_1*A_1_2*A_2_3; %Esto es la matriz de Denavit 
DHr= A_0_3((1:3),1:3);
disp('Matriz de rotacion')
disp(DHr)

% Matriz RT
RT = transpose(DHr);
RT=simplify(RT);
disp('R traspuesta')
disp(RT)

% Derivadas respecto a q1 y q2
dRS_dq1 = diff(DHr, q1)*q1p;
dRS_dq2 = diff(DHr, q2)*q2p;
dRS_dq3 = diff(DHr, q3)*q3p;

Rs=dRS_dq1+dRS_dq2+dRS_dq3;
Rs=simplify(Rs);
disp('Suma de derivadas parciales')
disp(Rs)

Omega  = Rs*RT;

Omega=simplify(Omega);
disp('Matriz Omega = R´* RT')
disp(Omega)
Vector(1,1) = Omega(3,2);
Vector(2,1) = Omega(1,3);
Vector(3,1) = Omega(2,1);

disp("Velocidades angulares");
disp(Vector);

function [A] = NMQ_LDGR(T,q,d,a,alpha)

    A = [cosd(T)*cos(q) - sind(T)*sin(q), (-sind(T)*cos(q) - sin(q)*cosd(T))*(cosd(alpha)), (sind(T)*cos(q) + sin(q)*cosd(T))*(sind(alpha)), a*(cosd(T)*cos(q) - sind(T)*sin(q)); 
         sind(T)*cos(q) + sin(q)*cosd(T), (cosd(T)*cos(q) - sind(T)*sin(q))*cosd(alpha), (-cosd(T)*cos(q) + sind(T)*sin(q))*sind(alpha), a*(sind(T)*cos(q) + sin(q)*cosd(T)); 
                           0,                                 sind(alpha),                                cosd(alpha),                                   d; 
                           0,                                     0,                                          0,                                       1];    
end