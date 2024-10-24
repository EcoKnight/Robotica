clc; close; clear;
syms q1 q2 q3 q4 q5 d1 d3 d5 L2 q1p q2p q3p

A_0_1 = NMQ_LDGR(0, q1, d1, 0, -90);
A_1_2 = NMQ_LDGR(0, q2, 0, 0, 90);
A_2_3 = NMQ_LDGR(0, 0, d3+q3, 0, 0);

A_0_3 = A_0_1*A_1_2*A_2_3; %Esto es la matriz de Denavit 

X = A_0_3(1,4);
Y = A_0_3(2,4);
Z = A_0_3(3,4);

Theta = q2; 
Phi = q1;
Psi = 0;

Xq1 = diff(X, q1);
Xq2 = diff(X, q2);
Xq3 = diff(X, q3);

Yq1 = diff(Y, q1);
Yq2 = diff(Y, q2);
Yq3 = diff(Y, q3);

Zq1 = diff(Z, q1);
Zq2 = diff(Z, q2);
Zq3 = diff(Z, q3);

Thetaq1 = diff(Theta, q1);
Thetaq2 = diff(Theta, q2);
Thetaq3 = diff(Theta, q3);

Phiq1 = diff(Phi, q1);
Phiq2 = diff(Phi, q2);
Phiq3 = diff(Phi, q3);

Psiq1 = diff(Psi, q1);
Psiq2 = diff(Psi, q2);
Psiq3 = diff(Psi, q3);

JA =([ Xq1 Xq2 Xq3;
          Yq1 Yq2 Yq3;
          Zq1 Zq2 Zq3;
          Thetaq1 Thetaq2 Thetaq3;
          Phiq1 Phiq2 Phiq3;
          Psiq1 Psiq2 Psiq3]);

VJA =JA* [q1p; q2p; q3p];

Xp =VJA(1,1);
Yp =VJA(2,1);
Zp = VJA(3,1);
Thetap = VJA(4,1);
Phip = VJA(5,1);
Psip= VJA(6,1);

M(1,1) = Xp;
M(2,1) = Yp;
M(3,1) = Zp;
M(4,1) = Thetap;
M(5,1) = Phip;
M(6,1) = Psip;
disp(M)

function [A] = NMQ_LDGR(T,q,d,a,alpha)

    A = [cosd(T)*cos(q) - sind(T)*sin(q), (-sind(T)*cos(q) - sin(q)*cosd(T))*(cosd(alpha)), (sind(T)*cos(q) + sin(q)*cosd(T))*(sind(alpha)), a*(cosd(T)*cos(q) - sind(T)*sin(q)); 
         sind(T)*cos(q) + sin(q)*cosd(T), (cosd(T)*cos(q) - sind(T)*sin(q))*cosd(alpha), (-cosd(T)*cos(q) + sind(T)*sin(q))*sind(alpha), a*(sind(T)*cos(q) + sin(q)*cosd(T)); 
                           0,                                 sind(alpha),                                cosd(alpha),                                   d; 
                           0,                                     0,                                          0,                                       1];    
end