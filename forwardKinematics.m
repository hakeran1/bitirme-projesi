function [T_0_EE,T] = forwardKinematics(th)
%forwardKinematics  Forward kinematics for specified RRRRRRR
%                   parameters and MDH parameters
%
%   [T_0_6, T] = forwardKinematics(th)
%   
%   T outputs a 3 dimensional array where index i is T_0_i

L1 = 0.1;
L2 = 0.3;
L3 = 0.2;
L4 = 0.06;

%    |-| alpha_i-1 |-| a_i-1 |-| d_i |-| th_i |-|
MDH = [0 0 L1 th(1); -pi/2 0 0 -pi/2 + th(2); 0 L2 0 th(3);
    -pi/2 0 L3 th(4); pi/2 0 0 th(5); -pi/2 0 0 th(6); 0 0 L4 0];

T_0_EE = eye(4);
T = zeros(4,4,7);

for i = 1:7
    T_0_EE = T_0_EE * [cos(MDH(i,4)) -sin(MDH(i,4)) 0 MDH(i,2);
        sin(MDH(i,4))*cos(MDH(i,1)) cos(MDH(i,4))*cos(MDH(i,1)) -sin(MDH(i,1)) -sin(MDH(i,1))*MDH(i,3);
        sin(MDH(i,4))*sin(MDH(i,1)) cos(MDH(i,4))*sin(MDH(i,1)) cos(MDH(i,1)) cos(MDH(i,1))*MDH(i,3);
        0 0 0 1];

    T(:,:,i) = T_0_EE;
end
T_0_EE(abs(T_0_EE)<10e-12) = 0;