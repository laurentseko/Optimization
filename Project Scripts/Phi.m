function [Y_point] = Phi(t,Y,theta,V_i,T_i)
C = 0.1; 
mu = 3.986*10^14; 
R_t = 6378137;
H = 7000;

R = Y(1:2);
V = Y(3:4);
M = Y(5);

r_norm = norm(R);
v_norm = norm(V);

rhi = 1.225*exp(-(r_norm-R_t)/H);
e_h = (1/r_norm)*[-R(2);R(1)];
e_r = (1/r_norm)*R;
delta = asind(R'*V/(r_norm*v_norm));
if delta<-90
    pause
end
if delta>90
    pause
end

U = e_h*cosd(delta + theta) + e_r*sind(delta + theta);

T = T_i*U;
W = -mu*(1/r_norm^3)*M*R;
D = -C*rhi*v_norm*V;

V_prime = (T+W+D)/M;
Y_point = [V;V_prime;-T_i/V_i];
end