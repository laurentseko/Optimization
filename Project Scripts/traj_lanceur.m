function [vit_fin,cont_value,vit_reelle,Mat_ret,temps_ret] = traj_lanceur(theta,masse_erg,masse_fus,H_c)
a = [15; 10; 10];
k = [0.1101 ; 0.1532 ; 0.2154];
V_i = [2647.2 ; 2922.4 ; 4344.3];
v0 = 100; R_t = 6378137;
% theta = theta*pi/180;
r_0 = [R_t ; 0]; v_0 = [v0*cosd(theta(1)) ; v0*sind(theta(1))];
y_0 = [r_0 ; v_0 ; masse_fus];

%fprintf("traject\n");
Mat_ret = [];
temps_ret = [];
M_i = masse_fus;
t_i = 0;
for i = 1:3
    T_i = a(i)*M_i; q_i = T_i/V_i(i); t_ci = masse_erg(i)/q_i;
    pas_t = linspace(t_i,t_ci + t_i , 100);
    n = length(pas_t);
    [t,Y] = ode45(@(t,Y) Phi(t,Y,theta(i+1),V_i(i),T_i),pas_t, y_0);
    
    Mat_ret = [Mat_ret; Y];
    temps_ret = [temps_ret ; t];
    
    y_0 = [Y(n,1); Y(n,2); Y(n,3); Y(n,4); Y(n,5)-k(i)*masse_erg(i)];
    M_i = y_0(5); 
    t_i = t_ci + t_i;
end
vit_reelle = norm(y_0(3:4));
vit_fin = 100/vit_reelle;
cont_value = [norm(y_0(1:2)) - (R_t + H_c) , y_0(1:2)'*y_0(3:4)];
end

