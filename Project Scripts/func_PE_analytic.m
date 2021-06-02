function [f,x_1,x_2] = func_PE_analytic(x,vit_ejec,Omega,vit_p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x_2=(vit_ejec(2)-vit_ejec(3)*(1-Omega(3)*x))/(Omega(2)*vit_ejec(2));
x_1=(vit_ejec(1)-vit_ejec(2)*(1-Omega(2)*x_2))/(Omega(1)*vit_ejec(1));
f=vit_ejec(1)*log(x_1)+vit_ejec(2)*log(x_2)+ vit_ejec(3)*log(x)-vit_p

end

