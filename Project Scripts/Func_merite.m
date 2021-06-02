function [val_F] = Func_merite(Problem,x,rho)
%fprintf("fonc_merite\n");
[fvalue,Contvalue] = Problem(x);
val_F = fvalue + rho*norm(Contvalue,1);
end

