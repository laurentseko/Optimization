function [G_L] = grd_lag(grad_f,grad_C,lamb)
G_L = grad_f' + grad_C*lamb;
end

