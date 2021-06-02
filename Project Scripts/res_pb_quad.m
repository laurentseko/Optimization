function [dep_x,dep_lamb] = res_pb_quad(grad_f,grad_cont,Hes,cont_value)
% A = (Hes\grad_cont);
% B = (Hes\grad_f');
% C = grad_cont'*(Hes\grad_cont);
dep_lamb = -linsolve(grad_cont'*(Hes\grad_cont), grad_cont'*(Hes\grad_f') - cont_value');

b = grad_cont*dep_lamb + grad_f';
dep_x = -linsolve(Hes,b);

end

