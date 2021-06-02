function [grad_f,grad_C] = Gradient(Problem,x,h)
n = length(x);
if(nargin<3)
h = 1e-3;
end
for i = 1:n
    H = zeros(n,1);
    H(i) = h;
    [fplus_h,Contplus_h] = Problem(x + H);
    [f,Cont] = Problem(x);
    grad_f(i) = (fplus_h - f)/h;
    grad_C(i,:) = (Contplus_h - Cont)*(1/h) ; 
end     
end