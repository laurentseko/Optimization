function [grad_f] = Gradient_newton(fun,x,h)
    [fplus_h] = fun(x + h);
    [f] = fun(x);
    grad_f = (fplus_h - f)/h;
   
end     

