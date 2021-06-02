function [x,f,k,recap_x] = Newton(func_PE_analytic, x_0, eps)
k=0;
x=x_0;
[f]=func_PE_analytic(x);
a=10;
recap_x=[x];
while abs(f)>eps && k<100
 x_plus1=x-(f/a);
 if x_plus1<2.3
     x_plus1=2.3
 else
 end
 a_plus1=(func_PE_analytic(x)-func_PE_analytic(x_plus1))/(x-x_plus1);
 [f]=func_PE_analytic(x_plus1);
 x=x_plus1;
 a=a_plus1;
 recap_x = [recap_x,x];
 k=k+1
end
end

