function [c] = test_cont(x,V,k,mu,vit_p)
% v = x.*k;
% d = cumsum(v(3:-1:1));
% d = d(3:-1:1);
I = ones(3,1);
% c = dot(V,log(I + x./(mu*I + d))) - 11527;

mf1 = dot([0;1;1]+k , x) + mu;
mf2 = dot([0;1]+k(2:3) , x(2:3)) + mu;
mf3 = k(3)*x(3) + mu;
Mf = [mf1;mf2;mf3];
c = dot(V , log(I + x./Mf)) - vit_p; 
end

