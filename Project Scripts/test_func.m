function [f] = test_func(x,k,mu)
I = ones(3,1);
f = dot((I+k) , x) + mu;
end

