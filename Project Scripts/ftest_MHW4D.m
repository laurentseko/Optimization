function [f] = ftest_MHW4D(X)
f = power(X(1)-1,2) + power(X(1)-X(2),2) + power(X(2)-X(3),3) + power(X(3)-X(4),4) + power(X(4)-X(5),4);
end