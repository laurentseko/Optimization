function [H] = Hessienne(H_old,y,d,choix)
switch choix
    case{"BFGS","bfgs"}
%if choix == "BFGS"
    if dot(y,d) > 0
        H = H_old + (y*y')*(1/dot(y,d)) - (H_old*d)*(d'*H_old)*(1/dot(d'*H_old,d));
    else
        H = H_old    ;
    end 
    case{"SRI","sri"}
%else
    c_k = y - H_old*d;
    if dot(d,c_k) ~= 0
        H = H_old + (c_k*c_k')*(1/dot(d,c_k));
    else
        H = H_old;
    end
end
if sum((eig(H)>0)) < length(eig(H))
    n = length(eig(H));
    H = H + (0.1 - min(eig(H)))*eye(n,n); 
    H = H + n*eye(n,n);
end
end      


