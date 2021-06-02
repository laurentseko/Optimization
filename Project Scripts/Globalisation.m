function[new_x,dep_x_valide,compteur]=Globalisation(F,grad_f,cont_value,dep_x,rho,x_k,born_min,born_max)
compteur=0;
c = 0.1; s = 1; max_iter = 100;
x_plus1 = x_k + s*dep_x;
for i = 1:length(x_plus1)
    if x_plus1(i) < born_min(i)
        x_plus1(i) =  born_min(i) ;
    elseif x_plus1(i) > born_max(i)
        x_plus1(i) = born_max(i);
    end
end
F_1 = F(x_plus1);
F_2 = F(x_k);
compteur=2;
if F_1 < F_2
        %fprintf("acceptable\n");
        new_x = x_plus1;
        dep_x_valide = s*dep_x;
else
    F_prime = dot(grad_f,dep_x) - rho*norm(cont_value,1) ;
    if F_prime < 0
        fprintf("recherche lineaire\n");
        F_1 = F(x_plus1);
        F_2 = F(x_k);
        iter = 0;
        compteur=4;
        while (F_1 >= (F_2 + c*s*F_prime)) && (iter < max_iter)
            s = s/2;
            x_plus1 = x_k + s*dep_x;
            for i = 1:length(x_plus1)
                if x_plus1(i) < born_min(i)
                    x_plus1(i) =  born_min(i) ;
                elseif x_plus1(i) > born_max(i)
                    x_plus1(i) = born_max(i);
                end
            end
            F_1 = F(x_plus1);
            iter=iter+1;
            compteur=compteur+1;
        end
        if iter == max_iter
            %fprintf("non acceptable\n");
            new_x = x_k;
            dep_x_valide = dep_x;   
        else
            new_x = x_plus1;
            dep_x_valide = s*dep_x;
            %fprintf("sortie\n");
        end
     else
       %fprintf("non acceptable\n");
        new_x = x_k;
        dep_x_valide = dep_x;
    end
end
compteur;
end