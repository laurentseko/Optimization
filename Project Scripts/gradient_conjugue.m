function [X,fX,gX,kX,hist]=gradient_conjugue(A,b,c,X0,eps,Kmax)
      
    %f(x) = 0.5<Ax,x> + <b,x> + c
    
    %argument d'appel
    %X0: vecteur de taille n
    %A: matrice de taille n*n
    %b: vecteur de taille n
    %c: nombre reel
    %eps: nombre reel, tolerance
    %Kmax: nombre de pas maximal
   
    %argument de sortie
    %X: valeur de x trouvee qui minimise f
    %fX: scalaire qui contient la valeur de la fonction en X
    %gX: vecteur de taille n qui contient le gradient de la fonction en X
    %kX: nombre de pas utilise
    %hist: un tableau qui contient les pas de X
    k = 0;
    x = X0;
    g = A*x + b;
    if nargout>4
    h = zeros(2,Kmax+1);
    h(:,1) = [x];
    end
    normg=norm(g);
    d=-g;
    while normg > eps && k < Kmax
        v = A*d;
        alpha = -dot(g,d)/dot(v,d);
        x = x + alpha*d;
        g = g + alpha*v;
        normg1=norm(g);
        beta=normg1^2/normg^2;
        d=-g+beta*d;
        normg=normg1;
        k = k+1;
        if nargout>4
        h = [h,x];
        end
    end
    X = x;
    if nargout > 1
        fX = Quad(X,A,b,c);
        if nargout > 2 
            gX = g;
            if nargout > 3
                kX = k;
                if nargout > 4
                    hist = h;
                end
            end
        end
    end
end
