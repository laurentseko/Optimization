function [x, f_opt,c_value,recap_x,recap_lamb,nb_app,fact_pen,k,val_f,val_c] = SQP(Problem, x_0, eps, lamb_0,born_min, born_max, choix, h,rho)
if nargin < 8
    h = 1e-5;
    if nargin < 7 
        choix = "bfgs";
        if nargin < 6
            born_max = +Inf*ones(1,length(x_0));
            if nargin  < 5
                born_min = -Inf*ones(1,length(x_0));
            end
        end
    end
end
 x = x_0;
 lamb = lamb_0;
 fact_pen=[rho];
 val_f=[];
 val_c=[];

 [fvalue,Contvalue] = Problem(x);
 val_f=[ val_f,fvalue];
 val_c=[ val_c,Contvalue];
 nb_app(1,1)=1;
 k = 0;
 [grad_f,grad_C] = Gradient(Problem,x,h);
 H_old = eye(length(x));
 F =@(x) Func_merite(Problem,x,rho);
 [dep_x,dep_lamb] = res_pb_quad(grad_f,grad_C,H_old,Contvalue);
 lamb = dep_lamb;
 nb_app(1,k+2)=(1+2*length(x));
 [x_1,dep_x_valide,compteur] = Globalisation(F,grad_f,Contvalue,dep_x,rho,x,born_min,born_max);
 [fvalue,Contvalue] = Problem(x);
 val_f=[ val_f,fvalue];
 val_c=[ val_c,Contvalue];
nb_app(1,k+2)=nb_app(1,k+2)+compteur+1;
 if norm(x_1-x)<eps
     x_1 = x.*x_1;
 end
 recap_x = [x,x_1];
 recap_lamb=[lamb_0,lamb];
 fact_pen=[fact_pen,rho];
%  [grad_f,grad_C] = Gradient(Problem,x_1,h);
%  g = grd_lag(grad_f,grad_C,lamb);
 while norm(x_1-x)>eps && k < 100
      [~,Contvalue] = Problem(x_1);
      [grad_f_1,grad_C_1] = Gradient(Problem,x_1,h);
      [grad_f,grad_C] = Gradient(Problem,x,h);
      
      y = grd_lag(grad_f_1,grad_C_1,lamb) - grd_lag(grad_f,grad_C,lamb);
      nb_app(1,k+3)=(1+4*length(x));
      H = Hessienne(H_old,y,dep_x_valide,choix);
      [dep_x,dep_lamb] = res_pb_quad(grad_f_1,grad_C_1,H,Contvalue);
      F =@(x) Func_merite(Problem,x,rho);
      [new_x,dep_x_valide,compteur] = Globalisation(F,grad_f_1,Contvalue,dep_x,rho,x_1,born_min,born_max);
      nb_app(1,k+3)=nb_app(1,k+3)+compteur;
      if norm(new_x-x_1)<1e-20
          fprintf("pas vrai\n");
          H_old = eye(length(x));
          rho = rho + max(abs(dep_lamb));
          fact_pen=[fact_pen,rho];
          [dep_x,dep_lamb] = res_pb_quad(grad_f_1,grad_C_1,H_old,Contvalue);
          F =@(x) Func_merite(Problem,x,rho);
          [new_x,dep_x_valide,compteur] = Globalisation(F,grad_f_1,Contvalue,dep_x,rho,x_1,born_min,born_max);
          nb_app(1,k+3)=nb_app(1,k+3)+compteur;
      else
         H_old = H; 
         fact_pen=[fact_pen,rho];
      end
      

      x = x_1;
      x_1 = new_x;
      lamb = dep_lamb;
     [fvalue,Contvalue] = Problem(x);
     val_f=[ val_f,fvalue];
     val_c=[ val_c,Contvalue];
     nb_app(1,k+3)=nb_app(1,k+3)+1;
%       [grad_f,grad_C] = Gradient(Problem,x_1,h);
%       g = grd_lag(grad_f,grad_C,lamb);
      k = k+1;
      recap_x = [recap_x , x_1];
      recap_lamb=[recap_lamb,lamb];
 end
 x = x_1;
 [f_opt,c_value] = Problem(x);
 recap_x;
 recap_lamb;
 fact_pen;
k=k+2;
end