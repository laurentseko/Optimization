clear
close all
clc

h=1e-5;
eps = 1e-7;
choix = "bfgs";


x_0 = [3.5];
mu = 2500;
vit_p  = 9271;
vit_cible = 7726;
vit_ejec = [2647.2 ; 2922.4 ; 4344.3];
k = [0.1101 ; 0.1532 ; 0.2154];
Omega=(k./([1;1;1]+k));


func_PE_analytic=@(x)func_PE_analytic(x,vit_ejec,Omega,vit_p)
[f] = func_PE_analytic(x_0);
[x_3,f,iter] = Newton(func_PE_analytic, x_0, eps)
[f,x_1,x_2] = func_PE_analytic(x_3);

me3 = ((mu)/(1-k(3)*(x_3-1)))*(x_3-1);
me2 = ((mu+(1+k(3))*me3)/(1-k(2)*(x_2-1)))*(x_2-1);
me1 = ((mu+(1+k(3))*me3+(1+k(2))*me2)/(1-k(1)*(x_1-1)))*(x_1-1);

masse_erg=[me1,me2,me3]
[masse_fus]= test_func(masse_erg,k,mu)



fprintf("---------------Probleme etagement------------\n");

s = sprintf(" %f",x_1);
fprintf("x_1 =[%s ]\n",s);
s = sprintf(" %f",x_2);
fprintf("x_2 =[%s ]\n",s);
s = sprintf(" %f",x_3);
fprintf("x_3 =[%s ]\n",s);
fprintf("f_value = %f\n", f);
s = sprintf(" %f",iter);
fprintf("iteration =%s \n",s);
s = sprintf(" %f",masse_erg);
fprintf("masses ergol =[%s ]\n",s);
fprintf("masse_fusee = %f\n", masse_fus);
