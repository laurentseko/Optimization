clear
close all
clc

h=1e-5;
eps = 1e-7;
choix = "bfgs";


%%
%-------------cas test x1^2+x2^2-------------------------------------------
rho=10;
x_0=[1;1]; 
lamb=-2;
fun=@ftest_castest1;
Cont=@fcont_castest1;
Prob = @(x) Problem(fun,Cont,x);
born_min = [-10,-10];
born_max = [10,10];

[x, f, c, recap_x,recap_lamb,nb_app,fact_pen,k,recap_f,recap_c] = SQP(Prob, x_0, eps, lamb,born_min,born_max, choix, h,rho);


nb_appel_fonction_prob=cumsum(nb_app);

fprintf("--------------- premier cas test ------------\n");
s = sprintf(" %f",x);
fprintf("pts min =[%s ]\n",s);
fprintf("f_value = %f\n", f);
s = sprintf(" %f",c);
fprintf("cont =[%s ]\n",s);
s = sprintf(" %f",k);
fprintf("iterations = %f\n", k);
%%
%-------------cas test MHW4D-----------------------------------------------


x_0=[-1;2;1;-2;-2];
lamb=[100;100;100];
rho=1000;
fun = @(x) ftest_MHW4D(x);
Cont = @(x) fcont_MHW4D(x);
Prob = @(x) Problem(fun,Cont,x);
born_min = [-1000,-1000,-1000,-1000,-1000];
born_max = [1000,1000,1000,1000,1000];


[x, f, c, recap_x,recap_lamb,nb_app,fact_pen,k,recap_f,recap_c] = SQP(Prob, x_0, eps, lamb,born_min,born_max, choix, h,rho);


nb_appel_fonction_prob=cumsum(nb_app);

fprintf("---------------Probleme MHW4D------------\n");
s = sprintf(" %f",x);
fprintf("pts min =[%s ]\n",s);
fprintf("f_value = %f\n", f);
s = sprintf(" %f",c);
fprintf("cont =[%s ]\n",s);
s = sprintf(" %f",k);
fprintf("iterations = %f\n", k);
% pause(3);

%%
%-------------cas test Ariane1 --------------------------------------------

lamb = 1000;
x_0 = [4000;3641;1263];
rho=1000;
mu = 1700;
vit_p  = 11527;
vit_ejec = [2647.2 ; 2922.4 ; 4344.3];
k = [0.1101 ; 0.1532 ; 0.2154];
born_min = [10000,1000,1000];
born_max = [1000000,900000,100000];

fun = @(x) test_func(x,k,mu);
Cont = @(x) test_cont(x,vit_ejec,k,mu,vit_p);
Prob = @(x) Problem(fun,Cont,x);
fun = @ftest;
Cont = @fcont;


[masse_erg, masse_fus, c,recap_x,recap_lamb,nb_app,fact_pen,iter,recap_f,recap_c] = SQP(Prob, x_0, eps, lamb,born_min,born_max, choix, h,rho);


nb_appel_fonction_prob=cumsum(nb_app);

fprintf("---------------Probleme etagement Ariane------------\n");
s = sprintf(" %f",masse_erg);
fprintf("pts min =[%s ]\n",s);
fprintf("f_value = %f\n", masse_fus);
s = sprintf(" %f",c);
fprintf("cont =[%s ]\n",s);
s = sprintf(" %f",iter);
fprintf("iterations = %f\n", iter);
pause(3);

%%
%-------------Cas d'etude Projet Lanceur-----------------------------------

lamb = 1000;
rho=1000;
x_0 = [4000;3641;1263];
rt = 6378137;
H_c=300000;
mu = 2500;
vit_p  = 9271;
vit_cible = 7726;
vit_ejec = [2647.2 ; 2922.4 ; 4344.3];
k = [0.1101 ; 0.1532 ; 0.2154];
born_min_masses = [45000,10000,3000];
born_max_masses = [750000,90000,25000];
born_min_theta = [0,-10,-50,-50];
born_max_theta = [20,30,50,50];
    
hist_masse_erg=[];
hist_cont_PE=[];
hist_masse_fus=[];
hist_theta=[];
hist_v_fin=[];
hist_cont_PT=[];
theta_0 = [10;5;7;-5];
iter=0;
delta_v=0.2*vit_p;

 while abs(delta_v)>0.01 && iter<100
    born_min = born_min_masses;
    born_max = born_max_masses;
    fun = @(x) test_func(x,k,mu);
    Cont = @(x) test_cont(x,vit_ejec,k,mu,vit_p);
    Prob = @(x) Problem(fun,Cont,x);
    [masse_erg, masse_fus, c] = SQP(Prob, x_0, eps, lamb,born_min,born_max, choix,h,rho);
    hist_masse_erg=[hist_masse_erg,masse_erg];
    hist_cont_PE=[hist_cont_PE,c];
    hist_masse_fus=[hist_masse_fus,masse_fus];
    born_min = born_min_theta;
    born_max = born_max_theta;
    fun = @(theta) traj_lanceur(theta,masse_erg,masse_fus,H_c);
    [theta, vit_fin_sqp] = SQP(fun, theta_0, eps,[10;10],born_min,born_max, choix,1e-4,rho);
    hist_theta=[hist_theta,theta];
    [~,cont_value,vitesse_finale,Mat_ret,temps_ret] = fun(theta);
    hist_cont_PT=[hist_cont_PT,cont_value'];
    hist_v_fin=[hist_v_fin,vitesse_finale];
    delta_v = vit_cible-vitesse_finale;
    vit_p=vit_p+delta_v;
    theta_0 = theta;
    x_0 = masse_erg;
    iter=iter+1
 end

fprintf("---------------Probleme etagement------------\n");
 
s = sprintf(" %f",masse_erg);
fprintf("pts min =[%s ]\n",s);
fprintf("masse_fusee = %f\n", masse_fus);
s = sprintf(" %f",c);
fprintf("cont =[%s ]\n",s);

fprintf("---------------Probleme trajectoire------------\n");

s = sprintf(" %f",theta);
fprintf("pts min =[%s ]\n",s);
fprintf("f_value = %f\n", vitesse_finale);
s = sprintf(" %f",cont_value);
fprintf("cont =[%s ]\n",s);

s = sprintf(" %f",iter);
fprintf("iteration =%s \n",s);

mat_vit = Mat_ret(:,3:4);
C = num2cell(mat_vit, 2);
vect_norm = cellfun(@norm,C);
mat_alt = Mat_ret(:,1:2);
D= num2cell(mat_alt, 2);
vect_alt = cellfun(@norm,D);

figure(1)
plot(temps_ret,Mat_ret(:,5)); 
title('Evolution de la masse du lanceur fonction du temps');
xlabel('temps (secondes)') ;
ylabel('masse (kilogrammes');
figure(2)
plot(temps_ret,vect_alt);
title('Evolution de l altitude fonction du temps ');
xlabel('temps(secondes)') ;
ylabel('Altitude (metres)');
figure(3)
plot(temps_ret,vect_norm);
title('Evolution de la vitesse fonction du temps');
xlabel('temps(secondes)') ;
ylabel('vitesse (m/s)');
figure(4)
plot(Mat_ret(:,1),Mat_ret(:,2));
title('trajectoire dans le repere cartesien');
xlabel('x (metres)') ;
ylabel('y (metres)');

figure(5)
XCentre = 0;
YCentre = 0;
Rayon = rt;
VThetaDeg = 0:1:360;
VTheta = VThetaDeg *pi / 180;
XCercle = XCentre + Rayon * cos(VTheta);
YCercle = YCentre + Rayon * sin(VTheta);
plot(XCercle, YCercle) 
hold on 
plot(Mat_ret(:,1),Mat_ret(:,2),'g');
title('Trajectoire du lanceur');
xlabel('x (metres)');
ylabel('y (metres)');
hold off

figure(6)
XCentre = 0;
YCentre = 0;
Rayon = rt;
VThetaDeg = 0:1:360;
VTheta = VThetaDeg *pi / 180;
XCercle = XCentre + Rayon * cos(VTheta);
YCercle = YCentre + Rayon * sin(VTheta);
XCentre_orb = 0;
YCentre_orb = 0;
Rayon_orb = rt+300000;
VThetaDeg = 0:1:360;
VTheta = VThetaDeg' *pi / 180;
XCercle_orb = XCentre + Rayon_orb * cos(VTheta');
YCercle_orb = YCentre + Rayon_orb * sin(VTheta');
plot(XCercle, YCercle,XCercle_orb, YCercle_orb,'m') 
hold on 
plot(Mat_ret(:,1),Mat_ret(:,2),'g');
title('Trajectoire du lanceur (Zoom)');
xlabel('x (metres)');
ylabel('y (metres)');
hold off


fun_PE = @(mu) test_func(masse_erg,k,mu);
fun_PT = @(masse_fus) traj_lanceur(theta,masse_erg,masse_fus,0);
[m_u,R_c] = bonus(fun_PT,fun_PE,rt);
plot(R_c,m_u)
xlabel('Rayon cible (m)');
ylabel('masse (kg)');
title('Rayon cible en fonction de la masse utile');
