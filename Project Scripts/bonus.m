function [m_u,R_c] = bonus(fun_PT,fun_PE,rt)
m_u=linspace(0,24000,400);
for i=1:length(m_u)
    [f] = fun_PE(m_u(i));
    [vit_fin,cont_value,vit_reelle,Mat_ret,temps_ret] =fun_PT(f);
    R_c(i)=cont_value(1)+rt;
end
end

