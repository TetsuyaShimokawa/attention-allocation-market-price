function figure_kappa_m_relation

beta=[0.5,0.5];
numb_S=length(beta);
tau_0=1;
eta=[0.01,0.01];
kappa_alfa_vec=[0:0.01:2];

for cs1=1:numb_S
    for cs2=1:numb_S
        if(cs1==cs2)
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau_0 + eta(cs1);
        else
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau_0;
        end
    end
end
for cs=1:numb_S
    sig_rs(1,cs)=beta(cs)*tau_0;
end
g=sig_rs*pinv(sig_ss);
beta_bar=tau_0*(g.*beta)';
Sig_G_t=beta_bar*(beta_bar');
for cs1=1:numb_S
    for cs2=1:numb_S
        Sig_G(cs1,cs2)=(tau_0*g(cs1)*beta(cs1))*(tau_0*g(cs2)*beta(cs2));
    end
end

I_s=eye(numb_S);
i=ones(1,numb_S);
for cka=1:length(kappa_alfa_vec)
    kappa_alfa=kappa_alfa_vec(cka);    
    m{cka}=pinv(Sig_G+2*(kappa_alfa)*I_s)*Sig_G*i';    
end

for cka=1:length(kappa_alfa_vec)
    m_plot(cka)=m{cka}(1);
end

figure
plot(kappa_alfa_vec,m_plot,'-')



