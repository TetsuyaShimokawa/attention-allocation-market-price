function figure_price_information_kappa_relation

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
    %m
    m{cka}=pinv(Sig_G+2*(kappa_alfa)*I_s)*Sig_G*i';   
    M=diag(m{cka});
    %H,I 0.01
    sita_onemsita=0.01/(1-0.01);    
    bunsi=(g*M*(sig_rs')*(g*M*(sig_rs'))');
    bunbo=(g*M*(sig_ss)*M*g' + (sita_onemsita^2));
    H(cka)=bunsi/bunbo;
    I(cka)=(1/2)*(log(tau_0)-log(tau_0-H(cka)));
    %H, I 0.99
    sita_onemsita=0.99/(1-0.99);
    bunsi=(g*M*(sig_rs')*(g*M*(sig_rs'))');
    bunbo=(g*M*(sig_ss)*M*g' + (sita_onemsita^2));
    H2(cka)=bunsi/bunbo;
    I2(cka)=(1/2)*(log(tau_0)-log(tau_0-H2(cka)));
end

% figure
% plot(kappa_alfa_vec,H,'-')

figure
plot(kappa_alfa_vec,I,'-')

% figure
% plot(kappa_alfa_vec,H2,'-')

figure
plot(kappa_alfa_vec,I2,'-')


