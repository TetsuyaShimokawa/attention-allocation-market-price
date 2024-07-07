function figure_information_loss_sita_relation
% information_loss and sita relation

load MP_SKF MP_SKF

tau2=MP_SKF.tau2;
tau02=MP_SKF.tau02;
sigma2=MP_SKF.sigma2;
beta=MP_SKF.beta;
eta2=MP_SKF.eta2;
g=MP_SKF.g;
m=MP_SKF.m;

numb_info=length(g);
for cs1=1:numb_info
    for cs2=1:numb_info
        if(cs1==cs2)
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau02 + eta2(cs1);
        else
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau02;
        end
    end
end
for cs=1:numb_info
    sig_rs(1,cs)=beta(cs)*tau02;
end

sita_vec=[0:0.05:0.95];
Entropy_r=(1/2)*log(2*pi*tau2);
for csita=1:length(sita_vec)
    sita=sita_vec(csita);    
    M=diag(m);
    dmy_tau=(g*M*sig_rs'*(g*M*sig_rs')') / (g*M*pinv(sig_ss)*M'*g' + ((sita/(1-sita))^2)*sigma2);    
    Entropy_s=(1/2)*log(2*pi*(tau02-dmy_tau));
    infoloss_vec(csita)=Entropy_s-Entropy_r;
end

figure
plot(sita_vec,infoloss_vec)