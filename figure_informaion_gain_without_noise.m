function figure_informaion_gain_without_noise

m=[0:0.025:1];

beta=[0.5,-0.5];
tau_0=1;
eta=[0.01,0.01];

alfa=100;
nt=0;

for cs1=1:2
    for cs2=1:2
        if(cs1==cs2)
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau_0 + eta(cs1);
        else
            sig_ss(cs1,cs2)=beta(cs1)*beta(cs2)*tau_0;
        end
    end
end
for cs=1:2
    sig_rs(1,cs)=beta(cs)*tau_0;
end
g=sig_rs*pinv(sig_ss);

for cm_1=1:length(m)
    for cm_2=1:length(m)
        m_1=m(cm_1);
        m_2=m(cm_2);
        M=[m_1,0;0,m_2];
        H(cm_1,cm_2)=((g*M*sig_rs')*(sig_rs*M'*g'))/(g*M*sig_ss*M'*g' +(nt/alfa));
    end
end

figure
surf(H)
