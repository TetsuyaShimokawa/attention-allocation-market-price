function MP_SKF = fitting_model_to_data
% Ref. market_price_fitting_linear_sparse_con.m
% Fitting an attention allocation linear sparse model to actual market prices

numb_signals=500;

% data import
Pdata=xlsread('market_prices.xlsx');%T×(S+1) [target price,signal prices]
MP=[];
for cl=1:length(Pdata(1,:))
    if(isempty(find(isnan(Pdata(:,cl))==1)))
        if(isempty(find(isinf(Pdata(:,cl))==1)))
            MP=[MP,Pdata(:,cl)];
        end
    end
end

dmy_index = randperm(length(MP(1,:))-1);%Randomly extract info stocks
info_stock_index = sort(dmy_index(1:numb_signals)+1);
MPD=MP(:,[1,info_stock_index]);
clear MP
MP=MPD;

% data_preprocessing 
target_price=MP(:,1);
signal_prices=MP(:,[2:end]);
MR=price2ret(MP);
for cs=1:length(MR(1,:))
    if(var(MR(:,cs))==0)
        cs
    else
        nMR(:,cs)=MR(:,cs)-nanmean(MR(:,cs));
    end
end


spn=length(nMR(:,1));
target_rtn=nMR([2:spn],1);
x_mat=nMR([1:spn-1],[2:end]);%Use 1-day prior return as the signals.
tau02=nanvar(target_rtn);

% capm beta
for cs=1:length(x_mat(1,:))
    yr=x_mat(:,cs);
    xr=target_rtn;
    [b,bint,reg] = regress(yr,xr);
    beta(cs)=b;
    eta2(cs)=nanvar(reg);
end

% sig_ss, sig_rs, g, sig_G
clear sig_ss sig_rs g sig_G
numb_info=length(signal_prices(1,:));
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
g=sig_rs*pinv(sig_ss);
for cs1=1:numb_info
    for cs2=1:numb_info
        sig_G(cs1,cs2)=tau02*g(cs1)*beta(cs1)*g(cs2)*beta(cs2);
    end
end


% rational belief
clear omega_r tau2
omega_r=g*x_mat';
omega_r=omega_r';%T×1
tau2=tau02-(sig_rs*pinv(sig_ss)*sig_rs');%1×1

% fitting ⇒ gamma, kappa,sita, sigma2
params0=[1,1,0.5,1];
optnew = optimset('Display','notify','MaxFunEvals',500,'TolFun',1e-6);  
func=@(params)llh_market_price_linear_sparse(target_price,x_mat,omega_r,tau2,g,sig_G,params);
lb = [0,0,0.001,0];
ub = [Inf,Inf,1,Inf];
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
[params_ml, opt_mllh, exitsig] = fmincon(func,params0,A,b,Aeq,beq,lb,ub,nonlcon,optnew);

%params
kappa=params_ml(1);
gamma=params_ml(2);
sita=params_ml(3);
sigma2=params_ml(4);
%aic
aic=2*opt_mllh+2*4;%[kappa,gamma,sita,sigma2]
%opt_m
alfa=(1/tau2)*exp(-(1/2)*(1/tau2)*(omega_r.^2));%T×1
m_mat=[];
for ct=1:length(target_rtn)
    m=pinv(sig_G+2*(kappa/alfa(ct,1))*eye(numb_info))*sig_G*ones(numb_info,1);
    m_mat(ct,:)=m';
end

% information loss
M=diag(nanmean(m_mat));
dmy_tau=(g*M*sig_rs'*(g*M*sig_rs')') / (g*M*pinv(sig_ss)*M'*g' + ((sita/(1-sita))^2)*sigma2);
Entropy_r=(1/2)*log(2*pi*tau2);
Entropy_s=(1/2)*log(2*pi*(tau02-dmy_tau));
infogain_r= (1/2)*log(2*pi*tau02)- Entropy_r;
infoloss=Entropy_s-Entropy_r;
infoloss_ratio = infoloss/abs(Entropy_r);

% results
MP_SKF.TopixData = MP(:,1);
MP_SKF.info_stock_index =info_stock_index;
MP_SKF.InfoStockPriceData = MP(:,[2:end]);
MP_SKF.tau02=tau02;
MP_SKF.tau2=tau2;
MP_SKF.kappa=kappa;
MP_SKF.gamma=gamma;
MP_SKF.sita=sita;
MP_SKF.sigma2=sigma2;
MP_SKF.llh=-opt_mllh;
MP_SKF.aic=aic;
MP_SKF.m=nanmean(m_mat);
MP_SKF.beta=beta;
MP_SKF.eta2=eta2;
MP_SKF.g=g;
MP_SKF.sig_G=sig_G;
MP_SKF.Entropy_r=Entropy_r;
MP_SKF.Entropy_s=Entropy_s;
MP_SKF.infomation_loss=infoloss;
MP_SKF.infomation_loss_ratio =infoloss_ratio;
end

function mllh = llh_market_price_linear_sparse(price,x_mat,omega_r,tau2,g,sig_G,params)
kappa=params(1);
gamma=params(2);
sita=params(3);
sigma2=params(4);

a_bar=1/tau2;

price=price-nanmean(price);

% alfa
alfa=(1/tau2)*exp(-(1/2)*(1/tau2)*(omega_r.^2));%T×1

% m llh
numb_s=length(x_mat(1,:));
llh=0;
for ct=1:length(x_mat(:,1))
    try
        m=pinv(sig_G+2*(kappa/alfa(ct,1))*eye(numb_s))*sig_G*ones(numb_s,1);
    catch
        m=ones(numb_s,1);
    end
    M=diag(m);
    omega_s(ct,1)=g*M*x_mat(ct,:)';   
    diff=(gamma*tau2*a_bar/sita)*price(ct,1) - ((1-sita)/sita)*omega_s(ct,1);
    llh=llh -(1/2)*log(2*pi*sigma2)-(1/2)*(sigma2)*(diff^2);
end

mllh=-llh;%out put minus log liklyhood
end