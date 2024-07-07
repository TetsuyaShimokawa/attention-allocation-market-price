function table_fitting_parameters
% fitting results

load MP_SKF MP_SKF

R=[];%[llh,aic,tau02,kappa,gamma,sita,sigma2]
R(1,1)=MP_SKF.llh;
R(1,2)=MP_SKF.aic;
R(1,3)=MP_SKF.tau02;
R(1,4)=MP_SKF.kappa;
R(1,5)=MP_SKF.gamma;
R(1,6)=MP_SKF.sita;
R(1,7)=MP_SKF.sigma2;

xlswrite('Results.xlsx',R,'fitting')
