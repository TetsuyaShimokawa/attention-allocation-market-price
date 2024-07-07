function main

% Data
% Time series data for Topix and information prices are summarized in price_data.xlsx.

% Model-fitting results are included in results.xlsx.

% The structure MP_SKF.mat contains time-series price data and model-fitting results.

% Functions

figure_kappa_m_relation
% Make Fig. in Section 3.2

figure_informaion_gain_without_noise
% Make Fig. in Section 4.1

figure_price_information_kappa_relation
% Make Fig. in Section 4.2

MP_SKF = fitting_model_to_data;
save MP_SKF MP_SKF
% The model fitting to market data to detect parameters.

figure_histgram_beta_Kalman_gain
% Make Fig. in Section 5.1

table_fitting_parameters
% Make Table in Section 5.1

figure_m_beta_relation
% Make Fig. in Section 5.1

figure_information_loss_sita_relation
% Make Fig. in Section 5.1







