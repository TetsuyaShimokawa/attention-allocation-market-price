function figure_m_beta_relation
% m and beta relation

load MP_SKF MP_SKF

% m and beta relation
figure
scatter(MP_SKF.beta,MP_SKF.m)

% m and Kalman_gain relation
figure
scatter(MP_SKF.g,MP_SKF.m)