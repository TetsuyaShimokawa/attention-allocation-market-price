function figure_histgram_beta_Kalman_gain
% histgrams

load MP_SKF MP_SKF

% beta
figure
histogram(MP_SKF.beta,20)
[nanmean(MP_SKF.beta),nanstd(MP_SKF.beta)]

% eta2
figure
edge_end=0.004;numb_bins=20;span_bins=edge_end/numb_bins;
edges = [0:span_bins:edge_end];
histogram(MP_SKF.eta2,edges)
[nanmean(MP_SKF.eta2),nanstd(MP_SKF.eta2)]

% kalman gain
figure
histogram(MP_SKF.g,20)
[nanmean(MP_SKF.g),nanstd(MP_SKF.g)]