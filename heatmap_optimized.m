clc;
clear all;
close all;

addpath('source');

% run set_gurobi.m
% run set_yalmip.m

load("data/ASL_spont_01.mat");
set(0, 'DefaultLineLineWidth', 0.8);

%% prepare waveforms
time = acq_table.time;
pressure = acq_table.pressure;
flow = acq_table.flow;
volume = acq_table.volume;
pmus_true = acq_table.pmus;

[ins_marks, exp_marks] = retrieve_parity_marks(volume * 10);
fprintf("retrieved %d ins/exp marks from volume parity bit\n", length(ins_marks));

offset = 50;
idx1 = 8099;
idx2 = 8100;
peep = 5;
interval = (ins_marks(idx1)-offset):(ins_marks(idx2)-offset);
interval_table = acq_table(interval, :);

interval_table.flow = fir_filter(8, 0.2, 100, interval_table.flow);
interval_table.pressure = fir_filter(8, 0.2, 100, interval_table.pressure);

exp_start = exp_marks(idx1) - ins_marks(idx1) + offset;

waveforms = table();
waveforms.time = interval_table.time;
waveforms.paw   = interval_table.pressure - peep;
waveforms.pressure = waveforms.paw; % not ideal
waveforms.flow  = interval_table.flow;
waveforms.volume = interval_table.volume;
waveforms.pmus  = interval_table.pmus;

insexp = ones(length(waveforms.paw), 1);
insexp(exp_start:end) = 0;
waveforms.insexp = insexp;

[f, t, linkplot] = plot_dataset(interval_table);

%% build parametric optimizer once
miqp_optimizer = build_pmus_miqp_optimizer(waveforms, false, true, 0, 1e-3);

%% grid search for R and C
R_values = linspace(5, 50, 5);   % Resistance [(cmH2O.s)/L]
C_values = linspace(10, 80, 5);  % Compliance [mL / cmH2O]

error_matrix = nan(length(C_values), length(R_values));
count = 0;
for iR = 1:length(R_values)
    R = R_values(iR);
    temp_error_col = nan(length(C_values), 1);

    for iC = 1:length(C_values)
        C = C_values(iC);
        count = count + 1;
        out = miqp_optimizer({R/1000, 1/C}); % parametric solve

        pmus_hat = out{1};
        paw_hat  = out{2};

        err = norm(waveforms.paw - paw_hat);
        temp_error_col(iC) = err;
        fprintf("count = %d\n", count);
    end
    error_matrix(:, iR) = temp_error_col;
end

%%
[best_error, linear_idx] = min(error_matrix(:));
[worst_error, ~] = max(error_matrix(:));
[best_iC, best_iR] = ind2sub(size(error_matrix), linear_idx);

best_R = R_values(best_iR);
best_C = C_values(best_iC);

fprintf('\nBest R: %.2f(cmH2O.s)/L\n', best_R);
fprintf('Best C: %.2f mL/cmH2O\n', best_C);
fprintf('best error: %.2f\n', best_error);
fprintf('worst error: %.2f\n', worst_error);

%% Heatmap
best_R = R_values(best_iR);
best_C = C_values(best_iC);

figure;
imagesc(R_values, C_values(1:end), log10(error_matrix(1:end, :)));
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Resistance R ((cmH_2O.s)/L)');
ylabel('Compliance C (mL/cmH_2O)');
title('Log10 norm(P_{aw} - P_{aw}^{est}) cost surface');

hold on;
plot(best_R, best_C, 'go', 'MarkerSize', 8, 'LineWidth', 2);
[~, waveforms_hat, params_true, params_hat, ~, ~] = ...
            pmus_miqp_fixed(waveforms, false, true, 0, 1e-3, best_R/1000, 1/best_C);
plot(params_true.resistance, 1000/params_true.elastance, 'r*', 'MarkerSize', 8, 'LineWidth', 2);
legend('optimized best', 'ASL real best');

pmus_optimized = waveforms_hat.pmus;
pmus_candidate = waveforms.pressure ...
            - best_R * waveforms.flow / 60 ...
            - waveforms.volume / best_C;
paw_est = waveforms_hat.paw;
[f, t, linkplot] = plot_dataset(waveforms);
plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_optimized);
%plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_candidate);
plot(linkplot(1), waveforms.time - waveforms.time(1), paw_est);
legend(linkplot(1), "paw ASL", "paw est");
legend(linkplot(3), "pmus ASL", "pmus est");
sgtitle(sprintf('(best) MIQP with (R = %.2f, C = %.2f)', best_R, best_C));
fprintf("\ncost best (R = %.2f, C = %.2f): %.2f\n", best_R, best_C, best_error);

true_R = params_true.resistance;
true_C = 1000/params_true.elastance;
[~, true_waveforms_hat, params_true, ~, ~, ~] = ...
            pmus_miqp_fixed(waveforms, false, true, 0, 1e-3, true_R/1000, 1/true_C);
fprintf("cost true (R = %.2f, C = %.2f): %.2f\n", ...
    true_R, true_C, norm(true_waveforms_hat.paw - waveforms.paw));
[f, t, linkplot] = plot_dataset(waveforms);
plot(linkplot(3), waveforms.time - waveforms.time(1), true_waveforms_hat.pmus);
%plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_candidate);
plot(linkplot(1), waveforms.time - waveforms.time(1), true_waveforms_hat.paw);
sgtitle(sprintf('(true) MIQP with (R = %.2f, C = %.2f)', true_R, true_C));
legend(linkplot(1), "paw ASL", "paw est");
legend(linkplot(3), "pmus ASL", "pmus est");


%%
test_R = 16.81;
test_C = 36.53;
[~, test_waveforms_hat, params_true, ~, ~, ~] = ...
            pmus_miqp_fixed(waveforms, false, true, 0, 1e-3, test_R/1000, 1/test_C);
fprintf("cost test (R = %.2f, C = %.2f): %.2f\n", ...
    test_R, test_C, norm(test_waveforms_hat.paw - waveforms.paw));
[f, t, linkplot] = plot_dataset(waveforms);
plot(linkplot(3), waveforms.time - waveforms.time(1), test_waveforms_hat.pmus);
%plot(linkplot(3), waveforms.time - waveforms.time(1), pmus_candidate);
plot(linkplot(1), waveforms.time - waveforms.time(1), test_waveforms_hat.paw);
sgtitle(sprintf('(test) MIQP with (R = %.2f, C = %.2f)', test_R, test_C));
legend(linkplot(1), "paw ASL", "paw est");
legend(linkplot(3), "pmus ASL", "pmus est");