clear all; clc;
addpath('baselines');
addpath('results');

Nleave = 10; 
Ntree_list = [100, 200, 500];
zeta_list = [0.4, 0.6, 0.8];

fprintf('\nECAI-2023 Paper #419 "Highly-Efficient Robinson-Foulds Distance Estimation with Matrix Correction"');
fprintf('\nDemo: correction on noisy tree distance in Section 4.2\n');

%% Demo of correction on noisy tree distance (ref. Section 4.2)
for i = 1:length(Ntree_list)
    Ntree = Ntree_list(i);
    % Dtrue = true distance matrix calculated from complete trees
    Dtrue = csvread(['Dtrue_L',num2str(Nleave),'_N',num2str(Ntree),'.csv']);
    for j = 1:length(zeta_list)
        zeta = zeta_list(j);
        fprintf(['  Nleave = ',num2str(Nleave),', Ntree = ',num2str(Ntree),': ']);
        fprintf('noise level = %1.1f\n', zeta);
        % D0 = noisy distance matrix by adding the Gaussian noise
        D0 = Dtrue + zeta * mean(Dtrue(:)) * normrnd(0,1,size(Dtrue));
        D0 = max(D0, 0);

        %% Matrix Correction
        D_dc = correct_dc(D0);
        D_trf = correct_trf(D0);
        D_emc = correct_emc(D0); % our method

        %% Evaluation
        D_list = {D0, D_dc, D_trf, D_emc};
        clear mse mae vcr
        for k = 1:length(D_list)
            D = D_list{k};
            [mse(1,k), mae(1,k), vcr(1,k)] = eval_distance(D, Dtrue);
        end
        results{j,i} = [mse; mae; vcr*100];

        statistic = array2table(roundn(results{j,i},-1), 'VariableNames', {'D^0','DC','TRF','EMC'}, 'RowNames', {'MSE','MAE','VCR/%'});
        disp(statistic);
    end
end
        