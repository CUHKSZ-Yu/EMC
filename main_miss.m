clear all; clc;
addpath('baselines');
addpath('results');

Nleave = 10; 
Ntree_list = [100, 200, 500];
ratio_list = [0.4, 0.6, 0.8];

fprintf('\nECAI-2023 Paper #419 "Highly-Efficient Robinson-Foulds Distance Estimation with Matrix Correction"');
fprintf('\nDemo: correction on incomplete trees in Section 4.1\n');

%% Demo of correction on incomplete tree sets (ref. Section 4.1)
for i = 1:length(Ntree_list)
    Ntree = Ntree_list(i);
    % Dtrue = true distance matrix calculated from complete trees
    Dtrue = csvread(['Dtrue_L',num2str(Nleave),'_N',num2str(Ntree),'.csv']);
    for j = 1:length(ratio_list)
        ratio = ratio_list(j);
        fprintf(['  Nleave = ',num2str(Nleave),', Ntree = ',num2str(Ntree),': ']);
        fprintf('missing ratio = %1.1f\n', ratio);
        % D0 = approximate distance matrix approximated from incomplete trees
        D0 = csvread(['Dmiss_L',num2str(Nleave),'_N',num2str(Ntree),'_R',num2str(10*ratio),'.csv']);

        %% Matrix Correction
        D_dc = correct_dc(D0);
        D_trf = correct_trf(D0);
        D_emc = correct_emc(D0); % our method

        %% Evaluation
        D_list = {D0, D_dc, D_trf, D_emc};
        clear mse mae vcr err
        for k = 1:length(D_list)
            D = D_list{k};
            [mse(1,k), mae(1,k), vcr(1,k)] = eval_distance(D, Dtrue);
        end
        results{j,i} = [mse; mae; vcr*100];

        statistic = array2table(roundn(results{j,i},-1), 'VariableNames', {'D^0','DC','TRF','EMC'}, 'RowNames', {'MSE','MAE','VCR/%'});
        disp(statistic);
    end
end
        