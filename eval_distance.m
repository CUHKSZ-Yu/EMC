function [mse, mae, vcr] = eval_distance(Dnew, Dtrue)
% function [mse, mae, vcr] = eval_distance(Dnew, Dtrue)
%
% @param Dnew    Calibrated distance matrix
% @param Dtrue   Gound truth of distance matrix from complete trees
% 
% @return mse    Mean squared error
% @return mae    Mean absolute error
% @return vcr    Violated constraint ratio

mse = immse(Dnew, Dtrue);
mae = mean(mean(abs(Dnew-Dtrue)));

%% check triangular inequalities
count = 0; % the number of violated triangular inequalities
n = size(Dnew, 1);
for i = 1 : n-1
    for j = i+1 : n
        num = nnz(Dnew(i,j) > Dnew(:,i)+Dnew(:,j));
        if num > 0
            count = count + num;
        end
    end
end
vcr = count / (nchoosek(n,3)*3); % Ratio of violated triangle inequalities

end

