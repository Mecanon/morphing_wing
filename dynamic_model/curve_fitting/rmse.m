function r=rmse(T_exp, eps_exp, T_num, eps_num)
% Function to calculate root mean square error from a data vector or matrix 
% and the corresponding estimates.
% Usage: r=rmse(data,estimate)
% Note: data and estimates have to be of same size
% Example: r=rmse(randn(100,100),randn(100,100));

data = interp1(T_exp, eps_exp, T_num);
estimate = eps_num;

r=sqrt(sum((data(:)-estimate(:)).^2)/numel(data));