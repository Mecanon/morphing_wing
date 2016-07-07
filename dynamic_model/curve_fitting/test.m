% Before running code, remember to check if csv files are in correct
% format, otherwise just open it in excel and save it.

data_50 = csvread('flexinol_isobaric_50MPa.csv', 5, 0);

T_exp = [0, 1, 2, 3, 4, 5, 6];
eps_exp = [0, 0, 0, 0, 0, 0, 0];
T_num = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5];
eps_num = [1, 1, 1, 1, 1, 1];

r=rmse(T_exp, eps_exp, T_num, eps_num);

disp(r)