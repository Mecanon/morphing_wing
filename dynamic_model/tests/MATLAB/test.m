function []=test(a_list)
a_list = cell2mat(a_list);
disp(a_list)
disp(a_list(1))
disp(a_list(1,1))
%size(a_list)
class(a_list(1))
a_list(1) == 1
a_list*2
end