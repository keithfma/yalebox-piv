function [] = test_piv_util_print_error(uerr, verr)
%
% Print a table of absolute error quantiles 
%
% Arguments:
% 
%   uerr, verr = Double, error matrix, difference between exact and approximate
%       solutions for displacement in the x- and y-directions
% %   

qnt = 0:0.25:1;
qu = quantile(abs(uerr(:)), qnt);
qv = quantile(abs(verr(:)), qnt);
qm = quantile(sqrt(uerr(:).^2+verr(:).^2), qnt);

fprintf('displacement vector error quantiles\n');
fprintf('qnt\t| uu\t\t| vv\t\t| mag\n'); 
for i = 1:length(qnt)
    fprintf('%.2f\t| %.2e\t| %.2e\t| %.2e\n', ...
        qnt(i), qu(i), qv(i), qm(i));
end
fprintf('\n');

end