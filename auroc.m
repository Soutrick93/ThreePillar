% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  auroc.m
%  Helper. Area under the ROC curve (tie-aware Mann-Whitney statistic).
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function a = auroc(y, s)
% AUROC  Area under ROC via the Mann-Whitney statistic (tie-aware).
%   y in {0,1}, s = scores. Returns NaN if only one class present.
    y = y(:); s = s(:);
    P = sum(y == 1); N = sum(y == 0);
    if P == 0 || N == 0, a = NaN; return; end
    [ss, ord] = sort(s); yy = y(ord);
    n = numel(ss); ranks = zeros(n, 1); i = 1;
    while i <= n
        j = i;
        while j < n && ss(j+1) == ss(i), j = j + 1; end
        ranks(i:j) = (i + j) / 2;    % average rank for ties
        i = j + 1;
    end
    a = (sum(ranks(yy == 1)) - P*(P+1)/2) / (P * N);
end
