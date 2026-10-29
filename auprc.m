% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  auprc.m
%  Helper. Area under the precision-recall curve (average precision).
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function ap = auprc(y, s)
% AUPRC  Average precision (area under precision-recall curve).
%   y in {0,1}, s = scores. Returns NaN if no positives.
    y = y(:); s = s(:);
    P = sum(y == 1);
    if P == 0, ap = NaN; return; end
    [~, ord] = sort(s, 'descend'); yy = y(ord);
    tp = cumsum(yy == 1); fp = cumsum(yy == 0);
    prec = tp ./ (tp + fp);
    rec  = tp / P;
    recprev = [0; rec(1:end-1)];
    ap = sum((rec - recprev) .* prec);
end
