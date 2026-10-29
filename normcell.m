% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  normcell.m
%  Helper. Normalise a cell-line name for cross-dataset matching.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function c = normcell(s)
% NORMCELL  Canonical form of a cell-line name for matching CellMiner to ALMANAC.
%   Drops the tissue prefix ("BR:"), parentheticals ("(TB)"), the "/ATCC"
%   suffix, and all non-alphanumerics, upper-cased.
%     "BR:MCF7"        -> "MCF7"
%     "LC:A549/ATCC"   -> "A549"
%     "LE:HL-60(TB)"   -> "HL60"
%     "RE:786-0"       -> "7860"
    c = upper(strtrim(s));
    k = strfind(c, ':'); if ~isempty(k), c = c(k(1)+1:end); end   % drop tissue prefix
    c = regexprep(c, '\(.*?\)', '');       % drop parentheticals
    c = regexprep(c, '/ATCC', '');         % drop ATCC subline tag
    c = regexprep(c, '[^A-Z0-9]', '');     % keep letters + digits only
end
