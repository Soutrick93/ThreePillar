% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  normalize_drug.m
%  Helper. Normalise a drug name (lowercase, strip salts/synonyms).
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function s = normalize_drug(name)
% NORMALIZE_DRUG  Canonicalise a drug name for matching across databases.
%   Lowercases, removes parenthetical/bracketed vendor or form notes, and
%   strips common salt / hydrate suffixes so that e.g.
%     "Erlotinib Hydrochloride"      -> "erlotinib"
%     "Navelbine Ditartrate (TN)"    -> "navelbine"
%     "Tamoxifen citrate"            -> "tamoxifen"
%   Chemical names with meaningful hyphens (e.g. 5-fluorouracil) are kept.
%   Works in MATLAB and Octave.

    if isempty(name), s = ''; return; end
    s = lower(strtrim(name));

    s = regexprep(s, '\([^)]*\)', ' ');   % drop (...) vendor/brand notes
    s = regexprep(s, '\[[^\]]*\]', ' ');   % drop [...]
    s = regexprep(s, '\s+', ' ');
    s = strtrim(s);

    salts = {' hydrochloride',' dihydrochloride',' hydrobromide',' hydroiodide', ...
             ' sulfate',' sulphate',' bisulfate',' hemisulfate',' citrate', ...
             ' dicitrate',' mesylate',' maleate',' malate',' tartrate', ...
             ' ditartrate',' bitartrate',' phosphate',' diphosphate', ...
             ' acetate',' succinate',' fumarate',' besylate',' tosylate', ...
             ' pamoate',' lactate',' gluconate',' isethionate',' sodium', ...
             ' potassium',' calcium',' chloride',' bromide',' nitrate', ...
             ' base',' anhydrous',' monohydrate',' dihydrate',' trihydrate', ...
             ' hcl',' hbr'};

    changed = true;
    while changed
        changed = false;
        for k = 1:numel(salts)
            L = numel(salts{k});
            if numel(s) > L && strcmp(s(end-L+1:end), salts{k})
                s = strtrim(s(1:end-L));
                changed = true;
            end
        end
    end

    s = regexprep(s, '\s+', ' ');
    s = strtrim(s);
end
