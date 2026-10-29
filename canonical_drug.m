% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  canonical_drug.m
%  Helper. Map a drug name to its canonical form.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function s = canonical_drug(name)
% CANONICAL_DRUG  normalize_drug plus a small, documented synonym map that
% rewrites known brand names / spelling variants in the ALMANAC drug list to
% the generic name DGIdb uses. Only high-confidence mappings are included;
% raw CAS / NSC / ChEMBL codes are deliberately left unmapped.
%
% Every entry here is a documented harmonisation to report in the paper.

    s = normalize_drug(name);

    persistent M
    if isempty(M)
        M = containers.Map('KeyType', 'char', 'ValueType', 'char');
        M('eloxatin')                 = 'oxaliplatin';
        M('carboplatinum')            = 'carboplatin';
        M('cis-platin')               = 'cisplatin';
        M('emcyt')                    = 'estramustine';
        M('trisenox')                 = 'arsenic trioxide';
        M('zanosar')                  = 'streptozocin';
        M('vepesid j')                = 'etoposide';
        M('vepesid')                  = 'etoposide';
        M('azacytidine')              = 'azacitidine';
        M('azacytidine, 5-')          = 'azacitidine';
        M('5-fluoro-2''-deoxyuridine')= 'floxuridine';
        M('antibiotic ad 32')         = 'valrubicin';
        M('antibiotic ay 22989')      = 'sirolimus';
        M('5-aminolevulinic acid')    = 'aminolevulinic acid';
        M('navelbine')                = 'vinorelbine';
        M('adm')                      = 'doxorubicin';
    end

    if isKey(M, s), s = M(s); end
end
