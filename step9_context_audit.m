% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  step9_context_audit.m
%  Mechanism family 2. Expression-gated, context-specific motif
%  audit (Table 2). Input: data/ ; output: results/context_audit.mat
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function step9_context_audit(dataDir, outDir)
% Reads inputs from dataDir and writes outputs to outDir. With NO arguments
% it uses the folder this script sits in for both, so the scripts, the data
% files and the outputs can all live together in one folder.
% (run_all_threepillar passes data\ and results\ explicitly.)
if nargin < 1 || isempty(dataDir), dataDir = fullfile(fileparts(mfilename('fullpath')),'data'); end
if nargin < 2 || isempty(outDir),  outDir  = fullfile(fileparts(mfilename('fullpath')),'results'); end
clc;

% ----------------------------- configuration -------------------------
cfg.almanac    = fullfile(dataDir,'almanac_synergy.csv');
cfg.dgidb      = fullfile(dataDir,'interactions.tsv');
cfg.signor_tsv = fullfile(dataDir,'all_data_15_09_26.tsv');  % raw SIGNOR
cfg.expr       = fullfile(dataDir,'nci60_expr.txt');
cfg.tau       = 0.5;     % a protein is "active" in a cell if expr > tau
cfg.pos_thr   = 10;      % ZIP >= pos_thr : synergistic row
cfg.neg_thr   = 0;       % ZIP <= neg_thr : non-synergistic row
cfg.nperm     = 1500;    % confound-matched permutations
cfg.nbin      = 10;      % quantile bins of target-count product
cfg.out       = fullfile(outDir,'context_audit.mat');
rng(0);

TAB = sprintf('\t');
fprintf('=== step9: context-specific motif audit (family 2) ===\n');

% ----------------------------- drug -> targets -----------------------
dg = read_named(cfg.dgidb, TAB, true);
has_claim = isfield(dg,'drug_claim_name');
drug2tg = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(dg.drug_name)
    g = strtrim(dg.gene_name{i});
    if isempty(g), continue; end
    ky = {};
    d1 = normalize_drug(dg.drug_name{i}); if ~isempty(d1), ky{end+1}=d1; end
    if has_claim
        d2 = normalize_drug(dg.drug_claim_name{i});
        if ~isempty(d2) && (isempty(ky) || ~strcmp(d2,ky{1})), ky{end+1}=d2; end
    end
    for k = 1:numel(ky)
        d = ky{k};
        if isKey(drug2tg,d), lst = drug2tg(d); else, lst = {}; end
        drug2tg(d) = unique([lst, {g}]);
    end
end
fprintf('DGIdb: %d drug keys\n', drug2tg.Count);

% ----------------------------- SIGNOR network ------------------------
[dsym, Ddir, Dsig] = local_build_signor(cfg.signor_tsv);   % from raw TSV
sig2idx = containers.Map(dsym, num2cell(1:numel(dsym)));
nnode = numel(dsym);
fprintf('SIGNOR: %d nodes, %d directed edges\n', nnode, nnz(Ddir));

% ----------------------------- expression ----------------------------
[cellNames, geneNames, expr] = read_cellminer(cfg.expr);
expr(isnan(expr)) = 0;
g2row = containers.Map(geneNames, num2cell(1:numel(geneNames)));
col_of = containers.Map('KeyType','char','ValueType','double');
for j = 1:numel(cellNames), col_of(normcell(cellNames{j})) = j; end

% active[node, cell] : gene expressed > tau, or (no expression) kept active
ncell = size(expr,2);
active = true(nnode, ncell);
for n = 1:nnode
    if isKey(g2row, dsym{n})
        active(n,:) = expr(g2row(dsym{n}), :) > cfg.tau;
    end
end
fprintf('expression: %d genes x %d cells; mean active frac %.2f\n', ...
        size(expr,1), ncell, mean(active(:)));

% ----------------------------- labelled rows -------------------------
S = read_named(cfg.almanac, ',');
dr = S.drug_row; dc = S.drug_col; cl = S.cell_line_name;
zip = cellfun(@str2double, S.synergy_zip);
N = numel(dr);

rowA = zeros(N,1); rowB = zeros(N,1); rcol = zeros(N,1); y = zeros(N,1);
% cache per-drug sig targets
dcanon = containers.Map('KeyType','char','ValueType','any');
m = 0;
for i = 1:N
    if isnan(zip(i)), continue; end
    a = canonical_drug(dr{i}); b = canonical_drug(dc{i});
    if strcmp(a,b), continue; end
    if ~isKey(dcanon,a), dcanon(a)=local_targets(a, drug2tg, sig2idx); end
    if ~isKey(dcanon,b), dcanon(b)=local_targets(b, drug2tg, sig2idx); end
    TA = dcanon(a); TB = dcanon(b);
    if isempty(TA) || isempty(TB), continue; end
    key = normcell(cl{i});
    if ~isKey(col_of,key), continue; end
    if zip(i) >= cfg.pos_thr, lab = 1; elseif zip(i) <= cfg.neg_thr, lab = 0; else, continue; end
    m = m+1; rowA(m)=i; rowB(m)=i; rcol(m)=col_of(key); y(m)=lab;
    adrug{m}=a; bdrug{m}=b;
end
rcol=rcol(1:m); y=y(1:m);
fprintf('labelled (pair,cell) rows: %d  pos=%.1f%%  cells=%d\n', ...
        m, 100*mean(y), numel(unique(rcol)));

% ----------------------------- motif features ------------------------
FN = {'coReg','convFFL','convCoh','convIncoh','crossReg','reach2'};
STAT = zeros(m,6); CTX = zeros(m,6);
nTA = zeros(m,1); nTB = zeros(m,1);

% static motif per unique pair (network is identical across cells)
pairkey = cell(m,1);
for i = 1:m
    sp = sort({adrug{i}, bdrug{i}}); pairkey{i} = [sp{1} '|' sp{2}];
end
staticMap = containers.Map('KeyType','char','ValueType','any');
for i = 1:m
    if ~isKey(staticMap, pairkey{i})
        TA = dcanon(adrug{i}); TB = dcanon(bdrug{i});
        [c1,c2,c3,c4,c5,c6] = ffl_pair_features(TA, TB, Ddir, Dsig);
        staticMap(pairkey{i}) = [c1 c2 c3 c4 c5 c6];
    end
end
for i = 1:m, STAT(i,:) = staticMap(pairkey{i}); end
for i = 1:m
    nTA(i) = numel(dcanon(adrug{i})); nTB(i) = numel(dcanon(bdrug{i}));
end

% context motif per (pair,cell): gate the network to that cell's active nodes
cells = unique(rcol);
for ci = 1:numel(cells)
    c = cells(ci);
    act = active(:,c);
    Dc = Ddir;  Dc(~act,:) = false; Dc(:,~act) = false;
    Sc = Dsig;  Sc(~act,:) = 0;     Sc(:,~act) = 0;
    idx = find(rcol == c);
    cache = containers.Map('KeyType','char','ValueType','any');
    for jj = 1:numel(idx)
        i = idx(jj);
        if isKey(cache, pairkey{i})
            CTX(i,:) = cache(pairkey{i}); continue;
        end
        TA = dcanon(adrug{i}); TB = dcanon(bdrug{i});
        TA = TA(act(TA)); TB = TB(act(TB));     % targets active in this cell
        [c1,c2,c3,c4,c5,c6] = ffl_pair_features(TA, TB, Dc, Sc);
        CTX(i,:) = [c1 c2 c3 c4 c5 c6]; cache(pairkey{i}) = CTX(i,:);
    end
    if mod(ci,10)==0, fprintf('  gated %d / %d cells\n', ci, numel(cells)); end
end

% ----------------------------- audit ---------------------------------
Z = nTA .* nTB;                          % promiscuity (target-count) product
tab = struct('feature',FN,'auc_static',[],'auc_context',[], ...
             'wc_static',[],'wc_context',[],'p',[],'q',[]);
ps = zeros(1,6);
res = struct();
for j = 1:6
    res.auc_static(j)  = auroc(y, STAT(:,j));
    res.auc_context(j) = auroc(y, CTX(:,j));
    res.wc_static(j)   = local_within_cell(y, STAT(:,j), rcol);
    res.wc_context(j)  = local_within_cell(y, CTX(:,j), rcol);
    ps(j) = local_cmperm(y, CTX(:,j), Z, cfg.nperm, cfg.nbin);
end
res.q = local_bh(ps); res.p = ps; res.feat = FN;
res.baseline_wc  = local_within_cell(y, Z, rcol);
res.baseline_auc = auroc(y, Z);
res.n = m; res.ncells = numel(cells);

fprintf('\n%-10s %8s %8s %9s %9s %7s %6s\n', ...
        'feature','AUCstat','AUCctx','wcStat','wcCtx','p','q');
for j = 1:6
    fprintf('%-10s %8.3f %8.3f %9.3f %9.3f %7.3f %6.3f\n', FN{j}, ...
        res.auc_static(j),res.auc_context(j),res.wc_static(j),res.wc_context(j),res.p(j),res.q(j));
end
fprintf('target-count baseline  within-cell AUROC = %.3f\n', res.baseline_wc);

CTX_RESULT = res;
if ~exist(outDir,'dir'), mkdir(outDir); end
save(cfg.out, 'CTX_RESULT');
fprintf('saved -> %s\n', cfg.out);
end

% ===================== local helper functions =======================
function idx = local_targets(d, drug2tg, sig2idx)
    idx = [];
    if isKey(drug2tg,d)
        g = drug2tg(d);
        for k = 1:numel(g)
            if isKey(sig2idx,g{k}), idx(end+1) = sig2idx(g{k}); end
        end
    end
    idx = unique(idx);
end

function a = local_within_cell(y, f, grp)
    u = unique(grp); acc = []; w = [];
    for k = 1:numel(u)
        sel = grp==u(k); yv = y(sel);
        if sum(yv)<3 || sum(1-yv)<3, continue; end
        a1 = auroc(yv, f(sel));
        if ~isnan(a1), acc(end+1)=a1; w(end+1)=sum(sel); end
    end
    a = sum(acc(:).*w(:))/sum(w);
end

function p = local_cmperm(y, f, Z, nperm, nbin)
    y = y(:); rr = tiedrank(f); n1 = sum(y==1); n0 = sum(y==0);
    obs = (sum(rr(y==1)) - n1*(n1+1)/2)/(n1*n0);
    e = quantile(Z, linspace(0,1,nbin+1)); e = unique(e); e(1) = e(1)-1e-9;
    nb = numel(e)-1; b = discretize(Z, e);
    binidx = cell(nb,1); bink1 = zeros(nb,1);
    for k = 1:nb, binidx{k} = find(b==k); bink1(k) = sum(y(binidx{k})==1); end
    cnt = 0;
    for p_ = 1:nperm
        pos = false(numel(y),1);
        for k = 1:nb
            if bink1(k) > 0
                ix = binidx{k};
                sel = ix(randperm(numel(ix), bink1(k)));
                pos(sel) = true;
            end
        end
        stat = (sum(rr(pos)) - n1*(n1+1)/2)/(n1*n0);
        if abs(stat-0.5) >= abs(obs-0.5), cnt = cnt + 1; end
    end
    p = (cnt+1)/(nperm+1);
end

function q = local_bh(p)
    p = p(:); mm = numel(p); [~,o] = sort(p); q = zeros(mm,1); prev = 1;
    for i = mm:-1:1
        prev = min(prev, p(o(i))*mm/i); q(o(i)) = min(prev,1);
    end
    q = q';
end

function [dsym, Ddir, Dsig] = local_build_signor(tsvpath)
% Build the directed, signed SIGNOR network from the raw complete-dataset
% TSV (human protein->protein causal edges). Mirrors the pipeline's step3a.
    TAB = sprintf('\t');
    R = read_named(tsvpath, TAB);
    a = R.ENTITYA; b = R.ENTITYB; eff = R.EFFECT;
    keep = strcmp(R.TAX_ID,'9606') & strcmp(R.TYPEA,'protein') & strcmp(R.TYPEB,'protein');
    for i = 1:numel(a)
        if keep(i) && (isempty(strtrim(a{i})) || isempty(strtrim(b{i}))), keep(i)=false; end
    end
    a = a(keep); b = b(keep); eff = eff(keep);
    sgn = zeros(numel(a),1);
    for i = 1:numel(a)
        e = lower(strtrim(eff{i}));
        if strncmp(e,'up-regulates',12),       sgn(i) =  1;
        elseif strncmp(e,'down-regulates',14),  sgn(i) = -1;
        end
    end
    dsym = unique([a; b]);
    sidx = containers.Map(dsym, num2cell(1:numel(dsym)));
    ei = cellfun(@(s) sidx(s), a); ej = cellfun(@(s) sidx(s), b);
    nnode = numel(dsym);
    sel = ei ~= ej; ei = ei(sel); ej = ej(sel); sgn = sgn(sel);
    P = sparse(ei, ej, 1,   nnode, nnode);
    S = sparse(ei, ej, sgn, nnode, nnode);
    Ddir = logical(P > 0);
    Dsig = spfun(@sign, S);
end
