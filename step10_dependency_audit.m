% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  step10_dependency_audit.m
%  Mechanism family 3. DepMap functional genetic dependency audit
%  (Table 3). Input: data/ ; output: results/dependency_audit.mat
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function step10_dependency_audit(dataDir, outDir)
% Reads inputs from dataDir and writes outputs to outDir. With NO arguments
% it uses the folder this script sits in for both, so scripts + data +
% outputs live together. (run_all_threepillar passes data\ and results\.)
if nargin < 1 || isempty(dataDir), dataDir = fullfile(fileparts(mfilename('fullpath')),'data'); end
if nargin < 2 || isempty(outDir),  outDir  = fullfile(fileparts(mfilename('fullpath')),'results'); end
clc;

cfg.almanac = fullfile(dataDir,'almanac_synergy.csv');
cfg.dgidb   = fullfile(dataDir,'interactions.tsv');
cfg.crispr  = fullfile(dataDir,'CRISPRGeneEffect.csv');
cfg.model   = fullfile(dataDir,'Model.csv');
cfg.omics   = local_find_omics(dataDir, outDir);   % omics_features.mat (for predictive)
cfg.pos_thr = 10; cfg.neg_thr = 0;
cfg.dep_mod = -0.5; cfg.dep_str = -1.0;
cfg.nperm   = 1500; cfg.nbin = 10;
cfg.out     = fullfile(outDir,'dependency_audit.mat');
rng(0); TAB = sprintf('\t');
fprintf('=== step10: functional dependency audit (family 3) ===\n');

% ----------------------------- drug -> targets -----------------------
dg = read_named(cfg.dgidb, TAB, true);
has_claim = isfield(dg,'drug_claim_name');
drug2tg = containers.Map('KeyType','char','ValueType','any');
for i = 1:numel(dg.drug_name)
    g = strtrim(dg.gene_name{i}); if isempty(g), continue; end
    ky = {};
    d1 = normalize_drug(dg.drug_name{i}); if ~isempty(d1), ky{end+1}=d1; end
    if has_claim
        d2 = normalize_drug(dg.drug_claim_name{i});
        if ~isempty(d2) && (isempty(ky)||~strcmp(d2,ky{1})), ky{end+1}=d2; end
    end
    for k = 1:numel(ky)
        d = ky{k};
        if isKey(drug2tg,d), lst = drug2tg(d); else, lst = {}; end
        drug2tg(d) = unique([lst, {g}]);
    end
end

% ----------------------------- DepMap matrix -------------------------
fprintf('reading DepMap Chronos matrix (slow) ...\n');
opts = detectImportOptions(cfg.crispr, 'VariableNamingRule','preserve');
vn   = opts.VariableNames;
opts = setvartype(opts, vn(1),     'char');     % first column = ModelID
opts = setvartype(opts, vn(2:end), 'double');   % gene columns numeric; blanks -> NaN
T       = readtable(cfg.crispr, opts);
modelID = string(T{:,1});
geneSym = regexprep(string(T.Properties.VariableNames(2:end)), '\s*\(.*\)$','');
depVals = T{:,2:end};  clear T
fprintf('  DepMap: %d models x %d genes\n', numel(modelID), numel(geneSym));
gcol    = containers.Map(cellstr(geneSym), num2cell(1:numel(geneSym)));
model2row = containers.Map(cellstr(modelID), num2cell(1:numel(modelID)));

% ----------------------------- cell -> ModelID -----------------------
MD = read_named(cfg.model, ',');
alias = containers.Map({'U251','HL60TB','MDAMB435','SR'}, ...
                       {'U251MG','HL60','MDAMB435S','SR786'});
norm2model = containers.Map('KeyType','char','ValueType','char');
flds = {'CellLineName','StrippedCellLineName'};
for i = 1:numel(MD.ModelID)
    for ff = 1:numel(flds)
        if isfield(MD, flds{ff})
            v = MD.(flds{ff}){i};
            if ~isempty(v)
                k = normcell(v);
                if ~isKey(norm2model,k), norm2model(k) = MD.ModelID{i}; end
            end
        end
    end
end

% ----------------------------- drug-target columns in DepMap ---------
drugtg = containers.Map('KeyType','char','ValueType','any');   % drug -> col indices
dk = keys(drug2tg);
for i = 1:numel(dk)
    g = drug2tg(dk{i}); cols = [];
    for j = 1:numel(g), if isKey(gcol,g{j}), cols(end+1)=gcol(g{j}); end; end
    drugtg(dk{i}) = cols;
end

% ----------------------------- labelled rows -------------------------
S = read_named(cfg.almanac, ',');
dr = S.drug_row; dc = S.drug_col; cl = S.cell_line_name;
zip = cellfun(@str2double, S.synergy_zip);
N = numel(dr);

% map each ALMANAC cell -> DepMap row (once)
ualmcell = unique(cl);
cell2row = containers.Map('KeyType','char','ValueType','double');
for i = 1:numel(ualmcell)
    k = normcell(ualmcell{i}); mm = '';
    if isKey(norm2model,k), mm = norm2model(k);
    elseif isKey(alias,k) && isKey(norm2model,alias(k)), mm = norm2model(alias(k)); end
    if ~isempty(mm) && isKey(model2row,mm), cell2row(ualmcell{i}) = model2row(mm); end
end

% precompute per-(drug,row) dependency aggregate for the rows we use
usedrows = unique(cell2vals(cell2row));
aggM = containers.Map('KeyType','char','ValueType','any');
    % aggregate = [min mean nMod nStr nTot]
for r = usedrows(:)'
    for i = 1:numel(dk)
        cols = drugtg(dk{i});
        if isempty(cols), continue; end
        s = depVals(r, cols); s = s(~isnan(s));
        if isempty(s), continue; end
        aggM(sprintf('%s@%d', dk{i}, r)) = ...
            [min(s), mean(s), sum(s<cfg.dep_mod), sum(s<cfg.dep_str), numel(s)];
    end
end

% assemble rows
rA={}; rB={}; grp=[]; y=[]; mrow=[]; X=[]; Z=[];
ug = containers.Map('KeyType','double','ValueType','double');
for i = 1:N
    if isnan(zip(i)), continue; end
    a = canonical_drug(dr{i}); b = canonical_drug(dc{i});
    if strcmp(a,b) || ~isKey(cell2row,cl{i}), continue; end
    r = cell2row(cl{i});
    ka = sprintf('%s@%d',a,r); kb = sprintf('%s@%d',b,r);
    if ~isKey(aggM,ka) || ~isKey(aggM,kb), continue; end
    if zip(i) >= cfg.pos_thr, lab=1; elseif zip(i) <= cfg.neg_thr, lab=0; else, continue; end
    A = aggM(ka); B = aggM(kb);
    X(end+1,:) = [A(1)+B(1), A(2)+B(2), A(3)*B(3), A(3)+B(3), A(4)*B(4), double(A(3)>0&&B(3)>0)];
    Z(end+1,1) = A(5)*B(5);
    if ~isKey(ug,r), ug(r)=ug.Count+1; end
    grp(end+1,1)=ug(r); y(end+1,1)=lab; mrow(end+1,1)=r;
end
m = numel(y);
fprintf('dependency-labelled (pair,cell) rows: %d  pos=%.1f%%  cells=%d\n', ...
        m, 100*mean(y), numel(unique(mrow)));

% ----------------------------- audit ---------------------------------
FN = {'dep_sum_min','dep_mean_sum','ndep_prod','ndep_sum','strong_prod','both_dep'};
res = struct(); res.feat = FN; res.p = zeros(1,6);
for j = 1:6
    res.auc(j) = auroc(y, X(:,j));
    res.wc(j)  = local_within_cell(y, X(:,j), grp);
    res.p(j)   = local_cmperm(y, X(:,j), Z, cfg.nperm, cfg.nbin);
end
res.q = local_bh(res.p);
res.baseline_auc = auroc(y, Z);
res.baseline_wc  = local_within_cell(y, Z, grp);
res.p_cellcount_dep_sum_min = local_cmperm_cell(y, X(:,1), Z, grp, 1000, 5);
res.n = m; res.ncells = numel(unique(mrow));

fprintf('\n%-14s %8s %10s %7s %6s\n','feature','AUROC','within-cell','p','q');
for j = 1:6
    fprintf('%-14s %8.3f %10.3f %7.3f %6.3f\n', FN{j}, res.auc(j), res.wc(j), res.p(j), res.q(j));
end
fprintf('target-count baseline within-cell AUROC = %.3f\n', res.baseline_wc);
fprintf('dep_sum_min cell+count-matched p = %.4f\n', res.p_cellcount_dep_sum_min);

% ----------------------------- predictive increment ------------------
if ~isempty(cfg.omics) && exist(cfg.omics,'file')
    fprintf('predictive test: using %s\n', cfg.omics);
    O = load(cfg.omics);
    OX=O.X; oy=double(O.y(:)); og=O.groups; oDA=O.rDA; oDB=O.rDB; oC=O.rC;
    omcols = find(og<=2); n2 = numel(oy);
    depB = nan(n2,6); keep = false(n2,1);
    for i = 1:n2
        if ~isKey(cell2row,oC{i}), continue; end
        r = cell2row(oC{i});
        ka=sprintf('%s@%d',canonical_drug(oDA{i}),r); kb=sprintf('%s@%d',canonical_drug(oDB{i}),r);
        if ~isKey(aggM,ka)||~isKey(aggM,kb), continue; end
        A=aggM(ka); B=aggM(kb);
        depB(i,:)=[A(1)+B(1),A(2)+B(2),A(3)*B(3),A(3)+B(3),A(4)*B(4),double(A(3)>0&&B(3)>0)];
        keep(i)=true;
    end
    Xom=OX(keep,omcols); Xod=[Xom, depB(keep,:)]; yy=oy(keep);
    da=oDA(keep); db=oDB(keep); cc=oC(keep);
    drugs=unique([da;db]); dfold=containers.Map(drugs,num2cell(mod(0:numel(drugs)-1,5)+1));
    cells=unique(cc);      cfold=containers.Map(cells,num2cell(mod(0:numel(cells)-1,5)+1));
    fA=cellfun(@(d)dfold(d),da); fB=cellfun(@(d)dfold(d),db); fC=cellfun(@(c)cfold(c),cc);
    fprintf('\npredictive rows with dependency: %d (pos %.1f%%)\n', numel(yy), 100*mean(yy==1));
    prot={'LDO','LCO'};
    for pm=1:2
        r_om=[]; r_od=[];
        for k=1:5
            if pm==1, te=(fA==k)|(fB==k); tr=~te; else, te=(fC==k); tr=~te; end
            if sum(tr)<50 || numel(unique(yy(tr)))<2 || numel(unique(yy(te)))<2, continue; end
            r_om(end+1)=local_rusboost(Xom(tr,:),yy(tr),Xom(te,:),yy(te));
            r_od(end+1)=local_rusboost(Xod(tr,:),yy(tr),Xod(te,:),yy(te));
        end
        d=r_od-r_om; ci=1.96*std(d)/sqrt(numel(d));
        res.(['pred_' prot{pm}])=[mean(r_om),mean(r_od),mean(d),ci];
        fprintf('%s omics %.3f  omics+dep %.3f  Delta %+.4f (95%% CI %+.4f, %+.4f)\n', ...
            prot{pm},mean(r_om),mean(r_od),mean(d),mean(d)-ci,mean(d)+ci);
    end
else
    warning('omics_features.mat not found; skipping predictive increment.');
end

DEP_RESULT = res;
if ~exist(outDir,'dir'), mkdir(outDir); end
save(cfg.out, 'DEP_RESULT');
fprintf('saved -> %s\n', cfg.out);
end

% ===================== local helper functions =======================
function v = cell2vals(mapObj)
    k = keys(mapObj); v = zeros(1,numel(k));
    for i=1:numel(k), v(i)=mapObj(k{i}); end
end
function p = local_find_omics(dataDir, outDir)
% Locate omics_features.mat (written by step2) for the predictive test:
% check outDir and dataDir, then a results\ subfolder, then results\ in up
% to three parent folders. Returns '' if not found (predictive is skipped).
    cand = {fullfile(outDir,'omics_features.mat'), ...
            fullfile(dataDir,'omics_features.mat'), ...
            fullfile(dataDir,'results','omics_features.mat')};
    u = dataDir;
    for k = 1:3, u = fileparts(u); cand{end+1} = fullfile(u,'results','omics_features.mat'); end
    p = '';
    for i = 1:numel(cand)
        if exist(cand{i},'file'), p = cand{i}; return; end
    end
end
function a = local_within_cell(y,f,grp)
    u=unique(grp); acc=[]; w=[];
    for k=1:numel(u)
        sel=grp==u(k); yv=y(sel);
        if sum(yv)<3 || sum(1-yv)<3, continue; end
        a1=auroc(yv,f(sel)); if ~isnan(a1), acc(end+1)=a1; w(end+1)=sum(sel); end
    end
    a=sum(acc(:).*w(:))/sum(w);
end
function p = local_cmperm(y,f,Z,nperm,nbin)
    y=y(:); rr=tiedrank(f); n1=sum(y==1); n0=sum(y==0);
    obs=(sum(rr(y==1))-n1*(n1+1)/2)/(n1*n0);
    e=quantile(Z,linspace(0,1,nbin+1)); e=unique(e); e(1)=e(1)-1e-9;
    nb=numel(e)-1; b=discretize(Z,e);
    bi=cell(nb,1); bk=zeros(nb,1);
    for k=1:nb, bi{k}=find(b==k); bk(k)=sum(y(bi{k})==1); end
    cnt=0;
    for it=1:nperm
        pos=false(numel(y),1);
        for k=1:nb, if bk(k)>0, ix=bi{k}; pos(ix(randperm(numel(ix),bk(k))))=true; end; end
        st=(sum(rr(pos))-n1*(n1+1)/2)/(n1*n0);
        if abs(st-0.5)>=abs(obs-0.5), cnt=cnt+1; end
    end
    p=(cnt+1)/(nperm+1);
end
function p = local_cmperm_cell(y,f,Z,grp,nperm,nb)
    y=y(:); obs=local_within_cell(y,f,grp); strata={}; sk=[]; ug=unique(grp);
    for g=1:numel(ug)
        ix=find(grp==ug(g)); Zc=Z(ix);
        e=quantile(Zc,linspace(0,1,nb+1)); e=unique(e);
        if numel(e)<2, e=[min(Zc)-1, max(Zc)+1]; end
        e(1)=e(1)-1e-9; bb=discretize(Zc,e);
        for k=1:numel(e)-1, s=ix(bb==k); if ~isempty(s), strata{end+1}=s; sk(end+1)=sum(y(s)==1); end; end
    end
    cnt=0;
    for it=1:nperm
        yp=zeros(numel(y),1);
        for s=1:numel(strata), if sk(s)>0, ix=strata{s}; yp(ix(randperm(numel(ix),sk(s))))=1; end; end
        if abs(local_within_cell(yp,f,grp)-0.5)>=abs(obs-0.5), cnt=cnt+1; end
    end
    p=(cnt+1)/(nperm+1);
end
function q = local_bh(p)
    p=p(:); mm=numel(p); [~,o]=sort(p); q=zeros(mm,1); prev=1;
    for i=mm:-1:1, prev=min(prev,p(o(i))*mm/i); q(o(i))=min(prev,1); end
    q=q';
end
function a = local_rusboost(Xtr,ytr,Xte,yte)
    t=templateTree('MaxNumSplits',8);
    Mdl=fitcensemble(Xtr,ytr,'Method','RUSBoost','NumLearningCycles',120,'Learners',t,'LearnRate',0.1);
    [~,sc]=predict(Mdl,Xte); pos=find(Mdl.ClassNames==1); a=auroc(yte,sc(:,pos));
end
