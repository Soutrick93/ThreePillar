% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  run_all_threepillar.m
%  Master script. Runs the full three-family benchmark end to end.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

clear; clc;

% -------------------- LOCATION (automatic: this module folder) -------
cfg.proj    = fileparts(mfilename('fullpath'));   % the folder this file is in
cfg.data    = fullfile(cfg.proj, 'data');         % <-- paste ALL inputs here
cfg.results = fullfile(cfg.proj, 'results');      % outputs appear here
cfg.figures = fullfile(cfg.results, 'figures');

% -------------------- INPUT FILE NAMES (all in data\) ----------------
cfg.f_almanac      = 'almanac_synergy.csv';
cfg.f_string_links = '9606.protein.links.v12.0.txt';
cfg.f_string_info  = '9606.protein.info.v12.0.txt';
cfg.f_dgidb        = 'interactions.tsv';
cfg.f_signor       = 'all_data_15_09_26.tsv';
cfg.f_expr         = 'nci60_expr.txt';
% step6 also needs data\oneil_synergy.csv
% step10 also needs data\CRISPRGeneEffect.csv and data\Model.csv (DepMap 26Q3)

% -------------------- PARAMETERS -------------------------------------
cfg.string_thr   = 700; cfg.pos_thr = 10; cfg.neg_thr = 0;
cfg.pos_minlines = 3;   cfg.neg_maxlines = 0; cfg.dist_cap = 12;
cfg.k_drug = 30; cfg.k_cell = 30; cfg.K = 5;

% -------------------- RUN --------------------------------------------
if ~exist(cfg.results,'dir'), mkdir(cfg.results); end
if ~exist(cfg.figures,'dir'), mkdir(cfg.figures); end
addpath(fileparts(mfilename('fullpath')));

fprintf('==== three-pillar full rerun ====\nproject: %s\n\n', cfg.proj);

% ---- existing pipeline: Family 1, predictor, Figures 2-6 ----
step1_features(cfg);
step2_omics(cfg);
step3_models(cfg);
step4_figures(cfg);
step5_baseline(cfg);
cfg.ext_pos_thr = 10; cfg.ext_neg_thr = 0;
step6_external(cfg);
step7_deep(cfg);          % DeepSynergy-style MLP row of the comparison table (Table 5)

% ---- new pillars: Families 2 and 3, and Figure 7 ----
step9_context_audit(cfg.data, cfg.results);
step10_dependency_audit(cfg.data, cfg.results);
plot_pillars(cfg.results, cfg.figures);

fprintf('\n==== done. results in %s ====\n', cfg.results);
fprintf('New outputs: results\\context_audit.mat, results\\dependency_audit.mat,\n');
fprintf('             results\\figures\\Fig7_pillars.pdf\n');
