% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  read_cellminer.m
%  Helper. Read the CellMiner NCI-60 expression matrix.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function [cellNames, geneNames, expr] = read_cellminer(path)
% READ_CELLMINER  Parse a CellMiner NCI-60 expression download.
%   Skips the metadata preamble, finds the header row (starts with
%   "Gene name"), and returns:
%     cellNames : 1 x C cell array of the 60 cell-line column names (original)
%     geneNames : G x 1 cell array of gene symbols
%     expr      : G x C matrix of expression values (NaN where blank)
%   Columns 1-6 are metadata (gene name, entrez id, chr, start, end, cytoband).

    TAB = char(9); META = 6;
    fid = fopen(path, 'r');
    if fid < 0, error('read_cellminer:open', 'Cannot open %s', path); end

    % --- find header ---
    cellNames = {}; ncell = 0;
    while true
        ln = fgetl(fid);
        if ~ischar(ln), fclose(fid); error('read_cellminer:nohdr', 'No header found'); end
        f = strsplit_keep(ln, TAB);
        if ~isempty(f) && strncmpi(strtrim(f{1}), 'Gene name', 9)
            last = numel(f); while last > META && isempty(strtrim(f{last})), last = last - 1; end
            cellNames = strtrim(f(META+1:last));
            ncell = numel(cellNames);
            break;
        end
    end

    % --- read gene rows ---
    geneNames = {}; rows = {};
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        f = strsplit_keep(ln, TAB);
        if numel(f) < META + 1, continue; end
        g = strtrim(f{1});
        if isempty(g), continue; end
        v = nan(1, ncell);
        for c = 1:ncell
            idx = META + c;
            if idx <= numel(f), v(c) = str2double(f{idx}); end
        end
        geneNames{end+1, 1} = g;
        rows{end+1, 1} = v;
    end
    fclose(fid);
    expr = cell2mat(rows);
end

function parts = strsplit_keep(ln, d)
    pos = [0, find(ln == d), numel(ln) + 1];
    parts = cell(1, numel(pos) - 1);
    for k = 1:numel(pos) - 1
        parts{k} = ln(pos(k)+1 : pos(k+1)-1);
    end
end
