% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  read_named.m
%  Helper. Read a delimited table into a struct of named columns.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function T = read_named(path, delim, skipComments)
% READ_NAMED  Read a delimited text file with a header row into a struct.
%   T = read_named(path, delim) returns a struct whose fields are the
%   (sanitised) column names, each holding a cell array of trimmed strings.
%   delim must be a single character (',' or sprintf('\t')).
%
%   T = read_named(path, delim, true) first skips any leading lines that
%   start with '#' (metadata/comment lines, e.g. DGIdb "# Data version:")
%   and uses the next line as the header. Default is false, because some
%   files (STRING protein.info) have a header that itself starts with '#'.
%
%   Works in both MATLAB and Octave. Preserves empty fields.

    if nargin < 3, skipComments = false; end

    fid = fopen(path, 'r');
    if fid < 0
        error('read_named:open', 'Cannot open file: %s', path);
    end

    % ---- find the header line ----
    hline = '';
    while true
        ln = fgetl(fid);
        if ~ischar(ln)
            fclose(fid); error('read_named:empty', 'No header found in: %s', path);
        end
        if isempty(strtrim(ln)), continue; end                 % skip blank lines
        if skipComments && ~isempty(ln) && ln(1) == '#', continue; end  % skip comments
        hline = ln; break;
    end

    heads = local_split(hline, delim);
    n = numel(heads);
    names = cell(1, n);
    for i = 1:n
        nm = regexprep(heads{i}, '[^A-Za-z0-9_]', '');   % strip # and punctuation
        if isempty(nm), nm = sprintf('col%d', i); end
        if any(nm(1) == '0':'9'), nm = ['x' nm]; end
        names{i} = nm;
    end

    % ---- accumulate data rows ----
    data = cell(0, n);
    row = 0;
    while true
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if isempty(strtrim(ln)), continue; end
        parts = local_split(ln, delim);
        if numel(parts) < n, parts(numel(parts)+1:n) = {''}; end
        row = row + 1;
        for i = 1:n
            data{row, i} = strtrim(parts{i});
        end
    end
    fclose(fid);

    T = struct();
    for i = 1:n
        if row == 0
            T.(names{i}) = {};
        else
            T.(names{i}) = data(:, i);
        end
    end
end

function parts = local_split(ln, d)
% Quote-aware split on a single-character delimiter.
%   - does not split on a delimiter that sits inside double quotes
%   - strips one layer of surrounding double quotes from each field
%   - handles doubled "" as an escaped quote inside a quoted field
%   - preserves empty fields
% For tab/space delimited files (no quotes) this behaves like a plain split.
    parts = {};
    cur = '';
    inq = false;
    i = 1; N = numel(ln);
    while i <= N
        c = ln(i);
        if inq
            if c == '"'
                if i < N && ln(i+1) == '"'      % escaped quote ""
                    cur(end+1) = '"'; i = i + 2;
                    continue;
                else
                    inq = false; i = i + 1; continue;
                end
            else
                cur(end+1) = c; i = i + 1;
                continue;
            end
        else
            if c == '"'
                inq = true; i = i + 1; continue;
            elseif c == d
                parts{end+1} = cur; cur = ''; i = i + 1;
                continue;
            else
                cur(end+1) = c; i = i + 1;
                continue;
            end
        end
    end
    parts{end+1} = cur;
end
