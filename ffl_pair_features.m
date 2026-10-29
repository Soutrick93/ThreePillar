% =====================================================================
%  ThreePillar  -  Anticancer drug-synergy mechanism benchmark
%  ffl_pair_features.m
%  Helper. Directed feed-forward-loop motif features for a drug pair.
% ---------------------------------------------------------------------
%  Author      : Dr. Soutrick Das
%  Affiliation : School of Computer Science and Artificial Intelligence,
%                SR University, Warangal, Telangana 506371, India
%  Email       : soutrickd@gmail.com
%  Repository  : https://github.com/Soutrick93/ThreePillar
% =====================================================================

function [coReg, convFFL, convCoh, convIncoh, crossReg, reach2] = ...
         ffl_pair_features(TA, TB, Ddir, Dsig)
% FFL_PAIR_FEATURES  Directed feed-forward-loop features for a drug pair.
%   TA, TB : index vectors of the two drugs' targets in the SIGNOR node space.
%   Ddir   : sparse directed adjacency (i acts on j).
%   Dsig   : sparse signed adjacency (+1/-1/0).
%
%   coReg     : # proteins directly downstream of BOTH drugs (convergence)
%   convFFL   : # feed-forward loops in which the two drug targets are the
%               loop's two regulators (X,Y) converging on a shared effector Z
%   convCoh   : of those, # coherent loops (signs consistent)
%   convIncoh : of those, # incoherent loops
%   crossReg  : # directed edges between the two target sets (one regulates
%               the other)
%   reach2    : # proteins reachable within 2 directed steps of BOTH drugs

    coReg = 0; convFFL = 0; convCoh = 0; convIncoh = 0; crossReg = 0; reach2 = 0;
    if isempty(TA) || isempty(TB), return; end
    TA = TA(:)'; TB = TB(:)';

    outA = any(Ddir(TA, :), 1);
    outB = any(Ddir(TB, :), 1);
    coReg    = full(sum(outA & outB));
    crossReg = full(nnz(Ddir(TA, TB)) + nnz(Ddir(TB, TA)));

    % 2-step downstream reachability, shared
    r2A = outA | any(Ddir(find(outA), :), 1);
    r2B = outB | any(Ddir(find(outB), :), 1);
    reach2 = full(sum(r2A & r2B));

    % convergent feed-forward loops: drug targets are the two regulators
    for x = TA
        for y = TB
            if x == y, continue; end
            if Ddir(x, y)                        % X=x, Y=y
                z = find(Ddir(x, :) & Ddir(y, :));
                z(z == x | z == y) = [];
                for z1 = z
                    convFFL = convFFL + 1;
                    [convCoh, convIncoh] = tally(Dsig(x,y), Dsig(y,z1), Dsig(x,z1), convCoh, convIncoh);
                end
            end
            if Ddir(y, x)                        % X=y, Y=x
                z = find(Ddir(y, :) & Ddir(x, :));
                z(z == x | z == y) = [];
                for z1 = z
                    convFFL = convFFL + 1;
                    [convCoh, convIncoh] = tally(Dsig(y,x), Dsig(x,z1), Dsig(y,z1), convCoh, convIncoh);
                end
            end
        end
    end
end

function [c, ic] = tally(sXY, sYZ, sXZ, c, ic)
% coherent if the direct sign equals the product of the two-step signs
    if sXY ~= 0 && sYZ ~= 0 && sXZ ~= 0
        if sign(sXY * sYZ) == sXZ, c = c + 1; else, ic = ic + 1; end
    end
end
