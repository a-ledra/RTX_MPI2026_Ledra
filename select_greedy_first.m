function finalIdx = select_greedy_first(feasible, nSelect)
%SELECT_GREEDY_FIRST Pick codebook rows with Dr. Edwards' "Greedy First" algorithm.
%   Starting from the largest feasible Walsh score and working downward,
%   repeatedly chooses the feasible entry that is at an odd distance below
%   the last pick and at least OPTDIS away from it, where
%   OPTDIS = last/(j-1/2) (the basin-of-attraction spacing). If no entry
%   qualifies, OPTDIS is reduced by 1 and the search is repeated. 
%
%   FINALIDX = SELECT_GREEDY_FIRST(FEASIBLE) selects 16 rows.
%
%   FINALIDX = SELECT_GREEDY_FIRST(FEASIBLE, NSELECT) selects NSELECT rows
%   (NSELECT = 2^w for a word length of w bits).
%
%   FEASIBLE is a vector of faithful row indices (Walsh scores), as found
%   by find(all(faithmat,1)). FINALIDX is an NSELECT-by-1 column of
%   selected entries of FEASIBLE, sorted ascending. If the algorithm runs
%   out of room, remaining entries repeat the previous pick, so the
%   minimum distance of the codebook will be 0.
%
%
%   Example:
%       feasible = find(all(vertcat(faithmat{:}) == 1, 1));
%       finalIdx = select_greedy_first(feasible, 16);
%
%   See also SELECT_GREEDY_LAST, SELECT_GREEDY_ROWS, SELECT_QUANTILE_ROWS.

    if nargin < 2 || isempty(nSelect)
        nSelect = 16;
    end

    feasible = sort(feasible(:));
    if isempty(feasible)
        error('select_greedy_first:emptyFeasible', 'FEASIBLE must not be empty.');
    end

    greedy = zeros(nSelect, 1);
    % The largest Walsh score anchors the codebook.
    greedy(nSelect) = feasible(end);

    for j = nSelect:-1:2
        % Optimum spacing for what is left, using the basin of attraction
        % distance.
        optdis = greedy(j) / (j - 1/2);

        % Keep only entries below the latest pick at an odd distance.
        testset = feasible - greedy(j);
        oddfeas = feasible((mod(testset, 2) == 1) & (testset < 0));

        % Highest entry that is optdis away or more.
        testfind = oddfeas(find(oddfeas <= (greedy(j) - optdis), 1, 'last'));

        % If nothing qualifies, shrink the distance by 1 and try again.
        while isempty(testfind)
            optdis = optdis - 1;
            testfind = oddfeas(find(oddfeas <= (greedy(j) - optdis), 1, 'last'));
            % Out of room before reaching the first entry: repeat the previous
            % pick, which makes the minimum distance 0.
            if optdis < 0
                testfind = greedy(j);
            end
        end
        greedy(j-1) = testfind;
    end

    finalIdx = sort(greedy);
end
