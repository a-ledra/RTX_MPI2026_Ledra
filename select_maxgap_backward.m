function [finalIdx, bestGap] = select_maxgap_backward(feasible, nSelect)
%SELECT_MAXGAP_BACKWARD Pick codebook rows with Lawan's Backward Max-Gap algorithm.
%   Chooses NSELECT faithful Walsh rows by maximizing the minimum sequency
%   gap between chosen rows. A binary search finds the largest integer gap
%   for which NSELECT rows can be picked; each tested gap is checked by
%   starting at the highest faithful sequency and moving backward.
%
%   FINALIDX = SELECT_MAXGAP_BACKWARD(FEASIBLE) selects 16 rows.
%
%   [FINALIDX, BESTGAP] = SELECT_MAXGAP_BACKWARD(FEASIBLE, NSELECT) selects
%   NSELECT rows (NSELECT = 2^w for a word length of w bits) and also
%   returns the largest minimum gap found.
%
%   FEASIBLE is a vector of faithful row indices (Walsh scores), as found
%   by find(all(vertcat(faithmat{:}) == 1, 1)). FINALIDX is an NSELECT-by-1
%   column of selected entries of FEASIBLE, sorted ascending. An error is
%   thrown if FEASIBLE has fewer than NSELECT entries.
%
%
%   Example:
%       feasible = find(all(vertcat(faithmat{:}) == 1, 1));
%       finalIdx = select_maxgap_backward(feasible, 16);
%
%   See also SELECT_MAXGAP_FORWARD, SELECT_GREEDY_ROWS, SELECT_QUANTILE_ROWS.

    if nargin < 2 || isempty(nSelect)
        nSelect = 16;
    end

    feasible = feasible(:);
    if length(feasible) < nSelect
        error('select_maxgap_backward:notEnoughRows', ...
            'Not enough faithful rows. Need %d rows, but only %d are faithful.', ...
            nSelect, length(feasible));
    end

    % In Walsh order, row index minus 1 is the sequency/sign-change count.
    faithfulFreq = feasible - 1;
    [freqSorted, order] = sort(faithfulFreq);
    idxSorted = feasible(order);

    % Binary-search bounds for the integer minimum gap.
    low = 0;
    high = max(freqSorted) - min(freqSorted);
    bestGap = 0;
    bestSelectedIdx = [];

    while low <= high
        gap = floor((low + high) / 2);

        chosenIdx = greedyChooseWithGapBackward(idxSorted, freqSorted, nSelect, gap);

        if length(chosenIdx) >= nSelect
            bestGap = gap;
            bestSelectedIdx = chosenIdx(1:nSelect);
            low = gap + 1;
        else
            high = gap - 1;
        end
    end

    % Sort the selected rows in increasing sequency order.
    finalIdx = sort(bestSelectedIdx(:));
end

function chosenIdx = greedyChooseWithGapBackward(idxSorted, freqSorted, K, gap)
% Test one gap by starting at the highest sequency and moving backward.

    chosenIdx = idxSorted(end);
    lastFreq = freqSorted(end);

    for j = length(idxSorted)-1:-1:1
        if lastFreq - freqSorted(j) >= gap
            chosenIdx(end+1) = idxSorted(j); %#ok<AGROW>
            lastFreq = freqSorted(j);

            if length(chosenIdx) == K
                return;
            end
        end
    end
end
