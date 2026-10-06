function codebook = select_codebook(feasible, nSelect, method)
%SELECT_CODEBOOK Choose a codebook from the faithful set with any of six algorithms.
%   CODEBOOK = SELECT_CODEBOOK(FEASIBLE, NSELECT, METHOD) picks NSELECT rows
%   from the vector FEASIBLE of faithful row indices (Walsh scores) and
%   returns them as a sorted NSELECT-by-1 column.
%
%   METHOD is one of:
%       'greedy_first'      Dr. Edwards, Greedy First
%       'greedy_last'       Dr. Edwards, Greedy Last
%       'maxgap_backward'   Lawan, Backward Max-Gap
%       'maxgap_forward'    Lawan, Forward Max-Gap
%       'greedy_rows'       Adrianna, greedy farthest-point sampling
%       'quantile_rows'     Adrianna, WalshScore quantile spacing
%
%   Returns NaN(NSELECT,1) if there are fewer than NSELECT faithful rows.
%
%   See also SELECT_GREEDY_FIRST, SELECT_GREEDY_LAST, SELECT_MAXGAP_BACKWARD,
%   SELECT_MAXGAP_FORWARD, SELECT_GREEDY_ROWS, SELECT_QUANTILE_ROWS.

    feasible = sort(feasible(:));

    if length(feasible) < nSelect
        warning('select_codebook:notEnoughRows', ...
            'Only %d feasible rows found; need %d.', length(feasible), nSelect);
        codebook = NaN(nSelect, 1);
        return;
    end

    switch lower(method)
        case 'greedy_first'
            codebook = select_greedy_first(feasible, nSelect);
        case 'greedy_last'
            codebook = select_greedy_last(feasible, nSelect);
        case 'maxgap_backward'
            codebook = select_maxgap_backward(feasible, nSelect);
        case 'maxgap_forward'
            codebook = select_maxgap_forward(feasible, nSelect);
        case {'greedy_rows', 'quantile_rows'}
            % These two work on a table. Every faithful row is equally
            % good, so a constant rank column leaves the choice to the
            % WalshScore spacing alone.
            T = table(feasible, feasible - 1, ones(size(feasible)), ...
                'VariableNames', {'OriginalRow', 'WalshScore', 'MeanAccuracy'});
            if strcmpi(method, 'greedy_rows')
                selTbl = select_greedy_rows(T, nSelect, height(T), 'MeanAccuracy');
            else
                selTbl = select_quantile_rows(T, nSelect, height(T), 'MeanAccuracy');
            end
            codebook = sort(selTbl.OriginalRow(:));
        otherwise
            error('select_codebook:unknownMethod', ...
                ['Unknown method ''%s''. Use greedy_first, greedy_last, ' ...
                 'maxgap_backward, maxgap_forward, greedy_rows or quantile_rows.'], method);
    end
end
