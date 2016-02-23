function cbc = piv_correlation_correction(xcr)
%
% Apply correlation-based correction (CBC) to the input correlation matrices.
% See [1] for method description.
%
% [1] Hart, D. P. (2000). PIV error correction. Experiments in Fluids, 29(1),
% 13–22. doi:10.1007/s003480050421

% NOTE: first attempt, correct using 4 neighbors if available. I expect that
% this will introduce edge artefacts since edge points are convolved with
% neighbors on one side only. If so, this could be avoided by requiring
% symmetric corrections, or skipping edge points, or could be accounted for by
% adjusting the point location.

[nr, nc] = size(xcr);
cbc = cell(nr, nc);

for ii = 1:nr
    for jj = 1:nc
        
        if ~isempty(xcr{ii,jj}); 
            % init correlation matrix
            cbc{ii,jj} = xcr{ii,jj};
        else
            % skip points outside the roi
            continue; 
        end
        
        % apply corrections for 4 neighbors, if available
        if ii>1 && ~isempty(xcr{ii-1,jj})
            cbc{ii,jj} = cbc{ii,jj}.*xcr{ii-1,jj};
        end        
        if ii<nr && ~isempty(xcr{ii+1,jj})
            cbc{ii,jj} = cbc{ii,jj}.*xcr{ii+1,jj};
        end        
        if jj>1 && ~isempty(xcr{ii,jj-1})
            cbc{ii,jj} = cbc{ii,jj}.*xcr{ii,jj-1};
        end        
        if jj<nc && ~isempty(xcr{ii,jj+1})
            cbc{ii,jj} = cbc{ii,jj}.*xcr{ii,jj+1};
        end        
        
    end
end