function [wr_u, raw_u] = dedup_wrapped(wr, raw, tol)
% Sorts the wrapped angles and drops those within tol of the previous one;
% raw carries the unwrapped values along.
    [wr_s, idx] = sort(wr(:));
    raw_s = raw(:); raw_s = raw_s(idx);

    wr_u = wr_s(1);
    raw_u = raw_s(1);

    for k = 2:numel(wr_s)
        if abs(wr_s(k) - wr_u(end)) > tol
            wr_u(end+1,1) = wr_s(k);
            raw_u(end+1,1) = raw_s(k);
        end
    end
end
