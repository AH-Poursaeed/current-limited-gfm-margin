function result = ias_system_margin(clusterArray)
% System margin: the smallest of the cluster margins, eq. (32).

H = numel(clusterArray);
M_vec = NaN(H,1);

for h = 1:H
    M_vec(h) = clusterArray(h).M_CL;
end

[M_sys, idx_crit] = min(M_vec);

result.M_sys    = M_sys;
result.M_vec    = M_vec;
result.idx_crit = idx_crit;

end
