function [delta1, delta2] = studyC_initial_sep(P)
% Pre-fault angles of the two converters, from Pe_k = Pref_k. fsolve is
% started at the angle of the equivalent single cluster; if it fails, the
% best point of a coarse grid is taken instead. Returns NaN when the power
% mismatch stays above 5e-4 pu, i.e. no operating point was found.

Eh_guess = (P.E1 + P.E2)/2;
Pmax_lin = Eh_guess * P.Vth / P.Xc;
Pref_h   = P.Pref1 + P.Pref2;
if Pmax_lin > abs(Pref_h)
    d0 = asin(Pref_h / Pmax_lin);
else
    d0 = pi/4;
end
x0 = [d0; d0];

function r = res(x)
    [Pe1, Pe2] = studyC_per_converter_pe(x(1), x(2), P, 'pre');
    r = [Pe1 - P.Pref1; Pe2 - P.Pref2];
end

opts = optimoptions('fsolve','Display','off','TolFun',1e-12,'TolX',1e-12, ...
    'MaxIterations',200,'MaxFunctionEvaluations',2000);
try
    [x, ~, exitflag] = fsolve(@res, x0, opts);
    if exitflag <= 0
        x = grid_fallback(P);
    end
catch
    x = grid_fallback(P);
end

r_final = res(x);
% a unit whose current limit is below what its reference needs has no
% operating point
if max(abs(r_final)) > 5e-4
    delta1 = NaN;
    delta2 = NaN;
    return;
end

% identical units: remove the last-digit asymmetry left by the solver
if abs(P.H1-P.H2)+abs(P.D1-P.D2)+abs(P.Pref1-P.Pref2)+abs(P.Imax1-P.Imax2)+abs(P.E1-P.E2)+abs(P.Xa-P.Xb) < 1e-12
    sep_avg = 0.5*(x(1)+x(2));
    x = [sep_avg; sep_avg];
end

delta1 = x(1); delta2 = x(2);
end

function x = grid_fallback(P)
n = 121;
g = linspace(-pi, pi, n);
best = inf; bx = [0;0];
for i = 1:n
    for j = 1:n
        [p1, p2] = studyC_per_converter_pe(g(i), g(j), P, 'pre');
        v = (p1-P.Pref1)^2 + (p2-P.Pref2)^2;
        if v < best
            best = v;
            bx = [g(i); g(j)];
        end
    end
end
x = bx;
end
