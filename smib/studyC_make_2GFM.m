function P = studyC_make_2GFM(varargin)
% Parameters of the two-converter hub as name/value pairs:
%
%   GFM1 --jXa--+
%               +-- hub --jXc-- Vth
%   GFM2 --jXb--+
%
% The defaults split case 1 of Table I into two equal units (H = 2.5 s,
% D = 5 pu, Pref = 0.25 pu, Imax = 0.6 pu each, Xc = 0.7 pu) placed at the
% hub. During the fault Vth is scaled by kV and Xc by kX; the post-fault
% network is the pre-fault one.

ip = inputParser;
ip.addParameter('H1', 2.5);
ip.addParameter('H2', 2.5);
ip.addParameter('D1', 5);
ip.addParameter('D2', 5);
ip.addParameter('Pref1', 0.25);
ip.addParameter('Pref2', 0.25);
ip.addParameter('Qref1', 0);
ip.addParameter('Qref2', 0);
ip.addParameter('Imax1', 0.6);
ip.addParameter('Imax2', 0.6);
ip.addParameter('E1', 1.0);
ip.addParameter('E2', 1.0);
ip.addParameter('Vth', 1.0);
ip.addParameter('theta_th', 0.0);
ip.addParameter('Xa', 0.0);
ip.addParameter('Xb', 0.0);
ip.addParameter('Xc', 0.7);
ip.addParameter('kV', 0.2);
ip.addParameter('kX', 1.0);
ip.addParameter('f0', 60);
ip.parse(varargin{:});
P = ip.Results;
P.ws = 2*pi*P.f0;

P.Vth_pre = P.Vth;
P.Xc_pre  = P.Xc;
P.Vth_fault = P.kV * P.Vth;
P.Xc_fault  = P.kX * P.Xc;
P.Vth_post = P.Vth;
P.Xc_post  = P.Xc;

end
