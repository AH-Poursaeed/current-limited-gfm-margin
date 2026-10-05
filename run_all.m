% Runs all studies and then compares the results with the numbers printed
% in the paper. MATPOWER has to be on the MATLAB path. Start it from any
% folder; every script writes its output next to itself, in ieee14/ or
% smib/. The fine mismatch sweep takes most of the time. Several scripts
% begin with clear and close all, so the workspace is emptied.

if exist('runpf', 'file') ~= 2
    error('run_all: MATPOWER is not on the path (runpf not found).');
end
if exist('fsolve', 'file') ~= 2
    error('run_all: the Optimization Toolbox is needed (fsolve not found).');
end
fprintf('MATLAB %s\n', version);

% modified IEEE 14-bus system, Sections IV.B and IV.C
cd(fullfile(fileparts(mfilename('fullpath')), 'ieee14'));
step1
step2
step3
step4
step5
studyA_validation_residuals
studyA_step4_timed
studyB_run
fig7_pdelta_stages
fig8_energy_barrier

% single converter and two-converter hub, Sections IV.A and IV.C
cd(fullfile('..', 'smib'));
ias_demo_ccpac_cases
studyD_run
studyE_frame_check
studyE_resistive_sweep
studyC_run
studyF_tolerance_sweep
studyF_emf_check

cd('..');
nfail = check_paper_values;
if nfail > 0
    error('run_all: %d checks failed', nfail);
end
