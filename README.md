# Transient stability margin and CCT for current-limited grid-forming converter clusters

MATLAB code for the numerical results of

> A. H. Poursaeed, F. Namdari, S. Shen and P. A. Crossley, "Current-Limited
> Grid-Forming Converter Clusters in Inverter-Dominated Grids: A Transient
> Stability Margin Framework," *IEEE Transactions on Industry Applications*,
> early access, 2026, doi: [10.1109/TIA.2026.3738689](https://doi.org/10.1109/TIA.2026.3738689).

The paper screens the transient stability of grid-forming (GFM) converter
clusters under current limiting. The network seen from each cluster is reduced
to a Thevenin equivalent for the pre-fault, fault-on and post-fault stage. On
these, a current-limited power-angle curve gives the post-fault stable and
unstable equilibria, a normalized margin from their separation, and a critical
clearing time (CCT) from the fault-on trajectory and the post-fault energy
barrier. The method is meant as a fast first screening step before detailed
EMT simulation, which it does not replace.

This repository holds the reduced-order computations: the single-converter
cases, the two-converter aggregation studies and the modified IEEE 14-bus
system.

## Requirements

- MATLAB. We ran the code with R2026a only, on Windows 11 and on Linux
  (University of Exeter ISCA cluster). `exportgraphics` and the append mode of
  `writetable` need R2020a at least; the releases in between were not tried.
- Optimization Toolbox, for `fsolve` in the two-converter studies.
- [MATPOWER](https://matpower.org) 8.1, for the power flow and the admittance
  matrices of the IEEE 14-bus part.
- For `pscad/run_sim.py` only: Python, PSCAD and its automation library
  `mhi.pscad`. The script times the PSCAD runs behind the EMT row of
  Table XI.

## Running everything

Install MATPOWER and put it on the MATLAB path (`install_matpower` in the
MATPOWER folder does this). Then, from the repository folder:

```matlab
run_all
```

This runs all studies and ends with `check_paper_values`, which compares the
results with the numbers printed in the paper and lists one line per check:

```
check                                                         deviation  tolerance
Table I    Delta_ref (deg)                                      0.00445      0.005  ok
Table I    Delta (deg)                                           0.0048      0.005  ok
...
69 checks, 0 failed
```

The deviation is the largest difference from the printed values. The
tolerance is half a unit of the last printed digit, except where the notes
below say otherwise. `run_all` ends with an error if a check fails.

The whole run takes about 20 minutes on one core of a cluster node and about
10 minutes on a laptop (Core i5-1135G7, on mains power). Most of that is the
fine mismatch sweep (`studyF_tolerance_sweep`), which appends one row per case
to `smib/studyF_tolerance_sweep.csv` and continues from there if it is
interrupted. For the same reason a second run does not repeat the sweep;
delete that file to compute it again. No display is needed, so
`matlab -batch run_all` works as well.

Three things to know before starting:

- `run_all` empties the MATLAB workspace and closes all figures, because
  several scripts begin with `clear` and `close all`. It ends in the
  repository folder.
- `data/` is also an input. The fine sweep reads
  `data/smib/studyC_aggregation_error.csv`, the results behind Table XII, to
  compare itself with them, so keep `data/` in place.
- `studyC_run`, `studyD_run`, `studyE_resistive_sweep` and
  `studyF_tolerance_sweep` first compare a few results with values stored in
  the code or in `data/` and stop with an error if they do not agree. For the
  last two the CCTs have to be identical, which means the same bisection
  steps. This was the case on both machines we used.

Every script writes its results into the folder it is in (`smib/` or
`ieee14/`). `data/` holds the csv files of our own run on the cluster and the
two timing logs behind Table XI, and

```matlab
check_paper_values data
```

checks those instead, without running anything.

## Which script produces what

Scripts are run from inside their folder. Within `ieee14/` the order matters:
`step1` to `step5` first, since each step reads what the previous one saved.

| In the paper | Script and the files it writes |
|---|---|
| Table I, Figs. 3-5: margins of the four single-converter cases | `smib/ias_demo_ccpac_cases.m` writes `ias_ccpac_margins.csv`, `fig3_power_angle.png`, `fig4_current.png`, `fig5_margin.png` |
| Table II: energy barrier and CCT of the four cases | `smib/ias_demo_ccpac_cases.m` writes `ias_ccpac_cct.csv` |
| Table III: CCT under circular, d-priority and q-priority limiting, and the sweep mentioned with it in Section IV.A (case 1 against the fault-on voltage factor kV) | `smib/studyD_run.m` writes `studyD_main_sweep.csv`, `studyD_kV_sweep.csv` |
| Table XIV: CCT with the priority split in the Thevenin-voltage frame and in the converter frame | `smib/studyE_frame_check.m` writes `studyE_frame_check.csv` |
| Section IV.A: limiter ordering against the R/X ratio of the coupling | `smib/studyE_resistive_sweep.m` writes `studyE_resistive_sweep.csv` |
| Table XII: aggregation error on the two-converter hub | `smib/studyC_run.m` writes `studyC_aggregation_error.csv` |
| Table XIII and the separation figures of Section II.A: admissible parameter spread within a cluster | `smib/studyF_tolerance_sweep.m` writes `studyF_tolerance_sweep.csv` and calls `smib/studyF_envelope.m`, which writes `studyF_admissibility_envelope.csv`; `smib/studyF_emf_check.m` writes `studyF_emf_check.csv` |
| Table VII: Thevenin equivalents per stage | `ieee14/step1.m`, then `ieee14/step2.m`, which writes `thevenin_equivalents.csv` |
| Table VIII: reduction against the full network | `ieee14/studyA_validation_residuals.m` writes `studyA_validation_residuals.csv` |
| Tables IX and X: equilibria, margin, energy barrier and CCT per cluster | `ieee14/step3.m`, `step4.m`, then `step5.m`, which writes `ieee14_results.csv` |
| Table XI, reduced-order row: time of the CCT computation | `ieee14/studyA_step4_timed.m` writes `studyA_step4_timing.csv` |
| Table XI, EMT row: time of a PSCAD run | `pscad/run_sim.py` writes `studyA_pscad_runtime.csv` |
| Section IV.C: margin and CCT under 22 changes of the external network | `ieee14/studyB_run.m` writes `studyB_robustness_table.csv` |
| Fig. 7: power-angle curves of the bus 2 cluster per stage | `ieee14/fig7_pdelta_stages.m` writes `fig7_pdelta_stages.png` |
| Fig. 8: energy barrier and clearing-time check | `ieee14/fig8_energy_barrier.m` writes `fig8_energy_barrier.png` |

## Layout

```
run_all.m              runs all studies, then check_paper_values
check_paper_values.m   compares results with the numbers in the paper
smib/                  single converter against a Thevenin source, and
                       the two-converter hub
  ias_*.m              power-angle curve with current limiting,
                       equilibria, margin, energy barrier, CCT
  studyC_*.m           two-converter hub and its aggregation
  studyD_*.m           priority current limiting
  studyE_*.m           priority limiting with resistance and with the
                       split taken in the converter frame
  studyF_*.m           fine mismatch sweep and admissible spread
ieee14/                modified IEEE 14-bus system
  step1.m ... step5.m  the computation, one step per script
  studyA_*.m           reduction check and timing
  studyB_*.m           changes of the external network
  fig7_*.m, fig8_*.m   figures
  other .m files       functions shared by the scripts above
pscad/run_sim.py       times the PSCAD runs through its automation
                       library
data/                  csv files of our run, same names as above
```

The parameters are set in `smib/ias_cluster_params_ccpac.m` for the four
single-converter cases, in `smib/studyC_make_2GFM.m` for the two-converter hub
and at the top of `ieee14/step1.m` for the IEEE 14-bus model. The four cases
are those of reference [12] of the paper; `ccpac` in the file names stands for
the current-constrained power-angle characterization of that reference.

## Notes on the numbers

- Table II gives V_crit with five significant digits. For case 4 it prints
  1.24330 where the code gives 1.243293, so the check allows 5e-5 for this
  column.
- The angles of Table I are grid values. The equilibria are taken at the
  nearest point of a grid of 2000 points on [-pi, pi], which has a step of
  0.18 deg.
- The +35.2 % and -30.9 % quoted with Table III are the change of the mean CCT
  over the four cases, relative to the mean CCT with the circular limiter. The
  mean of the column `dCCT_rel_pct` in `studyD_main_sweep.csv`, +37.1 % and
  -32.0 %, is the mean of the four relative changes, which is another number.
- The 64 points with R/X up to 0.25 in `studyE_resistive_sweep.csv` hold the
  four cases at R/X = 0 twice, once under each of the two ways in which the
  resistance is varied. These are the same points, and they are the cases of
  Table III.
- Tables II and III use a search window of [0, 2] s and a tolerance of
  1e-3 s, so their CCTs are multiples of 1/1024 s (0.40625 = 416/1024).
  Table XIV uses [0, 20] s and 1e-4 s, and differs from Table III by up to
  0.8 ms.
- The residuals of Table VIII are at round-off level (1e-16 to 1e-15 pu). We
  obtained the same digits on both machines we tried, but they need not agree
  digit for digit elsewhere; the check only asks for values below 1e-14 pu.
- Tables IX and X print four decimals. Most entries are cut rather than
  rounded (0.722270 is printed as 0.7222), and two V_cr values, 1.0473 and
  1.0663, are rounded. The check therefore allows one unit of the last digit.
- Fig. 8 marks the CCT read from the plotted curve by interpolation,
  0.23389 s. Table X gives the result of the bisection, 0.233869 s.
- Wall-clock times depend on the machine. `data/ieee14/studyA_step4_timing.csv`
  and `data/pscad/studyA_pscad_runtime.csv` are the logs of the runs quoted in
  Table XI and are kept as they are; a new run writes its own times to
  `ieee14/`. The first log has one row per cluster, and the 0.173 s of
  Table XI is its bus 8 row. `sim_time_step` in the second log is in
  microseconds.
- The figures are saved as the scripts draw them. Those in the paper were
  resized and relabelled for the column width, so they look different while
  showing the same curves.

## Notes on the result files

- CCTs are in s, also where the column name has no unit. `CCT_agg_EAC` is the
  CCT of the aggregated cluster by the energy criterion; `CCT_agg_TD` and
  `CCT_per_TD` are the time-domain CCTs of the aggregated cluster and of the
  two converters. `dCCT_rel` is their difference relative to `CCT_agg_TD`, and
  it is negative when the aggregated cluster gives the longer CCT.
- `dCCT_rel`, `Msys_relDev` and `CCT_relDev` are fractions; `err_at_max`,
  `err_at_next` and the columns ending in `_pct` are in per cent. Angles are in
  rad unless the name ends in `_deg`. `rho` is R/X.
- The cluster margin is called `M_CL` or `MhCL`. In `ieee14_results.csv`,
  `Ecrit` is the energy barrier, `Vcr_pu` is the terminal voltage at which the
  current limit is reached, and `alpha_u` equals `Delta_rad`.
- A CCT equal to the end of the search window means that the case was still
  stable there: 2 s in `studyD_kV_sweep.csv` for kV from 0.5, 20 s in
  `studyE_resistive_sweep.csv` (column `censored`), 10 s for bus 3 in
  `ieee14_results.csv` (column `CCT_is_gt`) and for the case `B2_pretrip_4_5`
  in `studyB_robustness_table.csv`, whose `CCT_relDev` has no meaning for that
  reason.
- In `studyC_aggregation_error.csv` and `studyF_tolerance_sweep.csv` some cases
  have no error value. In these no pre-fault operating point was found for the
  two converters, or for their aggregate; the `note` column of the second file
  says which. The point is searched with `fsolve` from one starting point and
  then on a coarse grid (`studyC_initial_sep`), so a missing value means that
  none was found. In three cases of the fine sweep, all at 0.20 pu separation,
  none was found for the aggregate written as two half-size units, although
  the aggregate has a CCT by the energy criterion. The cases without a value
  are left out of every maximum.
- In `studyF_admissibility_envelope.csv` an admissible spread ends at the last
  level before the first one that has no error value or an error above the
  bound. A `max_admissible` of 0 with `err_at_max` 0 means that not even the
  smallest tested level is within the bound.
- Apart from the two timing logs, `data/` was written on Linux. A run on
  Windows gave the same numbers in every csv file except the timings and 23
  values of the column `A_acc_faulton` in `studyE_resistive_sweep.csv`, which
  differ in the 15th digit.

## Citation

```bibtex
@article{poursaeed2026current,
  author  = {Poursaeed, Amir Hossein and Namdari, Farhad and
             Shen, Shuhang and Crossley, Peter A.},
  title   = {Current-Limited Grid-Forming Converter Clusters in
             Inverter-Dominated Grids: A Transient Stability Margin
             Framework},
  journal = {IEEE Transactions on Industry Applications},
  year    = {2026},
  note    = {Early access},
  doi     = {10.1109/TIA.2026.3738689}
}
```

## License

MIT, see [LICENSE](LICENSE). MATPOWER is distributed under its own license.
