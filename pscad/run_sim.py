"""Wall-clock time of the EMT run (Table XI, PSCAD/EMTDC row).

Attaches to a running PSCAD with the IEEE 14-bus case loaded, runs it once
cold (the first run also compiles the Fortran) and then N_WARM times, and
writes the times to a csv. Needs the PSCAD automation library (mhi.pscad).
"""
import csv
import platform
import time

import mhi.pscad

PROJECT = "ieee_14_bus_IAS"
N_WARM = 5
OUT_CSV = "studyA_pscad_runtime.csv"

pscad = mhi.pscad.application()
project = pscad.project(PROJECT)

# PSCAD gives the duration in s and the solution time step in microseconds
sim_duration = sim_step = None
try:
    p = project.parameters()
    sim_duration = p.get("time_duration")
    sim_step = p.get("time_step")
except Exception as e:
    print(f"could not read the simulation parameters: {e}")
print(f"project {PROJECT}, duration {sim_duration} s, time step {sim_step} us")

t0 = time.perf_counter()
project.run()
t_cold = time.perf_counter() - t0
print(f"cold run: {t_cold:.2f} s")

warm = []
for i in range(N_WARM):
    t0 = time.perf_counter()
    project.run()
    dt = time.perf_counter() - t0
    warm.append(dt)
    print(f"warm run {i + 1}/{N_WARM}: {dt:.2f} s")

mean = sum(warm) / len(warm)
mn = min(warm)
# population standard deviation (divide by N)
std = (sum((x - mean) ** 2 for x in warm) / len(warm)) ** 0.5
compile_est = t_cold - mean

print(f"warm: mean {mean:.3f} s, min {mn:.3f} s, std {std:.3f} s (n = {N_WARM})")
print(f"cold: {t_cold:.3f} s, compile about {compile_est:.2f} s")

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["field", "value"])
    w.writerow(["machine", platform.platform()])
    w.writerow(["processor", platform.processor() or "n/a"])
    w.writerow(["python", platform.python_version()])
    w.writerow(["project", PROJECT])
    w.writerow(["sim_duration_s", sim_duration])
    w.writerow(["sim_time_step", sim_step])
    w.writerow(["cold_run_s", round(t_cold, 3)])
    for i, dt in enumerate(warm):
        w.writerow([f"warm_run_{i + 1}_s", round(dt, 3)])
    w.writerow(["warm_mean_s", round(mean, 3)])
    w.writerow(["warm_min_s", round(mn, 3)])
    w.writerow(["warm_std_s", round(std, 3)])
    w.writerow(["compile_estimate_s", round(compile_est, 3)])

print(f"written {OUT_CSV}")
