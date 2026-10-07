# Requirements and verification (draft)

There is no customer for this project. The requirement values below are my own engineering
assumptions for a generic Earth-observation smallsat. I wrote them after seeing the first
mission runs and a 10-run Monte Carlo, not before (the table now uses a 30-run Monte Carlo), so a pass here shows the design meets the
stated numbers, not that the numbers were fixed blind. They are not taken from a real mission. Every result comes from
the simulation in this repository, under the assumptions listed in the README.

| ID | Requirement (assumed) | Verification | Evidence | Result | Status |
|----|-----------------------|--------------|----------|--------|--------|
| R1 | Detumble from a tip-off of up to 6 deg/s to below 0.5 deg/s within 2 orbits (11 600 s) | Simulation, Monte Carlo | `run_monte_carlo` | 28 of 30 runs pass; runs 24 and 27 take 12 120 s and 12 380 s (worst); mean 8 315 s | **fail in 2 of 30** |
| R2 | Capture nadir pointing (error below 0.1 deg) within 300 s of handover | Simulation, Monte Carlo | `run_monte_carlo` | worst 180 s (30 runs) | pass |
| R3 | Pointing error after capture: mean at most 0.01 deg, maximum at most 0.03 deg | Simulation, Monte Carlo, estimate in the loop | `run_monte_carlo` | mean 0.0048 deg (worst run 0.0061), worst max 0.0176 deg (30 runs, last 2000 s) | pass |
| R4 | Attitude knowledge error at most 30 arcsec rms | Simulation, estimator vs truth | `run_monte_carlo`, `run_estimator_robustness` | 14.3 arcsec rms (worst 17.1, 30 runs); 21 arcsec in the combined stress case | pass |
| R5 | Wheel momentum at most 10 % of capacity (0.2 N m s) after dumping | Simulation | `run_mission`, `run_monte_carlo` | 30-run table: worst 0.84 N m s, but this window (last orbit of a 17 400 s run) contains the capture transient of the two late handovers. Rerun of the 3 slowest runs (9, 24, 27) to 24 000 s: 0.013 to 0.014 N m s. The other 27 runs were not rerun | pass (checked on the 3 slowest runs) |
| R6 | Wheel torque at its limit in at most 1 % of samples after handover | Simulation | `run_mission` | 0.1 % | pass |
| R7 | Torquer dipole within the 30 A m^2 limit at all times | Simulation | `run_monte_carlo` | worst 28.6 A m^2 after handover, margin 5 % (logged every 10 s) | pass, thin margin |
| R8 | Estimator covariance consistent with its error (NEES between 1 and 6) under inertia, torque and drag errors | Analysis, Monte Carlo of the filter alone | `mekf_check`, `run_estimator_robustness` | NEES 2.4 to 3.6 for single errors; not met with unmodelled torque noise (13) or a spun-up body | partial |
| R9 | Simulink model agrees with the independent plain-MATLAB reference within stated limits | Test | `simulink_crosscheck(true)` | 40 of 40 checks | pass |
| R10 | Magnetic momentum dumping law produces the intended torque and never raises wheel momentum | Test | `dump_check` | 6 of 6 checks | pass |

## Notes on the Monte Carlo

- R1 fails for 2 of 30 runs against the limit I wrote down. The 11 600 s value has no mission driver
  (it is two orbits, chosen after seeing 10 runs with a worst of 10 630 s), so I am reporting the failure
  and not moving the limit. Detumble time depends on tumble direction relative to the field as much as
  on rate (run 3: 5.9 deg/s in 6 750 s; run 21: 3.7 deg/s in 10 280 s).
- The 30-run table uses a 17 400 s run (3 orbits) and statistics over the last 2000 s. For late handovers
  that window is not fully settled. The three-run rerun to 24 000 s shows the wheel momentum settles
  to about 0.013 N m s.

## Known gaps

- The Monte Carlo varies tip-off rate, initial attitude, gyro bias and noise seeds. It does not
  vary inertia, actuator error, disturbance size or sensor error models.
- R8 is partial: the filter is overconfident with white torque noise of 2e-4 N m per step and loses
  consistency when a constant torque spins the body to about 1 deg/s without a controller.
- All wheels, sensors and disturbances are idealised: no wheel jitter or friction model, no
  flexible modes, no tracker latency, bias or blinding, white gyro noise only.
- Logged peaks are sampled every 10 s.
