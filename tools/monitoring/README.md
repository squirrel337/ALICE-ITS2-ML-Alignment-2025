# ITS2 alignment monitoring

Everything here reads the **`TrkVtxer`** tree, written by `TrackVertexQualityEstimator` from
`TrackerFit` — the fit that uses the clusters alone. The other tree, `ResMonitor`, comes from the
refit that carries the primary vertex as an eighth measured point to drive the alignment gradient;
its residuals are pulled off the clusters by that constraint and its "DCA" is the residual of a
constrained point, so it is not used for physics. `compare_track_fits.C` reads both, on purpose,
to show the difference.

`s2` is the transverse residual and `s1` the longitudinal one: at the vertex point the module
writes `Residual_s1 = proj_GZc - meas_GZc` and `Residual_s2 = ±√(ΔGX² + ΔGY²)`. So s2 pairs with
`DCA_y` and s1 with `DCA_z`. Transverse is plotted first throughout.

Not to be confused with the repository's own `monitor/` directory, which holds the original
post-processing tooling and is untouched by any of this.

## Running

Everything runs from the repository root.

    ./tools/monitoring/make_monitoring_plots.sh -1     # baseline
    ./tools/monitoring/make_monitoring_plots.sh  4     # after five epochs

    root -l -b -q tools/monitoring/plot_epoch_trends.C
    root -l -b -q tools/monitoring/compare_track_fits.C

The ROOT macros take `(residual-dir, out-dir)`, defaulting to `MLPTrain_Step901/Residual` and
`tools/monitoring/plots`. The shell driver takes `<epoch> [outdir] [residual-dir]` and finds the
repository from its own location, so it works from anywhere. Plot directories are gitignored.

## What each file does

| file | produces |
| --- | --- |
| `make_monitoring_plots.sh` | runs the supplied `user_macros/` against one epoch → `plots/` |
| `plot_epoch_trends.C` | residual and DCA against epoch, all epochs on one axis |
| `compare_track_fits.C` | cluster-only vs vertex-constrained fit, side by side |
| `compare_dca_reference.C` | DCA measured against `v_est` vs against `v_reco` (diagnostic) |
| `quicklook.C` | first-pass cost / DCA / vertex / residual glance |
| `check_determinism.sh` | runs one composed job twice and shows where the two logs first differ |
| `run_to_run.sh` | runs one composed job N times from N pristine copies → `costs.tsv` |
| `plot_run_to_run.py` | turns one or two `costs.tsv` into statistics and a figure |
| `weight_spread.py` | per-parameter spread of the epoch-0 weights across repetitions |
| `compare_weight_files.py` | two weight files, per column / per layer / per sensor |

## Measuring run-to-run spread

The module is not bit-reproducible: two byte-identical job trees, same machine, same
binary, give different numbers. Any comparison between two versions of the module has
to clear that band first, so the band has to be measured.

    ./config/runctl.sh set GEOM_BACKEND=cache JOB_TAG=rt-cache JOB_NDATA=4000 JOB_NEPOCH=1
    ./config/runctl.sh compose
    ./tools/monitoring/run_to_run.sh runs/rt-cache cache 30 1

    # on a machine with O2 loaded, the same thing with the other backend
    ./config/runctl.sh set GEOM_BACKEND=o2 JOB_TAG=rt-o2 JOB_NDATA=4000 JOB_NEPOCH=1
    ./config/runctl.sh compose
    ./tools/monitoring/run_to_run.sh runs/rt-o2 o2 30 1

    ./tools/monitoring/plot_run_to_run.py spread.png \
        o2=runs/rt-o2.rtr-o2/costs.tsv cache=runs/rt-cache.rtr-cache/costs.tsv

    ./tools/monitoring/weight_spread.py weights.png \
        o2=runs/rt-o2.rtr-o2/costs.tsv cache=runs/rt-cache.rtr-cache/costs.tsv

`plot_run_to_run.py` compares the cost, one scalar per run. `weight_spread.py` reads
the `weights_Epoch_At_0.txt` that each repetition wrote — the fifth column of
`costs.tsv` — and measures how far each of the 24120 x 17 sensor parameters moves
between repetitions of the same configuration. That band is what any comparison
between two versions of the module has to beat to mean anything.

Three things that are easy to get wrong:

- **`nEPOCH` must be at least 1.** The epoch loop is `for (iepoch = 0; iepoch < nEpoch; ...)`
  at `YMultiLayerPerceptron.cxx:1748`, so `nEPOCH=0` never enters it — the job builds the
  network, prints no cost, writes no weights and exits. `nEPOCH=1` runs epoch 0 and writes
  `weights_Epoch_At_0.txt`, which is the parameter set to compare.
- **Every repetition needs its own copy.** A run consumes its directory: `run_train_circle.C`
  moves `UpdateSensorsList.txt` and `TrendingNetwork/` into `MLPTrain/`. Re-running in place
  is a different experiment. `run_to_run.sh` handles this; a hand-rolled loop usually does not.
- **One run at a time unless you have checked the memory budget.** A run cling-JITs the whole
  module and then holds the event sample: 7.8 GB resident at `nDATA=4000`. Two do not fit under
  a 13 GB cgroup. An OOM-killed run does not announce itself — it stops at exactly the line an
  `nEPOCH=0` run stops at, with no error, no cost and no weights, in about a third of the time.
  Check `/sys/fs/cgroup/memory/.../memory.limit_in_bytes` against `ps -o rss` before raising it.
- **The seed directory must survive.** `run_train_circle.C` reads `SetPrevUSL`, `SetPrevWeight`
  and `SetPrevWeightDU` from `MLPTrain_Step<FIRST_STEP-1>`. `run_to_run.sh` clears only
  `MLPTrain/` and the output steps at or above `FIRST_STEP`, and refuses to start if the seed
  is missing.

## Plots

From `make_monitoring_plots.sh`, into `plots/`, suffixed `_ep<N>`:

| file | content |
| --- | --- |
| `dca_vs_pT`, `dca_vs_eta`, `dca_vs_z`, `dca_vs_phi` | `DCA_y` and `DCA_z`, width and mean |
| `residualXY_vs_phi`, `residualXY_vs_z` | transverse residual per layer, HB0/HB1 split |
| `residualZ_vs_phi`, `residualZ_vs_z` | longitudinal residual per layer, HB0/HB1 split |

From `plot_epoch_trends.C`:

| file | content |
| --- | --- |
| `residual_trend_vs_epoch` | width and mean against epoch, per layer, transverse then longitudinal |
| `residual_XY_by_epoch`, `residual_Z_by_epoch` | distributions per layer, every epoch overlaid |
| `dca_trend_vs_epoch` | DCA mean and width against epoch |
| `dca_width_vs_pT_by_epoch` | DCA width against pT, one curve per epoch |

From `compare_track_fits.C`: `residual_XY_distributions`, `residual_Z_distributions`,
`dca_clusterfit_vs_{pT,eta,zvtx,phi}`.

## Supplied macros

`user_macros/` holds the macros as received, unmodified. They address the tree by the global names
`TrkVtxer` and `ResMonitor`; the driver points both at `TrkVtxer`.
`check_vertex_pT_plots_trkvtxer_color.C` is the one addition — `check_vertex_p_default_...` with
the draw variable changed from `p` to `pT` and the axis titles relabelled.

Two things the driver does rather than editing those files: it gives the HB0/HB1 overlays separate
colours and a key, and it blanks bins whose parent projection holds fewer than 200 tracks — the
fitter runs on every x bin including empty ones past the detector edge, and those come back with
error bars larger than the pad.

## Report builders

`build_report_monitoring.py`, `build_report_epoch_trends.py`, `build_report_cluster_fit.py` and
`build_report_dca_reference.py` embed the PNGs into standalone `report_*.html` pages. They expect
to run from this directory with the plots already produced, and they read the run logs that the
macros print, so redirect those to `<macro-name>.log` first.
