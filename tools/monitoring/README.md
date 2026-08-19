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
