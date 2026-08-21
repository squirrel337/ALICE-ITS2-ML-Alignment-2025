# ALICE ITS2 — ML-based Alignment (2025)

Self-supervised neural alignment of the ALICE Inner Tracking System 2, with **adaptive
vertex estimation**.

Each of the **24,120 ALPIDE sensors** carries its own small correction network. There is
no external alignment reference: the training signal is the track residual itself — a
circle fit in *r*–φ together with a *z*–β line fit — so the detector is aligned against
its own tracks.

Per sensor the module fits **11 network parameters** (5 neuron biases and 6 synapse
weights), from which **6 rigid-body degrees of freedom** are derived: three rotations
`Ralpha, Rbeta, Rgamma` and three translations `T1, T2, T3`.

| | |
|---|---|
| Layers | 7 — chip boundaries `0, 108, 252, 432, 3120, 6480, 14712, 24120` |
| Staves per layer | `12, 16, 20, 24, 30, 42, 48` |
| Sensors | 24,120 |
| Runtime dependency | **ROOT only** — see below |

---

## What is new since 2024

**Adaptive vertex estimation.** The 2024 module trains against whatever vertex the
reconstruction supplied. This one re-estimates it from the prongs of each event, and can
reject the event on what it finds. `UpdateVertexByAlignment` fits the `nTrackMax` prongs
pairwise, scores each as `dev = p × |dca| / 40 µm` in *y* and *z* separately, drops prongs
above `QUALITY_VERTEXING`, refits from what remains, and discards the event if too many
prongs were bad or the worst remaining `dev` exceeds `QUALITY_TRACKVERTEX`.

Those thresholds are momentum-weighted, so they are not windows: at `QUALITY_VERTEXING = 5`
a 1 GeV/c prong may miss by 200 µm and a 10 GeV/c prong by 20 µm. The console's Vertex tab
translates whatever you type into microns at both ends.

**The layer-0 defect is fixed.** In the 2024 tree a track with no layer-0 hit read `vxyz`
from three stack slots that were never written; 47 % of that sample was selected by
leftover stack contents rather than by physics. Here the track is projected to the vertex β
and the innermost layer that actually has a hit is searched for.

Full description: **[`docs/workflow.html`](docs/workflow.html)**.

---

## Two geometry backends

The module needed O<sup>2</sup> at *run* time, not just to build: `batch_train` loads
`YDetectorGeometry.cxx` through cling, so O<sup>2</sup> headers and libraries had to
resolve while the job ran. They no longer have to.

| `GEOM_BACKEND` | Reads geometry via | Needs O<sup>2</sup> at run time |
|---|---|---|
| **`o2`** *(default)* | `o2::its::GeometryTGeo`, as the module always has | yes — load it with `alienv` first |
| `cache` | a per-chip file built once by `tools/export_geometry_cache.C` | no |

The public interface is identical either way, so none of the 350+ `yGEOM->` call sites
change. Build the cache once:

```sh
# with O2 loaded, this is all you need -- O2 supplies the AlignParam dictionary
root -l -b -q 'tools/export_geometry_cache.C'

# without O2, build a dictionary first from the alignment file's own StreamerInfo
root -l -b -q 'tools/make_alignlib.C'
root -l -b -q 'tools/export_geometry_cache.C'
```

The cache records a fingerprint of the alignment it was built from, and `doctor` refuses a
cache that does not match the `ALIGN_FILE` the job will run with.

---

## Running a job

Everything one run needs lives in `config/runconsole.conf`. Drive it from a terminal or
from a ROOT window — the window shells out to the same script, so anything it can do
works over ssh with no display.

```sh
eval `alienv load -w $O2_DIR/sw O2/latest`   # default backend is o2, so load it

./config/runctl.sh ui         # set everything in one window
./config/runctl.sh doctor     # check this machine has what the run needs
./config/runctl.sh compose    # build the job directory
./config/runctl.sh run        # launch it
./config/runctl.sh log -f     # watch
```

`compose` builds a self-contained job directory and patches the knobs into **that job's**
copies of the headers and sources. The module checkout is only ever read.

Full reference: **[`config/README.md`](config/README.md)** and
**[`docs/run-console.html`](docs/run-console.html)**.

### Current configuration

Every default in `config/runconsole.conf` is what this module actually ships, so composing
without editing anything reproduces the tree as checked out. `doctor` warns on any knob
where the configuration and the module disagree.

| Knob | Value | File |
|---|---|---|
| `nDATA` | 50000 | `YMLPParallel.h` |
| `nEPOCH` | 5 | `YMLPParallel.h` |
| `nTrackMax` | 8 | `Ymlp/inc/DetectorConstant.h` |
| `DET_MAG` | −0.500673 T | `Ymlp/inc/DetectorConstant.h` |
| `FITMODEL` | 2 (circle) | `Ymlp/inc/DetectorConstant.h` |
| `QUALITY_VERTEXING` | 5 | `Ymlp/src/YMultiLayerPerceptron.cxx` |
| `QUALITY_TRACKVERTEX` | 10 | `Ymlp/src/YMultiLayerPerceptron.cxx` |
| `GEOM_BACKEND` | `o2` | `config/runconsole.conf` |

---

## Layout

| Path | What |
|---|---|
| `Ymlp/` | The module. `YMultiLayerPerceptron.cxx` is the core (~15.1k lines), `YAlignment.cxx` the driver |
| `run_train_circle.C` | Entry point — `batch_train` loads the geometry, then runs this for one step |
| `config/` | Run console: one configuration file, a library, a dispatcher |
| `tools/ConfigUI/` | The ROOT window |
| `tools/` | Geometry cache export, the O<sup>2</sup> comparison probes, the alignment fingerprint |
| `tools/monitoring/` | Residual and DCA monitoring, read from the cluster-only track fit |
| `monitor/` | The original post-training macros |
| `docs/` | [Workflow](docs/workflow.html) · [Run console](docs/run-console.html) |
| `NetworkParameters/` | Deployed weights |

**There is no compiled build.** `Ymlp/CMakeLists.txt` is 534 null bytes; the module is
cling-JIT only, loaded from source on every run. That is also what lets the console patch
a job's `.cxx` without rebuilding anything.

---

## Which fit a number came from

The module fits every track **twice**, and the two are not interchangeable.

| Fit | Tree | Use |
|---|---|---|
| `TrackerFit` | `TrkVtxer` | clusters alone — the physical trajectory |
| `GetCost_Beam_CircleFit` | `ResMonitor` | clusters **+ the primary vertex** as an eighth point, to drive the gradient |

Impact parameters and residuals belong to the first. In the second the vertex is part of
the fit, so its "DCA" is the residual of a constrained point pinned to
`Sigma_MEAS[7] = 4.74 µm`: measured that way the transverse width collapses to a flat 3 µm,
where the cluster-only fit gives 82 µm at 0.3–0.5 GeV/c falling to 14 µm above 4 GeV/c —
the expected multiple-scattering behaviour. Everything under `tools/monitoring/` reads
`TrkVtxer`; `compare_track_fits.C` reads both, deliberately, to show the difference.

`s2` is the transverse residual and `s1` the longitudinal one, from the vertex-point
definitions `Residual_s1 = proj_GZc − meas_GZc` and `Residual_s2 = ±√(ΔGX² + ΔGY²)`. So s2
pairs with `DCA_y`, s1 with `DCA_z`.

---

## Known limitations

These are measured, not hypothetical, and any number produced from this tree carries them.

**`Angle2Alpha` and `kB2C` are reconstructed.** `Ymlp/inc/YO2Compat.h` stands in for a
handful of O<sup>2</sup> header-only helpers, and these two were rebuilt from the
conventions rather than copied from O<sup>2</sup> source. In this module `kB2C` reaches
further than it did in 2024: besides the impact parameter it sets the seed radius
`R = pT/q/(kB2C·B)` that `TrackerFit` fixes rather than fits, which feeds
`UpdateVertexByAlignment` and the track–vertex quality gate. A wrong value moves the
estimated vertex and changes which events survive, not merely which impact parameter is
reported. Seven call sites depend on it. Both are isolated in that one header, so a
correction is a one-line change. **Verify before a training campaign.**

**The reconstructed vertex and the χ² vertex disagree by 1.81 mm in z.** An independent
check from the stored cluster positions alone — least squares on *z* against *r*,
extrapolated to *r* = 0, touching none of the module's fitting — puts it at 1.783 mm, so
the clusters agree with the χ² vertex rather than with `tv3`. That looks like a question
about the input, not a module defect, and it caps what the alignment can converge to.

**`VERTEXFIT` is inert.** `DetectorConstant.h` defines it and nothing in `Ymlp/` reads it.
The console's Job tab carries it because the 2024 console does.

**Geometry equivalence is settled; physics equivalence is not.** The cache reproduces
O<sup>2</sup>'s `getMatrixL2G` for all 24,120 chips with zero difference, and the composed
transforms — `LToG`, `GToS`, `NormalVector` over five pixel corners per chip — agree to
0.038 µm, inside one float ULP. What has *not* been shown is that a full training run gives
the same cost and the same residuals under both backends. Run one configuration under each
`GEOM_BACKEND` and compare; that is what the option exists for.

**Fitted residual widths are stable to about 0.5 %, not to the decimals they print.** The
distributions are heavy-tailed and the iterative ±2σ Gaussian fit moves under small
perturbations — a 0.038 µm geometry change once shifted a fitted L3 width by 2.4 µm while
the RMS and IQR of the same distribution moved 0.4 %. Prefer the robust statistics when the
difference matters.

---

## License

Apache 2.0 — see [LICENSE](LICENSE).
