# Run-to-run reproducibility, cache backend

The module used not to be reproducible: two byte-identical job trees, same
machine, same binary, gave different alignment parameters. The cause was an
out-of-bounds read in `BetaLinearization`, fixed in `87b8844` (PR #5). With
that fix repeated runs are **bit-identical**.

Both halves below were measured on the session container, cache backend,
`nEPOCH=1`, `DATA_FILE=XXXXinput.root`, `PARAMS_ARCHIVE=MLPTrain_Step900`,
each repetition from its own pristine copy of one composed job.

## Before and after, at nDATA=4000

The "before" column is the original ten-repetition measurement; "after" is
five repetitions of the same configuration with the fix.

| | before (n=10) | after (n=5) |
| --- | --- | --- |
| `weights_Epoch_At_0.txt` MD5 | 10 distinct | **1 distinct** |
| parameters differing between repetitions | 151600 / 410040 (37.0 %) | **0 / 410040** |
| RMS difference between two runs, median over pairs | 5.131e-09 | **0** |
| per-run parameter RMS, sd | 1.902e-10 (1.28e-03 relative) | **0** |
| worst single parameter, sd | 5.248e-07 | **0** |

Parameter scale for reference: RMS over all parameters is 1.4902e-07 before
and 1.4938e-07 after, so the old spread was **3.4 % of the parameter scale**,
with 87 sensors scattering by more than their own mean value.

Confirmed independently at `nDATA=300`, four repetitions each: the unfixed
build produced four distinct MD5s and 50257 / 410040 differing parameters,
the fixed build four identical MD5s and zero.

## What the cause was

`BetaLinearization` walked `nLAYER+1` entries unconditionally. Seven of its
nine call sites pass that many, because the cost path appends the primary
vertex as an eighth point. `TrackerFit` fits the clusters alone and passes
`nLAYER`, with `beta` and `dirXc` sized to match, so `hitUpdate[nLAYER]` read
past the end of the vector. `std::vector<bool>` packs bits into a word, so
that read returned an uninitialised bit rather than faulting, and its value
changed between runs. When it read true, `beta[nLAYER]` -- one past the end
of the caller's array -- became the anchor the whole chain was aligned to,
shifting the array by 2π or 4π.

Instrumented over 400 tracks, the shift distribution was 0: 326, ±2π: 30,
±4π: 44, and the bit read true in 8/10 calls in one run and 3/10 in the next.

The shift is absorbed everywhere it is consumed -- `circle3Dfit_Z` refits from
the shifted beta so `fparZR` carries the same offset (z unchanged to
7.7e-06 cm), and the ±4π candidate search in `TrackVertexQualityEstimator`
re-matches the vertex, with zero failures over 400 tracks. What does not
cancel is conditioning: `DetATA = nhits*beta2sum - betasum²` is a difference
of large numbers, and with beta near 11.4 instead of -1.1 the fitted slope
moved by up to 3e-7 relative, differently from run to run. The
`QUALITY_VERTEXING` and `Num_Of_Bad_Tracks > 2` gates turned that into
dropped events and a different epoch.

That 1e-7 perturbation looked negligible in isolation. It was not, because
the gates are discrete: measuring the end-to-end effect is what settled it,
not reasoning about the size of the intermediate.

## Reproducing

    ./config/runctl.sh set GEOM_BACKEND=cache JOB_TAG=repro JOB_NDATA=4000 \
        JOB_NEPOCH=1 DATA_FILE=XXXXinput.root PARAMS_ARCHIVE=MLPTrain_Step900
    ./config/runctl.sh compose
    ./tools/monitoring/run_to_run.sh runs/repro repro 5 1

    md5sum runs/repro.rtr-repro/r*/MLPTrain_Step901/weights/weights_Epoch_At_0.txt
    ./tools/monitoring/weight_spread.py out.png repro=runs/repro.rtr-repro/costs.tsv

One repetition is about 950 s at `nDATA=4000` and holds 7.8 GB resident, so
run them one at a time unless the memory budget has been checked.

Bit-identity is now the acceptance criterion: any difference between two runs
of one configuration is a regression, not noise. There is no longer a
tolerance band to clear when comparing two versions of the module.
