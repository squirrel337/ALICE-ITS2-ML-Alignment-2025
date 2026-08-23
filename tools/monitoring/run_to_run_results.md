# Run-to-run spread, cache backend

Measured 2026-08-23 on the session container. Ten repetitions of one identical
configuration, each from its own pristine copy of the same composed job.

    GEOM_BACKEND=cache  JOB_NDATA=4000  JOB_NEPOCH=1  FIRST_STEP=901
    DATA_FILE=XXXXinput.root  PARAMS_ARCHIVE=MLPTrain_Step900
    3096 train / 1032 test events, 24768 / 8256 tracks
    10/10 repetitions produced weights, ~850 s each, one at a time

Nothing differs between the repetitions. Same binary, same events, same seed,
same machine, no training beyond epoch 0. Everything below is therefore the
module's own non-determinism and nothing else.

## The band

| quantity | value |
| --- | --- |
| epoch-0 fit cost, mean | 0.994572 |
| cost sd | 3.197e-04  (321 ppm relative) |
| cost range over 10 runs | 1.158e-03  (1164 ppm) |
| cost sd, 95% CI | [2.20e-04, 5.84e-04] |
| parameter RMS, mean over runs | 1.4902e-07 |
| **RMS difference between two runs, median over 45 pairs** | **5.131e-09** |
| **as a fraction of the parameter scale** | **3.44 %** |

Two runs of the identical configuration land 3.4 % apart in parameter space.
That is the floor for any comparison between two versions of the module: a
difference smaller than this is not evidence of anything.

## Where the instability sits

It is not spread evenly. Of 410040 parameters (24120 sensors x 17), 162895
(39.7 %) differ between repetitions at all; the rest are identical every run.

| relative spread `sd/|mean|` | parameters | share of the ones that move |
| --- | --- | --- |
| > 1 %   | 42330 | 25.99 % |
| > 10 %  |  8541 |  5.24 % |
| > 50 %  |  2130 |  1.31 % |
| > 100 % |  1043 |  0.64 % |

Median `sd/|mean|` is 3.96e-03, p95 is 1.07e-01.

Per sensor, over the 13118 that move at all: median 5.4e-03, p95 1.3e-01,
worst 1.75e+02 -- that sensor's run-to-run scatter is 175 times its own mean
value. 824 sensors exceed 10 %, 186 exceed 50 %, and 87 have a scatter larger
than the value they are scattering around.

The variance is concentrated: the worst 50 sensors carry 55 % of it, the worst
500 carry 95 %. So most of the detector is reproducible and a small tail is not.

## By column

| col | sensors that move | median sd | p95 sd | max sd | median \|mean\| |
| --- | --- | --- | --- | --- | --- |
| 0, 1, 11, 12 | 0 % | identical in every repetition | | | |
| 2  | 54.1 % | 4.161e-13 | 2.608e-11 | 2.348e-08 | 1.466e-10 |
| 3  | 54.3 % | 1.420e-09 | 2.261e-08 | 5.248e-07 | 2.300e-07 |
| 4  | 54.1 % | 4.840e-11 | 1.546e-09 | 7.984e-08 | 7.108e-09 |
| 5  | 32.3 % | 8.315e-31 | 3.326e-30 | 2.341e-16 | 4.885e-15 |
| 6  | 54.2 % | 1.843e-14 | 7.767e-13 | 3.380e-10 | 2.864e-12 |
| 7  | 54.2 % | 4.013e-14 | 1.692e-12 | 7.362e-10 | 6.237e-12 |
| 8  | 47.5 % | 9.362e-17 | 2.341e-16 | 4.173e-13 | 3.286e-14 |
| 9  | 54.2 % | 2.859e-14 | 1.507e-12 | 8.795e-10 | 4.651e-12 |
| 10 | 54.2 % | 5.933e-11 | 1.905e-09 | 1.556e-07 | 9.751e-09 |
| 13 | 54.1 % | 4.329e-11 | 1.388e-09 | 1.133e-07 | 7.162e-09 |
| 14 | 54.1 % | 4.178e-13 | 2.494e-11 | 1.540e-08 | 9.170e-11 |
| 15 | 54.1 % | 4.052e-13 | 2.522e-11 | 1.749e-08 | 9.646e-11 |
| 16 | 54.1 % | 1.268e-13 | 2.329e-11 | 2.389e-08 | 1.192e-10 |

Column 3 carries both the largest values and the largest scatter, and its
worst-case sd (5.2e-07) is twice the median value of the column itself.
Column 5 moves only at the denormal level and is effectively constant.

## What is not measured here

The o2 arm. O2 is not installed in this container -- no `O2_ROOT`, no `alienv`,
no CVMFS -- and `GEOM_BACKEND=o2` resolves `o2::its::GeometryTGeo` through cling
at run time, so it cannot be run here at all. Run the same ten repetitions on a
machine with O2 loaded and pass both `costs.tsv` to `weight_spread.py` for the
comparison.

At n=10 per arm the 95 % CI on a ratio of standard deviations spans roughly
[0.5, 2.0] even when the two arms are identical, so that comparison can only
resolve a difference larger than about a factor of two. n=30 would bring it to
about [0.7, 1.5].

## Reproducing

    ./config/runctl.sh set GEOM_BACKEND=cache JOB_TAG=rt-cache JOB_NDATA=4000 \
        JOB_NEPOCH=1 DATA_FILE=XXXXinput.root PARAMS_ARCHIVE=MLPTrain_Step900
    ./config/runctl.sh compose
    ./tools/monitoring/run_to_run.sh runs/rt-cache cache 10 1
    ./tools/monitoring/weight_spread.py   plots/run_to_run_weights.png cache=runs/rt-cache.rtr-cache/costs.tsv
    ./tools/monitoring/plot_run_to_run.py plots/run_to_run_cost.png    cache=runs/rt-cache.rtr-cache/costs.tsv
