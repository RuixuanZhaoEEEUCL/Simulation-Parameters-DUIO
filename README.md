# Supplementary Documents for DUIO

Supplementary materials for the journal paper [**Distributed Unknown Input Observer Design: A Geometric Approach**](https://arxiv.org/abs/2603.18335).

## Files

- [Simulation_Parameters.pdf](Simulation_Parameters.pdf) — Numerical parameters used in the paper's simulations.
- [duio_algorithm_1.m](duio_algorithm_1.m) — Algorithm 1: subspace bases and canonical projections.
- [duio_algorithm_2.m](duio_algorithm_2.m) — Algorithm 2: good/bad spectral separation.
- [duio_algorithm_3.m](duio_algorithm_3.m) — Algorithm 3: output-injection gain design.

## MATLAB Usage

Run the three functions in order for each sensor node, using its state matrix `A`, unknown-input matrix `Bbar`, and measurement matrix `C`. Place all three `.m` files on your MATLAB path. The functions implement local offline geometric observer design; network consensus is not included.

**Requirements:** MATLAB R2022a or later. Control System Toolbox is required for Algorithm 3 when free poles must be assigned.

## Algorithm 1 — Subspace bases and projections

### Syntax

```matlab
[Wi,PWi,Si,PSi] = duio_algorithm_1(A,Bbar,C)
[Wi,PWi,Si,PSi,info] = duio_algorithm_1(A,Bbar,C,opts)
```

### Description

Computes the infimal conditioned-invariant subspace `W*` containing the unknown-input image, the associated infimal unobservability subspace `S*`, and their quotient projections. The same function applies to continuous- and discrete-time models.

### Input Arguments

| Argument | Description |
|---|---|
| `A` | Real `n`-by-`n` state matrix. |
| `Bbar` | Real `n`-by-`mu` matrix of inputs unknown to this node. Use `[]` when there are no unknown inputs. |
| `C` | Real `p`-by-`n` local output matrix. Use `[]` when there are no measurements. |
| `opts` | Optional options structure; see **Options** below. |

### Output Arguments

Here, `w = size(Wi,2)` and `s = size(Si,2)`.

| Output | Size | Description |
|---|---|---|
| `Wi` | `n`-by-`w` | Orthonormal basis of `W*`. |
| `PWi` | `(n-w)`-by-`n` | Projection with orthonormal rows and `Ker(PWi) = Im(Wi)`. |
| `Si` | `n`-by-`s` | Orthonormal basis of `S*`. |
| `PSi` | `(n-s)`-by-`n` | Projection with orthonormal rows and `Ker(PSi) = Im(Si)`. |
| `info` | Structure | Subspace dimensions, iteration histories, and numerical residuals. |

## Algorithm 2 — Good/bad spectral separation

### Syntax

```matlab
[Wgi,PWgi] = duio_algorithm_2(A,C,Wi,Si,PWi,mode)
[Wgi,PWgi,info] = duio_algorithm_2(A,C,Wi,Si,PWi,mode,opts)
```

### Description

Uses ordered real Schur decomposition to separate the fixed modes on `S*/W*` into good and bad regions. Adds the bad directions to `W*` and returns the basis and projection for `Wg*`.

### Input Arguments

| Argument | Description |
|---|---|
| `A`, `C` | The same state and output matrices used in Algorithm 1. |
| `Wi`, `Si`, `PWi` | Outputs from Algorithm 1, passed without modification. |
| `mode` | `'ct'` for continuous-time Hurwitz stability; `'dt'` for discrete-time Schur stability. |
| `opts` | Optional structure specifying the stability margin and numerical tolerances. |

### Output Arguments

Here, `wg = size(Wgi,2)`.

| Output | Size | Description |
|---|---|---|
| `Wgi` | `n`-by-`wg` | Orthonormal basis of `Wg*`. |
| `PWgi` | `(n-wg)`-by-`n` | Projection with orthonormal rows and `Ker(PWgi) = Im(Wgi)`. |
| `info` | Structure | Fixed, good, and bad eigenvalues; Schur factors; dimensions; and numerical residuals. |

Use `info.GoodEigenvalues` and `info.BadEigenvalues` to inspect the classification. If `wg = n`, `PWgi` has zero rows.

## Algorithm 3 — Output-injection gain

### Syntax

```matlab
Li = duio_algorithm_3(A,C,Wi,PWi,Si,PSi,desiredPoles,mode)
[Li,info] = duio_algorithm_3(A,C,Wi,PWi,Si,PSi,desiredPoles,mode,opts)
```

### Description

Assigns the free poles on `X/S*` and returns the output-injection gain `Li`, using the convention **`A + Li*C`**. The good fixed poles selected by Algorithm 2 are retained.

### Input Arguments

| Argument | Description |
|---|---|
| `A`, `C` | The same state and output matrices used in Algorithm 1. |
| `Wi`, `PWi`, `Si`, `PSi` | Outputs from Algorithm 1, passed without modification. |
| `desiredPoles` | Vector of exactly `n-size(Si,2)` free poles. Use `[]` if this number is zero. |
| `mode`, `opts` | The same time domain and stability-region settings used in Algorithm 2. |

Desired poles must lie strictly inside the selected good region. Nonreal poles must occur in complex-conjugate pairs. Repeated-pole multiplicity must not exceed the rank of the quotient output matrix `C0`.

Pole placement requires Control System Toolbox when there are free poles.

### Output Arguments

| Output | Size | Description |
|---|---|---|
| `Li` | `n`-by-`p` | Output-injection gain. |
| `info` | Structure | Free quotient matrices `A0`, `C0`, gain `L0`, assigned poles `ActualFreePoles`, and numerical residuals. |

## Options

Pass a structure as the last argument. Omitted fields use their defaults.

| Field | Default | Description |
|---|---|---|
| `StabilityMargin` | `0` | Design margin for the good/bad region. |
| `BoundaryTolerance` | `1e-9` | Additional numerical guard at the spectral boundary. |
| `RankTolerance` | `1e-10` | Relative coefficient for SVD rank decisions. |
| `CheckTolerance` | `1e-7` | Numerical residual acceptance threshold. |
| `PoleTolerance` | `1e-6` | Relative tolerance for matching assigned poles. |
| `SpectralMethod` | `'schur'` | Spectral separation method; only `'schur'` is supported. |

With `m = opts.StabilityMargin` and `b = opts.BoundaryTolerance`, the implemented good regions are:

| Mode | Good region |
|---|---|
| `'ct'` | `real(lambda) < -m-b` |
| `'dt'` | `abs(lambda) < 1-m-b` |

All remaining modes, including boundary modes, are bad. Both `m` and `b` must be nonnegative; discrete-time mode also requires `m+b < 1`. Setting both to zero gives strict Hurwitz or Schur stability. A positive margin excludes slowly decaying stable modes.

Use the same `mode` and `opts` in Algorithms 2 and 3. Selecting `'dt'` does not discretize `A`; supply a discrete-time model.

## Examples

### Continuous-time design

```matlab
A = [0  1  2  0;
     0 -2  0  0;
     0  0  1  0;
     0  0  0 -3];
Bbar = [1; 0; 0; 0];
C = [1 0 0 0;
     0 0 0 1];

mode = 'ct';
opts = struct('StabilityMargin',0);

% Algorithm 1: compute the geometric bases and projections.
[Wi,PWi,Si,PSi] = duio_algorithm_1(A,Bbar,C,opts);

% Algorithm 2: select the good and bad fixed modes.
[Wgi,PWgi,info2] = duio_algorithm_2(A,C,Wi,Si,PWi,mode,opts);

% Algorithm 3: assign the n-size(Si,2) free poles.
desiredPoles = -4;
[Li,info3] = duio_algorithm_3( ...
    A,C,Wi,PWi,Si,PSi,desiredPoles,mode,opts);

% Reduced-order error dynamics.
Abar = PWgi*(A+Li*C)*PWgi.';
disp([size(Wi,2),size(Si,2),size(Wgi,2)]);  % [1 3 2]
disp(eig(Abar));                           % -4 and -2; order may vary
```

The fixed good pole `-2` is retained, the fixed bad pole `1` is excluded from the reduced-order error dynamics, and the free pole is assigned to `-4`.

### Discrete-time design with a stability margin

```matlab
A = [0.2 1   2      0;
     0   0.5 0      0;
     0   0   0.9999 0;
     0   0   0      0.3];
Bbar = [1; 0; 0; 0];
C = [1 0 0 0;
     0 0 0 1];

mode = 'dt';
opts = struct('StabilityMargin',1e-3);

[Wi,PWi,Si,PSi] = duio_algorithm_1(A,Bbar,C,opts);
[Wgi,PWgi,info2] = duio_algorithm_2(A,C,Wi,Si,PWi,mode,opts);

desiredPoles = 0.2;
[Li,info3] = duio_algorithm_3( ...
    A,C,Wi,PWi,Si,PSi,desiredPoles,mode,opts);

Abar = PWgi*(A+Li*C)*PWgi.';
disp([size(Wi,2),size(Si,2),size(Wgi,2)]);  % [1 3 2]
disp(eig(Abar));                           % 0.2 and 0.5; order may vary
disp(info2.BadEigenvalues);                % 0.9999
```

The margin classifies the slowly decaying fixed pole `0.9999` as bad. The reduced-order error dynamics retain `0.5` and use the assigned free pole `0.2`.
