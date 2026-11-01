function [Wgi,PWgi,info] = duio_algorithm_2(A,C,Wi,Si,PWi,mode,opts)
%DUIO_ALGORITHM2 Good/bad fixed-spectrum separation using ordered real Schur.
%   [Wgi,PWgi] = DUIO_ALGORITHM2(A,C,Wi,Si,PWi,mode)
%   [Wgi,PWgi,info] = DUIO_ALGORITHM2(A,C,Wi,Si,PWi,mode,opts)
%
%   Syntax uses the SAME A and C supplied to DUIO_ALGORITHM1. Wi, Si and
%   PWi are its orthonormal basis/projection outputs. Inputs must be finite,
%   real numeric matrices; internal calculations use full double precision.
%
%   INPUT ARGUMENTS
%     A     n-by-n state matrix (nonempty, square).
%     C     p-by-n local output matrix; [] means no output.
%     Wi    n-by-w orthonormal basis of W*.
%     Si    n-by-s orthonormal basis of S*, with Im(Wi) contained in Im(Si).
%     PWi   (n-w)-by-n matrix with orthonormal rows and Ker(PWi)=Im(Wi).
%     mode  'ct': continuous time, Hurwitz good region.
%           'dt': discrete time, Schur-stable good region.
%     opts  Scalar structure. Defaults:
%           StabilityMargin    = 0     (good/bad REGION DESIGN parameter)
%           BoundaryTolerance  = 1e-9  (absolute spectral boundary guard)
%           RankTolerance      = 1e-10 (SVD cutoff coefficient)
%           CheckTolerance     = 1e-7  (residual acceptance threshold)
%           PoleTolerance      = 1e-6  (shared with Algorithm 3; unused here)
%           SpectralMethod     = 'schur' (the only supported backend in v2)
%
%   GOOD/BAD CLASSIFICATION
%     With m=StabilityMargin and b=BoundaryTolerance, a fixed pole is good if
%       'ct': real(lambda) < -m-b,
%       'dt': abs(lambda)  <  1-m-b.
%     All other poles are bad, including boundary poles. m,b must be finite
%     and nonnegative; 'dt' additionally requires m+b<1. m has inverse-time
%     units in continuous time and is dimensionless in discrete time.
%     Setting m=b=0 gives strict Hurwitz/Schur stability. A positive margin
%     can classify otherwise stable but slowly decaying modes as bad.
%     BoundaryTolerance is NOT the matrix-rank tolerance.
%
%   OUTPUT ARGUMENTS
%     Wgi   n-by-wg orthonormal basis of W_g*.
%     PWgi  (n-wg)-by-n matrix with orthonormal rows and Ker(PWgi)=Im(Wgi).
%           If wg=n, PWgi has size 0-by-n, NOT a row of zeros.
%     info  Structure with V, AV, La, fixed/good/bad eigenvalues, bad
%           coordinate basis, reordered Schur factors, dimensions, region
%           settings, near-boundary eigenvalues and numerical residuals.
%
%   ALGORITHM
%     Paper Algorithm 2, lines 1--3, 7--8, are implemented directly. Lines
%     4--6 (minimal polynomial, factorization and bad-factor kernel) are
%     evaluated by their equivalent invariant-subspace operation using
%     SCHUR, ORDEIG and ORDSCHUR. No numerical minimal polynomial is formed.
%     In exact arithmetic, leading reordered Schur vectors span the same
%     generalized spectral subspace as Ker(beta_b(AV)). The lift V*Xbad is
%     essential: Xbad is in quotient coordinates, not state coordinates.
%
%   EXAMPLE
%     opts = struct('StabilityMargin',0.01);
%     [Wi,PWi,Si,PSi] = duio_algorithm1(A,Bbar,C,opts);
%     [Wgi,PWgi,info] = duio_algorithm2(A,C,Wi,Si,PWi,'dt',opts);
%     % Choose n-size(Si,2) free poles in the SAME good region, then call
%     % duio_algorithm3. No free poles: desiredPoles = [].
%
%   LIMITATIONS
%     Numerical rank and spectral separation are tolerance-dependent.
%     Inputs inconsistent with the required geometry, or failed residual
%     checks, raise an error. Closely clustered spectra near the selected
%     boundary may be ill-conditioned; inspect info.NearBoundaryEigenvalues.
%     This routine does not assign free poles or design network consensus.
%     MATLAB R2022a+; no Control System or Symbolic Math Toolbox required.
%

narginchk(6,7);
if nargin < 7, opts = struct(); end
o = local_options(opts);
mode = local_mode(mode,o);
[A,C,Wi,Si,PWi,n,~,w,s] = local_geometry(A,C,Wi,Si,PWi,o);

% Algorithm 2, line 1: orthogonal complement of W* inside S*.
V = local_orth(Si-Wi*(Wi.'*Si),o);
q = s-w;
if size(V,2) ~= q
    error('duio:DimensionMismatch','V must have s-w columns. Check the bases/tolerance.');
end
local_check(norm(Wi.'*V,'fro'),o,'W^T*V = 0');
local_check(norm(V.'*V-eye(q),'fro'),o,'V orthonormality');

% Algorithm 2, line 2: one friend satisfying PWi*(A+La*C)*Wi = 0.
[La,friendResidual] = local_friend(A,C,Wi,PWi,o);
ALa = A+La*C;
PS = local_null(Si.',o).';
sResidual = local_norm(PS*ALa*Si)/max(1,local_norm(ALa));
local_check(sResidual,o,'S invariant under the chosen friend');

% Algorithm 2, line 3: dynamics on S*/W* in V-coordinates.
AV = V.'*ALa*V;

% Algorithm 2, lines 4--6: equivalent ordered real Schur computation.
% ORDEIG, not EIG(T), supplies the ordering required by ORDSCHUR.
if q == 0
    Us = zeros(0,0); Ts = zeros(0,0);
    lambda = zeros(0,1); orderedLambda = zeros(0,1);
    isBad = false(0,1); nb = 0;
else
    [U,T] = schur(AV,'real');
    lambda = ordeig(T);
    isBad = ~local_good(lambda,mode,o);
    % Keep each real 2-by-2 Schur block intact. Both members of a conjugate
    % pair have the same real part/modulus under the supported regions.
    j = 1;
    while j < q
        if T(j+1,j) ~= 0
            isBad(j:j+1) = any(isBad(j:j+1));
            j = j+2;
        else
            j = j+1;
        end
    end
    nb = sum(isBad);
    if nb == 0 || nb == q
        Us = U; Ts = T;
    else
        try
            [Us,Ts] = ordschur(U,T,isBad);
        catch ME
            error('duio:SchurReorderingFailed', ...
                'Schur reordering failed; inspect scaling/boundary separation: %s',ME.message);
        end
    end
    orderedLambda = ordeig(Ts);
    if any(local_good(orderedLambda(1:nb),mode,o)) || ...
            any(~local_good(orderedLambda(nb+1:end),mode,o))
        error('duio:SpectralSeparation', ...
            ['Classification changed during Schur reordering. Inspect poles near ' ...
             'the boundary, improve scaling, or use a justified boundary guard.']);
    end
end
Xbad = Us(:,1:nb);

% Algorithm 2, line 7: lift quotient coordinates and join the W* basis.
Wgi = local_orth([Wi,V*Xbad],o);
if size(Wgi,2) ~= w+nb
    error('duio:DimensionMismatch','The Wg lift lost rank. Check the bases/tolerance.');
end

% Algorithm 2, line 8: canonical quotient projection.
PWgi = local_null(Wgi.',o).';

% Diagnostics and consistency checks. These do not alter the algorithm.
info.V = V;
info.AV = AV;
info.La = La;
info.Method = 'schur';
info.SchurVectors = Us;
info.OrderedSchurMatrix = Ts;
info.OrderedEigenvalues = orderedLambda;
info.BadCoordinateBasis = Xbad;
info.FixedEigenvalues = lambda;
info.BadMask = isBad;
info.GoodEigenvalues = lambda(~isBad);
info.BadEigenvalues = lambda(isBad);
info.QuotientDimension = q;
info.BadDimension = nb;
info.WgDimension = size(Wgi,2);
info.FriendResidual = friendResidual;
info.SInvarianceResidual = sResidual;
info.SchurResidual = local_norm(AV-Us*Ts*Us.')/max(1,local_norm(AV));
info.SchurOrthonormality = norm(Us.'*Us-eye(q),'fro');
info.BadSubspaceResidual = local_norm(AV*Xbad-Xbad*(Xbad.'*AV*Xbad)) / ...
    max(1,local_norm(AV));
info.WgInvarianceResidual = local_norm(PWgi*ALa*Wgi)/max(1,local_norm(ALa));
info.WgOrthonormality = norm(Wgi.'*Wgi-eye(size(Wgi,2)),'fro');
info.ProjectionOrthonormality = norm(PWgi*PWgi.'-eye(n-size(Wgi,2)),'fro');
info.WContainmentResidual = local_norm(Wi-Wgi*(Wgi.'*Wi));
info.NestingResidual = local_norm(PS*Wgi);
info.ProjectionResidual = local_norm(PWgi*Wgi);
info.Mode = mode;
info.Options = o;
if strcmp(mode,'ct')
    info.DesignBoundary = -o.StabilityMargin;
    info.EffectiveBoundary = -o.StabilityMargin-o.BoundaryTolerance;
    info.GoodRegion = 'real(lambda) < -StabilityMargin-BoundaryTolerance';
    distanceToBoundary = abs(real(lambda)+o.StabilityMargin);
else
    info.DesignBoundary = 1-o.StabilityMargin;
    info.EffectiveBoundary = 1-o.StabilityMargin-o.BoundaryTolerance;
    info.GoodRegion = 'abs(lambda) < 1-StabilityMargin-BoundaryTolerance';
    distanceToBoundary = abs(abs(lambda)-(1-o.StabilityMargin));
end
info.NearBoundaryEigenvalues = lambda(distanceToBoundary <= o.BoundaryTolerance);
local_check(info.SchurResidual,o,'Schur factorization');
local_check(info.SchurOrthonormality,o,'Schur orthonormality');
local_check(info.BadSubspaceResidual,o,'Bad spectral subspace invariance');
local_check(info.WgInvarianceResidual,o,'Wg invariant under the chosen friend');
local_check(info.WgOrthonormality,o,'Wg orthonormality');
local_check(info.ProjectionOrthonormality,o,'PWg row orthonormality');
local_check(info.WContainmentResidual,o,'W contained in Wg');
local_check(info.NestingResidual,o,'Wg contained in S');
local_check(info.ProjectionResidual,o,'PWg*Wg = 0');
end

% =========================================================================
% Local numerical utilities (no extra helper file is required).
% =========================================================================
function o = local_options(o)
if nargin == 0 || isempty(o), o = struct(); end
if ~isstruct(o) || ~isscalar(o)
    error('duio:InvalidOptions', 'opts must be a scalar structure.');
end
def = struct('RankTolerance',1e-10, 'CheckTolerance',1e-7, ...
    'PoleTolerance',1e-6, 'BoundaryTolerance',1e-9, ...
    'StabilityMargin',0, 'SpectralMethod','schur');
names = fieldnames(o);
for k = 1:numel(names)
    if ~isfield(def,names{k})
        error('duio:InvalidOptions', 'Unknown option: %s.', names{k});
    end
    def.(names{k}) = o.(names{k});
end
o = def;
positive = {'RankTolerance','CheckTolerance','PoleTolerance'};
for k = 1:numel(positive)
    v = o.(positive{k});
    if ~isnumeric(v) || ~isreal(v) || ~isscalar(v) || ...
            isnan(v) || isinf(v) || v <= 0 || v >= 1
        error('duio:InvalidOptions', '%s must be between 0 and 1.',positive{k});
    end
end
nonnegative = {'BoundaryTolerance','StabilityMargin'};
for k = 1:numel(nonnegative)
    v = o.(nonnegative{k});
    if ~isnumeric(v) || ~isreal(v) || ~isscalar(v) || ...
            isnan(v) || isinf(v) || v < 0
        error('duio:InvalidOptions', '%s must be nonnegative.',nonnegative{k});
    end
end
% Normalize validated scalar numeric options. MATLAB does not allow an
% integer StabilityMargin to be added to a nonscalar double eigenvalue array.
scalarOptions = [positive, nonnegative];
for k = 1:numel(scalarOptions)
    name = scalarOptions{k};
    o.(name) = double(o.(name));
end
if ~(ischar(o.SpectralMethod) || ...
        (isstring(o.SpectralMethod) && isscalar(o.SpectralMethod)))
    error('duio:InvalidOptions', 'SpectralMethod must be ''schur'' in version 2.');
end
o.SpectralMethod = lower(char(o.SpectralMethod));
if ~strcmp(o.SpectralMethod,'schur')
    error('duio:InvalidOptions', 'SpectralMethod must be ''schur'' in version 2.');
end
end

function M = local_matrix(M,name)
if ~isnumeric(M) || ~ismatrix(M) || ~isreal(M) || ...
        any(isnan(M(:))) || any(isinf(M(:)))
    error('duio:InvalidMatrix','%s must be a finite real numeric matrix.',name);
end
M = full(double(M));
end

function s = local_norm(M)
if isempty(M), s = 0; else, s = norm(M,2); end
end

function t = local_cutoff(M,o,reference)
if nargin < 3, reference = 1; end
t = o.RankTolerance * max([1, reference, local_norm(M)]);
end

function Q = local_orth(M,o,reference)
if nargin < 3, reference = 1; end
if isempty(M)
    Q = zeros(size(M,1),0);
else
    Q = orth(M,local_cutoff(M,o,reference));
end
end

function Z = local_null(M,o,reference)
if nargin < 3, reference = 1; end
if size(M,2) == 0
    Z = zeros(0,0);
elseif size(M,1) == 0
    Z = eye(size(M,2));
else
    Z = null(M,local_cutoff(M,o,reference));
end
end

function D = local_pinv(M,o,reference)
if nargin < 3, reference = 1; end
if isempty(M)
    D = zeros(size(M,2),size(M,1));
else
    D = pinv(M,local_cutoff(M,o,reference));
end
end

function local_check(value,o,what)
if (isnan(value) || isinf(value)) || value > o.CheckTolerance
    error('duio:ResidualCheck', ...
        '%s: residual %.3g exceeds %.3g. Check scaling and rank tolerance.', ...
        what,value,o.CheckTolerance);
end
end

function [A,C,W,S,PW,n,p,w,s] = local_geometry(A,C,W,S,PW,o)
A = local_matrix(A,'A');
n = size(A,1);
if n == 0 || size(A,2) ~= n
    error('duio:DimensionMismatch','A must be nonempty and square.');
end
C = local_matrix(C,'C');
if isequal(size(C),[0,0]), C = zeros(0,n); end
if size(C,2) ~= n
    error('duio:DimensionMismatch','C must have n columns.');
end
p = size(C,1);
W = local_matrix(W,'W'); S = local_matrix(S,'S');
if isequal(size(W),[0,0]), W = zeros(n,0); end
if isequal(size(S),[0,0]), S = zeros(n,0); end
PW = local_matrix(PW,'PW');
w = size(W,2); s = size(S,2);
if size(W,1) ~= n || size(S,1) ~= n || w > s || s > n || ...
        ~isequal(size(PW),[n-w,n])
    error('duio:DimensionMismatch','Inconsistent W, S, or PW dimensions.');
end
local_check(norm(W.'*W-eye(w),'fro'),o,'W orthonormality');
local_check(norm(S.'*S-eye(s),'fro'),o,'S orthonormality');
local_check(norm(PW*PW.'-eye(n-w),'fro'),o,'PW row orthonormality');
local_check(norm(PW*W,'fro'),o,'PW*W = 0');
local_check(norm(W-S*(S.'*W),'fro'),o,'W contained in S');
end

function [La,res] = local_friend(A,C,W,PW,o)
% Solve PW*La*(C*W) = -PW*A*W. PW has orthonormal rows.
D = C*W;
La = -PW.'*(PW*A*W)*local_pinv(D,o,max(1,local_norm(C)));
if any(isnan(La(:))) || any(isinf(La(:)))
    error('duio:ResidualCheck','Nonfinite output-injection friend.');
end
res = local_norm(PW*(A+La*C)*W) / ...
    max([1,local_norm(A),local_norm(La)*local_norm(C)]);
local_check(res,o,'Conditioned-invariance friend equation');
end

function mode = local_mode(mode,o)
if ~(ischar(mode) || (isstring(mode) && isscalar(mode)))
    error('duio:InvalidMode','mode must be ''ct'' or ''dt''.');
end
mode = lower(char(mode));
if ~any(strcmp(mode,{'ct','dt'}))
    error('duio:InvalidMode','mode must be ''ct'' or ''dt''.');
end
if isinf(o.StabilityMargin+o.BoundaryTolerance)
    error('duio:InvalidOptions','StabilityMargin + BoundaryTolerance must be finite.');
end
if strcmp(mode,'dt') && o.StabilityMargin+o.BoundaryTolerance >= 1
    error('duio:InvalidOptions', ...
        'In dt mode, StabilityMargin + BoundaryTolerance must be below 1.');
end
end

function yes = local_good(lambda,mode,o)
if strcmp(mode,'ct')
    yes = real(lambda) < -o.StabilityMargin-o.BoundaryTolerance;
else
    yes = abs(lambda) < 1-o.StabilityMargin-o.BoundaryTolerance;
end
end
