function [Wi,PWi,Si,PSi,info] = duio_algorithm_1(A,Bbar,C,opts)
%DUIO_ALGORITHM1 Bases and canonical projections (paper Algorithm 1).
%   [Wi,PWi,Si,PSi,info] = duio_algorithm1(A,Bbar,C)
%   [...] = duio_algorithm1(A,Bbar,C,opts)
%
% INPUTS (one sensor node):
%   A     n-by-n real state matrix.
%   Bbar  n-by-mu locally UNKNOWN input matrix; [] means no unknown input.
%   C     p-by-n local output matrix; [] means no output.
%   opts  Optional shared options structure; see README.md.
%
% OUTPUTS:
%   Wi    n-by-w orthonormal basis of the infimal conditioned-invariant W*.
%   PWi   (n-w)-by-n projection with Ker(PWi) = Im(Wi), orthonormal rows.
%   Si    n-by-s orthonormal basis of the infimal unobservability S*.
%   PSi   (n-s)-by-n projection with Ker(PSi) = Im(Si), orthonormal rows.
%   info  Dimension histories and numerical residuals.
%
% Each numbered block refers to the supplied, corrected pseudocode.
% In the pseudocode Ker(M) means a null-space BASIS MATRIX: null(M).
% orth/null use an explicit SVD cutoff; the iteration is not otherwise changed.
% The n iterations are retained exactly (no early exit).
%

narginchk(3,4);
if nargin < 4, opts = struct(); end
o = local_options(opts);
A = local_matrix(A,'A');
n = size(A,1);
if n == 0 || size(A,2) ~= n
    error('duio:DimensionMismatch','A must be nonempty and square.');
end
Bbar = local_matrix(Bbar,'Bbar'); C = local_matrix(C,'C');
if isequal(size(Bbar),[0,0]), Bbar = zeros(n,0); end
if isequal(size(C),[0,0]), C = zeros(0,n); end
if size(Bbar,1) ~= n || size(C,2) ~= n
    error('duio:DimensionMismatch','Bbar must have n rows; C must have n columns.');
end

% Algorithm 1, line 1: W_i^0 <- 0_(n-by-n).
Wk = zeros(n,n);
wHistory = zeros(n+1,1);

% Algorithm 1, lines 2--6: conditioned-invariant iteration.
for k = 1:n
    % Line 3: P_i^k <- (Ker((W_i^(k-1))^T))^T.
    Pk = local_null(Wk.',o).';
    % Line 4: T_i^k <- Ker([P_i^k; C_i]).
    Tk = local_null([Pk;C],o);
    % Line 5: W_i^k <- orth([Bbar_i, A*T_i^k]).
    Wk = local_orth([Bbar,A*Tk],o);
    wHistory(k+1) = size(Wk,2);
end

% Algorithm 1, lines 7--8: recompute the FINAL projection, then return W*.
PWi = local_null(Wk.',o).';
Wi = Wk;

% Algorithm 1, line 9: S_i^0 <- I_n.
Sk = eye(n);
sHistory = zeros(n+1,1); sHistory(1) = n;

% Algorithm 1, lines 10--14: unobservability-subspace iteration.
for k = 1:n
    % Line 11: P_i^k <- (Ker((S_i^(k-1))^T))^T.
    Pk = local_null(Sk.',o).';
    % Line 12: T_i^k <- Ker([P_i^k*A; C_i]).
    Tk = local_null([Pk*A;C],o);
    % Line 13: S_i^k <- orth([W_i^*, T_i^k]).
    Sk = local_orth([Wi,Tk],o);
    sHistory(k+1) = size(Sk,2);
end

% Algorithm 1, lines 15--16: FINAL projection and S*.
PSi = local_null(Sk.',o).';
Si = Sk;

% Numerical verification of the computed fixed points (not extra iterations).
Wnext = local_orth([Bbar,A*local_null([PWi;C],o)],o);
Snext = local_orth([Wi,local_null([PSi*A;C],o)],o);
info.WDimension = size(Wi,2);
info.SDimension = size(Si,2);
info.WDimensionHistory = wHistory;
info.SDimensionHistory = sHistory;
info.WFixedPointResidual = norm(Wi*Wi.'-Wnext*Wnext.','fro');
info.SFixedPointResidual = norm(Si*Si.'-Snext*Snext.','fro');
info.UnknownInputResidual = local_norm(PWi*Bbar)/max(1,local_norm(Bbar));
info.NestingResidual = local_norm(PSi*Wi);
info.WOrthonormality = norm(Wi.'*Wi-eye(size(Wi,2)),'fro');
info.SOrthonormality = norm(Si.'*Si-eye(size(Si,2)),'fro');
info.Options = o;
local_check(info.WFixedPointResidual,o,'W iteration fixed point');
local_check(info.SFixedPointResidual,o,'S iteration fixed point');
local_check(info.UnknownInputResidual,o,'Im(Bbar) contained in W');
local_check(info.NestingResidual,o,'W contained in S');
local_check(info.WOrthonormality,o,'W orthonormality');
local_check(info.SOrthonormality,o,'S orthonormality');
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
% Convert validated numeric options to double.
scalarOptions = {'RankTolerance', 'CheckTolerance', ...
    'PoleTolerance', 'BoundaryTolerance', 'StabilityMargin'};

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
