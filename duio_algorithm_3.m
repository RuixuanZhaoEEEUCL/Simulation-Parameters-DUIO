function [Li,info] = duio_algorithm_3(A,C,Wi,PWi,Si,PSi,desiredPoles,mode,opts)
%DUIO_ALGORITHM3 Output injection by free-pole assignment (paper Algorithm 3).
%   [Li,info] = duio_algorithm3(A,C,Wi,PWi,Si,PSi,desiredPoles,mode)
%   [...] = duio_algorithm3(A,C,Wi,PWi,Si,PSi,desiredPoles,mode,opts)
%
% Wi,PWi,Si,PSi come from Algorithm 1 with the SAME A and C.
% desiredPoles: exactly n-size(Si,2) stable poles, counted with multiplicity.
%               Pass [] when n-size(Si,2)==0. Conjugate pairs are required.
% mode: 'ct' or 'dt'; use the same mode and margins as Algorithm 2.
%
% Li is n-by-p and uses the paper's PLUS convention A+Li*C.
% Only poles on X/S* are assigned. The good fixed poles on S*/W_g*
% are retained, NOT overwritten by desiredPoles.
%
% MATLAB Control System Toolbox is required only when there are free poles.
% The call to place includes a MINUS sign because place returns K for A-B*K.
% Redundant output rows are compressed only inside line 6; this parameterizes
% a valid L0 without changing A0+L0*C0 or the supplied output coordinates.
%

narginchk(8,9);
if nargin < 9, opts = struct(); end
o = local_options(opts);
mode = local_mode(mode,o);
[A,C,Wi,Si,PWi,n,p,~,s] = local_geometry(A,C,Wi,Si,PWi,o);
PSi = local_matrix(PSi,'PSi');
nFree = n-s;
if ~isequal(size(PSi),[nFree,n])
    error('duio:DimensionMismatch','PSi must have size (n-s)-by-n.');
end
local_check(norm(PSi*PSi.'-eye(nFree),'fro'),o,'PSi row orthonormality');
local_check(norm(PSi*Si,'fro'),o,'PSi*Si = 0');

if ~isnumeric(desiredPoles) || ...
        (~isvector(desiredPoles) && ~isempty(desiredPoles)) || ...
        any(isnan(desiredPoles(:))) || any(isinf(desiredPoles(:)))
    error('duio:InvalidPoles','desiredPoles must be a finite numeric vector.');
end
desiredPoles = double(desiredPoles(:));
if numel(desiredPoles) ~= nFree
    error('duio:PoleCount','Exactly n-s = %d desired poles are required.',nFree);
end
if ~all(local_good(desiredPoles,mode,o))
    error('duio:UnstablePoles', ...
        'All desired poles must be strictly inside the specified good region.');
end
[paired,~] = local_match_spectra(desiredPoles,conj(desiredPoles),1e-10);
if ~paired
    error('duio:ConjugatePoles','Nonreal desired poles must occur in conjugate pairs.');
end

% Algorithm 3, line 1: one friend of W*.
[La,friendResidual] = local_friend(A,C,Wi,PWi,o);
ALa = A+La*C;
sResidual = local_norm(PSi*ALa*Si)/max(1,local_norm(ALa));
local_check(sResidual,o,'S invariant under the chosen friend');

% Algorithm 3, line 2: quotient dynamics on X/S*.
PSdag = local_pinv(PSi,o);
A0 = PSi*ALa*PSdag;

% Algorithm 3, line 3: orthonormal basis of C*S*; r need not equal s.
SY = local_orth(C*Si,o,max(1,local_norm(C)));
r = size(SY,2);

% Algorithm 3, line 4: projection eliminating C*S* from the output space.
PY = local_null(SY.',o).';

% Algorithm 3, line 5: quotient output map, size (p-r)-by-(n-s).
C0 = PY*C*PSdag;
if ~isequal(size(C0),[p-r,nFree])
    error('duio:DimensionMismatch','Unexpected quotient-output dimension.');
end

% Algorithm 3, line 6: choose L0 with the prescribed free poles.
if nFree == 0
    L0 = zeros(0,p-r);  % Empty quotient: no pole-placement call is needed.
    outputRank = 0;
else
    % Remove redundant output directions only for the numerical place call.
    Uout = local_orth(C0,o,max(1,local_norm(C)));
    outputRank = size(Uout,2);
    Ceff = Uout.'*C0;
    if outputRank == 0 || local_observability_rank(A0,Ceff,o) ~= nFree
        error('duio:NotObservable', ...
            'The free quotient pair is not numerically observable. Check inputs/tolerance.');
    end
    % MATLAB place restriction, NOT a new theoretical existence condition.
    for k = 1:nFree
        if sum(desiredPoles == desiredPoles(k)) > outputRank
            error('duio:RepeatedPoles', ...
                ['A repeated pole has multiplicity greater than rank(C0)=%d. ' ...
                 'Choose distinct poles or use another pole-placement solver.'],outputRank);
        end
    end
    if exist('place','file') == 0 && exist('place','builtin') == 0
        error('duio:MissingPlace', ...
            'Algorithm 3 requires place from MATLAB Control System Toolbox.');
    end
    try
        K = place(A0.',Ceff.',desiredPoles.');
    catch ME
        error('duio:PolePlacementFailed','place failed: %s',ME.message);
    end
    L0 = -K.'*Uout.';  % A0 + L0*C0 = A0 - K^T*Ceff.
    if norm(imag(L0),'fro') > o.CheckTolerance*max(1,norm(L0,'fro'))
        error('duio:ConjugatePoles','Pole placement did not produce a real gain.');
    end
    L0 = real(L0);
end

% Algorithm 3, line 7: lift through P_S*, NOT through P_Wg*.
Li = La+PSdag*L0*PY;
if ~isequal(size(Li),[n,p]) || any(isnan(Li(:))) || any(isinf(Li(:)))
    error('duio:DimensionMismatch','Li must be a finite n-by-p matrix.');
end

% Verify the implemented free-pole assignment and the lifting identity.
Aclosed = A0+L0*C0;
actualPoles = eig(Aclosed);
[matched,poleError] = local_match_spectra(actualPoles,desiredPoles,o.PoleTolerance);
if ~matched
    error('duio:PoleVerification', ...
        'The computed free poles do not match the requested poles to tolerance %.3g.', ...
        o.PoleTolerance);
end
if ~all(local_good(actualPoles,mode,o))
    error('duio:PoleVerification', ...
        'Computed free poles are outside the requested good region.');
end
info.La = La;
info.A0 = A0;
info.SY = SY;
info.PY = PY;
info.C0 = C0;
info.L0 = L0;
info.ImageOutputDimension = r;
info.EffectiveOutputRank = outputRank;
info.FreeDimension = nFree;
info.DesiredPoles = desiredPoles;
info.ActualFreePoles = actualPoles;
info.FreePoleRelativeError = poleError;
info.FriendResidual = friendResidual;
info.SInvarianceResidual = sResidual;
info.OutputAnnihilationResidual = local_norm(PY*C*Si)/max(1,local_norm(C));
info.LiftingResidual = local_norm(PSi*(A+Li*C)*PSdag-Aclosed)/max(1,local_norm(Aclosed));
info.FixedBlockPreservationResidual = local_norm((Li-La)*C*Si) / ...
    max([1,local_norm(Li-La)*local_norm(C)]);
info.FinalWInvarianceResidual = local_norm(PWi*(A+Li*C)*Wi)/max(1,local_norm(A+Li*C));
info.Mode = mode;
info.Options = o;
local_check(info.OutputAnnihilationResidual,o,'PY*C*Si = 0');
local_check(info.LiftingResidual,o,'Lifted quotient matrix');
local_check(info.FixedBlockPreservationResidual,o,'Fixed S* dynamics preserved');
local_check(info.FinalWInvarianceResidual,o,'Final W* invariance');
end

function r = local_observability_rank(A,C,o)
% Span of C^T,A^T*C^T,... with normalization to avoid growing powers.
n = size(A,1);
Q = local_orth(C.',o);
As = A/max(1,local_norm(A));
for k = 1:n
    if size(Q,2) == n, break; end
    Q = local_orth([Q,As.'*Q],o);
end
r = size(Q,2);
end

function [ok,err] = local_match_spectra(actual,target,tol)
% Multiset comparison by bipartite matching; respects repeated eigenvalues.
% No Statistics/Optimization Toolbox is needed.
actual = actual(:); target = target(:);
n = numel(actual);
if n ~= numel(target), ok = false; err = Inf; return; end
if n == 0, ok = true; err = 0; return; end
D = abs(bsxfun(@minus,actual,target.'));
D = bsxfun(@rdivide,D,max(1,abs(target.')));
adj = D <= tol;
match = zeros(1,n); % match(j) is the actual pole assigned to target j.
for i = 1:n
    seen = false(1,n);
    [found,match,~] = local_augment(i,adj,match,seen);
    if ~found, ok = false; err = Inf; return; end
end
err = 0;
for j = 1:n, err = max(err,D(match(j),j)); end
ok = true;
end

function [found,match,seen] = local_augment(i,adj,match,seen)
found = false;
for j = 1:size(adj,2)
    if adj(i,j) && ~seen(j)
        seen(j) = true;
        if match(j) == 0
            match(j) = i; found = true; return;
        end
        [augmentOK,match,seen] = local_augment(match(j),adj,match,seen);
        if augmentOK
            match(j) = i; found = true; return;
        end
    end
end
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
