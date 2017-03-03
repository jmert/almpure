function alms=gen_delta_alms(ell,lmax,mmax,delta,seed)
% alms=gen_delta_alms(ell,lmax,mmax,delta,seed)
%
% Draws random alms for C_l = 1 (or l(l+1)C_l = 1).
%
% INPUTS
%
%   ell     Delta ell location. Can be an array of values in which case a
%           corresponding number of sets of alms are generated.
%           (Note that ell == 0 is ignored.)
%
%   lmax    Maximum ell. Controls dimensions of alms array.
%
%   mmax    Maximum m mode, where mmax <= lmax. If ell > mmax, then only
%           the first mmax alm values are filled.
%
%   delta   Optional. Selects whether delta response is in terms of C_l
%           (the 'cl' option) or D_l (the 'dl' option; i.e. l(l+1)*Cl).
%           Defaults to 'cl'.
%
%   seed    Optional. Defaults to 1336+(1:length(ell)), or if a scalar,
%           replicated to length(ell).
%
% OUTPUTS
%
%   alms    Array of alm coefficients. If length(ell)==1, then alms is
%           2D with output shape (lmax+1 x mmax+1). If length(ell)>1, then
%           alms is 3D with output shape (length(ell) x lmax+1 x mmax+1).
%           (3D form is compatible with input for alm2map().)
%
% EXAMPLE
%
%   % Generate alms for a E-mode only, ell=80 spectrum
%   alms = gen_delta_alms([0,80,0], 700, 700, [1,2,3]);
%

  if ~exist('delta','var') || isempty(delta)
    delta = 'cl';
  end

  if ~exist('seed','var') || isempty(seed)
    seed = 1336 + (1:length(ell));
  end
  if isnumeric(seed)
    if length(seed) == 1
      seed = repmat(seed, length(ell), 1);
    end
    for ii=1:length(seed)
      rngstr{ii} = RandStream.create('mlfg6331_64', 'Seed',seed(ii), ...
          'NumStreams',length(seed), 'NormalTransform','Ziggurat');
    end
  else
    if ~iscell(seed)
      rngstr = repmat({seed}, length(ell), 1);
    end
  end

  alms = zeros(length(ell), lmax+1, mmax+1);
  for ii=1:length(ell)
    if ell(ii) == 0
      continue
    end
    alms(ii,:,:) = gen_one(ell(ii), lmax, mmax, delta, rngstr{ii});
  end

  alms = squeeze(alms);
end

function alms=gen_one(ell,lmax,mmax,delta,rngstr)
  if strcmp(delta,'dl')
    scale = sqrt(2*pi/(ell*(ell+1)));
  else
    scale = 1;
  end
  rt2 = scale/sqrt(2);

  alms = zeros(lmax+1,mmax+1);
  % m == 0 must be real
  alms(ell+1,1) = scale*rngstr.randn();
  % Other m's are complex
  ii = min(ell+1,mmax+1);
  alms(ell+1,2:ii) = complex(rt2*rngstr.randn(1,ii-1), rt2*rngstr.randn(1,ii-1));
end
