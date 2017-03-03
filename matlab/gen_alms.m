function alms=gen_alms(cl,mmax,seed)
% alms=gen_alms(cl,mmax,seed)
%
% Draws random alms with variance cl
%
% INPUTS
%
%   cl      C_l spectra, each as a column vector.
%
%   mmax    Maximum m mode, where mmax <= lmax. If ell > mmax, then only
%           the first mmax alm values are filled.
%
%   seed    Optional. Defaults to 1336+(1:size(cl,2)). May also be a cell array
%           of random number generators for each spectrum, initialized to have
%           lmax+1-many substreams.
%
% OUTPUTS
%
%   alms    Array of alm coefficients. If size(cl,2)==1, then alms is
%           2D with output shape (lmax+1 x mmax+1). If size(cl,2)>1, then
%           alms is 3D with output shape (size(cl,2) x lmax+1 x mmax+1).
%           (3D form is compatible with input for alm2map().)
%
% EXAMPLE
%
%   cl   = [0, 1./(1:700).^2]';
%   alms = gen_alms([cl,cl,cl], [], [1,2,3]);
%

  lmax  = size(cl, 1);
  nspec = size(cl, 2);

  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax;
  end

  if ~exist('seed','var') || isempty(seed)
    seed = 1336 + (1:nspec);
  end
  if length(seed) ~= nspec
    error('%d seeds or RNGs are required when given %d spectra', nspec, nspec)
  end
  if isnumeric(seed)
    for ii=1:nspec
      try
        rngstr{ii} = RandStream.create('mlfg6331_64', 'Seed',seed(ii), ...
            'NumStreams',lmax+1, 'NormalTransform','Ziggurat');
      catch
        rngstr{ii} = RandStream.create('mlfg6331_64', 'Seed',seed(ii), ...
            'NumStreams',lmax+1, 'RandnAlg','Ziggurat');
      end
    end
  end

  alms = zeros(nspec, lmax+1, mmax+1);
  for ii=1:nspec
    alms(ii,:,:) = gen_one(cl(:,ii), lmax, mmax, rngstr{ii});
  end

  alms = squeeze(alms);
end

function alms=gen_one(cl,lmax,mmax,rngstr)
  rt2 = 1/sqrt(2);

  alms = zeros(lmax+1,mmax+1);
  for ell=0:(lmax-1)
    scale = sqrt(cl(ell+1));

    % Initialize substream in a predictable way (so that both mmax==lmax and
    % mmax<lmax cases give the same initial random deviates).
    rngstr.Substream = ell+1;
    % m == 0 must be real
    alms(ell+1,1) = scale*rngstr.randn();
    % Other m's are complex
    ii = min(ell+1,mmax+1);
    alms(ell+1,2:ii) = complex(scale*rt2*rngstr.randn(1,ii-1), ...
        scale*rt2*rngstr.randn(1,ii-1));
  end
end
