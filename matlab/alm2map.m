function map=alm2map(alms,nside)
% map=alm2map(alms,nside)
%
% Generates maps given a set of alms using the S2HAT library.
%
% INPUTS
%
%   alms     A 4-D array of complex alms for the non-negative m values. The
%            dimensions correspond to (nstokes, nlmax, nmmax, nmaps) where
%            each of the values means:
%
%              nstokes  Number of Stokes parameters to generate. Should be
%                       either 1 or 3.
%
%              nlmax    l-mode maximum, where 1 < nlmax
%
%              nmmax    m-mode maximum, where 1 < nmmax <= nlmax
%
%              nmaps    Number of distinct maps
%                       Note! nmaps ~= 1 is not yet supported.
%
%   nside    The nside of the output HEALPix map to render from the
%            given alms. If empty, defaults to nlmax/2.
%
% OUTPUTS
%
%   map      Synthesized map(s). Dimensions will be (12*nside*nside, nstokes)
%            where each of T,Q,U maps are stored in ring pixel ordering.
%
% EXAMPLE
%
%   nstokes = 3; nlmax = 80; nmmax = nlmax; nmaps = 1;
%   nside = 128;
%   alms = zeros(nstokes, nlmax+1, nmmax+1, nmaps);
%   % Make a random delta-80 map
%   alms(:,nmmax+1,:,:) = complex(...
%       randn(nstokes,nmmax+1,nmaps), randn(nstokes,nmmax+1,nmaps));
%
%   hmap.nside = nside;
%   hmap.ordering = 'ring';
%   hmap.pixel = 1:nside2npix(nside);
%   hmap.map = alm2map(alms, nside);
%

  if ~exist('nside','var') || isempty(nside)
    nside = 0.5*size(alms,2); % 0.5 * nlmax
  end

  if isreal(alms)
    alms = complex(alms, zeros(size(alms)));
  end

  if size(alms,1) ~= 1 && size(alms,1) ~= 3
    error('expected alms for 1 or 3 Stokes parameters')
  end

  if size(alms,4) ~= 1
    error('nmaps > 1 is not yet supported')
  end

  map = alm2map_c(alms, int32(nside));
end

