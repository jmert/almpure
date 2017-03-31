function alms=s2hat_map2alm(map,apmask,lmax,mmax,qwghts)
% alms=s2hat_map2alm(map,apmask,lmax,mmax,qwghts)
%
% Decomposes the given map into alms using the S2HAT library.
%
% INPUTS
%
%   map      A 3-D array of full-sky map pixels (stored in ring ordering) for
%            up to three Stokes parameters. The dimensions are expected to be
%            size (npix, nstokes, nmaps) where each of the values means:
%
%              npix     Number of pixels in a full-sky map. This is used to
%                       automatically determine the HEALPix nside.
%
%              nstokes  Number of Stokes parameters available. Should be
%                       either 1 or 3.
%
%              nmaps    Number of distinct maps
%                       Note! nmaps ~= 1 is not yet supported.
%
%   apmask   Apodization masks of size (npix, 1, nmaps) in which case the mask
%            is applied to all stokes parameters for each map or size
%            (npix, nstokes, nmaps) for independent masks. If empty, no extra
%            apodization is applied.
%
%   lmax     Maximum l-mode to decompose.
%
%   mmax     Maximum m-mode to decompose, where 0 <= nmmax <= nlmax. If not
%            given or empty, then nmmax = nlmax.
%
%   qwghts   Optional, defaults to all ones. Column vector of ring quadrature
%            weights. The length must be 2*nside (the number of rings in the
%            Norhtern Hemisphere, including the equator). Southern hemisphere
%            quadrature weights are assumed by symmetry.
%
% OUTPUTS
%
%   alms     alms for the decomposed maps. Dimensions will be
%            (nstokes,lmax+1,mmax+1,nmaps).
%
% EXAMPLE
%

  if ~exist('apmask','var') || isempty(apmask)
    apmask = ones(size(map));
  end
  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax;
  end

  if ~any(ndims(map) == [2 3])
    error('map must have 2 or 3 dimensions')
  end

  npix    = size(map, 1);
  nstokes = size(map, 2);
  nmaps   = size(map, 3);

  if nstokes ~= 1 && nstokes ~= 3
    error('expected map for 1 or 3 Stokes parameters')
  end

  if nmaps ~= 1
    error('nmaps > 1 is not yet supported')
  end

  if npix ~= size(apmask,1)
    error('map and apmask must have same number of pixels')
  end
  if nstokes ~= size(apmask,2) && size(apmask,2) ~= 1
    error('map and apmask have incompatible number of stokes parameters')
  end
  if nmaps ~= size(apmask,3)
    error('map and apmask must have the same number of maps')
  end

  nside   = sqrt(npix/12);

  if ~exist('qwghts','var') || isempty(qwghts)
    if 12*nside^2 ~= size(map,1)
      error('Could not determine the NSIDE of input map.')
    end
    nringsN = 2*nside;
    qwghts = ones(nringsN, 1);
  end
  % If we have nrings instead of (nrings+1)/2 entries, silently accept this
  % and truncate only if the vector is symmetric about the equator.
  if size(qwghts,1) == 4*nside-1
    nringsN = 2*nside;
    if ~all(qwghts(1:nringsN-1,1) == flipud(qwghts(nringsN+[1:nringsN-1],1)))
      error('Ring weights are not symmetric about the equator.')
    end
    qwghts = qwghts(1:nringsN,1);
  end

  if 2*nside ~= size(qwghts,1)
    error('qwghts must cover %d northern hemisphere rings', 2*nside)
  end
  if 1 ~= size(qwghts,2)
    error('qwghts must have dimension 2 of length 1')
  end

  if ~all(cvec(apmask == 1))
    % Apply the apodization mask before sending to map2alm_c, if one was
    % provided. This is mainly for symmetry with map2almpure which requires
    % an apodization mask be provided.
    map = bsxfun(@times, map, apmask);
  end

  alms = s2hat_map2alm_c(map, int32(lmax), int32(mmax), qwghts);
end

