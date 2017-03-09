function alms=map2almpure(map,apmask,lmax,mmax,qwghts)
% alms=map2almpure(map,apmask,lmax,mmax)
%
% Decomposes the given map into Smith-style "pure" alms using the PS2HAT
% library
%
% INPUTS
%
%   map      A 3-D array of full-sky map pixels (stored in ring ordering) for
%            two Stokes (Q and U) parameters. The dimensions are expected to
%            be size (npix, 2, nmaps) where each of the values means:
%
%              npix     Number of pixels in a full-sky map. This is used to
%                       automatically determine the HEALPix nside.
%
%              nmaps    Number of distinct maps
%                       Note! nmaps ~= 1 is not yet supported.
%
%   apmask   Apodization masks of size (npix, nmaps) corresponding to the
%            mask to be applied to map.
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
%   alms     E and B alms for the decomposed maps. Dimensions will be
%            (2,lmax+1,mmax+1,nmaps).
%
% EXAMPLE
%

  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax;
  end

  if ~any(ndims(map) == [2 3])
    error('map must have 2 or 3 dimensions')
  end
  if ndims(apmask) ~= 2
    error('apmask must have 2 dimensions')
  end

  npix    = size(map, 1);
  nstokes = size(map, 2);
  nmaps   = size(map, 3);

  if nstokes ~= 2
    error('expected map for 2 Stokes parameters')
  end
  if nmaps ~= 1
    error('nmaps > 1 is not yet supported')
  end

  if npix ~= size(apmask,1)
    error('map and apmask must have same number of pixels')
  end
  if nmaps ~= size(apmask,2)
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

  % Do not apodize the map. map2almpure_c internally does this
  % itself.
  alms = map2almpure_c(map, apmask, int32(lmax), int32(mmax), qwghts);
end

