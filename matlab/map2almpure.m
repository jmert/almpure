function alms=map2almpure(map,apmask,lmax,mmax)
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
% OUTPUTS
%
%   alms     E and B alms for the decomposed maps. Dimensions will be
%            (2,lmax+1,mmax+1,nmaps).
%
% EXAMPLE
%

  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax
  end

  if ~any(ndims(map) == [2 3])
    error('map must have 2 or 3 dimensions')
  end
  if ~any(ndims(apmask) == [1 2])
    error('apmask must have 1 or 2 dimensions')
  end

  if size(map,2) ~= 2
    error('expected map for 2 Stokes parameters')
  end
  if size(map,3) ~= 1 || size(apmask,2) ~= 1
    error('nmaps > 1 is not yet supported')
  end

  if size(map,1) ~= size(apmask,1)
    error('map and apmask must have same number of pixels')
  end

  % Do not apodize the map. map2almpure_c internally does this
  % itself.
  alms = map2almpure_c(map, apmask, int32(lmax), int32(mmax));
end

