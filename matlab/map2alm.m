function alms=map2alm(map,lmax,mmax)
% alms=map2alm(map,lmax,mmax)
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
%   lmax     Maximum l-mode to decompose.
%
%   mmax     Maximum m-mode to decompose, where 0 <= nmmax <= nlmax. If not
%            given or empty, then nmmax = nlmax.
%
% OUTPUTS
%
%   alms     alms for the decomposed maps. Dimensions will be
%            (nstokes,lmax+1,mmax+1,nmaps).
%
% EXAMPLE
%

  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax
  end

  if ~any(ndims(map) == [2 3])
    error('map must have 2 or 3 dimensions')
  end

  if size(map,2) ~= 1 && size(map,2) ~= 3
    error('expected map for 1 or 3 Stokes parameters')
  end

  if size(map,3) ~= 1
    error('nmaps > 1 is not yet supported')
  end

  alms = map2alm_c(map, int32(lmax), int32(mmax));
end

