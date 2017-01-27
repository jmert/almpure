function [cl,dl]=powspec(map,apmask,lmax,mmax)
% [cl,dl]=alm2cl(map,apmask,lmax,mmax)
%
% Computes the 6 (or 9) auto- and cross-spectra of the given T, Q, and U maps
%
% INPUTS
%
%   map      A 3-D array of full-sky map pixels (stored in ring ordering) for
%            three Stokes (T, Q, and U) parameters. The dimensions are expected
%            to be size (npix, 3, nmaps) where each of the values means:
%
%              npix     Number of pixels in a full-sky map. This is used to
%                       automatically determine the HEALPix nside.
%
%              nmaps    Number of distinct maps.
%
%   apmask   Apodization masks corresponding to the masks to be applied to the
%            maps. Size is expected to be (npix, 2, nmaps) where the Q and U
%            maps share a common apodization.
%
%   lmax     Maximum l-mode to decompose.
%
%   mmax     Maximum m-mode to decompose, where 0 <= nmmax <= nlmax. If not
%            given or empty, then nmmax = nlmax.
%
% OUTPUTS
%
%   cl
%
%   dl
%
% EXAMPLE
%

  if ~exist('mmax','var') || isempty(mmax)
    mmax = lmax;
  end

  if ~any(ndims(map) == [2 3])
    error('map must have 2 or 3 dimensions')
  end
  if ~any(ndims(apmask) == [2 3])
    error('apmask must have 2 or 3 dimensions')
  end

  if size(map,2) ~= 3
    error('expected map for 3 Stokes parameters')
  end

  if size(map,1) ~= size(apmask,1)
    error('map and apmask must have same number of pixels')
  end
  if size(apmask,3) ~= 1 && size(map,3)~=size(apmask,3)
    error('apmask must be shared with all or specified for each set of maps')
  end

  npix  = size(map,1);
  nmaps = size(map,3);
  nautos = nmaps;
  ncross = nmaps*(nmaps-1) / 2;
  nspecs = nautos + ncross;

  % Preallocate for speed
  cl = zeros(lmax+1, 9, nspecs);
  cltmp = zeros(lmax+1, 9);

  map(isnan(map)) = 0.0;
  apmask(isnan(apmask)) = 0.0;

  % Compute the map auto spectra
  for ii=1:nautos
    masksel = min(size(apmask,3), ii);

    almsT{ii} = map2alm(    map(:,1,  ii), apmask(:,1,masksel), lmax, mmax);
    almsP{ii} = map2almpure(map(:,2:3,ii), apmask(:,2,masksel), lmax, mmax);

    aT = squeeze(almsT{ii});
    aE = squeeze(almsP{ii}(1,:,:));
    aB = squeeze(almsP{ii}(2,:,:));

    % Store in diagonal + row order (same as cmbfast??)
    %
    % TT EE BB TE TB EB (ET BT BE)

    % TT
    cl(:,1,ii) = alm2cl( aT );
    % EE
    cl(:,2,ii) = alm2cl( aE );
    % BB
    cl(:,3,ii) = alm2cl( aB );
    % TE
    cl(:,4,ii) = alm2cl( aT, aE );
    % TB
    cl(:,5,ii) = alm2cl( aT, aB );
    % EB
    cl(:,6,ii) = alm2cl( aE, aB );
  end

  if ncross == 0
    cl = cl(:,1:6,1);
  end

  off = nautos + 1;
  for ii=1:nautos
    aT = squeeze(almsT{ii});
    aE = squeeze(almsP{ii}(1,:,:));
    aB = squeeze(almsP{ii}(2,:,:));

    for jj=(ii+1):nautos
      aT2 = squeeze(almsT{jj});
      aE2 = squeeze(almsP{jj}(1,:,:));
      aB2 = squeeze(almsP{jj}(2,:,:));

      % TT
      cl(:,1,off) = alm2cl( aT, aT2 );
      % EE
      cl(:,2,off) = alm2cl( aE, aE2 );
      % BB
      cl(:,3,off) = alm2cl( aB, aB2 );
      % TE
      cl(:,4,off) = alm2cl( aT, aE2 );
      % TB
      cl(:,5,off) = alm2cl( aT, aB2 );
      % EB
      cl(:,6,off) = alm2cl( aE, aB2 );
      % ET
      cl(:,7,off) = alm2cl( aT2, aE );
      % BT
      cl(:,8,off) = alm2cl( aT2, aB );
      % BE
      cl(:,9,off) = alm2cl( aE2, aB );

      off = off + 1;
    end
  end

  if nargout > 1
    dl = bsxfun(@(l,cl) l.*(l+1)/(2*pi).*cl, cvec(0:lmax), cl);
  end
end
