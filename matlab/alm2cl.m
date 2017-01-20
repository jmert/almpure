function cl=alm2cl(alms1,alms2)
% cl=alm2cl(alms1,alms2)
%
% Computes the Cl given one or two sets of alms.
%
% INPUTS
%   alms1    A set of alms. If 2D, dimensions should correspond to
%            (lmax+1 x mmax+1). If 3D, then dimensions are assumed to be
%            (N x lmax+1 x mmax+1) (where N is some number of spectra,
%            consistent with output of map2alm()).
%
%   alms2    Optional. If specified, the cross spectrum between alms1 and
%            alms2 is computed, otherwise, defaults to copying alms1.
%
% OUTPUTS
%
%   cl       Power spectra of given alms. Output has shape (Cl,spec).
%
% EXAMPLE
%
%   alms = map2alm(map, lmax, mmax);
%   cl   = alm2cl(alms);
%

  if ~exist('alms2','var') || isempty(alms2)
    alms2 = alms1;
  end

  if ndims(alms1) == 2
    alms1 = reshape(alms1, [1 size(alms1)]);
  end
  if ndims(alms2) == 2
    alms2 = reshape(alms2, [1 size(alms2)]);
  end

  cl = zeros(size(alms1,2), size(alms1,1));
  for ii=1:size(alms1,1)
    cl(:,ii) = alm2cl_one(squeeze(alms1(ii,:,:)), squeeze(alms2(ii,:,:)));
  end
end

function cl=alm2cl_one(alms1,alms2)
  cl = zeros(size(alms1,1),1);
  for ll=1:size(alms1,1)
    ii = 1:min(ll, size(alms1,2));
    cl(ll) = mean(alms1(ll,ii).*conj(alms2(ll,ii)));
  end
end

