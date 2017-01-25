function [cl,dl]=alm2cl(alms1,alms2)
% [cl,dl]=alm2cl(alms1,alms2)
%
% Computes the C_l (and D_l) spectra given one or two sets of alms.
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
%   dl       Scaled version of cl, given by l(l+1)/2pi * cl
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

  if nargout > 1
    lmax = size(cl,1)-1;
    dl = bsxfun(@(l,cl) l.*(l+1)/(2*pi).*cl, cvec(0:lmax), cl);
  end
end

function cl=alm2cl_one(alms1,alms2)
  cl = zeros(size(alms1,1),1);
  for ll=1:size(alms1,1)
    ii = 1:min(ll, size(alms1,2));
    % Take just the real part. For auto-spectra, the imaginary components
    % cancel out anyway. For cross-spectra, the cross is defined by
    %   Cl = 1/(2l+1) sum_m { (alm1*conj(alm2) + conj(alm1)*alm2) / 2 }
    % which just gives the real part.
    cl(ll) = mean( real(alms1(ll,ii).*conj(alms2(ll,ii))) );
  end
end

