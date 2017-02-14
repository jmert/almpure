function [bc,bp]=cl2bandpowers(cl,binedges)
% [bc,bp]=cl2bandpowers(cl,binedges)
%
% Averages the given C_ls according to the given bin edges (where bins are
% the left-inclusive half-open interval).
%
% INPUTS
%   cl          Input C_l spectra, where dimension 1 runs from 0 to some lmax
%
%   binedges    The boundaries of each bandpower bin.
%
% OUTPUTS
%   bc          Band centers for each bin
%
%   bp          The binned C_l, i.e. the bandpowers
%
% EXAMPLE
%
%   be = [0,20:35:580]; % BICEP binning
%   [~,dl] = powspec(map, apmask, lmax, lmax, true);
%   [bc,bp] = cl2bandpowers(dl, be);
%

  nbins = length(binedges) - 1;
  l = 0:size(cl,1);

  bc = zeros(nbins, 1);
  bp = zeros(nbins, size(cl,2));

  for ii=1:nbins
    bs = binedges(ii);
    be = binedges(ii+1) - 1;
    bc(ii) = 0.5 * (bs + be);
    bp(ii,:) = mean(cl(bs<=l&l<be,:), 1);
  end
end
