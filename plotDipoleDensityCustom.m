% plotDipoleDensity(xyzList, FWHM_mm, sliceDimension)
%
% Uses 'MNI' coordinates as default (as EEG.dipfit.coordformat is).
% sliceDimension options are, 1==Axial, 2==Sagittal, 3==Coronal.
%
% History
% 12/19/2021 Makoto. Modified.
% 07/18/2020 Makoto. Created.

function plotDipoleDensityCustom(xyzList, FWHM_mm, sliceDimension)

% Prepare Fieldtrip input data structure.
source = [];
for dipIdx = 1:size(xyzList,1)
    source(dipIdx).posxyz = xyzList(dipIdx,1:3);
    source(dipIdx).momxyz = [0 0 0];
    source(dipIdx).rv     = 0;
end

% Prepare custom color.
customJet = jet(128);
customJet(1,:) = 0.3; % this is to make the background black
cmin = 0;
cmax = [];

% Prepare visualization option.
switch sliceDimension
    case 1
        plotargs = {'mriview', 'top',  'mrislices', -40:20:60, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax}; % For PowPowCAT.
    case 2
        plotargs = {'mriview', 'side', 'mrislices', -70:10:70, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax};
    case 3
        plotargs = {'mriview', 'rear', 'mrislices', -90:10:60, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax};
end

% Plot slices.
[dens3d, mri] = dipoledensity(source, 'coordformat', 'MNI', ...
                   'methodparam', FWHM_mm/2.355, 'plot', 'on', 'norm2JointProb', 'on', 'plotargs', plotargs);

% Display centroid xyz (SD).
figTitleCtString = sprintf('MNI [%.0f %.0f %.0f] (SD [%.0f %.0f %.0f])',...
                            mean(xyzList(:,1)), mean(xyzList(:,2)), mean(xyzList(:,3)),...
                            std(xyzList(:,1)),  std(xyzList(:,2)),  std(xyzList(:,3)));
set(gcf, 'Name', figTitleCtString, 'NumberTitle', 'off')