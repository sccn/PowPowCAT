% calc_PowPowCAT()--For batch-processing PowPowCAT.
%
% Usage: EEG = calc_PowPowCAT(EEG, upperFreqLimit, inputDataType, methodType, numIterations)
%
% Input: EEG -- EEGLAB data structure.
%        upperFreqLimit -- Frequency in Hz.
%        inputDataType -- 1, electrode data; 2, ICA time series.
%        methodType -- 1, Pearson's correlation; 2, Speaman's correlation.
%        numIterations -- Number of iteration for permutaion test to build
%                         the distribution of surrogate statistics.
%
% Output: EEG--EEGLAB data structure. Calculated results are stored under
%         EEG.etc.PowPowCAT.
%
% History
% 12/25/2020 Makoto. Batch version created for Pal Gunnar Larsson.

% Copyright (C) 2020, Makoto Miyakoshi (mmiyakoshi@ucsd.edu) , SCCN,INC,UCSD
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function EEG = calc_PowPowCAT(EEG, upperFreqLimit, inputDataType, methodType, numIterations)

% Compute PSD with 1-sec window with 50% overlap
% Natty ----------- Disabling and Changing to more linear log transform.
deviationFromLog = 5;
freqBins = logspace(log10(1+deviationFromLog), log10(upperFreqLimit+deviationFromLog), 100)-deviationFromLog; %5th
% -----------------

% Select input data.
switch inputDataType
    case 1
        inputData = EEG.data;
        inputLabel = 'channel';
    case 2
        inputData = EEG.icaact;
        inputLabel = 'IC';
end

% Display message if epoched data are detected. (12/04/2020 Makoto)
if size(inputData,3) == 1 % Continuous data.
    disp('Continuous data detected. Sliding window width is set to 1 s.')
else % Epoched data.
    warning('Epoched data detected. Sliding window width is set to the epoch length.')
    warning('Epoched data analysis is severely insensitive compared with continuous data analysis.')
end

progress('init','STFT started.');
for icIdx = 1:size(inputData,1)
    
    % Compute short-term Fourier Transform
    % Obtain the PSD length from the output of the first iteration (Thanks Pal Gunnar Lasson! 1/10/2017)
    if icIdx == 1;
        if size(inputData,3) == 1 % Continuous data. (12/04/2020 Makoto)
            [~, freqs, times, firstPSD] = spectrogram(inputData(icIdx,:), EEG.srate, floor(EEG.srate/2), freqBins, EEG.srate);
        else % Epoched data.
            [~, freqs, times, firstPSD] = spectrogram(inputData(icIdx,:), EEG.pnts,  0, freqBins, EEG.srate);
        end
        
        PSD = zeros(100, size(firstPSD,2), size(inputData,1));
        PSD(:,:,icIdx) = firstPSD;
    else
        if size(inputData,3) == 1 % Continuous data. (12/04/2020 Makoto)
            [~, ~, ~, PSD(:,:,icIdx)] = spectrogram(inputData(icIdx,:), EEG.srate, floor(EEG.srate/2), freqBins, EEG.srate);
        else % Epoched data.
            [~, ~, ~, PSD(:,:,icIdx)] = spectrogram(inputData(icIdx,:), EEG.pnts,  0, freqBins, EEG.srate);
        end
    end
    progress(icIdx/size(inputData,1));
end

% Remove the chunks that contain 'boundary' for the case of continuous data.
if size(inputData,3) == 1
    if isempty(EEG.event)
        boundaryIdx = [];
    else
        boundaryIdx = find(strcmp({EEG.event.type}, 'boundary'));
    end
    if any(boundaryIdx)
        boundaryLatencyInSec = [EEG.event(boundaryIdx).latency]*1/EEG.srate;
        [~, removeIdx] = histc(boundaryLatencyInSec, times-0.5); % 'times' are bin centers, so converting to bin edges.
    else
        removeIdx = [];
    end
    goodTimesIdx = setdiff(1:length(times), removeIdx);
    times = times(goodTimesIdx);
    PSD   = PSD(:,goodTimesIdx,:);
end

% Clean data by robust SD. (12/04/2020 Makoto)
switch methodType
    
    case 1
    disp('Pearson''s correlation is selected. 20% of data points are rejected for cleaning.')
    residualVariance2 = bsxfun(@minus, PSD, median(PSD,2));
    residualVariance_norm = bsxfun(@rdivide, residualVariance2, var(residualVariance2,0,2));
    residualVariance_norm_abs = abs(residualVariance_norm);
    residualVarianceSum = squeeze(sum(sum(residualVariance_norm_abs,3),1));
    goodChunkIdx2 = find(residualVarianceSum<=prctile(residualVarianceSum, 80));
    cleanPSD = PSD(:,goodChunkIdx2,:);
    clear PSD residualVariance
    % Compute cross-frequency coupling for each ICs (spectrum covariance across moving windows).
    covMatrix  = zeros(size(cleanPSD,1), size(cleanPSD,1), size(cleanPSD,3));
    for icIdx = 1:size(covMatrix,3)
        covMatrix(:,:,icIdx)  = corrcoef(cleanPSD(:,:,icIdx)');
    end
    
    case 2
    disp('Spearman''s correlation is selected. No data point rejected.')
    cleanPSD = PSD;
    % Compute cross-frequency coupling for each ICs (spectrum covariance across moving windows).
    covMatrix  = zeros(size(cleanPSD,1), size(cleanPSD,1), size(cleanPSD,3));
    for icIdx = 1:size(covMatrix,3)
        covMatrix(:,:,icIdx)  = corr(cleanPSD(:,:,icIdx)', 'type', 'Spearman');
    end
end

% Natty ----------- Adding
% Subsample the PSD time series into 200 points
cleanPsdPerm   = permute(cleanPSD, [2 1 3]);
cleanPsdPerm2D = cleanPsdPerm(:,:);
resampledPSD   = resample(cleanPsdPerm2D, 200, size(cleanPsdPerm2D,1));
resampledPsd3D = permute(reshape(resampledPSD, [200, 100, size(cleanPSD,3)]), [2 1 3]);
psdTimeSeriesNorm   = bsxfun(@rdivide, resampledPsd3D, std( resampledPsd3D,0,2));
psdTimeSeriesDemean = bsxfun(@minus,   psdTimeSeriesNorm, mean(psdTimeSeriesNorm,2));
clear cleanPsdPerm cleanPsdPerm2D resampledPSD resampledPsd3D psdTimeSeriesNorm
% -----------------

% Perform permutation test for each IC to save memory (12/09/2020 Makoto.)
processTimeList = zeros(size(covMatrix,3),1);
pvalMatrix      = zeros(size(covMatrix));
disp('Permutation test started. Please wait...')
progress('init','Permutation test started.');
for icIdx = 1:size(covMatrix,3)
    surroMatrix = single(zeros(size(cleanPSD,1), size(cleanPSD,1), numIterations));
    for iterationIdx = 1:numIterations
        
        permIdxMatrix = zeros(size(cleanPSD,1), size(cleanPSD,2));
        for freqIdx = 1:size(cleanPSD,1)
            permIdxMatrix(freqIdx,:) = randperm(size(cleanPSD,2));
        end
        surroPSD_currentIc = cleanPSD(:,:,icIdx);
        surroPSD = surroPSD_currentIc(permIdxMatrix);
        if     methodType==1
            surroMatrix(:,:,iterationIdx) = corrcoef(surroPSD');
        elseif methodType==2
            surroMatrix(:,:,iterationIdx) = corr(surroPSD', 'type', 'Spearman');
        end
    end
    
    % Perform Tim's non-parametric statistics.
    tmpCov   = covMatrix(:,:,icIdx);
    pvalMatrix(:,:,icIdx) = stat_surrogate_pvals(surroMatrix, tmpCov, 'both');
    
    progress(icIdx/size(covMatrix,3));
end

% FDR correction using all pixels from all ICs for global correction (this one is a very easy stats...)
pvalMask = pvalMatrix;
pvalMatrix_noDiagonal = pvalMatrix;
pvalMatrix_noDiagonal = bsxfun(@plus, pvalMatrix_noDiagonal, 2*eye(length(freqBins)));
nonDiagonalPvalues    = pvalMatrix_noDiagonal(pvalMatrix_noDiagonal<=1); % Must be integer multiple of 100^2-100
fdr005 = fdr(nonDiagonalPvalues, 0.05);
fdr001 = fdr(nonDiagonalPvalues, 0.01);

% Store the result.
EEG.etc.PowPowCAT.covMatrix  = covMatrix;
EEG.etc.PowPowCAT.freqs      = freqs;
EEG.etc.PowPowCAT.meanPSD    = squeeze(mean(cleanPSD,2));  % dB-converted
EEG.etc.PowPowCAT.stdPSD     = squeeze(std(cleanPSD,0,2)); % dB-converted
EEG.etc.PowPowCAT.pvalMatrix = pvalMatrix;
EEG.etc.PowPowCAT.fdr005 = fdr005;
EEG.etc.PowPowCAT.fdr001 = fdr001;
EEG.etc.PowPowCAT.inputData       = inputDataType;
EEG.etc.PowPowCAT.inputDataLabel  = {'EEG.data' 'EEG.icaact'};
EEG.etc.PowPowCAT.surroIter       = numIterations;
EEG.etc.PowPowCAT.methodType      = methodType;
EEG.etc.PowPowCAT.methodTypeLabel = {'Pearson' 'Spearman'};
EEG.etc.PowPowCAT.upperFreqLimit  = upperFreqLimit;
% Natty -----------
EEG.etc.PowPowCAT.timeSeriesPSD = psdTimeSeriesDemean;
% -----------------