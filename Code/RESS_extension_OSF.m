%% This script accompanies the manuscript:
%    "Attention differentially modulates the amplitude 
%     of resonance frequencies in the visual cortex" 2019
%     from Gulbinaite, Roozendaal, and VanRullen
% 
% Following additional files are needed to run this code: 
%             - filterFGx.m
%             - eeglab toolbox (for plotting RESS filter weights using EEGLAB topoplot function)
%%
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

% Load the file
load('...\Subj1.mat');
outfilename = '...\Subj1_RESS.mat';

%% Conditions
condition_labels = {'LEFT', 'RIGHT'};
freq_titles = {'3','4','5','6','7','7_5','8','8_5','9','9_5','10','10_5','11','11_5','12','12_5','13','15','17','19','21','23','25','27','29','31','33','35','37','39','41','43','45','47','53','55','60','65','70','75','80'};

frex = [3:1:7 7.5:0.5:13 15:2:47 53 55 60:5:80];    % flicker frequencies
freqlist = 201:1:length(frex)+200;                  % EEG markers 

% condition EEG markers (cue / LED state / validity):
% LEFT/ FLICKER/ VALID     - 193
% LEFT/ FLICKER/ INVALID   - 194
% LEFT/ STATIC/  VALID     - 195
% LEFT/ STATIC/  INVALID   - 196

% RIGHT/ FLICKER/ VALID    - 197
% RIGHT/ FLICKER/ INVALID  - 198
% RIGHT/ STATIC/  VALID    - 199
% RIGHT/ STATIC/  INVALID  - 200

condlist = [193 194 199 200; 195:198];              % for creating RESS filters
attcondL = [193 194; 199 200];                      % attend left; ignore left;
attcondR = [197 198; 195 196];                      % attend right; ignore right

% Matrix of harmonics - used to select trials that do not contain flicker
% frequency of interest and its harmonics
harm_mat = zeros(length(frex),30);
for freqi = 1: length(frex)
    tmp = frex(freqi):frex(freqi):80;
    harm_mat(freqi,1:length(tmp)) = tmp;
end
%% RESS filter settings 
peakwidt  = 0.6; % FWHM at peak frequency
neighfreq = 1;   % distance of neighboring frequencies away from peak frequency, +/- in Hz
neighwidt = 1;   % FWHM of the neighboring frequencies
regu      = .01; % 1% shrinkage regularization of R covariance matrices, 
                 % reduces the influence of noise on the resulting eigendecomposition   
time4ress = [1500 7000];
tidx      = dsearchn(EEG.times',time4ress');
       
% FFT settings
nFFT  = ceil( EEG.srate/.1 );
hz = EEG.srate/2*linspace(0,1,nFFT/2+1);

% Parameters for computing SNR
numbins  = 10; % +-1 Hz
skipbins = 5; % skip 0.5 Hz
%% Initialize matrices to store the results
[ress ress_noise] = deal(zeros(EEG.pnts,EEG.trials));
[resstopos resstopos_noise] = deal(zeros(size(condlist,1),length(freqlist),EEG.nbchan));

%% Get the markers
EEGcueMarks = zeros(1,EEG.trials);
EEGstimMarks =  zeros(1,EEG.trials);
for ei = 1:length(EEG.epoch)
    [~,zeroloc] = min(abs(cell2mat(EEG.epoch(ei).eventlatency)-0));
    EEGcueMarks(ei) = EEG.epoch(ei).eventtype{zeroloc};
    EEGstimMarks(ei) = EEG.epoch(ei).eventtype{zeroloc+1};
end

% Go to double precision - important for eigval decomposition
EEG.data = double(EEG.data);

%% MAIN LOOP: 
% RESS filters are created separately for LEFT flicker (condi = 1) and RIGHT flicker (condi= 2)

for condi = 1:2
    for freqi = 1:length(frex)
        
        % Select trials based on codition (i.e. flicker frequency)
        freqtrials = find(EEGstimMarks == freqlist(freqi) &  ismember(EEGcueMarks,condlist(condi,:)));
        subdata = EEG.data(:,:,freqtrials);
        
        %% get covariance for flicker frequency
        
        % filter at flicker frequency
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi),peakwidt,[]);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covAt = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% get covariance for below flicker frequency
        
        % filter at neighboring frequency -1 Hz
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi)-neighfreq,neighwidt);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covLo = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% get covariance for above flicker frequency
        
        % filter at neighboring frequency +1 Hz
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi)+neighfreq,neighwidt);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covHi = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% Do generalized eigen value decomposition
        covR = (covHi+covLo)/2;
        [evecs,evals] = eig(covAt,covR + regu*diag(diag(covR)));
        
        % next three lines sorts vectors/values ascending (max component is last column)
        [~,sidx] = sort(diag(evals),1,'descend');
        evecs = real(evecs(:,sidx)); % take only real part in case complex
        evals = diag(real(evals));
        evals = evals(sidx);
        evecs_norm = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors
        
        %% Find component with highest SNR at the flicker frequency
        
        % reconstruct top 6
        tmpresstx = zeros(6,EEG.pnts,length(freqtrials));
        for triali=1:length(freqtrials)
            tmpresstx(:,:,triali) = ( subdata(:,:,triali)' * evecs_norm(:,1:6) )';
        end
        
        % power spectra
        tmpresspwr = mean(abs(fft(tmpresstx(:,tidx(1):tidx(2),:),nFFT,2)).^2,3);
        
        % compute SNR at flicker frequency
        freqidx = dsearchn(hz',frex(freqi));
        denom = squeeze(mean(tmpresspwr(:,[freqidx-numbins-skipbins:freqidx-skipbins freqidx+skipbins:freqidx+skipbins+numbins]),2));
        numer = squeeze(tmpresspwr(:,freqidx));
        tmpsnr = numer./denom;
        [~,comp2plot] = max(tmpsnr);

        %% Apply spatial filter to data
        evecs_norm = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors
        
        % get maps
        topos = zscore( covAt * evecs_norm );
        
        % force sign of components
        for ci=1:EEG.nbchan
            [~,idx] = max(abs(topos(:,ci)));                 % find strongest weight
            topos(:,ci) = topos(:,ci) * sign(topos(idx,ci)); % force to positive sign
        end
        
        resstopos(condi,freqi,:) = topos(:,comp2plot);
        
        % apply filter to data
        for triali=1:length(freqtrials)
            ress(:,freqtrials(triali)) = ( subdata(:,:,triali)' * evecs_norm(:,comp2plot) )';
        end
        
        %% Plot
        figure
        subplot(221)
        plot(evals,'-o')
        title([num2str(frex(freqi)) 'Hz' condition_labels{condi}])
        
        subplot(222)
        topoplot(double(squeeze( real( topos(:,comp2plot)))),EEG.chanlocs,'plotrad',.75,'numcontour',0,'electrodes','off');
        title(['Comp:' num2str(comp2plot)])
        
        %% Plot RESS component power spectrum 
        
        % Do FFT & compute SNR
        fft_tmp_av = mean(abs(fft(ress(tidx(1):tidx(2),freqtrials),nFFT)/diff(tidx)).^2,2);
        fft_tmp_av = squeeze(fft_tmp_av(floor(1:nFFT/2+1)));
        
        for hzi=numbins+1:length(hz)-numbins-1
            numer = fft_tmp_av(hzi);
            denom = mean( fft_tmp_av([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
            fft_tmp_av(hzi) = numer./denom;
        end
        
        % Plot
        subplot(2,2,[3,4])
        plot(hz,fft_tmp_av,'ro-','Linewidth',1.5);
        hold on
        plot(hz,fft_tmp_av,'bo-','Linewidth',1.5);
        set(gca,'xlim',[frex(freqi)-4 frex(freqi)+4])
        ylabel('SNR')

        %% =============================================================
        % Building FREQUENCY TAGGING null hypothesis distribution:
        %
        % To determine SNR values that can be expected due to overfitting
        % the noise and to obtain a statistical estimate of robustness of
        % flicker-induced EEG activity (null hypothesis distribution of
        % frequency tagging), for each flicker frequency and hemifield
        % we create RESS spatial filters using trials, in which neither
        % flicker frequency nor harmonically related responses were
        % present.
        %==============================================================
        
        freqs2use = frex;
        freqs2use(ismember(freqs2use,harm_mat(freqi,:) )) = [];
        freqs2use_idx = dsearchn(frex',freqs2use');
        
        freqtrials_noise = find(ismember(EEGstimMarks,freqlist(freqs2use_idx)) &  ismember(EEGcueMarks,condlist(condi,:)));
        subdata = EEG.data(:,:,freqtrials_noise);
        
        %% get covariance for flicker frequency
        
        % filter at frequency
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi),peakwidt);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covAt = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% get covariance for below flicker frequency
        
        % filter at frequency
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi)-neighfreq,neighwidt);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covLo = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% get covariance for above flicker frequency
        
        % filter at frequency
        filtTmp = filterFGx(subdata,EEG.srate,frex(freqi)+neighfreq,neighwidt);
        
        % get flicker time periods and compute covariance
        filtTmp = reshape(filtTmp(:,tidx(1):tidx(2),:),EEG.nbchan,[]);
        filtTmp = bsxfun(@minus,filtTmp,mean(filtTmp,2));
        covHi = (filtTmp*filtTmp') / (diff(tidx)-1);
        
        %% Do generalized eigen value decomposition
        covR = (covHi+covLo)/2;
        [evecs,evals] = eig(covAt,covR + regu*diag(diag(covR)));
        
        % next three lines sorts vectors/values ascending (max component is last column)
        [~,sidx] = sort(diag(evals),1,'descend');
        evecs = real(evecs(:,sidx)); % take only real part in case complex
        evals = diag(real(evals));
        evals = evals(sidx);
        evecs_norm = bsxfun(@rdivide,evecs,sqrt(sum(evecs.^2,1))); % normalize vectors
        
        [~,comp2plot] = max(evals);
        
        % get maps
        topos = zscore( covAt * evecs_norm ) ;
        for ci=1:EEG.nbchan
            [~,idx] = max(abs(topos(:,ci)));                 % find strongest weight
            topos(:,ci) = topos(:,ci) * sign(topos(idx,ci)); % force to positive sign
        end
        
        resstopos_noise(condi,freqi,:) = topos(:,comp2plot);
        
        % apply filter to data
        for triali=1:length(freqtrials_noise)
            ress_noise(:,freqtrials_noise(triali)) = (squeeze(EEG.data(:,:,freqtrials_noise(triali)))'*evecs_norm(:,comp2plot))';
        end
        
    end
end

%% save RESS results
save(outfilename,'resstopos','ress','ress_noise','resstopos_noise','hz','nFFT','tidx','time4ress','EEGstimMarks','EEGcueMarks','condition_labels','peakwidt','neighfreq','neighwidt')
