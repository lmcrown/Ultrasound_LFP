function [OUT]=fus_EEG_postprocess(csc_name,PLOT_IT,downsample,Clean_LFP)
%saves artifact threshold- you could probably make a globals file to save
%this in
%Goal is to create a starter function that can be used to itterate and
%generate an average hippocampal PSD for saline and MK801 animals

%will want to downsample and detrend
%then remove artifacts
%then plot the PSD in an OUT file along with rat number, drug conditon
% clear all
OUT.aborted=false;
PLOT_IT=0;

if nargin<1
    
    parts=strsplit(pwd,'\');
    OUT.animal=parts{9}(1:4);
    OUT.drug=parts{7};
    OUT.day=parts{8};
    % head = Read_nlx_header('HIPP-CSC7_4586135492_End.ncs');
    % if csc=='hip'
    
    file=find_files('*HIPP*.ncs');
    if isempty(file) && str2double(OUT.animal)>=1042
        file=find_files('CSC8*.ncs'); %nancy says for all animals 1042 and on that csc 8 is hippocampus, 16 is mpfc
    elseif isempty(file) && str2double(OUT.animal)<=1042
        OUT.aborted=true;
        disp('no cortical file found')
        return
    end

else
    file={csc_name};
end
%CSC 4 thalamus, CSC 12 MSN, CSC 8 HIPP, 19 MPFC

%     num=extractAfter(file{1},'HIPP-CSC');  %seems like not all have the
%     spilit file stuff so I"m just taking the end of it anyways, so I
%     might as well just take the whole file
%     chan=num(1);
%     cd('Split-relabeled files\')
%     split_file=find_files(['CSC' chan '*']);
% end

if downsample==1
    [LFP,sFreq]=convert_dwnspl_detrend(file{1},downsample_fq); %second input is desired downsampled freq -NOTE THIS CHANGES TO SECONDS
else
[LFP,~,sFreq] = nlx2matCSC_Matrix(file{1});
end
% Pull out time stamps
[Timestamps, EventIDs, TTLs, Extras, EventStrings, Header] = Nlx2MatEV('Events.nev', [1 1 1 1 1], 1, 1, [] );
EVT={EventStrings, Timestamps'};
EVT=table(EventStrings, Timestamps');

% NOW IS WHERE YOU NEED STANDARDIZED NAMES

ix_Us=LFP(:,1)>EVT.Var2(8) & LFP(:,1)<EVT.Var2(13); %will need to match strings when this is standardized
inj=EVT.Var2(9);
stim_start=EVT.Var2(11);
stim_end=EVT.Var2(12);
stim_ix=LFP(:,1)>stim_start & LFP(:,1)< stim_end;
baseline_ix=LFP(:,1)>EVT.Var2(8) & LFP(:,1)<inj;


if PLOT_IT==1
figure;
plot(LFP(ix_Us,1),LFP(ix_Us,2))
vline(inj) %injection
vline(stim_start)
vline(stim_end)
end

notch=Notch_filter_cowen(sFreq,59,60);

notched_LFP=filtfilt(notch,LFP(:,2));
% LFP_i=LFP(ix,:);
% figure
% plot(LFP(ix,1),LFP(ix,2))
% [~,y]=ginput(1);
% OUT.art_thresh=y;
if Clean_LFP==1;
[BIX,artifact_times_usec] = LD_Clean_LFP(LFP(ix,:),[],6e4,downsample_fq);
perc_bad = sum(BIX)/length(BIX);
fprintf('BAD percent: %2.2f\n',perc_bad*100)
if perc_bad > .3 %if its grabbing 10 minutes this means that ill get 6
    disp('Too much bad data')
    OUT.aborted=true;
    return
end
LFP=LFP(~BIX,:) ;
end


%save the LFP snip that you have so that you don't need to download again


Theta_filt = designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',5, 'HalfPowerFrequency2',12, ...
    'SampleRate',sFreq,'designmethod', 'butter');
HighGamma = designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',130, 'HalfPowerFrequency2',180, ...
    'SampleRate',sFreq,'designmethod', 'butter'); %but he stims at 100hz i think so why not look at that?
LowGamma = designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',30, 'HalfPowerFrequency2',50, ...
    'SampleRate',sFreq,'designmethod', 'butter');
Delta=designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',0.5, 'HalfPowerFrequency2',3, ...
    'SampleRate',sFreq,'designmethod', 'butter');
Broad_gamma=designfilt('bandpassiir','FilterOrder',12, ...
    'HalfPowerFrequency1',30, 'HalfPowerFrequency2',120, ...
    'SampleRate',sFreq,'designmethod', 'butter');


% se=pentropy(newLFP(:,2),downsample_fq,'FrequencyLimits',[30 150]);
hgpow(:,1)=LFP(:,1);
lgpow(:,1)=LFP(:,1);
thetapow(:,1)=LFP(:,1);
deltapow(:,1)=LFP(:,1);
broadpow(:,1)=LFP(:,1);

hgpow(:,2)=abs(hilbert(filtfilt(HighGamma,LFP(:,2))));
lgpow(:,2)=abs(hilbert(filtfilt(LowGamma,LFP(:,2))));
thetapow(:,2)=abs(hilbert(filtfilt(Theta_filt,LFP(:,2))));
deltapow(:,2)=abs(hilbert(filtfilt(Delta,LFP(:,2))));
broadpow(:,2)=abs(hilbert(filtfilt(Broad_gamma,LFP(:,2))));

%smooth it
smooth_theta(:,1)=thetapow(ix_Us,1);
smooth_theta(:,2)=smoothdata(thetapow(ix_Us,2),'movmedian',10e3); 

smooth_hgamma(:,1)=hgpow(ix_Us,1);
smooth_hgamma(:,2)=smoothdata(hgpow(ix_Us,2),'movmedian',10e3);

%find some power values for gamma and theta percent change from basline
 OUT.base_theta=nanmean(smooth_theta(baseline_ix,2));
 OUT.base_hgamma=nanmean(smooth_hgamma(baseline_ix,2));

% stim_theta=nanmean(thetapow(stim_ix));
bl_norm_theta(:,1)=smooth_theta(:,1);
bl_norm_theta(:,2)=(100*((smooth_theta(:,2)-OUT.base_theta)./OUT.base_theta));
 
bl_norm_hgamma(:,1)=smooth_hgamma(:,1);
bl_norm_hgamma(:,2)=(100*((smooth_hgamma(:,2)-OUT.base_hgamma)./OUT.base_hgamma));
 
% figure; 
% plot(bl_norm_theta(:,1),bl_norm_theta(:,2))
%  figure; 
%  plot(bl_norm_hgamma(:,1),bl_norm_hgamma(:,2))



if PLOT_IT==1
figure; 
subplot 211
plot(((bl_norm_theta(:,1)-bl_norm_theta(1,1))./10e6),bl_norm_theta(:,2)) %This is baseline normalized and smoothed across milliseconds
vline((inj-bl_norm_theta(1,1))/10e6) %injection
vline((stim_start-bl_norm_theta(1,1))/10e6)
vline((stim_end-bl_norm_theta(1,1))/10e6)
xlabel('Seconds')
ylabel('Theta Power (% change from BL)')
pubify_figure_axis_robust

subplot 212
plot(((bl_norm_hgamma(:,1)-bl_norm_hgamma(1,1))./10e6),bl_norm_hgamma(:,2)) %This is baseline normalized and smoothed across milliseconds
vline((inj-bl_norm_hgamma(1,1))/10e6) %injection
vline((stim_start-bl_norm_hgamma(1,1))/10e6)
vline((stim_end-bl_norm_hgamma(1,1))/10e6)
xlabel('Seconds')
ylabel('High Gamma Power (% change from BL)')
pubify_figure_axis_robust
end
%Make plot of time series of theta power?
OUT.broad_pow=nanmean(broadpow);
OUT.hg_delt=nanmean(hgpow./deltapow);
OUT.lg_delt=nanmean(lgpow./deltapow);
OUT.theta_delt=nanmean(thetapow./deltapow);
OUT.broad_delt=nanmean(broadpow./deltapow);
OUT.raw_theta=nanmean(thetapow);
% LFP(:,2)=LFP_mv;
[pxx,f] =pmtm(LFP(ix,2),5,[2:0.5:120],sFreq);
[pxx_noart,f] =pmtm(newLFP(:,2),5,[2:0.5:120],sFreq);

new_num=mean([pxx_noart(f==59),pxx_noart(f==61)]);
pxx(f==60)=new_num;
new_num2=mean([pxx_noart(f==119),pxx_noart(f==121)]);
pxx_noart(f==120)=new_num2;
dbpxx=10*log10(pxx_noart);


b_coeffs=robustfit(f,dbpxx);

b_coeffs_lg=robustfit(f(f<55),dbpxx(f<55));
b_coeffs_hg=robustfit(f(f>65 &f<90),dbpxx(f>65 &f<90));

y_ep_lg=b_coeffs_lg(1)+(b_coeffs_lg(2)*f(f<55));
y_ep_hg=b_coeffs_lg(1)+(b_coeffs_hg(2)*f(f>65 &f<90));

OUT.resid_lg(:)=dbpxx(f<55)-y_ep_lg;
OUT.resid_hg(:)=dbpxx(f>65 &f<90)-y_ep_hg;
OUT.slope_lg=b_coeffs_lg(2);
OUT.slope_hg=b_coeffs_hg(2);


y_ep=b_coeffs(1)+(b_coeffs(2)*f);

OUT.resid=dbpxx-y_ep;

OUT.psd=pxx;
OUT.dbpsd=dbpxx;
OUT.freqs=f;
OUT.model=y_ep;
OUT.slope=b_coeffs(2);

%find theta frequency with instafreq
[thetaFreq,t]=instfreq(newLFP(:,2),sFreq,'FrequencyLimits',[5 10]); %instafreq seems to want to work in seconds
thetafreq=median(thetaFreq);
OUT.thetafreq=thetafreq;

[lowgammafreq,t]=instfreq(newLFP(:,2),sFreq,'FrequencyLimits',[30 50]); %instafreq seems to want to work in seconds
lowgammaf=median(lowgammafreq);
OUT.lowgammafrex=lowgammaf;
%low gamma= 30-50
%high gamma= 130-180
[highgammafreq,t]=instfreq(newLFP(:,2),sFreq,'FrequencyLimits',[130 180]); %instafreq seems to want to work in seconds
highgammaf=median(highgammafreq);
OUT.highgammafrex=highgammaf;

figure
plot(f,OUT.dbpsd)
hold on
plot(f,OUT.model)

figure;
plot(f,10*log10(pxx))
hold on
plot(f,10*log10(pxx_noart))
legend('original','no art');
if nargin==0
title(sprintf('Animal %s %s %s',OUT.animal, OUT.drug, OUT.day))
else
    title(sprintf('%s',csc_name))
end


%this would need to be changed depending on your path
OUT.pxx=pxx;
OUT.pxx_noart=pxx_noart;
OUT.frex=f;
LFPdata.lfp=newLFP;
% LFPdata.out=OUT;
%  place2save='E:\Darrin\Reduced_Files\PSD7_10min_pre'; %changing this for location
%  save([place2save '\Rat' OUT.animal],'LFPdata')

