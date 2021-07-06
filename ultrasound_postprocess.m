%Have a look at ultrasound data

chan=find_files('CSC10*.ncs');
ref=find_files('CSC15*.ncs');
% [LFP_chan,~,sFreq] = nlx2matCSC_Matrix(chan{1});
% [LFP_ref,~,~] = nlx2matCSC_Matrix(ref{1});
downsample_fq=1000;

[full_timeLFP,~,sFreq] = nlx2matCSC_Matrix(chan{1}); %can you downsample after the fact in cheetah to fix the huge file issue?


% rereference

 [LFP_chan,~]=convert_dwnspl_detrend(chan{1},downsample_fq); 
 [LFP_ref,sFreq]=convert_dwnspl_detrend(ref{1},downsample_fq); %time will be in seconds

reref(:,1)=LFP_chan(:,1);
reref(:,2)=LFP_chan(:,2)-LFP_ref(:,2);
% LD_Load_Events
FieldSelection = [1 1 1 1 1];
ExtractHeader=1;
ExtractMode=1;
[eT, EventIDs, Ttls, Extras, EventStrings, Header] = Nlx2MatEV( 'Events.nev', FieldSelection, ExtractHeader, ExtractMode);

figure; plot(reref(:,1),reref(:,2))
figure; plot(LFP_chan(:,1),LFP_chan(:,2))

%try to see if this helps...
nfilt=Notch_filter_cowen(sFreq,59,61);
LFP(:,1)=full_timeLFP(:,1);
LFP(:,2)=filtfilt(nfilt,full_timeLFP(:,2));


figure; plot(full_timeLFP(:,1),full_timeLFP(:,2))

figure; plot(LFP(:,1),LFP(:,2))

%created event IDs so that you can put them in a matrix
evID=nan(Rows(EventStrings),2);

%Start Recording= 1
%End Recording =2
%TTL output =3 for trig
%injection =4

Index = find(contains(EventStrings,'Starting Recording'));
evID(Index,1)=1;
evID(Index,2)=eT(Index);

Index2 = find(contains(EventStrings,'Stopping Recording'));
evID(Index2,1)=2;
evID(Index2,2)=eT(Index2);

Index3 = find(contains(EventStrings,'TTL Output on AcqSystem1_0 board 0 port 2 value (0x0001)')); %there are some events that have the 0x0040 that i dont know what it is
evID(Index3,1)=3;
evID(Index3,2)=eT(Index3);

Index4 = find(contains(EventStrings,'injection'));
evID(Index4,1)=4;
evID(Index4,2)=eT(Index4);


Index5 = find(contains(EventStrings,'zip kick applied'));
evID(Index5,1)=5;
evID(Index5,2)=eT(Index5);

figure;
plot(full_timeLFP(:,1),full_timeLFP(:,2))
vline(evID(evID==1,2))
vline(evID(evID==3,2))

cutLFP=full_timeLFP( full_timeLFP(:,1)>evID(evID==3,2),:);

%% if you want to ginput some points and look at those
[x,~]=ginput(2);
ix=full_timeLFP>x(1) & full_timeLFP<x(2);

 [pxx,f] =pmtm(full_timeLFP(ix,2),5,[2:1:100],sFreq);
 figure;
 plot(f,pxx)

 
 
 d=Notch_filter_cowen(sFreq,59,61);
 LFP(:,1)=full_timeLFP(:,1);
 LFP(:,2)=filtfilt(d,full_timeLFP(:,2));
 
  [pxx2,f] =pmtm(LFP(ix,2),5,[2:1:100],sFreq);
 figure;
 plot(f,pxx2)
 
 %% try some other time points on the notched version
 [x,~]=ginput(2);
ix=LFP>x(1) & LFP<x(2);

 [pxx,f] =pmtm(LFP(ix,2),5,[2:1:100],sFreq);
 figure;
 plot(f,pxx)
 
 %% seems to work... time units are diff 
 [s,f,t]=spectrogram(LFP(:,2),500,300,[1:1:100],sFreq);
 figure;
 imagesc(t,f,10*log10(abs(s)))
 set(gca,'YDir','normal')

 %% wavelet
 
 [W,f] = cwt(LFP(ix,2),'bump',sFreq);
[wvt] = cwt_fix(W,f,[1:1:30]);
figure;
imagesc(LFP(ix,1),[1:1:30],wvt')
set(gca,'YDir','normal')
%%

endtime=LFP(end,1)-(2*60); %2 min before end
starttime=LFP(end,1)-(12*60);  %make 10 min range
ix=LFP(:,1)>starttime & LFP(:,1)<endtime;

 [pxx,f] =pmtm(cutLFP(ix,2),5,[2:1:100],sFreq);
 figure;plot(f,pxx)
 

 [pxx,f] =pmtm(LFP_chan(:,2),5,[1:1:100],sFreq);
 figure;
 plot(f,pxx)
 
 
 
[W,f] = cwt(LFP_chan(:,2),'bump',sFreq);
[wvt] = cwt_fix(W,f,[1:1:100]);
figure
imagesc(LFP_chan(:,1),[1:1:100],wvt)


ix=cutLFP(:,1)>1.619093e15 & cutLFP(:,1)< 1.619095e15;
figure; plot(cutLFP(ix,1),cutLFP(ix,2))


             [wvt] = cwt_fix(W,f,[1:1:100]);


% event_times=eT./10e6
%% trying to fuck around with before and after power shit
time_MK=205*60*2000;
time=length(full_timeLFP)/2000;
full_timeLFP(:,1)=linspace(0,time,length(full_timeLFP));
[pxx,f] =pmtm(full_timeLFP(5339397:end,2),5,[1:1:100],sFreq);
figure;
plot(f,pxx)
hold on
[pxx2,f] =pmtm(full_timeLFP(1:100000,2),5,[1:1:100],sFreq);
plot(f,pxx2)
