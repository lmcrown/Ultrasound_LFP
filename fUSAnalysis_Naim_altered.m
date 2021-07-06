function [Time, CleanData, Mean, MeanError] = fUSAnalysis_Naim (Hemisphere,Concentration,Mice,Baseline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function to analyze the fUS data exported from Icostudio
% The data must be in txt format with headers [Time, Right side, Time, Left
% side]
% Output:
%       Time: Time series in seconds
%       Concentration: Matrix of mean, filterd data, every column is a
%       mouse. Order is same as order in "Mice" input
%       Mean: Avergage of all the mice signal
%       MeanError: Error with the mean signal
% Input:
%       Hemisphere: takes two values (R,L) for selectin the side
%       Concentration: takes three values (0,0.1,1) for injection cons.
%       Mice: takes an array of the number of mice to be analzyed [1 2 3 4
%       10 11]
%       Baseline: time of baseline in minutes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','all')
disp('Reading .txt Data...')
% Read data from .txt file
for Index = 1:length(Mice)
    FileName = strcat(num2str(Concentration),'_Mouse_',num2str(Mice(Index)),'.txt');
    Temp =  table2array(readtable(FileName));
    Temp(:,3) = [];
    if Hemisphere == 'R'
        Time(:,Index) = Temp(:,1);
        Data(:,Index) = Temp(:,2);
    elseif Hemisphere == 'L'
        Time(:,Index) = Temp(:,1);
        Data(:,Index) = Temp(:,3);
    else
        disp('Wronge Hemisphere input')
        return
    end
end
disp('Calculating Baseline and Change...')
% Get baseline and precentage change
Baselinetime = Baseline*60;
for Index = 1:length(Mice)
    Baseline = mean(Data(1:Baselinetime,Index));
    ChangeData(:,Index) = ((Data(:,Index)-Baseline)./Baseline)*100;
end
disp('Cleaning the output data from spikes...')
% Clean Data using mean filter and wavelet
wname = 'sym4';
level = 4;
sorh  = 's';
thr = 19.5;
NumbeOfPoints = 75;
FilterCoeff = ones(1, NumbeOfPoints)/NumbeOfPoints;
for Index = 1:length(Mice)
     ChangeData(:,Index) = medfilt1( ChangeData(:,Index),25);
    CleanData(:,Index) = filter(FilterCoeff, 1, ChangeData(:,Index));
    CleanData(:,Index) = CleanData(:,Index);
    [CleanData(:,Index),~,~,~,~] = wdencmp('gbl',CleanData(:,Index),wname,level,thr,sorh,1);
end
disp('Calculating mean and mean error...')
% Get mean and Error
Mean = mean(CleanData,2);
MeanError = std(CleanData')'/sqrt(length(Mice));
disp('DONE')
end