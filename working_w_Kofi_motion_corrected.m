load('F:\Keck Medicine of USC\Mice_MotionCorrected_Hippocampi_Data.mat')

figure; plot(hippoL_m3_1_0mg)


data=[hippoL_m1_1_0mg, hippoL_m2_1_0mg, hippoL_m3_1_0mg, hippoL_m4_1_0mg, hippoL_m5_1_0mg, hippoL_m7_1_0mg, hippoL_m8_1_0mg, hippoL_m9_1_0mg,...
    hippoL_m1_saline, hippoL_m2_saline, hippoL_m3_saline, hippoL_m4_saline, hippoL_m5_saline, hippoL_m7_saline, hippoL_m8_saline, hippoL_m9_saline];
MG1IX=logical([ones(8,1);zeros(8,1)]');
MG0IX=logical([zeros(8,1);ones(8,1)]');

d1=smoothdata(data,1,'movmedian');
std_thresh=2.5*std(d1,1);



baseix=5*60;
for imouse=1:Cols(data)
    badix=d1(:,imouse)>mean(d1(:,imouse))+std_thresh(imouse) | d1(:,imouse) < mean(d1(:,imouse))-std_thresh(imouse);
    d1(badix,imouse)=nan;
    d1(1:20,:)=nan;
    d1(3500:end,:)=nan;
base(imouse)=nanmedian(d1(20:baseix,imouse));
bl_norm(:,imouse)=((d1(20:3500,imouse)-base(imouse))./base(imouse)).*100;
end

%bl_norm col 7 is nuts- removing...
bl_norm(:,7)=nan;
%so is 12
bl_norm(:,12)=nan;
figure; plot_confidence_intervals(bl_norm(:,MG0IX)',[],[],[0 0 1]); hold on; plot_confidence_intervals(bl_norm(:,MG1IX)',[],[],[1 0 0]);
ylabel('% Change from BL')
xlabel('Time (s)')
pubify_figure_axis_robust

drugix=5*60+20*60;
changeix=2250;
bar_plot_with_p(nanmedian(bl_norm(drugix:end,:)),MG0IX,MG1IX, 'Av. change from basline, post 20 min injection')
xticklabels({'Saline' 'MK-801'})
ylabel('Median Change from BL 20min post injection')

bar_plot_with_p(nanmedian(bl_norm(changeix:end,:)),MG0IX,MG1IX, 'Av. change from basline')
xticklabels({'Saline' 'MK-801'})
ylabel('Median Change from BL 20min post injection')
%% try the ones that had no iso changes
%mouse 7,8,9 of 1 mg/kg had isochanges
%mouse 7 and 8 of saline had changes

data=[hippoL_m1_1_0mg, hippoL_m2_1_0mg, hippoL_m3_1_0mg, hippoL_m4_1_0mg, hippoL_m5_1_0mg,...
    hippoL_m1_saline, hippoL_m2_saline, hippoL_m3_saline, hippoL_m5_saline, hippoL_m9_saline];
MG1IX=logical([ones(5,1);zeros(5,1)]');
MG0IX=logical([zeros(5,1);ones(5,1)]');

%  d1=data;
 d1=smoothdata(data,1,'movmedian');
% std_thresh=2.5*nanstd(d1,1);


baseix=5*60;
for imouse=1:Cols(data)
%     badix=d1(:,imouse)>mean(d1(:,imouse))+std_thresh(imouse) | d1(:,imouse) < mean(d1(:,imouse))-std_thresh(imouse);
%     d1(badix,imouse)=nan;
    d1(1:20,:)=nan;
    d1(3500:end,:)=nan;
base(imouse)=nanmedian(d1(20:baseix,imouse));
bl_norm(:,imouse)=((d1(20:3500,imouse)-base(imouse))./base(imouse)).*100;
end

%bl_norm col 7 is nuts- removing...
% bl_norm(:,7)=nan;
figure; plot_confidence_intervals(bl_norm(:,MG0IX)',[],[],[0 0 1]); hold on; plot_confidence_intervals(bl_norm(:,MG1IX)',[],[],[1 0 0]);
ylabel('% Change from BL')
xlabel('Time (s)')
pubify_figure_axis_robust

drugix=5*60+20*60;
changeix=2250;
bar_plot_with_p(nanmean(bl_norm(end-(60*5):end,:)),MG0IX,MG1IX, 'Av. change from basline')
xticklabels({'Saline' 'MK-801'})
ylabel('Median Change from BL 20min post injection')

%% now trying a low pass filter
%mouse 7,8,9 of 1 mg/kg had isochanges
%mouse 7 and 8 of saline had changes

data=[hippoL_m1_1_0mg, hippoL_m2_1_0mg, hippoL_m3_1_0mg, hippoL_m4_1_0mg, hippoL_m5_1_0mg,...
    hippoL_m1_saline, hippoL_m2_saline, hippoL_m3_saline, hippoL_m5_saline, hippoL_m9_saline];
MG1IX=logical([ones(5,1);zeros(5,1)]');
MG0IX=logical([zeros(5,1);ones(5,1)]');

%  d1=data;
 d1=smoothdata(data,1,'movmedian');
% std_thresh=2.5*nanstd(d1,1);


baseix=5*60;
for imouse=1:Cols(data)
%     badix=d1(:,imouse)>mean(d1(:,imouse))+std_thresh(imouse) | d1(:,imouse) < mean(d1(:,imouse))-std_thresh(imouse);
%     d1(badix,imouse)=nan;
    d1(1:20,:)=nan;
    d1(3500:end,:)=nan;
base(imouse)=nanmedian(d1(20:baseix,imouse));
bl_norm(:,imouse)=((d1(20:3500,imouse)-base(imouse))./base(imouse)).*100;
end

%bl_norm col 7 is nuts- removing...
% bl_norm(:,7)=nan;
figure; plot_confidence_intervals(bl_norm(:,MG0IX)',[],[],[0 0 1]); hold on; plot_confidence_intervals(bl_norm(:,MG1IX)',[],[],[1 0 0]);
ylabel('% Change from BL')
xlabel('Time (s)')
pubify_figure_axis_robust

drugix=5*60+20*60;
changeix=2250;
bar_plot_with_p(nanmean(bl_norm(end-(60*5):end,:)),MG0IX,MG1IX, 'Av. change from basline')
xticklabels({'Saline' 'MK-801'})
ylabel('Median Change from BL 20min post injection')
