%% correction by cluster-based permutation test (EEG or Eye)
%(1)find clusters from EEG-IEM or eye classification timecourse
%(2)correct significance based on null distribution of cluster stat
%(3)plot EEG/Eye timecourse with corrected sig. bars

clear all
dataset = 'EEG'; %or 'eye'

%load data for IEM reconstruction of EEG, or eye results
load E:/2backrr_results/cvDR.mat
%load E:/2backrr_eyelink/acc_by_sub_2B_UMI_hampel.mat %for eye data

t = 4; %which task
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
rand('twister',sum(100*clock));

if strcmp(dataset,'EEG')
    bf=[0.0156    0.4219    1.0000    0.4219    0.0156  0];% basis function
    pnts = [60 188 376 148];% number of timepoints for each task
    
    %if the reconstruction is of N-1 in 2back -- only one ISI
    %pnts(3) = 188;
    
    nb_win =80; %window for average data for each timepoint/smoothing window
    nbi_t =1:pnts(t); %timepoints of task
    
    %correlate channel output with basis function across timepoints
    for nbi = 1:length(nbi_t)-nb_win/20
        for isub=1:length(subjects)
            bfcorr(nbi,isub)=corr(bf', squeeze(mean(Orichan_outp(isub,nbi_t(nbi):nbi_t(nbi)+nb_win/20,:),2)));
        end
    end
    %Fisher-z transform bfcorr and conduct 2-tailed t test
    for nbi = 1:length(nbi_t)-nb_win/20
        bfcorri = 0.5*(log(1+bfcorr(nbi,:))-log(1-bfcorr(nbi,:)));
        [h(nbi),p(nbi),ci{nbi},stats{nbi}] = ttest(bfcorri,0,'dim',2);
        tstat(nbi) = stats{nbi}.tstat;
    end
    
elseif strcmp(dataset,'eye')
    
    tpnum = size(acc_allsub,2);
    for tp = 1:tpnum
        acc_perm  = acc_allsub(:,tp)-1/6;
        [h(tp),p(tp),ci{tp},stats{tp}] = ttest(acc_perm',0,'dim',2);
        tstat(tp) = stats{tp}.tstat;
    end
end

%%%%%%%%%%%%%%%% find clusters, separate clusters for pos. and neg. t's %%%%%%%%%%%%%%%%
%first delay
%i = 1: +1 -> +new | -1 -> -new | 0 -> nothing
%i > 1:  0/+1, -1 -> -new | -1, -1 -> -add
%        0/-1, +1 -> +new | +1, +1 -> +add
%        0/1, 0 -> nothing

%(see end of script if 2-back with 2-epoch trials)
k_pos = 0; %cluster index
k_neg = 0;
clst_pos = [];
clst_neg = [];
if h(1) == 1
    if tstat(1) > 0
        k_pos = k_pos+1;
        clst_pos(k_pos,1) = 1; %1st column: cluster size
        clst_pos(k_pos,2) = 1;%2nd column: starting timepoint of cluster
        clst_pos(k_pos,3) = tstat(1); %3rd column: sum of t-values in cluster
    else
        k_neg = k_neg+1;
        clst_neg(k_neg,1) = 1; %1st column: cluster size
        clst_neg(k_neg,2) = 1;%2nd column: starting timepoint of cluster
        clst_neg(k_neg,3) = tstat(1); %3rd column: sum of t-values in cluster
    end
end


for i = 2:length(h) %length(h) is end of time
    if h(i) == 1
        if tstat(i) > 0
            if h(i-1) == 1 && tstat(i-1)>1
                clst_pos(k_pos,1) = clst_pos(k_pos,1) + 1;
                clst_pos(k_pos,3) = clst_pos(k_pos,3) + tstat(i);
            else
                k_pos = k_pos+1;
                clst_pos(k_pos,1) = 1;
                clst_pos(k_pos,2) = i;
                clst_pos(k_pos,3) = tstat(i);
            end
        else %tstat(i) < 0
            if h(i-1) == 1 && tstat(i-1)<0
                clst_neg(k_neg,1) = clst_neg(k_neg,1) + 1;
                clst_neg(k_neg,3) = clst_neg(k_neg,3) + tstat(i);
            else
                k_neg = k_neg+1;
                clst_neg(k_neg,1) = 1;
                clst_neg(k_neg,2) = i;
                clst_neg(k_neg,3) = tstat(i);
            end
        end
    end
end

%%%%%%%%%%%%%%%% load up the cluster stat distributions %%%%%%%%%%%%%%%%
cd E:/2backrr_results
load('eeg_clst_stat_cvDR_fin.mat');
load('cvDR.mat');
%for eye, for example
%load E:/2backrr_eyelink/clst_stat_umi.mat
%load E:/2backrr_eyelink/acc_by_sub_2B_UMI_hampel.mat


%%%%%%%%%%%%%%%% calculate significance %%%%%%%%%%%%%%%%
clst_stat=abs(clst_stat);

%for positive clusters
for i = 1: size(clst_pos,1) %loop though clusters
    clst_pos(i,3) = sum(tstat(clst_pos(i,2): clst_pos(i,2)+clst_pos(i,1)-1));%sum of t-values in each cluster
    clst_pos(i,4) = numel(find(clst_stat > clst_pos(i,3)))/10000; %proportion of clusters above observed
    if clst_pos(i,4) <.025
        clst_pos(i,5) = 1; %significance of cluster
    else
        clst_pos(i,5) = 0;
    end
end

%for negative clusters
for i = 1: size(clst_neg,1)
    clst_neg(i,3) = sum(tstat(clst_neg(i,2): clst_neg(i,2)+clst_neg(i,1)-1));
    clst_neg(i,4) = numel(find(clst_stat > abs(clst_neg(i,3))))/10000;
    if clst_neg(i,4) <.025
        clst_neg(i,5) = 1;
    else
        clst_neg(i,5) = 0;
    end
end

%%%%%%%%%%%%%%%% loop through all clusters and adjust h's %%%%%%%%%%%%%%%%
for i = 1: size(clst_pos,1)
    if clst_pos(i,5) == 0
        h(clst_pos(i,2):clst_pos(i,2)+clst_pos(i,1)-1) = 0;
    end
end
for i = 1: size(clst_neg,1)
    if clst_neg(i,5) == 0
        h(clst_neg(i,2):clst_neg(i,2)+clst_neg(i,1)-1) = 0;
    end
end

%% plot EEG timecourse with corrected sig bars

h_plot = [0 0 h 0 0]; %center in the smoothed window
end_tp = [1000 3550 7300 2750];

% %if the reconstruction is of N-1 in 2back -- only one ISI
%     end_tp(3) = 3550;

load cvDR.mat

%duplicate one channel to make curve symmetrical
Orichan_mat = nan(length(subjects),pnts(t),7);
Orichan_mat(:,:,2:7) = Orichan_outp;
Orichan_mat(:,:,1) = Orichan_outp(:,:,6);

figure;
imagesc(linspace(-200,end_tp(t),pnts(t)),1:7,squeeze(mean(Orichan_mat,1))');
colorbar

hold on;
plot(linspace(-200,end_tp(t),pnts(t)),(h_plot*-1)+8,'sw','markerfacecolor','w','markersize',6) %Red squares if neg. reconstruction
title(['DR cross validation'])
set(gca,'yticklabel',{'-90','-60','-30','0','30','60','90'});
ylabel('Channel')
xlabel('Time(ms)')
set(gca,'fontsize',12)

%% plot Eye timecourse with corrected sig bars
job = 1; %1 = 2B-UMI, 2 = 2B-PMI, 3 = 1B, 4 = DR
t = 3;
fignum = 1;

trial_length = [1000 3550 3550 2750];
times =-200:4:trial_length(t);
if job == 4
    times = times(1:end-1); %only for DR
end
acc_avgsub = mean(acc_allsub,1);
tpnum = size(acc_allsub,2);

options.handle = figure(fignum);
if job == 1
    options.color_area = [128 193 219]./255;    % Blue theme
    options.color_line = [ 52 148 186]./255;
else
    options.color_area = [243 169 114]./255;    % Orange theme
    options.color_line = [236 112  22]./255;
end
options.x_axis = times;
options.alpha      = 0.5;
options.line_width = 2;
options.error      = 'sem';
fig = plot_areaerrorbar(acc_allsub,options);

hold on
plot(times,ones(1,tpnum)*0.16667,'k--','LineWidth',1.3);
plot(zeros(101,1),0.164:0.0001:0.174,'b--','LineWidth',1);
ylim([0.164 0.174]);
xlim([-200 trial_length(t)]);

hplot = h;
hplot(h==0) = NaN;
if job == 2
    plot(times,hplot*0.1655,'Color',options.color_line,'LineWidth', 5);
else
    plot(times,hplot*0.165,'Color',options.color_line,'LineWidth', 5);
end
title('2-Back')
xlabel('Time (ms)');
ylabel('Posterior Probability');
%legend(H(1:2),'UMI','PMI'); %H keeps figure handles (fig) for UMI and PMI

if job == 4 %plot vertical line for probe onset in DR
    plot(ones(101,1)*2250,0.164:0.0001:0.174,'b--','LineWidth',1);
end

%%% Color the sig. neg. clusters in 1-back black
% if job == 3
%     hplot_neg = hplot;
%     hplot_neg(51:end) = NaN;
%     plot(times,hplot_neg*0.165,'Color','k','LineWidth', 5);
% end
%% if 2back task (with discontinuity between the two epochs)
% identify clusters separately for each epoch (i.e. after n and after n+1)

%     k_pos = 0; %cluster index
%     k_neg = 0;
%     clst_pos = [];
%     clst_neg = [];
%
%     % first epoch
%     if h(1) == 1
%         if tstat(1) > 0
%             k_pos = k_pos+1;
%             clst_pos(k_pos,1) = 1; %1st column: cluster size
%             clst_pos(k_pos,2) = 1;%2nd column: starting timepoint of cluster
%             clst_pos(k_pos,3) = tstat(1); %3rd column: sum of t-values in cluster
%         else
%             k_neg = k_neg+1;
%             clst_neg(k_neg,1) = 1; %1st column: cluster size
%             clst_neg(k_neg,2) = 1;%2nd column: starting timepoint of cluster
%             clst_neg(k_neg,3) = tstat(1); %3rd column: sum of t-values in cluster
%         end
%     end
%
%
%     for i = 2:188  %length(h) is end of time
%         if h(i) == 1
%             if tstat(i) > 0
%                 if h(i-1) == 1 && tstat(i-1)>0
%                     clst_pos(k_pos,1) = clst_pos(k_pos,1) + 1;
%                     clst_pos(k_pos,3) = clst_pos(k_pos,3) + tstat(i);
%                 else
%                     k_pos = k_pos+1;
%                     clst_pos(k_pos,1) = 1;
%                     clst_pos(k_pos,2) = i;
%                     clst_pos(k_pos,3) = tstat(i);
%                 end
%             else %tstat(i) < 0
%                 if h(i-1) == 1 && tstat(i-1)<0
%                     clst_neg(k_neg,1) = clst_neg(k_neg,1) + 1;
%                     clst_neg(k_neg,3) = clst_neg(k_neg,3) + tstat(i);
%                 else
%                     k_neg = k_neg+1;
%                     clst_neg(k_neg,1) = 1;
%                     clst_neg(k_neg,2) = i;
%                     clst_neg(k_neg,3) = tstat(i);
%                 end
%             end
%         end
%     end
%
%     %second epoch
%      if h(189) == 1
%         if tstat(189) > 0
%             k_pos = k_pos+1;
%             clst_pos(k_pos,1) = 1;
%             clst_pos(k_pos,2) = 189;
%             clst_pos(k_pos,3) = tstat(189);
%         else
%             k_neg = k_neg+1;
%             clst_neg(k_neg,1) = 1;
%             clst_neg(k_neg,2) = 189;
%             clst_neg(k_neg,3) = tstat(189);
%         end
%     end
%
%
%     for i = 190:length(h)  %length(h) is end of time
%         if h(i) == 1
%             if tstat(i) > 0
%                 if h(i-1) == 1 && tstat(i-1)>0
%                     clst_pos(k_pos,1) = clst_pos(k_pos,1) + 1;
%                     clst_pos(k_pos,3) = clst_pos(k_pos,3) + tstat(i);
%                 else
%                     k_pos = k_pos+1;
%                     clst_pos(k_pos,1) = 1;
%                     clst_pos(k_pos,2) = i;
%                     clst_pos(k_pos,3) = tstat(i);
%                 end
%             else %tstat(i) < 0
%                 if h(i-1) == 1 && tstat(i-1)<0
%                     clst_neg(k_neg,1) = clst_neg(k_neg,1) + 1;
%                     clst_neg(k_neg,3) = clst_neg(k_neg,3) + tstat(i);
%                 else
%                     k_neg = k_neg+1;
%                     clst_neg(k_neg,1) = 1;
%                     clst_neg(k_neg,2) = i;
%                     clst_neg(k_neg,3) = tstat(i);
%                 end
%             end
%         end
%     end