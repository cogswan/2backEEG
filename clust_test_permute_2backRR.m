%% cluster-based permutation (generate null distribution)
%for EEG or eye data
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
    %     pnts(3) = 188;

    nb_win =80; %window for average data for each timepoint/smoothing window
    nbi_t =1:pnts(t); %timepoints of task
    
    %correlate channel output with basis function across timepoints
    for nbi = 1:length(nbi_t)-nb_win/20
        for isub=1:length(subjects)
            bfcorr(nbi,isub)=corr(bf', squeeze(mean(Orichan_outp(isub,nbi_t(nbi):nbi_t(nbi)+nb_win/20,:),2)));
        end
    end
end

%start 10000 permutations
for iperm = 1:10000
    if strcmp(dataset,'EEG')
        %Fisher-z transform bfcorr, and randomly flip the sign to create
        %null distribution
        for nbi = 1:length(nbi_t)-nb_win/20
            bfcorri = 0.5*(log(1+bfcorr(nbi,:))-log(1-bfcorr(nbi,:))).*(round(rand(1,length(subjects)))*-2+1);
            [h(nbi),p(nbi),ci{nbi},stats{nbi}] = ttest(bfcorri,0,'dim',2);
            tstat(nbi) = stats{nbi}.tstat;
        end
        
    elseif strcmp(dataset,'eye')
        tpnum = size(acc_allsub,2);
        %randomly flip sign to create null distribution
        for tp = 1:tpnum
            acc_perm = (acc_allsub(:,tp)-1/6).*(round(rand(1,length(subjects)))*-2+1)';
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
    
    
    %calculate and save largest (in absolute value) cluster-level stat for both pos.
    %and neg. clusters in each iteration
    
    if  isempty(clst_pos) && isempty(clst_neg)
        clst_stat(iperm,1) = 0;
    elseif isempty(clst_pos) %if only neg. clusters
        clst_stat(iperm,1) = min(clst_neg(:,3));
    elseif isempty(clst_neg) %if only pos. clusters
        clst_stat(iperm,1) = max(clst_pos(:,3));
    else                    %if both neg and pos clusters
        clst_stat_neg = min(clst_neg(:,3));
        clst_stat_pos = max(clst_pos(:,3));
        if abs(clst_stat_neg) > clst_stat_pos
            clst_stat(iperm,1) = clst_stat_neg;
        else
            clst_stat(iperm,1) = clst_stat_pos;
        end
    end
    
    iperm
end

%save cluster stat distribution
%if EEG data, for example, DR cross validation
cd E:/2backrr_results/
save('eeg_clst_stat_cvDR_fin.mat','clst_stat_pos','clst_stat_neg');

%if eye data, for example
%cd E:/2backrr_eyelink/
%save('clst_stat_umi.mat','clst_stat_pos','clst_stat_neg');

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