%% Perform probabilistic classifier on processed eye data 

clear all

label = {'PMI'}; %or 'UMI'; if 2-back
t = 3; %which task

trial_length = [1000 3550 3550 2750]; %end time for each task
blocknum = [6 4 4 8]; %number of blocks for each task
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};

for s = 1:length(subjects)
    
    load(['F:\2backrr\2backRR_eye_processed\' subjects{s} '_eye_' tasks{t} '.mat']);
    
    %%%%%%%%%% start setting up data for cross-validation  %%%%%%%%%%
    if t == 3 %if 2back, adjust orimat on PMI/UMI
        switch label
            case 'PMI'
                orimat_all = orimat_all_pmi;
            case 'UMI'
                orimat_all = orimat_all_umi;
        end
    end
    
    %adjust block list for special cases
    if s == 12 && t == 3 %missing block, failed to transfer
        blocklist = 1:3;
    elseif s == 3 && t == 4 %block only has 5 stimOnset markers
        blocklist = [1:3 5:8];
    elseif s == 10 && t == 4 %04.17.20 added, first 2 DRs overwritten
        blocklist = 3:8;
    else
        blocklist = 1:blocknum(t);
    end
    
    acc_avg = zeros(1, size(datmat_all,3));
    
    for itest = blocklist
        trnblock = itest;
        tstblock = setdiff(blocklist,trnblock);
        
        tstdata = datmat_all(:,blockmat_all == itest,:); %1 is X, 2 is Y
        trndata = datmat_all(:,blockmat_all ~= itest,:);
        %%datmat_all is F * (trialnum*iblock) * time
        %tstdata/trndata: F * num_trial (N: tst/trn) * timepoints
        
        %orientation vector for training and testing data
        tstOri = orimat_all(blockmat_all == itest);
        trnOri = orimat_all(blockmat_all ~= itest);
        
        %%%%%%% start probabilistic classifer %%%%%%%%%%%%
        %train
        cfg0.demean = 'yes';
        cfg0.discardNan = 'yes';
        cfg0.gamma = 0.01;
        
        tpnum = size(tstdata,3); % this is the whole trial timecourse
        
        for tp = 1:tpnum
            trnX = trnOri; %N (labels)
            trnY = squeeze(trndata(:,:,tp)); %F * N, 2 * num_trial (N - trn)
            decoder(tp) = train_probClass(cfg0, trnX, trnY);
        end
        
        %test
        cfg0.demean = 'trainData';
        
        for tp = 1:tpnum
            tstY = squeeze(tstdata(:,:,tp));
            pPost(tp,:,:) = decode_probClass(cfg0, decoder(tp), tstY); %pPost = tp * C * N
        end
        
        % calculate correct classification averaged across trials and
        % classes
        for ori = 1:6
            acc(ori,:) = squeeze(nanmean(pPost(:,ori,tstOri == ori),3)); %acc = ori * tp
        end
        
        %over all folds
        acc_avg = acc_avg + nanmean(acc,1)/numel(blocklist); %number of folds
    end
    
    acc_allsub(s,:) = acc_avg; %accuracy over all timepoints for each subject
    
    acc_avgsub = nanmean(acc_allsub,1);
    
    if t == 3 %if 2-back
        save(['acc_by_sub_' tasks{t} '_' label '_hampel.mat'],'acc_allsub','acc_avgsub');
    else
        save(['acc_by_sub_' tasks{t} '_hampel.mat'],'acc_allsub','acc_avgsub');
    end
end
