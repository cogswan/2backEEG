%% IEM training and testing

%% Initialization
clear all;

% Run EEGLAB
run('E:/eeglab2019_1/eeglab.m');

%% Build basis function (BF)
x_full = linspace(0,pi-(pi/180),180);
x = x_full(1:30:end);
x_deg = rad2deg(x);
bf = sin(x).^6;

bf = bf./max(bf);         % norm to have unit height
peak = find(bf==1);           % find the max point(s)
peak = peak(1);                 % only keep the index of the first point that =='s the max
bf = wshift('1D', bf, peak-1);    % shift it so that it peaks at x==0
%  figure;plot(x_deg,wshift('1D',bf,0)) %For testing purposes

%%  Downsample the data, train the early delay model (timepoint 58; 940-1040ms)
bf = [1.0000    0.4219    0.0156     0    0.0156    0.4219]; %basis function derived from above
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
t0 = 4; %training task

for sub = 1:length(subjects)
    
    EEG_train=pop_loadset('filename',[subjects{sub}, '_rmeye_' tasks{t0} '.set'], 'filepath',['E:/runica/rmeye_fin/']); %delayed recog.
    EEG_train = pop_resample(EEG_train,50);
    
    
    nChans = 6; % 6 orientation channels
    bin_Ori_trn= [EEG_train.event.ori]'; %orientation labels for training trials
    Oriblock = [EEG_train.event.blocknum]';
    
    %Remove first two DR blocks of S09 due to missing behav data
    if t0 == 4 && sub == 10
        EEG_train.data= EEG_train.data(:,:,Oriblock > 2);
        bin_Ori_trn = bin_Ori_trn(Oriblock > 2);
        Oriblock = Oriblock(Oriblock > 2);
    end
    
    trntrialinds = 1:size(EEG_train.data,3);
    
    %deriving channel weights for early delay model (940-1040ms; timepoint
    %58-63)
    if t0 == 1 %FL
        train_epoch = 25:34; %288-471ms
    elseif t0 == 4 %DR
        train_epoch = 58:63; %940-1040ms
    end
    
    for trn_tr = train_epoch
        trnOridata = squeeze(EEG_train.data(:,trn_tr,:));
        trn_Ori = trnOridata;               % training data
        Oritrng = bin_Ori_trn;              % vector of trial labels for training data
        OriX_pre = zeros(size(trn_Ori,2), nChans);
        
        %circularly shift tuning function to align at center (3rd channel)
        for ii=1:size(trn_Ori,2)
            OriX_pre(ii,:) = wshift('1D', bf, -(Oritrng(ii)-1));
        end
        Oriw_pre = OriX_pre\trn_Ori';
        Oriw_pre_all(trn_tr-train_epoch(1)+1,:,:) = Oriw_pre; %train_epoch(1) is 58 for DR
    end
    
    
    %save testing data and trained channel weights
    if t0 == 1
        cd E:/2backrr_results/FL_TrnW_288-471_fin %DR_TrnW_200-500_fin
        save([subjects{sub} '_FL_TrnW_288-471_fin'],'Oriw_pre_all');
    elseif t0 == 4
        cd E:/2backrr_results/DR_TrnW_940-1040_fin %DR_TrnW_200-500_fin
        save([subjects{sub} '_DR_TrnW_940-1040_fin'],'Oriw_pre_all');
    end
    clearvars -except sub subjects bf t0 t tasks
end

%%  Set up 2back test trials (concatenate adjacent 2 epochs)
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
t = 3;

for sub = 1:length(subjects)
    
    EEG_test=pop_loadset('filename',[subjects{sub}, '_rmeye_2B.set'], 'filepath',['E:/runica/rmeye_fin/']); %2-back
    EEG_test = pop_resample(EEG_test,50);
    
    %build the proper test trial, concatenate every 2 adjacent epochs
    for i=1:size(EEG_test.event,2) %remove first trial of each block
        if mod(EEG_test.event(i).urevent-1,128)==0
            idx(i)=0;
        else
            idx(i)=1;
        end
    end
    
    EEG_test_tmp.event=EEG_test.event(idx==1);
    EEG_test_tmp.data=EEG_test.data(:,:,idx==1);
    cnt = 0;
    for i=1:size(EEG_test_tmp.event,2)-1
        if EEG_test_tmp.event(i+1).urevent-EEG_test_tmp.event(i).urevent==1
            cnt=cnt+1;
            EEG_test_new.data(:,:,cnt)=cat(2,EEG_test_tmp.data(:,:,i), EEG_test_tmp.data(:,:,i+1));
            EEG_test_new.event(cnt)=EEG_test_tmp.event(i);
        end
    end
    EEG_test_data = EEG_test_new.data;
    bin_Ori_tst =  [EEG_test_new.event.ori]';
    Oriblock = [EEG_test_new.event.blocknum]';
    
    
    cd E:/2backrr_results/2b_2tr_dat_fin
    save([subjects{sub} '_2b_2tr_fin'], 'EEG_test_data','bin_Ori_tst','Oriblock');
    clearvars -except sub subjects bf t0 t tasks
end
%% Testing DR/FL trained IEMs on tasks
clear all
bf = [1.0000    0.4219    0.0156     0    0.0156    0.4219];
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
filepath ='E:/2backrr_results/';

t0 = 4;
if t0 == 1 %FL
    train_epoch = 25:34; %288-471ms
elseif t0 == 4 %DR
    train_epoch =58:63; %940-1040ms
end

t = 3; %test on 2-back
for sub = 1:length(subjects)
    
    if t0 == 4
        load (['E:/2backrr_results/DR_TrnW_940-1040_fin/' subjects{sub} '_DR_TrnW_940-1040_fin.mat']);
    elseif t0 == 1
        load (['E:/2backrr_results/FL_TrnW_288-471_fin/' subjects{sub} '_FL_TrnW_288-471_fin']);
    end
    
    %reconstruct 2-back from n-1
    %   load (['E:/2backrr_results/2b_next_tr_dat_fin/' subjects{sub} '_2b_next_tr_dat_fin']);
    
    %reconstruct from n
    if t == 3
        load (['E:/2backrr_results/2b_2tr_dat_fin/' subjects{sub} '_2b_2tr_fin']);
    else
        EEG_test = pop_loadset('filename',[subjects{sub}, '_rmeye_' tasks{t} '.set'], 'filepath',['E:/runica/rmeye_fin/']);
        EEG_test = pop_resample(EEG_test,50);
        EEG_test_data = EEG_test.data;
        bin_Ori_tst =  [EEG_test.event.ori]';
    end
    
    nChans = 6;
    
    tsttrialinds = 1:size(EEG_test_data,3);
    
    Orichan_outp_avg = zeros(size(EEG_test_data,2),6);
    
    for trn_ind = 1:numel(train_epoch)
        for tst_tr = 1:size(EEG_test_data,2)
            
            tstOridata = squeeze(EEG_test_data(:,tst_tr,:));
            
            tst_Ori = tstOridata;         % testing data
            
            Oritstg = bin_Ori_tst;        % trial labels for tst data.
            
            Oriw_pre = squeeze(Oriw_pre_all(trn_ind,:,:));
            
            Ori_chan_pre = (Oriw_pre'\tst_Ori)';
            
            for ii=1:size(bin_Ori_tst,1)
                Ori_chan_pre(ii,:) = wshift('1D', Ori_chan_pre(ii,:),  bin_Ori_tst(ii)-ceil(nChans/2));
            end
            
            Orichan_outp_avg(tst_tr,:) = Orichan_outp_avg(tst_tr,:) + (nanmean(Ori_chan_pre,1)/numel(train_epoch));
        end
    end
    
    Orichan_outp(sub,:,:) = Orichan_outp_avg;
    
    clearvars -except subjects sub bf Orichan_outp t tasks t0 train_epoch filepath
end
save ([filepath tasks{t0} tasks{t} '_288-471_fin.mat'],'Orichan_outp'); %if t0 == 1

%% train & test on same task/2B (leave one run out/cross validation)
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
t = 3; %2B
bf = [1.0000    0.4219    0.0156     0    0.0156    0.4219];
numblock = [6 4 4 8]; %number of blocks in each task
pnts = [60 188 376 148]; %number of timepoints in each task
Orichan_outp_fin = nan(length(subjects),pnts(t),pnts(t),6);

for sub = 1:length(subjects)
    if t == 3
        load (['E:/2backrr_results/2b_2tr_dat_fin/' subjects{sub} '_2b_2tr_fin']);
        EEG.data = EEG_test_data;
    else
        EEG_test = pop_loadset('filename',[subjects{sub}, '_rmeye_' tasks{t} '.set'], 'filepath',['E:/runica/rmeye_fin/']);
        EEG_test = pop_resample(EEG_test,50);
        EEG.data = EEG_test.data;
        bin_Ori_tst =  [EEG_test.event.ori]';
        Oriblock = [EEG_test.event.blocknum]'; %array of block number
    end
    
    %Remove first two DR blocks of S09 due to missing behav data
    if t == 4 && sub == 10
        EEG.data = EEG.data(:,:,Oriblock > 2);
        bin_Ori_tst = bin_Ori_tst(Oriblock > 2);
        Oriblock = Oriblock(Oriblock > 2);
        runs = 3:numblock(t);
    else
        runs = 1:numblock(t);
    end
    
    nChans = 6;
    
    bin_Ori_trn = bin_Ori_tst;
    
    trntrialinds = 1:size(EEG.data,3);
    tsttrialinds = 1:size(EEG.data,3);
    
    Orichan_outp_avg = nan(pnts(t),pnts(t),6);
    
    for trn_tr = 1:size(EEG.data,2)
        tst_tr = trn_tr; %this is Cross Validation
        %  for  tst_tr = 1:size(EEG.data,2) %this is Temporal Generalization
        
        trnOridata = squeeze(EEG.data(:,trn_tr,:));
        tstOridata = squeeze(EEG.data(:,tst_tr,:));
        
        Ori_chan_pre = nan(size(trntrialinds,2), nChans);
        
        for rr = runs
            fprintf('Computing iteration %d out of %d\n', rr, size(runs,2));
            looOriblock = Oriblock;
            
            trn_Ori = trnOridata(:,looOriblock ~= rr);         % data from training runs
            tst_Ori = tstOridata(:,looOriblock == rr);         % data from testing run
            
            Oritrng = bin_Ori_trn(looOriblock ~= rr);             % vector of trial labels for training data
            Oritstg = bin_Ori_tst(looOriblock == rr);             % trial labels for testing data
            
            OriX_pre = zeros(size(trn_Ori,2), nChans);           % to store design matrix: rows are observations (trials), and columns store the predictors for each orientation channel
            
            for ii=1:size(trn_Ori,2)
                OriX_pre(ii,:) = wshift('1D', bf, -(Oritrng(ii)-1));
            end
            
            Oriw_pre = OriX_pre\trn_Ori';
            
            Orix_pre = (Oriw_pre'\tst_Ori)';
            
            Ori_chan_pre(looOriblock==rr,:) = Orix_pre;
        end
        
        for ii=1:size(bin_Ori_trn,1)
            Ori_chan_pre(ii,:) = wshift('1D', Ori_chan_pre(ii,:), bin_Ori_trn(ii)-ceil(nChans/2));
        end
        
        Orichan_outp_avg(trn_tr,tst_tr,:) = nanmean(Ori_chan_pre,1);
        
    end
    
    Orichan_outp_fin(sub,:,:,:)= Orichan_outp_avg;
end

%for cross validation
for i = 1:size(Orichan_outp_fin,2)
    Orichan_outp(:,i,:) = Orichan_outp_fin(:,i,i,:);
end

save (['E:/2backrr_results/cv' tasks{t} '.mat'],'Orichan_outp');
%save (['E:/2backrr_results/tempGen_' tasks{t} '.mat','Orichan_outp']);
