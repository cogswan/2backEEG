%% Analysis of behavioral data

clear all;
subjs = {'2','3','4','5','6','8','10','12','7','9','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};

%initialize matrices to record data
acc = cell(4,1);
RT =  cell(4,1);
dp =  cell(4,1);

for t = 1:length(tasks) %loop through tasks
    blknum = [6 4 4 8]; %number of blocks for each task
    RTcol = [8 8 8 10]; %the column number for each task's recorded behav data
    
    for sub = 1:length(subjects)
        cd (['E:/2backrr_behav/S' subjects{sub} '/Behavioral/']);
        for blk = 1: blknum(t)
            
            load (['S' subjs{sub} '_' tasks{t} '_B'  num2str(blk) '.mat']);
            
            acc{t}(sub,blk) = accuracy;
            
            if t > 1
                dp{t}(sub,blk) = dprime;
            end
            
            RT{t}(sub,blk) =  nanmean(record(record(:,RTcol(t))~=99, RTcol(t)));
            
            %account for the first 2 DR blocks where behavioral data were
            %lost (S09)
            if (sub == 10 && t == 4 && blk < 3)
                acc{t}(sub,blk) = NaN;
                dp{t}(sub,blk) = NaN;
                RT{t}(sub,blk) = NaN;
            end
        end
    end
end

%% descriptive stats for behav data for eacdh task
acc_avg  = [];
acc_std = [];
dp_avg =[];
dp_std = [];
RT_avg =[];
RT_std = [];

for t = 1:4
    
    acc_avg_sub{t} = nanmean(acc{t},2);
    acc_avg(t) = nanmean(acc_avg_sub{t});
    acc_std(t) = nanstd(acc_avg_sub{t});
    
    
    dp_avg_sub{t}= nanmean(dp{t},2);
    dp_avg(t) = nanmean(dp_avg_sub{t});
    dp_std(t) = nanstd(dp_avg_sub{t});
    
    RT_avg_sub{t}= nanmean(RT{t},2);
    RT_avg(t) = nanmean(RT_avg_sub{t});
    RT_std(t) = nanstd(RT_avg_sub{t});
end


%% find if there are outliers in performance for each task
t = 4;
low_cut = nanmean(acc_avg_sub{t}) - 3*nanstd(acc_avg_sub{t});
find(acc_avg_sub{t}<low_cut)
