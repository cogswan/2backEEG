%% significance tests on hypothesized windows for eye data

clear all
clc
%for 1B, 2B: delay is 1150-3150 ms
%for DR: delay is 1000-2000 ms
for job = 1:4
    for tf = 1:2 %timeframe:1 is stimulus, 2 is delay
        clearvars -except job tf p_all t_all d_all
        switch job
            case 1
                load E:/2backrr_eyelink/acc_by_sub_2B_UMI_hampel.mat
                t = 3;
                win0 = [0 500; 1150 3150];
            case 2
                load E:/2backrr_eyelink/acc_by_sub_2B_PMI_hampel.mat
                t = 3;
                win0 = [0 500; 1150 3150];
            case 3
                load E:/2backrr_eyelink/acc_by_sub_1B_hampel.mat
                t = 2;
                win0 = [0 500; 1150 3150];
            case 4
                load E:/2backrr_eyelink/acc_by_sub_DR_hampel.mat
                t = 4;
                win0 = [0 500; 1000 2000];
        end
        
        subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
        tasks = {'FL','1B','2B','DR'};
        
        task_dur  = [1000, 3550, 7300, 2750];
        times = linspace(-200,task_dur(t),size(acc_allsub,2));
        win = dsearchn(times',[win0(tf,1) win0(tf,2)]');
        sub_avg = mean(acc_allsub(:,win(1):win(2)),2);
        [h,p,ci,stats] = ttest(sub_avg' - 1/6, 0, 'dim', 2);
        
        ttest_mean = mean(sub_avg' - 1/6);
        d = (ttest_mean - 0)/stats.sd;
        
        if tf == 1
            save(['E:/2backrr_eyelink/sig_test/' tasks{t} '_stim.mat'])
        else
            save(['E:/2backrr_eyelink/sig_test/' tasks{t} '_delay.mat'])
        end
        p_all(job,tf) = p; %every row is a job, column are differen timeframes
        t_all(job,tf) = stats.tstat;
        d_all(job,tf) = d;
    end
end

%% FDR correction within task
p_2b = [p_all(1,1) p_all(1,2) p_all(2,2)];
p_1b = [p_all(3,1) p_all(3,2)];
p_dr = [p_all(4,1) p_all(4,2)];

[h1_2b, crit_p_2b, adj_ci_cvrg_2b, p_corr_fdr_2b]=fdr_bh(p_2b);
[h1_1b, crit_p_1b, adj_ci_cvrg_1b, p_corr_fdr_1b]=fdr_bh(p_1b);
[h1_dr, crit_p_dr, adj_ci_cvrg_dr, p_corr_fdr_dr]=fdr_bh(p_dr);
