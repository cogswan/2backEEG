%% Significant tests of IEM reconstruction of certain windows
clear all; clc
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};

t = 4;
load E:/2backrr_results/cvDR.mat

bf=[0.0156    0.4219    1.0000    0.4219    0.0156 0]; %basis function
task_dur  = [1000, 3550, 7300, 2750]; %end time of trial of each task
times = linspace(-200,task_dur(t),size(Orichan_outp,2)); %vector of timecourse

%% sig. test for single window
win = dsearchn(times',[940 1040]');

%correlate reconstructed tuning curve with basis function
for isub=1:length(subjects)
    bfcorr(isub)=corr(bf', squeeze(mean(mean(Orichan_outp(isub,win(1):win(2),:),1),2)));
end

%conduct t test with Fischer z-transformed correlation coeff.
[h,p,ci,stats] = ttest(0.5*(log(1+bfcorr)-log(1-bfcorr)),0);

%calculate cohen's d
bfcorri= bfcorr;
ttest_mean = mean(0.5*(log(1+bfcorri)-log(1-bfcorri)));
d = (ttest_mean - 0)/stats.sd;

%% sig. test for two windows
%test principal hypothesis (PH) 1 & 2: delay - neg and pos reconstructions

d1 = dsearchn(times',[1150 3150]'); %first delay
d2 = dsearchn(times',[4900 6900]'); %second delay

%correlate reconstructed tuning curve with basis function
for isub=1:length(subjects)
    bfcorr(1,isub)=corr(bf', squeeze(mean(mean(Orichan_outp(isub,d1(1):d1(2),:),1),2)));
    bfcorr(2,isub)=corr(bf', squeeze(mean(mean(Orichan_outp(isub,d2(1):d2(2),:),1),2)));
end

%conduct t test with Fischer z-transformed correlation coeff.
for nbi = 1:2
    bfcorri= squeeze(bfcorr(nbi,:));
    [h(nbi),p(nbi),ci{nbi},stats{nbi}] = ttest(0.5*(log(1+bfcorri)-log(1-bfcorri)),0);
end

%FDR correction
[h1, crit_p, adj_ci_cvrg, p_corr_fdr]=fdr_bh(p);

%calculate cohen's d
for nbi = 1:2
    bfcorri= squeeze(bfcorr(nbi,:));
    ttest_mean(nbi) = mean(0.5*(log(1+bfcorri)-log(1-bfcorri)));
    d(nbi) = (ttest_mean(nbi) - 0)/stats{nbi}.sd;
end

