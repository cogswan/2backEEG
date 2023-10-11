%% plot/visualize IEM reconstruction: heat map or curve over a window
clear all
subjects = {'02','03','04','05','06','08','10','12','07','09','14','11','13','18','20','22','15','24','17','26','28','21','30','23','25','27','29','32','33','35'};
tasks = {'FL','1B','2B','DR'};
t0 = 4; %training task
t = 3; %testing task

%load in data
load('DR2B_940-1040_fin.mat')

%duplicate one channel to make curve symmetrical
Orimatrx = Orichan_outp;
Orimatrx_new=nan(length(subjects),size(Orichan_outp,2),7);%length(subjects));
Orimatrx_new(:,:,1)=Orimatrx(:,:,6);
Orimatrx_new(:,:,2:7)=Orimatrx(:,:,:);

task_dur  = [1000, 3550, 7300, 2750]; %end time of the trial for each task
pnts = [60 188 376 148]; %number of timepoints in each task
times0  = linspace(-200, task_dur(t0),pnts(t0)); %timecourse of training task
times = linspace(-200, task_dur(t), size(Orichan_outp,2)); %timecourse of testing task

%% Plot reconstruction heat map (ONE training window)
figure;
imagesc(1:7, linspace(-200,task_dur(t),size(Orichan_outp,2)), squeeze(mean(Orimatrx_new,1)))
colorbar
xlabel('Channel')
ylabel('Time(ms)')
set(gca,'xticklabel',{'-90','-60','-30','0','30','60','90'})
set(gca,'fontsize',12)

%% Plot reconstructed tuning curves of single windows
tstarts = 1150;
tends = 3150;
win=dsearchn(times',[tstarts tends]');

figure;
a = min(squeeze(mean(mean(Orimatrx_new(:,win(1):win(2),:,:),1),2)));
plot(squeeze(mean(mean(Orimatrx_new(:,win(1):win(2),:,:),1),2))-a, 'LineWidth', 3)
set(gca,'xticklabel',{'-90','-60','-30','0','30','60','90'})
xlabel('Channel')
ylabel('Response(a.u)')
set(gca,'fontsize',12)
