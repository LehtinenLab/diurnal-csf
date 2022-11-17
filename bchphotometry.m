cd '/Users/pkalugin/Desktop/full_runs'
filename = '220314_3xmice_run12022-06-05T17_49_47.5269760-04_00';
a = tdfread(filename);
b = a.(string(fieldnames(a)));
savefile = zeros(floor(size(b,1)/3),13);
% order of columns in savefile:
% tvals, mouse 1/470m, mouse 1/415m, mouse 1/470m clean, mouse 1/415m clean, etc
%%
mouse = 3;
tvals = b(1:3:end,1); % milliseconds counting off from midnight - lights ON 2.88E7, lights OFF 7.20E7
bgvals = b(2:3:end,mouse+1);
yvals470 = b(1:3:end,mouse+1);
yvals415 = b(3:3:end,mouse+1);
xvals = (1:min(length(yvals415),length(yvals470)))';
tvals = tvals(1:size(savefile,1));
bgvals = bgvals(1:size(savefile,1));
yvals470 = yvals470(1:size(savefile,1));
yvals415 = yvals415(1:size(savefile,1));
xvals = xvals(1:size(savefile,1));
trim470 = 1:150000; 
trim415 = 1:200000;
% 470nm cutoffs (assuming order is 470-560-415):
% 45000/32000/180000 for first, (27000)/40000:180000/15000 for second, (110000)/(95000)/(250000) for third
% 150000/150000/350000 for fourth, (140000)/140000/(140000) for fifth
% 4/26 run: 225000/225000/225000
% 5/3 run: 225000/225000/225000
% dark-dark runs:
% 5/12 run: 200000/200000/200000
% 5/19 run: 200000/200000/30000
% 5/27 run: ()/()/()
% WT runs:
% 6/1 run: ()/190000/100000
% 6/5 run: 20000/200000/200000
% 415nm cutoffs:
% 45000/45000/180000 for first, (45000)/34000/16000 for second,
% (100000)/100000/45000 for third, (100000)/150000/150000 for fourth,
% 10000/100000/(100000) for fifth
% 4/26 run: 30000/150000/150000
% 5/3 run: 100000/90000/100000
% dark-dark runs:
% 5/12 run: 50000/200000/200000
% 5/19 run: 190000/190000/190000
% 5/27 run: ()/()/()
% WT runs:
% 6/1 run: ()/200000/60000
% 6/5 run: 20000/60000/60000
xvals470t = xvals(trim470);
xvals415t = xvals(trim415);
%xvalst = (1:length(xvals(trim)))';
%yvalst415 = yvals415(trim);
%yvalst470 = yvals470(trim);
%bgvalst = bgvals(trim);
%% mixed fits
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
%% linear fit
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-100000,-100000],...
               'Upper',[100000,100000],...
               'StartPoint',[1 1]);
g = fittype('a+b*x','options',fo);
%%
figure,plot(yvals470-bgvals);
yvals470m = medfilt1(yvals470-bgvals,1001);
figure,plot(yvals470m);
figure,plot(yvals415-bgvals);
yvals415m = medfilt1(yvals415-bgvals,1001);
figure,plot(yvals415m);

yvals470mt = yvals470m(trim470);
yvals415mt = yvals415m(trim415);

%% stretched exponential
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
f1 = fit(xvals470t,yvals470mt,g);
f2 = fit(xvals415t,yvals415mt,g);
figure,plot(f1,xvals,yvals470m);
figure,plot(f2,xvals,yvals415m);
%figure,plot(f0);
%cleantr = yvals2-bgvals2-f1(xvals);
%cleantr = vals2-f1(xvals);
%cleantr = (yvals470-bgvals-f1(xvals)+f1.a)./(yvals415-bgvals-f2(xvals)+f2.a);
%cleantr = yvals470-bgvals-f1(xvals)+f1.a;
cleantr = [yvals470m yvals415m yvals470m-f1(xvals) yvals415m-f2(xvals)];
figure,plot(xvals,cleantr);
figure,plot(tvals,cleantr);
hold on
xline(2.88E7);
xline(7.20E7);
hold off
%% record
savefile(:,1) = tvals(1:size(savefile,1));
savefile(:,4*mouse-2:4*mouse+1) = cleantr(1:size(savefile,1),:);
%% save file
writematrix(savefile,strcat(filename,"_mfilt_raw_clean_all.csv"));
%% hyperbolic
fo2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-1000000,-1000000,-1000000],...
               'Upper',[1000000,1000000,1000000],...
               'StartPoint',[0 0 0]);
g2 = fittype('a+(b/(c-x))','options',fo2);
f2 = fit(xvalst,yvalst2,g2);
figure,plot(f2,xvals,yvals2);
figure,plot(xvals,yvals2-f2(xvals));


%% figures
%lights off at 7.20E7, so 3 hrs before to 5 hrs after is 6.12E7 to 0.36E7
%(next day) - take 48000 readings after 5pm
%lights on at 2.88E7, so 3 hrs before to 5 hrs after is 1.80E7 to 4.68E7
cd '/Users/pkalugin/Desktop/220411photometryfits'
a1 = readmatrix('220314_3xmice_run12022-03-17T11_31_50.5075456-04_00_clean415_run1.csv');
a2 = readmatrix('220314_3xmice_run12022-03-22T15_51_52_Correct.0959488-04_00_clean415_run2.csv');
a3 = readmatrix('220314_3xmice_run12022-03-24T13_48_37.3833600-04_00_clean415_run3.csv');
a4 = readmatrix('220314_3xmice_run12022-03-28T16_11_35.9201280-04_00_clean415_run4.csv');
a5 = readmatrix('220314_3xmice_run12022-04-05T11_04_24.4745856-04_00_clean415_run5.csv');
a6 = readmatrix('220314_3xmice_run12022-04-26T15_12_12.4783104-04_00_clean415_run6.csv');
a7 = readmatrix('220314_3xmice_run12022-05-03T14_19_50.6039680-04_00_clean415_run7.csv');
%% 470nm curves
cd '/Users/pkalugin/Desktop/220411photometryfits'
a1 = readmatrix('220314_3xmice_run12022-03-17T11_31_50.5075456-04_00_clean470_run1.csv');
a2 = readmatrix('220314_3xmice_run12022-03-22T15_51_52_Correct.0959488-04_00_clean470_run2.csv');
a3 = readmatrix('220314_3xmice_run12022-03-24T13_48_37.3833600-04_00_clean470_run3.csv');
a4 = readmatrix('220314_3xmice_run12022-03-28T16_11_35.9201280-04_00_clean470_run4.csv');
a5 = readmatrix('220314_3xmice_run12022-04-05T11_04_24.4745856-04_00_clean470_run5.csv');
a6 = readmatrix('220314_3xmice_run12022-04-26T15_12_12.4783104-04_00_clean470_run6.csv');
a7 = readmatrix('220314_3xmice_run12022-05-03T14_19_50.6039680-04_00_clean470_run7.csv');
%% fresh curves, raw vs bleach corrected (no linear pretrend)
cd '/Users/pkalugin/Desktop/full_runs'
a1 = readmatrix('220314_3xmice_run12022-03-17T11_31_50.5075456-04_00_mfilt_raw_clean_all.csv');
a2 = readmatrix('220314_3xmice_run12022-03-22T15_51_52_Correct.0959488-04_00_mfilt_raw_clean_all.csv');
a3 = readmatrix('220314_3xmice_run12022-03-24T13_48_37.3833600-04_00_mfilt_raw_clean_all.csv');
a4 = readmatrix('220314_3xmice_run12022-03-28T16_11_35.9201280-04_00_mfilt_raw_clean_all.csv');
a5 = readmatrix('220314_3xmice_run12022-04-05T11_04_24.4745856-04_00_mfilt_raw_clean_all.csv');
a6 = readmatrix('220314_3xmice_run12022-04-26T15_12_12.4783104-04_00_mfilt_raw_clean_all.csv');
a7 = readmatrix('220314_3xmice_run12022-05-03T14_19_50.6039680-04_00_mfilt_raw_clean_all.csv');
a8 = readmatrix('220314_3xmice_run12022-05-12T12_12_20.2895104-04_00_mfilt_raw_clean_all.csv');
a9 = readmatrix('220314_3xmice_run12022-05-19T12_05_44.4585216-04_00_mfilt_raw_clean_all.csv');
a10 = readmatrix('220314_3xmice_run12022-05-27T11_35_49.0945536-04_00_mfilt_raw_clean_all.csv'); %recording restarted on second day
a11 = readmatrix('220314_3xmice_run12022-06-01T18_09_12.6388224-04_00_mfilt_raw_clean_all.csv');
a12 = readmatrix('220314_3xmice_run12022-06-05T17_49_47.5269760-04_00_mfilt_raw_clean_all.csv');

%% POWERS
photopowers470 = [69 95 82;...
    round(0.9*279) round(1.1*279) 279;...
    round(0.9*219) round(1.1*219) 219;...
    round(0.9*198) round(1.1*198) 198;...
    173 224 191;...
    205 225 211;...
    256 287 255;...
    190 240 215;...
    215 250 215;...
    %185 225 205;...
    201 238 223;...
    260 245 200;...
    ];



%% FINAL FIGURES

% Fig 4c (!? linear detrend WT trace)
% 24 hrs
figure
ha = fill([17989 89945 89945 17989],[0.09 0.09 0.10 0.10],[0 0 0]);
hold on
%stdshade(excl470d0bleachttrdff(:,excl470d0bleachttrdffvolwt(13:end,1))',0.6,map(230,:))
stdshade(excl470d0bleachttrdff',0.6,[0 188/255 47/255]) %green TTR-mN
stdshade(excl470d0bleachwtdff',0.6,[165 166 169]/255) %gray WT
xline(17989)
xline(89945)
hold off

% 8 hrs (final Fig 4b - no subsetting responders)
figure
%hb = fill(floor([17989 47971 47971 17989])/49,[0.08 0.08 0.09 0.09],[0 0 0]);
hold on
%stdshade(excl470d0bleachttrdff(1:47971,excl470d0bleachttrdffvolwt(13:end,1))',0.6,map(230,:))
stdshade2(excl470d0bleachttrdff(1:49:47971,:)',1,[0 188/255 47/255])
stdshade2(excl470d0bleachwtdff(1:49:47971,:)',1,[165 166 169]/255)
xline(floor(17989)/49)
xticks(floor([0 5996 11993 17989 23985 29982 35978 41974 47971]/49))
xticklabels({'-3','-2','-1','0','1','2','3','4','5'})
xlim([0 47971]/49)
xlabel('Time to lights off (hrs)')
ylabel('dF/F')
%yticklabels({'-1','0','1','2','3','4','5','6','7','8','9'})
hold off
print -painters -depsc f4d.eps


% Fig S4d (work out esthetics - new boxplot functions from Matlab 2020?)
figure,boxplot(excl470d0bleachalldffvolm(:,1),excl470d0bleachalldffvolm(:,2))
hold on
scatter(excl470d0bleachalldffvolm(:,2),excl470d0bleachalldffvolm(:,1))
yline(min(excl470d0bleachalldffvolm(26:30,1)))
yline(max(excl470d0bleachalldffvolm(26:30,1)))
yticks([-0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45]);
hAx=gca;
hAx.Position(3)=hAx.Position(3)*15/25;
ylim([-0.15 0.25])
hold off

% Fig 4d (final Fig 4c)
% heatmap
map = customcolormap_preset('red-white-blue');
figure
hc = fill([17989 47971 47971 17989],[-1 -1 -2 -2],[0 0 0]);
hold on
hd = fill([-3700 -2200 -2200 -3700],[0.5 0.5 14.5 14.5],map(230,:),'LineStyle','none');
he = fill([-3700 -2200 -2200 -3700],[14.5 14.5 21.5 21.5],map(150,:),'LineStyle','none');
hf = fill([-3700 -2200 -2200 -3700],[21.5 21.5 25.5 25.5],map(26,:),'LineStyle','none');
imagesc(excl470d0bleachttrdff(1:47972,idxttrdff)',[-0.1 0.1])
green = cat(3, (102/255)*ones(size(excl470d0bleachttrdffpn2s(1:47972,:)')),(204/255)*ones(size(excl470d0bleachttrdffpn2s(1:47972,:)')),(0/255)*ones(size(excl470d0bleachttrdffpn2s(1:47972,:)'))); 
hd = imagesc(green);
mask = zeros(size(green,1),size(green,2));
for i = 1:14
    mask(i,excl470d0bleachttrdffpn2srises(i)-100:excl470d0bleachttrdffpn2srises(i)+100) = 1;
end
set(hd, 'AlphaData', mask)
axis([-3700 47971 -2 25.5])
set(gca, 'YDir','reverse')
colormap(gca,map)
colorbar
xline(17989,'LineWidth',1)
yline(14.5,'LineWidth',1)
yline(21.5,'LineWidth',1)
xticks([0 5996 11993 17989 23985 29982 35978 41974 47971])
xticklabels({'-3','-2','-1','0','1','2','3','4','5'})
hold off
print -painters -depsc f4dhfinal.eps

% second part of final Fig 4c
% 8 hrs
figure
%hb = fill([17989 47971 47971 17989],[0.15 0.15 0.16 0.16],[0 0 0]);
hold on
stdshade(excl470d0bleachttrdff(1:49:47971,excl470d0bleachttrdffvolwt(1:4,1))',1,map(26,:))
stdshade(excl470d0bleachttrdff(1:49:47971,excl470d0bleachttrdffvolwt(5:11,1))',1,map(150,:))
stdshade(excl470d0bleachttrdff(1:49:47971,excl470d0bleachttrdffvolwt(12:end,1))',1,map(230,:))
xline(floor(17989)/49)
xticks(floor([0 5996 11993 17989 23985 29982 35978 41974 47971])/49)
xticklabels({'-3','-2','-1','0','1','2','3','4','5'})
xlim([0 47971]/49)
ylim([-0.11 0.15])
%yticks([-0.12 -0.08 -0.04 0 0.04 0.08 0.12 0.16])
%yticklabels({'-12','-8','-4','0','4','8','12','16'})
ylabel('dF/F')
hold off
print -painters -depsc f4dpn.eps

% boxplot (final Fig 4d)
figure,boxplot(excl470d0bleachttrdffpn2srises-17989)
hold on
scatter(ones(1,length(excl470d0bleachttrdffpn2srises)),(excl470d0bleachttrdffpn2srises-17989))
yline(0)
xticks([])
yticks([-11993 -5996 0 5996 11993])
yticklabels({'-2','-1','0','1','2'})
ylim([-12000 12000])
xlim([0.8 1.2])
hold off


% Fig S4e+ (final Fig 4!)
% light-dark ttr vs wt 60 hrs starting day 1 (post-bleach, linear detrend)
figure,plot(1)
%he = fill([17989 89945 89945 17989],[-0.035 -0.035 -0.03 -0.03],[0 0 0]);
hold on
%hf = fill([17989+143912 215868 215868 17989+143912],[-0.035 -0.035 -0.03 -0.03],[0 0 0]);
stdshade(excl470ldrawttrdffd(1:49:end,:)',1,[0 188/255 47/255])
stdshade(excl470ldrawwtdffd(1:49:end,:)',1,[165 166 169]/255)
xline(floor(17989)/49)
xline(floor(89945)/49)
xline(floor(17989+143912)/49)
xticks(floor([11993 29982 47971 65960 83949 101938 119927 137916 11993+143912 29982+143912 47971+143912 65960+143912])/49)
xticklabels({'6pm', '9pm', 'midnight', '3am', '6am', '9am', 'noon', '3pm', '6pm', '9pm', 'midnight', '3am'})
xtickangle(90)
xlim(floor([0 215868])/49)
yticks([-0.03 -0.02 -0.01 0 0.01 0.02 0.03])
%yticklabels({'-3', '-2', '-1', '0', '1', '2', '3'})
ylim([-0.03 0.03])
hold off
print -painters -depsc f4sd.eps
% hard to claim that ttr-mneon isn't fully bleached out

figure
hold on
plot(mean(excl470ldrawttrdffd,2),'color',[0 188/255 47/255])
plot(mean(excl470ldrawwtdffd,2),'color',[165 166 169]/255)
plot(mean(excl470ldrawttrdffd,2)-mean(excl470ldrawwtdffd,2),'color',map(230,:))
xline(17989)
xline(89945)
xline(17989+143912)
xticks([11993 29982 47971 65960 83949 101938 119927 137916 11993+143912 29982+143912 47971+143912 65960+143912])
xticklabels({'6pm', '9pm', 'midnight', '3am', '6am', '9am', 'noon', '3pm', '6pm', '9pm', 'midnight', '3am'})
xtickangle(90)
xlim([0 215868])
yticks([-0.03 -0.02 -0.01 0 0.01 0.02 0.03])
%yticklabels({'-3', '-2', '-1', '0', '1', '2', '3'})
ylim([-0.03 0.03])
hold off
print -painters -depsc f4se.eps

% dark-dark ttr 60 hrs starting day 1 (post-bleach, linear detrend)
figure,plot(1)
%hg = fill([17989 359780 359780 17989],[-0.0592 -0.0592 -0.05 -0.05],[0 0 0]);
hold on
stdshade(all470ddrawttrdffd(1:49:end,:)',1,[0 188/255 47/255])
xline(floor(17989)/49)
xline(floor(89945)/49)
xline(floor(17989+143912)/49)
xline(floor(89945+143912)/49)
xline(floor(17989+143912+143912)/49)
xticks(floor([11993 29982 47971 65960 83949 101938 119927 137916 11993+143912 29982+143912 47971+143912 65960+143912 83949+143912 101938+143912 119927+143912 137916+143912 11993+143912+143912 29982+143912+143912 47971+143912+143912 65960+143912+143912]/49))
xticklabels({'6pm', '9pm', 'midnight', '3am', '6am', '9am', 'noon', '3pm', '6pm', '9pm', 'midnight', '3am', '6am', '9am', 'noon', '3pm', '6pm', '9pm', 'midnight', '3am'})
xtickangle(90)
xlim(floor([0 359780]/49))
yticks([-0.05 -0.04 -0.03 -0.02 -0.01 0 0.01 0.02 0.03 0.04 0.05 0.06])
%yticklabels({'-5', '-4', '-3', '-2', '-1', '0', '1', '2', '3', '4', '5', '6'})
ylim([-0.05 0.06])
hold off
print -painters -depsc f4h.eps


% Fig 4f
% heatmap d0
map = customcolormap_preset('red-white-blue');
figure
hc = fill([17989 47971 47971 17989],[-1 -1 -2 -2],[0 0 0]);
hold on
%hd = fill([0 47971 47971 0],[0.5 0.5 25.5 25.5],[0 0 0]);
imagesc(all470dddetrttrdffd0pn2s(1:47972,:)',[-1 1])
green = cat(3, zeros(size(all470dddetrttrdffd0pn2s(1:47972,:)')),ones(size(all470dddetrttrdffd0pn2s(1:47972,:)')),zeros(size(all470dddetrttrdffd0pn2s(1:47972,:)'))); 
hd = imagesc(green);
mask = zeros(size(green,1),size(green,2));
for i = 1:5
    mask(i,all470dddetrttrdffd0pn2srises(i)-50:all470dddetrttrdffd0pn2srises(i)+50) = 1;
end
set(hd, 'AlphaData', mask)
axis([-3000 47971 -2 9.5])
set(gca, 'YDir','reverse')
colormap(gca,map)
colorbar
xline(17989,'LineWidth',1)
%yline(13.5,'LineWidth',1)
%yline(21.5,'LineWidth',1)
hold off

% heatmap d1
map = customcolormap_preset('red-white-blue');
figure
hc = fill([0 47971 47971 0],[-1 -1 -2 -2],[0 0 0]);
hold on
%hd = fill([0 47971 47971 0],[0.5 0.5 25.5 25.5],[0 0 0]);
imagesc(all470dddetrttrdffd1pn2s(1+143912:47972+143912,:)',[-1 1])
green = cat(3, zeros(size(all470dddetrttrdffd1pn2s(1+143912:47972+143912,:)')),ones(size(all470dddetrttrdffd1pn2s(1+143912:47972+143912,:)')),zeros(size(all470dddetrttrdffd1pn2s(1+143912:47972+143912,:)'))); 
hd = imagesc(green);
mask = zeros(size(green,1),size(green,2));
for i = 1:6
    mask(i,(all470dddetrttrdffd1pn2srises(i)-50-143912:all470dddetrttrdffd1pn2srises(i)+50-143912)) = 1;
end
set(hd, 'AlphaData', mask)
axis([-3000 47971 -2 9.5])
set(gca, 'YDir','reverse')
colormap(gca,map)
colorbar
xline(17989,'LineWidth',1)
%yline(13.5,'LineWidth',1)
%yline(21.5,'LineWidth',1)
hold off

% heatmap d2
map = customcolormap_preset('red-white-blue');
figure
hc = fill([0 47971 47971 0],[-1 -1 -2 -2],[0 0 0]);
hold on
%hd = fill([0 47971 47971 0],[0.5 0.5 25.5 25.5],[0 0 0]);
imagesc(all470dddetrttrdffd2pn2s(1+143912+143912:47972+143912+143912,:)',[-1 1])
green = cat(3, zeros(size(all470dddetrttrdffd2pn2s(1+143912+143912:47972+143912+143912,:)')),ones(size(all470dddetrttrdffd2pn2s(1+143912+143912:47972+143912+143912,:)')),zeros(size(all470dddetrttrdffd2pn2s(1+143912+143912:47972+143912+143912,:)'))); 
hd = imagesc(green);
mask = zeros(size(green,1),size(green,2));
for i = 1:6
    mask(i,all470dddetrttrdffd2pn2srises(i)-50-143912-143912:all470dddetrttrdffd2pn2srises(i)+50-143912-143912) = 1;
end
set(hd, 'AlphaData', mask)
axis([-3000 47971 -2 9.5])
set(gca, 'YDir','reverse')
colormap(gca,map)
colorbar
xline(17989,'LineWidth',1)
%yline(13.5,'LineWidth',1)
%yline(21.5,'LineWidth',1)
hold off

% boxplot
figure,boxplot([all470dddetrttrdffd0pn2srises-17989 all470dddetrttrdffd1pn2srises-17989-143912 all470dddetrttrdffd2pn2srises-17989-143912-143912],[ones(1,length(all470dddetrttrdffd0pn2srises)) ones(1,length(all470dddetrttrdffd1pn2srises))+1 ones(1,length(all470dddetrttrdffd2pn2srises))+2])
hold on
scatter([ones(1,length(all470dddetrttrdffd0pn2srises)) ones(1,length(all470dddetrttrdffd1pn2srises))+1 ones(1,length(all470dddetrttrdffd2pn2srises))+2],[all470dddetrttrdffd0pn2srises-17989 all470dddetrttrdffd1pn2srises-17989-143912 all470dddetrttrdffd2pn2srises-17989-143912-143912])
yline(0)
hold off


% Fig S4f
figure
boxplot(([mean(all470dddetrttrdffd0(29982:35978,:),1); mean(all470dddetrttrdffd1(29982+143912:35978+143912,:),1); mean(all470dddetrttrdffd2(29982+143912+143912:35978+143912+143912,:),1)])')
hold on
scatter([ones(1,9) ones(1,9)+1 ones(1,9)+2],([mean(all470dddetrttrdffd0(29982:35978,:),1) mean(all470dddetrttrdffd1(29982+143912:35978+143912,:),1) mean(all470dddetrttrdffd2(29982+143912+143912:35978+143912+143912,:),1)]))
yline(0.01)


% Fig 4b/supplementary recording summary
figure
hold on
for i = 11:12
    fill([-hrspreloff(i)-3 hrstotal(i)-hrspreloff(i)-3 hrstotal(i)-hrspreloff(i)-3 -hrspreloff(i)-3],-[i-1-3 i-1-3 i-0.5-3 i-0.5-3],[165 166 169]/255,'LineStyle','none');
end
%hold off

%figure
%hold on
for i = 1:7
    fill([-hrspreloff(i)-3 hrstotal(i)-hrspreloff(i)-3 hrstotal(i)-hrspreloff(i)-3 -hrspreloff(i)-3],-[i-1 i-1 i-0.5 i-0.5],[0 188/255 47/255],'LineStyle','none');
end
%hold off

%figure
%hold on
for i = 8:10
    j = 0;
    if i == 10
        j = 1;
    end
    fill([-hrspreloff(i)-3+j*24 hrstotal(i)-hrspreloff(i)-3+j*24 hrstotal(i)-hrspreloff(i)-3+j*24 -hrspreloff(i)-3+j*24],-[i-1+2 i-1+2 i-0.5+2 i-0.5+2],[0 188/255 47/255],'LineStyle','none');
end
xticks([-12 0 12 24 36 48 60 72 84 96])
hold off


% FigS4f-g period/amplitude summary
perldts = readmatrix('ldttr3000drain.csv');
perldws = readmatrix('ldwt3000drain.csv');
perddts = readmatrix('ddttr3000drain.csv');
amps = [max(ldavs(1:143912,:))-min(ldavs(1:143912,:)), max(ddavs(1:143912,:))-min(ddavs(1:143912,:))];
% S4f pers from RAIN
figure
b = bar([mean(perldts(:,4)) mean(perldws(:,4)) mean(perddts(:,4))],0.6,'k');
hold on
errhigh = [std(perldts(:,4))/sqrt(length(perldts(:,4))) std(perldws(:,4))/sqrt(length(perldws(:,4))) std(perddts(:,4))/sqrt(length(perddts(:,4)))];
er = errorbar(1:3,[mean(perldts(:,4)) mean(perldws(:,4)) mean(perddts(:,4))],errhigh);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylim([22 25])
hold off

b.FaceColor = 'flat';
b.CData(1,:) = [0 188/255 47/255];
b.CData(2,:) = [165 166 169]/255;
b.CData(3,:) = map(230,:);
b.CData(5,:) = [0 188/255 47/255];
ylim([22 25])
% S4g ranges
figure
boxplot([ranldts ranldws randdts],[ones(1,length(ranldts)) ones(1,length(ranldws))+1 ones(1,length(randdts))+2])
hold on
scatter([ones(1,length(ranldts)) ones(1,length(ranldws))+1 ones(1,length(randdts))+2],[ranldts ranldws randdts])
hold off
% p-values:
[h1,p1] = ttest2(ranldts, ranldws);
[h2,p2] = ttest2(ranldts, randdts);


%% Fit sinusoids to multi-day traces
% averages
ldavs = [mean(excl470ldrawttrdffd,2),mean(excl470ldrawwtdffd,2),mean(excl470ldrawttrdffd,2)-mean(excl470ldrawwtdffd,2)];
ddavs = mean(all470ddrawttrdffd,2);

for i = 1:3
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000,0.3,-10000000],...
        'Upper',[1000,10000000,0.6,10000000],...
        'StartPoint',[0 1 0.43 0]);
    gd = fittype('(a+b*sin((c*x/10000)+(d/100000)))/1000','options',fod);
    f1 = fit((1:215868)',ldavs(:,i),gd);
    figure,plot(f1,(1:215868)',ldavs(:,i));
    %figure,plot(ldavs(:,i)-f1(1:215868));
    (2*pi/f1.c*10000)*24/143912
    f1.b/1000
    %a4trimdd36hd(:,i) = a4trimdd36h(:,i)-f1(1:215868);
end


% mean amplitudes over first 24 hours
max(ldavs(1:143912,:))-min(ldavs(1:143912,:))
max(ddavs(1:143912,:))-min(ddavs(1:143912,:))


fod = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[-1000,-10000000,0.3,-10000000],...
    'Upper',[1000,10000000,0.6,10000000],...
    'StartPoint',[0 1 0.43 0]);
gd = fittype('(a+b*sin((c*x/10000)+(d/100000)))/1000','options',fod);
f1 = fit((1:359780)',ddavs,gd);
figure,plot(f1,(1:359780)',ddavs);
%figure,plot(ldavs(:,i)-f1(1:215868));
(2*pi/f1.c*10000)*24/143912
f1.b/1000

%% individual traces
% ld ttr-mn
perldts = zeros(1,10);
ampldts = zeros(1,10);
for i = 1:10
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000,0.3,-10000000],...
        'Upper',[1000,10000000,0.6,10000000],...
        'StartPoint',[0 1 0.43 0]);
    gd = fittype('(a+b*sin((c*x/10000)+(d/100000)))/1000','options',fod);
    f1 = fit((1:215868)',excl470ldrawttrdffd(:,i),gd);
    figure,plot(f1,(1:215868)',excl470ldrawttrdffd(:,i));
    %figure,plot(ldavs(:,i)-f1(1:215868));
    perldts(i) = (2*pi/f1.c*10000)*24/143912;
    ampldts(i) = f1.b/1000;
    %a4trimdd36hd(:,i) = a4trimdd36h(:,i)-f1(1:215868);
end

% ld wt
perldws = zeros(1,5);
ampldws = zeros(1,5);
for i = 1:5
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000,0.3,-10000000],...
        'Upper',[1000,10000000,0.6,10000000],...
        'StartPoint',[0 1 0.43 0]);
    gd = fittype('(a+b*sin((c*x/10000)+(d/100000)))/1000','options',fod);
    f1 = fit((1:215868)',excl470ldrawwtdffd(:,i),gd);
    figure,plot(f1,(1:215868)',excl470ldrawwtdffd(:,i));
    %figure,plot(ldavs(:,i)-f1(1:215868));
    perldws(i) = (2*pi/f1.c*10000)*24/143912;
    ampldws(i) = f1.b/1000;
    %a4trimdd36hd(:,i) = a4trimdd36h(:,i)-f1(1:215868);
end

% dd ttr-mn
perddts = zeros(1,9);
ampddts = zeros(1,9);
for i = 1:9
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000,0.3,-10000000],...
        'Upper',[1000,10000000,0.6,10000000],...
        'StartPoint',[0 1 0.43 0]);
    gd = fittype('(a+b*sin((c*x/10000)+(d/100000)))/1000','options',fod);
    f1 = fit((1:359780)',all470ddrawttrdffd(:,i),gd);
    figure,plot(f1,(1:359780)',all470ddrawttrdffd(:,i));
    %figure,plot(ldavs(:,i)-f1(1:215868));
    perddts(i) = (2*pi/f1.c*10000)*24/143912;
    ampddts(i) = f1.b/1000;
    %a4trimdd36hd(:,i) = a4trimdd36h(:,i)-f1(1:215868);
end


% ranges over first 24 hours
ranldts = max(excl470ldrawttrdffd(1:143912,:),[],1)-min(excl470ldrawttrdffd(1:143912,:),[],1);
ranldws = max(excl470ldrawwtdffd(1:143912,:),[],1)-min(excl470ldrawwtdffd(1:143912,:),[],1);
randdts = max(all470ddrawttrdffd(1:143912,:),[],1)-min(all470ddrawttrdffd(1:143912,:),[],1);

%% Each recording's start/stop times
hrspreloff = [idx01(1)...
idx02(1)...
idx03(1)...
idx04(1)...
idx05(1)...
idx06(1)...
idx07(1)...
idx08(1)...
idx09(1)...
idx010(1)...
idx011(1)...
idx012(1)]*24/143912;

hrstotal = [size(a1,1)...
size(a2,1)...
size(a3,1)...
size(a4,1)...
size(a5,1)...
size(a6,1)...
size(a7,1)...
size(a8,1)...
size(a9,1)...
size(a10,1)...
size(a11,1)...
size(a12,1)]*24/143912;

% need to stretch from 8 hours pre to 70 (ld + wt)/94 (dd) hours post
% 8, 9, 10 are dd
% recording 10 was started 24 hours late






%% RAIN outputs
% light-dark: ttr, wt, ttr-wt (df/f then detrend, means)
ldrain = [((0:215867)*36/215868)',mean(excl470ldrawttrdffd,2),mean(excl470ldrawwtdffd,2),mean(excl470ldrawttrdffd,2)-mean(excl470ldrawwtdffd,2)];
ldrainds30 = zeros(floor(size(ldrain,1)/30),size(ldrain,2));
for i = 1:floor(size(ldrain,1)/30)
    ldrainds30(i,1) = ldrain((i-1)*30+1,1);
    ldrainds30(i,2:end) = mean(ldrain(((i-1)*30+1):i*30,2:end),1);
end
ldrainds300 = zeros(floor(size(ldrain,1)/300),size(ldrain,2));
for i = 1:floor(size(ldrain,1)/300)
    ldrainds300(i,1) = ldrain((i-1)*300+1,1);
    ldrainds300(i,2:end) = mean(ldrain(((i-1)*300+1):i*300,2:end),1);
end
ldrainds3000 = zeros(floor(size(ldrain,1)/3000),size(ldrain,2));
for i = 1:floor(size(ldrain,1)/3000)
    ldrainds3000(i,1) = ldrain((i-1)*3000+1,1);
    ldrainds3000(i,2:end) = mean(ldrain(((i-1)*3000+1):i*3000,2:end),1);
end
ldrainds1500 = zeros(floor(size(ldrain,1)/1500),size(ldrain,2));
for i = 1:floor(size(ldrain,1)/1500)
    ldrainds1500(i,1) = ldrain((i-1)*1500+1,1);
    ldrainds1500(i,2:end) = mean(ldrain(((i-1)*1500+1):i*1500,2:end),1);
end
writematrix(ldrain,'220818ldrain.csv');
writematrix(ldrainds30,'220818ldrainds30.csv');
writematrix(ldrainds300,'220818ldrainds300.csv');
writematrix(ldrainds3000,'220818ldrainds3000.csv');
writematrix(ldrainds1500,'220818ldrainds1500.csv');

% dark-dark (df/f then detrend, means)
ddrain = [((0:359779)*60/359780)',mean(all470ddrawttrdffd,2)];
ddrainds30 = zeros(floor(size(ddrain,1)/30),size(ddrain,2));
for i = 1:floor(size(ddrain,1)/30)
    ddrainds30(i,1) = ddrain((i-1)*30+1,1);
    ddrainds30(i,2:end) = mean(ddrain(((i-1)*30+1):i*30,2:end),1);
end
ddrainds300 = zeros(floor(size(ddrain,1)/300),size(ddrain,2));
for i = 1:floor(size(ddrain,1)/300)
    ddrainds300(i,1) = ddrain((i-1)*300+1,1);
    ddrainds300(i,2:end) = mean(ddrain(((i-1)*300+1):i*300,2:end),1);
end
ddrainds3000 = zeros(floor(size(ddrain,1)/3000),size(ddrain,2));
for i = 1:floor(size(ddrain,1)/3000)
    ddrainds3000(i,1) = ddrain((i-1)*3000+1,1);
    ddrainds3000(i,2:end) = mean(ddrain(((i-1)*3000+1):i*3000,2:end),1);
end
ddrainds1500 = zeros(floor(size(ddrain,1)/1500),size(ddrain,2));
for i = 1:floor(size(ddrain,1)/1500)
    ddrainds1500(i,1) = ddrain((i-1)*1500+1,1);
    ddrainds1500(i,2:end) = mean(ddrain(((i-1)*1500+1):i*1500,2:end),1);
end
writematrix(ddrain,'220818ddrain.csv');
writematrix(ddrainds30,'220818ddrainds30.csv');
writematrix(ddrainds300,'220818ddrainds300.csv');
writematrix(ddrainds3000,'220818ddrainds3000.csv');
writematrix(ddrainds1500,'220818ddrainds1500.csv');


% NEW: RAIN outputs for all traces downsampled by 3000
ldrainttr = [((0:215867)*36/215868)',excl470ldrawttrdffd];
ldrainwt = [((0:215867)*36/215868)',excl470ldrawwtdffd];
ddrainttr = [((0:359779)*60/359780)',all470ddrawttrdffd];
ldrainttr3000 = zeros(floor(size(ldrainttr,1)/3000),size(ldrainttr,2));
for i = 1:floor(size(ldrainttr,1)/3000)
    ldrainttr3000(i,1) = ldrainttr((i-1)*3000+1,1);
    ldrainttr3000(i,2:end) = mean(ldrainttr(((i-1)*3000+1):i*3000,2:end),1);
end
ldrainwt3000 = zeros(floor(size(ldrainwt,1)/3000),size(ldrainwt,2));
for i = 1:floor(size(ldrainwt,1)/3000)
    ldrainwt3000(i,1) = ldrainwt((i-1)*3000+1,1);
    ldrainwt3000(i,2:end) = mean(ldrainwt(((i-1)*3000+1):i*3000,2:end),1);
end
ddrainttr3000 = zeros(floor(size(ddrainttr,1)/3000),size(ddrainttr,2));
for i = 1:floor(size(ddrainttr,1)/3000)
    ddrainttr3000(i,1) = ddrainttr((i-1)*3000+1,1);
    ddrainttr3000(i,2:end) = mean(ddrainttr(((i-1)*3000+1):i*3000,2:end),1);
end
writematrix(ldrainttr,'220820ldrainttr.csv');
writematrix(ldrainwt,'220820ldrainwt.csv');
writematrix(ddrainttr,'220820ddrainttr.csv');
writematrix(ldrainttr3000,'220820ldrainttr3000.csv');
writematrix(ldrainwt3000,'220820ldrainwt3000.csv');
writematrix(ddrainttr3000,'220820ddrainttr3000.csv');

%% light-dark long-term
% 4-7 (TTR) vs 11-12 (WT) - should have enough for two lights off from d1
a4trimdd36h = a4((idx04(1)+143912):(idx04(1)+359779),:);
a5trimdd36h = a5((idx05(1)+143912):(idx05(1)+359779),:);
a6trimdd36h = a6((idx06(1)+143912):(idx06(1)+359779),:);
a7trimdd36h = a7((idx07(1)+143912):(idx07(1)+359779),:);

a11trimdd36h = a11((idx011(1)+143912):(idx011(1)+359779),:);
a12trimdd36h = a12((idx012(1)+143912):(idx012(1)+359779),:);

%% detrending
a4trimdd36hd = zeros(size(a4trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a4trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a4trimdd36h(:,i));
    %figure,plot(a4trimdd36h(:,i)-f1(1:215868));
    a4trimdd36hd(:,i) = a4trimdd36h(:,i)-f1(1:215868);
end

a5trimdd36hd = zeros(size(a5trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a5trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a5trimdd36h(:,i));
    %figure,plot(a5trimdd36h(:,i)-f1(1:215868));
    a5trimdd36hd(:,i) = a5trimdd36h(:,i)-f1(1:215868);
end

a6trimdd36hd = zeros(size(a6trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a6trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a6trimdd36h(:,i));
    %figure,plot(a6trimdd36h(:,i)-f1(1:215868));
    a6trimdd36hd(:,i) = a6trimdd36h(:,i)-f1(1:215868);
end

a7trimdd36hd = zeros(size(a7trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a7trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a7trimdd36h(:,i));
    %figure,plot(a7trimdd36h(:,i)-f1(1:215868));
    a7trimdd36hd(:,i) = a7trimdd36h(:,i)-f1(1:215868);
end

a11trimdd36hd = zeros(size(a11trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a11trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a11trimdd36h(:,i));
    %figure,plot(a11trimdd36h(:,i)-f1(1:215868));
    a11trimdd36hd(:,i) = a11trimdd36h(:,i)-f1(1:215868);
end

a12trimdd36hd = zeros(size(a12trimdd36h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:215868)',a12trimdd36h(:,i),gd);
    %figure,plot(f1,(1:215868)',a12trimdd36h(:,i));
    %figure,plot(a12trimdd36h(:,i)-f1(1:215868));
    a12trimdd36hd(:,i) = a12trimdd36h(:,i)-f1(1:215868);
end

%% df/f
% with bleach
for i = 1:3
    a4trimdd36hdff(:,i) = (a4trimdd36h(:,i*4)-mean(a4trimdd36h(1:5996,i*4),1))/mean(a4trimdd36h(1:5996,(i*4-2)),1);
    a5trimdd36hdff(:,i) = (a5trimdd36h(:,i*4)-mean(a5trimdd36h(1:5996,i*4),1))/mean(a5trimdd36h(1:5996,(i*4-2)),1);
    a6trimdd36hdff(:,i) = (a6trimdd36h(:,i*4)-mean(a6trimdd36h(1:5996,i*4),1))/mean(a6trimdd36h(1:5996,(i*4-2)),1);
    a7trimdd36hdff(:,i) = (a7trimdd36h(:,i*4)-mean(a7trimdd36h(1:5996,i*4),1))/mean(a7trimdd36h(1:5996,(i*4-2)),1);
    
    a11trimdd36hdff(:,i) = (a11trimdd36h(:,i*4)-mean(a11trimdd36h(1:5996,i*4),1))/mean(a11trimdd36h(1:5996,(i*4-2)),1);
    a12trimdd36hdff(:,i) = (a12trimdd36h(:,i*4)-mean(a12trimdd36h(1:5996,i*4),1))/mean(a12trimdd36h(1:5996,(i*4-2)),1);
end

excl470ldbleachttrdff = [a4trimdd36hdff a5trimdd36hdff(:,2) a6trimdd36hdff a7trimdd36hdff];
excl470ldbleachwtdff = [a11trimdd36hdff a12trimdd36hdff(:,[2 3])];

% no bleach (raw)!
for i = 1:3
    a4trimdd36hdffr(:,i) = (a4trimdd36h(:,i*4-2)-mean(a4trimdd36h(1:5996,i*4-2),1))/mean(a4trimdd36h(1:5996,(i*4-2)),1);
    a5trimdd36hdffr(:,i) = (a5trimdd36h(:,i*4-2)-mean(a5trimdd36h(1:5996,i*4-2),1))/mean(a5trimdd36h(1:5996,(i*4-2)),1);
    a6trimdd36hdffr(:,i) = (a6trimdd36h(:,i*4-2)-mean(a6trimdd36h(1:5996,i*4-2),1))/mean(a6trimdd36h(1:5996,(i*4-2)),1);
    a7trimdd36hdffr(:,i) = (a7trimdd36h(:,i*4-2)-mean(a7trimdd36h(1:5996,i*4-2),1))/mean(a7trimdd36h(1:5996,(i*4-2)),1);
    
    a11trimdd36hdffr(:,i) = (a11trimdd36h(:,i*4-2)-mean(a11trimdd36h(1:5996,i*4-2),1))/mean(a11trimdd36h(1:5996,(i*4-2)),1);
    a12trimdd36hdffr(:,i) = (a12trimdd36h(:,i*4-2)-mean(a12trimdd36h(1:5996,i*4-2),1))/mean(a12trimdd36h(1:5996,(i*4-2)),1);
end

excl470ldrawttrdff = [a4trimdd36hdffr a5trimdd36hdffr(:,2) a6trimdd36hdffr a7trimdd36hdffr];
excl470ldrawwtdff = [a11trimdd36hdffr a12trimdd36hdffr(:,[2 3])];

% raw plus linear detrend!
for i = 1:3
    a4trimdd36hdffd(:,i) = (a4trimdd36hd(:,i*4-2)-mean(a4trimdd36hd(1:5996,i*4-2),1))/mean(a4trimdd36h(1:5996,(i*4-2)),1);
    a5trimdd36hdffd(:,i) = (a5trimdd36hd(:,i*4-2)-mean(a5trimdd36hd(1:5996,i*4-2),1))/mean(a5trimdd36h(1:5996,(i*4-2)),1);
    a6trimdd36hdffd(:,i) = (a6trimdd36hd(:,i*4-2)-mean(a6trimdd36hd(1:5996,i*4-2),1))/mean(a6trimdd36h(1:5996,(i*4-2)),1);
    a7trimdd36hdffd(:,i) = (a7trimdd36hd(:,i*4-2)-mean(a7trimdd36hd(1:5996,i*4-2),1))/mean(a7trimdd36h(1:5996,(i*4-2)),1);
    
    a11trimdd36hdffd(:,i) = (a11trimdd36hd(:,i*4-2)-mean(a11trimdd36hd(1:5996,i*4-2),1))/mean(a11trimdd36h(1:5996,(i*4-2)),1);
    a12trimdd36hdffd(:,i) = (a12trimdd36hd(:,i*4-2)-mean(a12trimdd36hd(1:5996,i*4-2),1))/mean(a12trimdd36h(1:5996,(i*4-2)),1);
end

excl470lddetrttrdff = [a4trimdd36hdffd a5trimdd36hdffd(:,2) a6trimdd36hdffd a7trimdd36hdffd];
excl470lddetrwtdff = [a11trimdd36hdffd a12trimdd36hdffd(:,[2 3])];

%% detrending the raw df or df/f
excl470ldrawttrdffd = zeros(size(excl470ldrawttrdff));
for i = 1:10
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-10000000],...
        'Upper',[10000000],...
        'StartPoint',[0]);
    gd = fittype('a*x/1000000','options',fod);
    f1 = fit((1:215868)',excl470ldrawttrdff(:,i),gd);
    %figure,plot(f1,(1:215868)',excl470ldrawttrdff(:,i));
    %figure,plot(excl470ldrawttrdff(:,i)-f1(1:215868));
    excl470ldrawttrdffd(:,i) = excl470ldrawttrdff(:,i)-f1(1:215868);
end

excl470ldrawttrdfd = zeros(size(excl470ldrawttrdf));
for i = 1:10
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-10000000],...
        'Upper',[10000000],...
        'StartPoint',[0]);
    gd = fittype('a*x/1000000','options',fod);
    f1 = fit((1:215868)',excl470ldrawttrdf(:,i),gd);
    %figure,plot(f1,(1:215868)',excl470ldrawttrdf(:,i));
    %figure,plot(excl470ldrawttrdf(:,i)-f1(1:215868));
    excl470ldrawttrdfd(:,i) = excl470ldrawttrdf(:,i)-f1(1:215868);
end

excl470ldrawwtdffd = zeros(size(excl470ldrawwtdff));
for i = 1:5
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-10000000],...
        'Upper',[10000000],...
        'StartPoint',[0]);
    gd = fittype('a*x/1000000','options',fod);
    f1 = fit((1:215868)',excl470ldrawwtdff(:,i),gd);
    %figure,plot(f1,(1:215868)',excl470ldrawwtdff(:,i));
    %figure,plot(excl470ldrawwtdff(:,i)-f1(1:215868));
    excl470ldrawwtdffd(:,i) = excl470ldrawwtdff(:,i)-f1(1:215868);
end

excl470ldrawwtdfd = zeros(size(excl470ldrawwtdf));
for i = 1:5
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-10000000],...
        'Upper',[10000000],...
        'StartPoint',[0]);
    gd = fittype('a*x/1000000','options',fod);
    f1 = fit((1:215868)',excl470ldrawwtdf(:,i),gd);
    %figure,plot(f1,(1:215868)',excl470ldrawwtdf(:,i));
    %figure,plot(excl470ldrawwtdf(:,i)-f1(1:215868));
    excl470ldrawwtdfd(:,i) = excl470ldrawwtdf(:,i)-f1(1:215868);
end

%% df
% with bleach
for i = 1:3
    a4trimdd36hdf(:,i) = (a4trimdd36h(:,i*4)-mean(a4trimdd36h(1:5996,i*4),1));
    a5trimdd36hdf(:,i) = (a5trimdd36h(:,i*4)-mean(a5trimdd36h(1:5996,i*4),1));
    a6trimdd36hdf(:,i) = (a6trimdd36h(:,i*4)-mean(a6trimdd36h(1:5996,i*4),1));
    a7trimdd36hdf(:,i) = (a7trimdd36h(:,i*4)-mean(a7trimdd36h(1:5996,i*4),1));
    
    a11trimdd36hdf(:,i) = (a11trimdd36h(:,i*4)-mean(a11trimdd36h(1:5996,i*4),1));
    a12trimdd36hdf(:,i) = (a12trimdd36h(:,i*4)-mean(a12trimdd36h(1:5996,i*4),1));
end

excl470ldbleachttrdf = [a4trimdd36hdf a5trimdd36hdf(:,2) a6trimdd36hdf a7trimdd36hdf];
excl470ldbleachwtdf = [a11trimdd36hdf a12trimdd36hdf(:,[2 3])];

% no bleach (raw)!
for i = 1:3
    a4trimdd36hdfr(:,i) = (a4trimdd36h(:,i*4-2)-mean(a4trimdd36h(1:5996,i*4-2),1));
    a5trimdd36hdfr(:,i) = (a5trimdd36h(:,i*4-2)-mean(a5trimdd36h(1:5996,i*4-2),1));
    a6trimdd36hdfr(:,i) = (a6trimdd36h(:,i*4-2)-mean(a6trimdd36h(1:5996,i*4-2),1));
    a7trimdd36hdfr(:,i) = (a7trimdd36h(:,i*4-2)-mean(a7trimdd36h(1:5996,i*4-2),1));
    
    a11trimdd36hdfr(:,i) = (a11trimdd36h(:,i*4-2)-mean(a11trimdd36h(1:5996,i*4-2),1));
    a12trimdd36hdfr(:,i) = (a12trimdd36h(:,i*4-2)-mean(a12trimdd36h(1:5996,i*4-2),1));
end

excl470ldrawttrdf = [a4trimdd36hdfr a5trimdd36hdfr(:,2) a6trimdd36hdfr a7trimdd36hdfr];
excl470ldrawwtdf = [a11trimdd36hdfr a12trimdd36hdfr(:,[2 3])];

% raw plus linear detrend!
for i = 1:3
    a4trimdd36hdfd(:,i) = (a4trimdd36hd(:,i*4-2)-mean(a4trimdd36hd(1:5996,i*4-2),1));
    a5trimdd36hdfd(:,i) = (a5trimdd36hd(:,i*4-2)-mean(a5trimdd36hd(1:5996,i*4-2),1));
    a6trimdd36hdfd(:,i) = (a6trimdd36hd(:,i*4-2)-mean(a6trimdd36hd(1:5996,i*4-2),1));
    a7trimdd36hdfd(:,i) = (a7trimdd36hd(:,i*4-2)-mean(a7trimdd36hd(1:5996,i*4-2),1));
    
    a11trimdd36hdfd(:,i) = (a11trimdd36hd(:,i*4-2)-mean(a11trimdd36hd(1:5996,i*4-2),1));
    a12trimdd36hdfd(:,i) = (a12trimdd36hd(:,i*4-2)-mean(a12trimdd36hd(1:5996,i*4-2),1));
end

excl470lddetrttrdf = [a4trimdd36hdfd a5trimdd36hdfd(:,2) a6trimdd36hdfd a7trimdd36hdfd];
excl470lddetrwtdf = [a11trimdd36hdfd a12trimdd36hdfd(:,[2 3])];

%% dark-dark long-term

a8trimdd60h = a8((idx08(1)+143912):(idx08(1)+503691),:);
a9trimdd60h = a9((idx09(1)+143912):(idx09(1)+503691),:);
a10trimdd60h = a10(idx010(1):(idx010(1)+359779),:);

%% detrending
a8trimdd60hd = zeros(size(a8trimdd60h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:359780)',a8trimdd60h(:,i),gd);
    %figure,plot(f1,(1:359780)',a8trimdd60h(:,i));
    %figure,plot(a8trimdd60h(:,i)-f1(1:359780));
    a8trimdd60hd(:,i) = a8trimdd60h(:,i)-f1(1:359780);
end

a9trimdd60hd = zeros(size(a9trimdd60h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:359780)',a9trimdd60h(:,i),gd);
    %figure,plot(f1,(1:359780)',a9trimdd60h(:,i));
    %figure,plot(a9trimdd60h(:,i)-f1(1:359780));
    a9trimdd60hd(:,i) = a9trimdd60h(:,i)-f1(1:359780);
end

a10trimdd60hd = zeros(size(a10trimdd60h));
for i = 2:13
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-1000,-10000000],...
        'Upper',[1000,10000000],...
        'StartPoint',[0 1]);
    gd = fittype('a+((b*x)/1000000)','options',fod);
    f1 = fit((1:359780)',a10trimdd60h(:,i),gd);
    %figure,plot(f1,(1:359780)',a10trimdd60h(:,i));
    %figure,plot(a10trimdd60h(:,i)-f1(1:359780));
    a10trimdd60hd(:,i) = a10trimdd60h(:,i)-f1(1:359780);
end

%% df/f
% with bleach
for i = 1:3
    a8trimdd60hdff(:,i) = (a8trimdd60h(:,i*4)-mean(a8trimdd60h(1:5996,i*4),1))/mean(a8trimdd60h(1:5996,(i*4-2)),1);
    a9trimdd60hdff(:,i) = (a9trimdd60h(:,i*4)-mean(a9trimdd60h(1:5996,i*4),1))/mean(a9trimdd60h(1:5996,(i*4-2)),1);
    a10trimdd60hdff(:,i) = (a10trimdd60h(:,i*4)-mean(a10trimdd60h(1:5996,i*4),1))/mean(a10trimdd60h(1:5996,(i*4-2)),1);
end

all470ddbleachttrdff = [a8trimdd60hdff a9trimdd60hdff a10trimdd60hdff];

% no bleach (raw)!
for i = 1:3
    a8trimdd60hdffr(:,i) = (a8trimdd60h(:,i*4-2)-mean(a8trimdd60h(1:5996,i*4-2),1))/mean(a8trimdd60h(1:5996,(i*4-2)),1);
    a9trimdd60hdffr(:,i) = (a9trimdd60h(:,i*4-2)-mean(a9trimdd60h(1:5996,i*4-2),1))/mean(a9trimdd60h(1:5996,(i*4-2)),1);
    a10trimdd60hdffr(:,i) = (a10trimdd60h(:,i*4-2)-mean(a10trimdd60h(1:5996,i*4-2),1))/mean(a10trimdd60h(1:5996,(i*4-2)),1);
end

all470ddrawttrdff = [a8trimdd60hdffr a9trimdd60hdffr a10trimdd60hdffr];

% raw plus linear detrend! d0
for i = 1:3
    a8trimdd60hdffd0(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1:5996,i*4-2),1))/mean(a8trimdd60h(1:5996,(i*4-2)),1);
    a9trimdd60hdffd0(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1:5996,i*4-2),1))/mean(a9trimdd60h(1:5996,(i*4-2)),1);
    a10trimdd60hdffd0(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1:5996,i*4-2),1))/mean(a10trimdd60h(1:5996,(i*4-2)),1);
end

all470dddetrttrdffd0 = [a8trimdd60hdffd0 a9trimdd60hdffd0 a10trimdd60hdffd0];
%excl470ddbleachttrdff = [a8trimd024hdff a9trimd024hdff];


% raw plus linear detrend! d1
for i = 1:3
    a8trimdd60hdffd1(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1+143912:5996+143912,i*4-2),1))/mean(a8trimdd60h(1+143912:5996+143912,(i*4-2)),1);
    a9trimdd60hdffd1(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1+143912:5996+143912,i*4-2),1))/mean(a9trimdd60h(1+143912:5996+143912,(i*4-2)),1);
    a10trimdd60hdffd1(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1+143912:5996+143912,i*4-2),1))/mean(a10trimdd60h(1+143912:5996+143912,(i*4-2)),1);
end

all470dddetrttrdffd1 = [a8trimdd60hdffd1 a9trimdd60hdffd1 a10trimdd60hdffd1];
%excl470ddbleachttrdff = [a8trimd024hdff a9trimd024hdff];

% raw plus linear detrend! d2
for i = 1:3
    a8trimdd60hdffd2(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1))/mean(a8trimdd60h(1+143912+143912:5996+143912+143912,(i*4-2)),1);
    a9trimdd60hdffd2(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1))/mean(a9trimdd60h(1+143912+143912:5996+143912+143912,(i*4-2)),1);
    a10trimdd60hdffd2(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1))/mean(a10trimdd60h(1+143912+143912:5996+143912+143912,(i*4-2)),1);
end

all470dddetrttrdffd2 = [a8trimdd60hdffd2 a9trimdd60hdffd2 a10trimdd60hdffd2];
%excl470ddbleachttrdff = [a8trimd024hdff a9trimd024hdff];


%% detrending the raw df/f
all470ddrawttrdffd = zeros(size(all470ddrawttrdff));
for i = 1:9
    fod = fitoptions('Method','NonlinearLeastSquares',...
        'Lower',[-10000000],...
        'Upper',[10000000],...
        'StartPoint',[0]);
    gd = fittype('a*x/1000000','options',fod);
    f1 = fit((1:359780)',all470ddrawttrdff(:,i),gd);
    %figure,plot(f1,(1:359780)',all470ddrawttrdff(:,i));
    %figure,plot(all470ddrawttrdff(:,i)-f1(1:359780));
    all470ddrawttrdffd(:,i) = all470ddrawttrdff(:,i)-f1(1:359780);
end

%% df
% with bleach
for i = 1:3
    a8trimdd60hdf(:,i) = (a8trimdd60h(:,i*4)-mean(a8trimdd60h(1:5996,i*4),1));
    a9trimdd60hdf(:,i) = (a9trimdd60h(:,i*4)-mean(a9trimdd60h(1:5996,i*4),1));
    a10trimdd60hdf(:,i) = (a10trimdd60h(:,i*4)-mean(a10trimdd60h(1:5996,i*4),1));
end

all470ddbleachttrdf = [a8trimdd60hdf a9trimdd60hdf a10trimdd60hdf];

% no bleach (raw)!
for i = 1:3
    a8trimdd60hdfr(:,i) = (a8trimdd60h(:,i*4-2)-mean(a8trimdd60h(1:5996,i*4-2),1));
    a9trimdd60hdfr(:,i) = (a9trimdd60h(:,i*4-2)-mean(a9trimdd60h(1:5996,i*4-2),1));
    a10trimdd60hdfr(:,i) = (a10trimdd60h(:,i*4-2)-mean(a10trimdd60h(1:5996,i*4-2),1));
end

all470ddrawttrdf = [a8trimdd60hdfr a9trimdd60hdfr a10trimdd60hdfr];

% raw plus linear detrend! d0
for i = 1:3
    a8trimdd60hdfd0(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1:5996,i*4-2),1));
    a9trimdd60hdfd0(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1:5996,i*4-2),1));
    a10trimdd60hdfd0(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1:5996,i*4-2),1));
end

all470dddetrttrdfd0 = [a8trimdd60hdfd0 a9trimdd60hdfd0 a10trimdd60hdfd0];


% raw plus linear detrend! d1
for i = 1:3
    a8trimdd60hdfd1(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1+143912:5996+143912,i*4-2),1));
    a9trimdd60hdfd1(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1+143912:5996+143912,i*4-2),1));
    a10trimdd60hdfd1(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1+143912:5996+143912,i*4-2),1));
end

all470dddetrttrdfd1 = [a8trimdd60hdfd1 a9trimdd60hdfd1 a10trimdd60hdfd1];

% raw plus linear detrend! d2
for i = 1:3
    a8trimdd60hdfd2(:,i) = (a8trimdd60hd(:,i*4-2)-mean(a8trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1));
    a9trimdd60hdfd2(:,i) = (a9trimdd60hd(:,i*4-2)-mean(a9trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1));
    a10trimdd60hdfd2(:,i) = (a10trimdd60hd(:,i*4-2)-mean(a10trimdd60hd(1+143912+143912:5996+143912+143912,i*4-2),1));
end

all470dddetrttrdfd2 = [a8trimdd60hdfd2 a9trimdd60hdfd2 a10trimdd60hdfd2];

%% peak-normalize to hours +2 to +3 (use absolute value of peak to keep directionality!)
all470dddetrttrdffd0pn2 = all470dddetrttrdffd0./abs(mean(all470dddetrttrdffd0(29982:35978,:),1));

[~,idxttrdffddd0] = sort(-mean(all470dddetrttrdffd0(29982:35978,:),1),2);
all470dddetrttrdffd0pn2s = all470dddetrttrdffd0pn2(:,idxttrdffddd0);


all470dddetrttrdffd1pn2 = all470dddetrttrdffd1./abs(mean(all470dddetrttrdffd1(29982+143912:35978+143912,:),1));

[~,idxttrdffddd1] = sort(-mean(all470dddetrttrdffd1(29982+143912:35978+143912,:),1),2);
all470dddetrttrdffd1pn2s = all470dddetrttrdffd1pn2(:,idxttrdffddd1);


all470dddetrttrdffd2pn2 = all470dddetrttrdffd2./abs(mean(all470dddetrttrdffd2(29982+143912+143912:35978+143912+143912,:),1));

[~,idxttrdffddd2] = sort(-mean(all470dddetrttrdffd2(29982+143912+143912:35978+143912+143912,:),1),2);
all470dddetrttrdffd2pn2s = all470dddetrttrdffd2pn2(:,idxttrdffddd2);

%% Calculate rise times, i.e. where signal is 50% of peak walking backwards from peak
all470dddetrttrdffd0pn2srises = zeros(1,5);
for i = 1:5
    af = find(all470dddetrttrdffd0pn2s(1:47972,i) > 1);
    ag = find(all470dddetrttrdffd0pn2s(1:af(1),i) < 0.5);
    all470dddetrttrdffd0pn2srises(i) = ag(length(ag))+1;
end

all470dddetrttrdffd1pn2srises = zeros(1,6);
for i = 1:6
    af = find(all470dddetrttrdffd1pn2s(1+143912:47972+143912,i) > 1);
    ag = find(all470dddetrttrdffd1pn2s(1+143912:af(1)+143912,i) < 0.5);
    all470dddetrttrdffd1pn2srises(i) = ag(length(ag))+1+143912;
end

all470dddetrttrdffd2pn2srises = zeros(1,6);
for i = 1:6
    af = find(all470dddetrttrdffd2pn2s(1+143912+143912:47972+143912+143912,i) > 1);
    ag = find(all470dddetrttrdffd2pn2s(1+143912+143912:af(1)+143912+143912,i) < 0.5);
    all470dddetrttrdffd2pn2srises(i) = ag(length(ag))+1+143912+143912;
end


%% day 0 lights off (10x TTR-mNeon vs 2x WT)
idx01 = find(a1(:,1)>6.12E7);
idx02 = find(a2(:,1)>6.12E7);
idx03 = find(a3(:,1)>6.12E7);
idx04 = find(a4(:,1)>6.12E7);
idx05 = find(a5(:,1)>6.12E7);
idx06 = find(a6(:,1)>6.12E7);
idx07 = find(a7(:,1)>6.12E7);
idx08 = find(a8(:,1)>6.12E7);
idx09 = find(a9(:,1)>6.12E7);
idx010 = find(a10(:,1)>6.12E7); % this is the one where first lights off is actually second!!
idx011 = find(a11(35052:end,1)>3.6E6)+35051;
idx012 = find(a12(36982:end,1)>3.6E6)+36981;

% looks like it's actualy 143912 readings per 24 hours, not 144000
a1trimd024h = a1(idx01(1):(idx01(1)+143911),:);
a2trimd024h = a2(idx02(1):(idx02(1)+143911),:);
a3trimd024h = a3(idx03(1):(idx03(1)+143911),:);
a4trimd024h = a4(idx04(1):(idx04(1)+143911),:);
a5trimd024h = a5(idx05(1):(idx05(1)+143911),:);
a6trimd024h = a6(idx06(1):(idx06(1)+143911),:);
a7trimd024h = a7(idx07(1):(idx07(1)+143911),:);
a8trimd024h = a8(idx08(1):(idx08(1)+143911),:);
a9trimd024h = a9(idx09(1):(idx09(1)+143911),:);
%a10trimd024h = a10(idx010(1):(idx010(1)+143911),:);
a11trimd024h = a11(idx011(1):(idx011(1)+143911),:);
a12trimd024h = a12(idx012(1):(idx012(1)+143911),:);

%% NEW 6/27/2022 calculate dF/F relative to first hour mean
% F is average of the RAW signal (not the fit - I think)
% idea is to pull out column X-mean(column X)/mean(column X-2) for X = 4,8,12
for i = 1:3
    a1trimd024hdff(:,i) = (a1trimd024h(:,i*4)-mean(a1trimd024h(1:5996,i*4),1))/mean(a1trimd024h(1:5996,(i*4-2)),1);
    a2trimd024hdff(:,i) = (a2trimd024h(:,i*4)-mean(a2trimd024h(1:5996,i*4),1))/mean(a2trimd024h(1:5996,(i*4-2)),1);
    a3trimd024hdff(:,i) = (a3trimd024h(:,i*4)-mean(a3trimd024h(1:5996,i*4),1))/mean(a3trimd024h(1:5996,(i*4-2)),1);
    a4trimd024hdff(:,i) = (a4trimd024h(:,i*4)-mean(a4trimd024h(1:5996,i*4),1))/mean(a4trimd024h(1:5996,(i*4-2)),1);
    a5trimd024hdff(:,i) = (a5trimd024h(:,i*4)-mean(a5trimd024h(1:5996,i*4),1))/mean(a5trimd024h(1:5996,(i*4-2)),1);
    a6trimd024hdff(:,i) = (a6trimd024h(:,i*4)-mean(a6trimd024h(1:5996,i*4),1))/mean(a6trimd024h(1:5996,(i*4-2)),1);
    a7trimd024hdff(:,i) = (a7trimd024h(:,i*4)-mean(a7trimd024h(1:5996,i*4),1))/mean(a7trimd024h(1:5996,(i*4-2)),1);
    a8trimd024hdff(:,i) = (a8trimd024h(:,i*4)-mean(a8trimd024h(1:5996,i*4),1))/mean(a8trimd024h(1:5996,(i*4-2)),1);
    a9trimd024hdff(:,i) = (a9trimd024h(:,i*4)-mean(a9trimd024h(1:5996,i*4),1))/mean(a9trimd024h(1:5996,(i*4-2)),1);
    %a10trimd024hdff(:,i) = (a10trimd024h(:,i*4)-mean(a10trimd024h(1:5996,i*4),1))/mean(a10trimd024h(1:5996,(i*4-2)),1);
    a11trimd024hdff(:,i) = (a11trimd024h(:,i*4)-mean(a11trimd024h(1:5996,i*4),1))/mean(a11trimd024h(1:5996,(i*4-2)),1);
    a12trimd024hdff(:,i) = (a12trimd024h(:,i*4)-mean(a12trimd024h(1:5996,i*4),1))/mean(a12trimd024h(1:5996,(i*4-2)),1);
end

%% just raw Fs mean over first 1 hr

for i = 1:3
    firsthrfs(:,i) = [mean(a1trimd024h(1:5996,(i*4-2)),1);...
    mean(a2trimd024h(1:5996,(i*4-2)),1);...
    mean(a3trimd024h(1:5996,(i*4-2)),1);...
    mean(a4trimd024h(1:5996,(i*4-2)),1);...
    mean(a5trimd024h(1:5996,(i*4-2)),1);...
    mean(a6trimd024h(1:5996,(i*4-2)),1);...
    mean(a7trimd024h(1:5996,(i*4-2)),1);...
    mean(a8trimd024h(1:5996,(i*4-2)),1);...
    mean(a9trimd024h(1:5996,(i*4-2)),1);...
    mean(a11trimd024h(1:5996,(i*4-2)),1);...
    mean(a12trimd024h(1:5996,(i*4-2)),1);...
    ];
end

% F/mW
fmwratios = firsthrfs./photopowers470;
fmwratios = reshape(fmwratios,1,33);
fmwratiosttr = fmwratios([1 2 3 4 5 6 7 8 9 10 11 12 14 16 17 18 19 20 21 22 23 24 25 26 27]);
fmwratioswt = fmwratios([28 29 30 32 33]);

%% NEW 6/27/2022 mean dF/F traces (bleach corrected)

all470d0bleachttrdff = [a1trimd024hdff a2trimd024hdff a3trimd024hdff...
    a4trimd024hdff a5trimd024hdff a6trimd024hdff a7trimd024hdff...
    a8trimd024hdff a9trimd024hdff];
all470d0bleachwtdff = [a11trimd024hdff a12trimd024hdff];
excl470d0bleachttrdff = [a1trimd024hdff a2trimd024hdff a3trimd024hdff...
    a4trimd024hdff a5trimd024hdff(:,2) a6trimd024hdff a7trimd024hdff...
    a8trimd024hdff a9trimd024hdff];
excl470d0bleachwtdff = [a11trimd024hdff a12trimd024hdff(:,[2 3])];

%% NEW 6/27/2022 peak-normalize to hours +2 to +3 (use absolute value of peak to keep directionality!)
% EDIT 9/13/2022 peak-normalize to hours +1.5 to +2.5
all470d0bleachttrdffpn2 = all470d0bleachttrdff./abs(mean(all470d0bleachttrdff(26984:32980,:),1));
all470d0bleachwtdffpn2 = all470d0bleachwtdff./abs(mean(all470d0bleachwtdff(26984:32980,:),1));
excl470d0bleachttrdffpn2 = excl470d0bleachttrdff./abs(mean(excl470d0bleachttrdff(26984:32980,:),1));
excl470d0bleachwtdffpn2 = excl470d0bleachwtdff./abs(mean(excl470d0bleachwtdff(26984:32980,:),1));

[~,idxttrdff] = sort(-mean(excl470d0bleachttrdff(26984:32980,:),1),2);
excl470d0bleachttrdffpn2s = excl470d0bleachttrdffpn2(:,idxttrdff);
excl470d0bleachttrdffpn2sd = excl470d0bleachttrdffpn2s(1:49:47972,:);

[~,idxwtdff] = sort(-mean(excl470d0bleachwtdff(26984:32980,:),1),2);
excl470d0bleachwtdffpn2s = excl470d0bleachwtdffpn2(:,idxwtdff);
excl470d0bleachwtdffpn2sd = excl470d0bleachwtdffpn2s(1:49:47972,:);

%% NEW 9/13/2022 z-score the peak-normalized traces!
excl470d0bleachttrdffpn2sz = zeros(size(excl470d0bleachttrdffpn2s));
for i=1:25
    excl470d0bleachttrdffpn2sz(:,i) = excl470d0bleachttrdffpn2s(:,i)/std(excl470d0bleachttrdffpn2s(1:5996,i));
end

excl470d0bleachwtdffpn2sz = zeros(size(excl470d0bleachwtdffpn2s));
for i=1:5
    excl470d0bleachwtdffpn2sz(:,i) = excl470d0bleachwtdffpn2s(:,i)/std(excl470d0bleachwtdffpn2s(1:5996,i));
end


%% NEW 9/4/2022 z-score relative to hours -3 to -1, then peak-normalize, then sort relative to wt
excl470d0bleachttrdffz = zeros(size(excl470d0bleachttrdff));
for i=1:25
    excl470d0bleachttrdffz(:,i) = excl470d0bleachttrdff(:,i)/std(excl470d0bleachttrdff(1:5996,i));
end

excl470d0bleachwtdffz = zeros(size(excl470d0bleachwtdff));
for i=1:5
    excl470d0bleachwtdffz(:,i) = excl470d0bleachwtdff(:,i)/std(excl470d0bleachwtdff(1:5996,i));
end


excl470d0bleachttrdffzpn2 = excl470d0bleachttrdffz./abs(mean(excl470d0bleachttrdffz(29982:35978,:),1));
excl470d0bleachwtdffzpn2 = excl470d0bleachwtdffz./abs(mean(excl470d0bleachwtdffz(29982:35978,:),1));

[~,idxttrdffz] = sort(-mean(excl470d0bleachttrdffz(29982:35978,:),1),2);
excl470d0bleachttrdffzpn2s = excl470d0bleachttrdffzpn2(:,idxttrdffz);
excl470d0bleachttrdffzpn2sd = excl470d0bleachttrdffzpn2s(1:49:47972,:);

[~,idxwtdffz] = sort(-mean(excl470d0bleachwtdffz(29982:35978,:),1),2);
excl470d0bleachwtdffzpn2s = excl470d0bleachwtdffzpn2(:,idxwtdff);
excl470d0bleachwtdffzpn2sd = excl470d0bleachwtdffzpn2s(1:49:47972,:);



excl470d0bleachttrdffzvol = zeros(size(excl470d0bleachttrdffz,2),2);
for j = 1:size(excl470d0bleachttrdffz,2)
    excl470d0bleachttrdffzvol(j,1) = mean(excl470d0bleachttrdffz(29982:35978,j))-mean(excl470d0bleachttrdffz(1:5996,j));
    [~,excl470d0bleachttrdffzvol(j,2)] = ttest2(excl470d0bleachttrdffz(29982:35978,j),excl470d0bleachttrdffz(1:5996,j));
end

excl470d0bleachwtdffzvol = zeros(size(excl470d0bleachwtdffz,2),2);
for j = 1:size(excl470d0bleachwtdffz,2)
    excl470d0bleachwtdffzvol(j,1) = mean(excl470d0bleachwtdffz(29982:35978,j))-mean(excl470d0bleachwtdffz(1:5996,j));
    [~,excl470d0bleachwtdffzvol(j,2)] = ttest2(excl470d0bleachwtdffz(29982:35978,j),excl470d0bleachwtdffz(1:5996,j));
end

excl470d0bleachalldffzvolm = [[excl470d0bleachttrdffzvol(:,1); excl470d0bleachwtdffzvol(:,1)] [ones(25,1)+1; zeros(5,1)+1]];



excl470d0bleachttrdffzvolwt = zeros(size(excl470d0bleachttrdffz,2),3);
for j = 1:size(excl470d0bleachttrdffz,2)
    excl470d0bleachttrdffzvolwt(j,1) = j;
    excl470d0bleachttrdffzvolwt(j,2) = mean(excl470d0bleachttrdffz(29982:35978,j))-mean(mean(excl470d0bleachwtdffz(29982:35978,:)));
    [~,excl470d0bleachttrdffzvolwt(j,3)] = ttest2(excl470d0bleachttrdffz(29982:35978,j),reshape(excl470d0bleachwtdffz(29982:35978,:),1,[]));
end

excl470d0bleachttrdffzvolwt = sortrows(excl470d0bleachttrdffzvolwt,2);


%% NEW 6/27/2022 auROC-style analysis on exclusion bleach corrected traces
% DOES NOT SEEM TO BE NECESSARY/THE WAY TO GO
% 10 min windows relative to first hour baseline (bins of 999 readings)
% this is a sort of noise-normalized peak normalization(?)
% steps are:
% 1. compute ecdf for both bins of interest
% 2. evenly sample from min of mins to max of maxes of ecdfs e.g. 10000x
% 3. plot ecdfs relative to each other over this range
% 4. calculae auROC

excl470d0bleachttrdffauroc = zeros(144,size(excl470d0bleachttrdff,2));
for j = 1:size(excl470d0bleachttrdff,2)
    for w = 1:144  %144 windows of 10 min in 24 hours
        [f1,x1] = ecdf(excl470d0bleachttrdff(1:5996,j));
        [f2,x2] = ecdf(excl470d0bleachttrdff(((w-1)*999+1):((w-1)*999+1000),j));
        xs = min(min(x1),min(x2));
        xf = max(max(x1),max(x2));
        xrange = xs:((xf-xs)/10000):xf;
        v1 = interp1(x1(2:end),f1(2:end),xrange);
        v2 = interp1(x2(2:end),f2(2:end),xrange);
        n1 = find(isnan(v1));
        n2 = find(isnan(v2));
        if isempty(n1) == 0
            if n1(1) == 1
                tr = 0;
                for i = 1:length(n1)
                    tr = tr + (n1(i) ~= i);
                end
                if tr > 0
                    v1(1:n1(end-tr)) = 0;
                    v1(n1(end-tr+1):end) = 1;
                else
                    v1(1:n1(end)) = 0;
                end
            else
                v1(n1(1):end) = 1;
            end
        end
        if isempty(n2) == 0
            if n2(1) == 1
                tr = 0;
                for i = 1:length(n2)
                    tr = tr + (n2(i) ~= i);
                end
                if tr > 0
                    v2(1:n2(end-tr)) = 0;
                    v2(n2(end-tr+1):end) = 1;
                else
                    v2(1:n2(end)) = 0;
                end
            else
                v2(n2(1):end) = 1;
            end
        end
        %plot(1-v1,1-v2);
        excl470d0bleachttrdffauroc(w,j) = trapz(v1,1-v2);
    end
end

%% ttest for each bin relative to baseline?
% ALSO NOT GREAT
excl470d0bleachttrdfftt = zeros(144,size(excl470d0bleachttrdff,2));
for j = 1:size(excl470d0bleachttrdff,2)
    for w = 1:144  %144 windows of 10 min in 24 hours
        [~,excl470d0bleachttrdfftt(w,j)] = ttest2(excl470d0bleachttrdff(1:5996,j),excl470d0bleachttrdff(((w-1)*999+1):((w-1)*999+1000),j));
    end
end

%% Volcano plot of all traces, comparing hours -3 to -2 to +2 to +3
% EDIT 9/13/2022 compare to hours +1.5 to +2.5

excl470d0bleachttrdffvol = zeros(size(excl470d0bleachttrdff,2),2);
for j = 1:size(excl470d0bleachttrdff,2)
    excl470d0bleachttrdffvol(j,1) = mean(excl470d0bleachttrdff(26984:32980,j))-mean(excl470d0bleachttrdff(1:5996,j));
    [~,excl470d0bleachttrdffvol(j,2)] = ttest2(excl470d0bleachttrdff(26984:32980,j),excl470d0bleachttrdff(1:5996,j));
end

excl470d0bleachwtdffvol = zeros(size(excl470d0bleachwtdff,2),2);
for j = 1:size(excl470d0bleachwtdff,2)
    excl470d0bleachwtdffvol(j,1) = mean(excl470d0bleachwtdff(26984:32980,j))-mean(excl470d0bleachwtdff(1:5996,j));
    [~,excl470d0bleachwtdffvol(j,2)] = ttest2(excl470d0bleachwtdff(26984:32980,j),excl470d0bleachwtdff(1:5996,j));
end

excl470d0bleachalldffvolm = [[excl470d0bleachttrdffvol(:,1); excl470d0bleachwtdffvol(:,1)] [ones(25,1)+1; zeros(5,1)+1]];

%% Volcano plot of all traces, comparing hours -3 to -2 to +2 to +3 to WT
% this seems to be the most useful!
% EDIT 9/13/2022 compare to hours +1.5 to +2.5

excl470d0bleachttrdffvolwt = zeros(size(excl470d0bleachttrdff,2),3);
for j = 1:size(excl470d0bleachttrdff,2)
    excl470d0bleachttrdffvolwt(j,1) = j;
    excl470d0bleachttrdffvolwt(j,2) = mean(excl470d0bleachttrdff(26984:32980,j))-mean(mean(excl470d0bleachwtdff(26984:32980,:)));
    [~,excl470d0bleachttrdffvolwt(j,3)] = ttest2(excl470d0bleachttrdff(26984:32980,j),reshape(excl470d0bleachwtdff(26984:32980,:),1,[]));
end

excl470d0bleachttrdffvolwt = sortrows(excl470d0bleachttrdffvolwt,2);

%% Calculate rise time, i.e. where signal is 20% of peak walking backwards from peak
excl470d0bleachttrdffpn2srises = zeros(1,14);
for i = 1:14
    af = find(excl470d0bleachttrdffpn2s(1:47972,i) > 1);
    ag = find(excl470d0bleachttrdffpn2s(1:af(1),i) < 0.2);
    excl470d0bleachttrdffpn2srises(i) = ag(length(ag))+1;
end



%% calculate df relative to first two hour mean
% lights off at 17990 (idx(1) + 17988)!
% lights on at 89946 (idx(1) + 89945)!
a1trimd024hdf = a1trimd024h-mean(a1trimd024h(1:11993,:),1);
a2trimd024hdf = a2trimd024h-mean(a2trimd024h(1:11993,:),1);
a3trimd024hdf = a3trimd024h-mean(a3trimd024h(1:11993,:),1);
a4trimd024hdf = a4trimd024h-mean(a4trimd024h(1:11993,:),1);
a5trimd024hdf = a5trimd024h-mean(a5trimd024h(1:11993,:),1);
a6trimd024hdf = a6trimd024h-mean(a6trimd024h(1:11993,:),1);
a7trimd024hdf = a7trimd024h-mean(a7trimd024h(1:11993,:),1);
a8trimd024hdf = a8trimd024h-mean(a8trimd024h(1:11993,:),1);
a9trimd024hdf = a9trimd024h-mean(a9trimd024h(1:11993,:),1);
%a10trimd024hdf = a10trimd024h-mean(a10trimd024h(1:11993,:),1);
a11trimd024hdf = a11trimd024h-mean(a11trimd024h(1:11993,:),1);
a12trimd024hdf = a12trimd024h-mean(a12trimd024h(1:11993,:),1);

%% peak normalize df to hrs 1-3 after lights off
a1trimd024hdfpn = a1trimd024hdf./mean(a1trimd024hdf(23986:35978,:),1);
a2trimd024hdfpn = a2trimd024hdf./mean(a2trimd024hdf(23986:35978,:),1);
a3trimd024hdfpn = a3trimd024hdf./mean(a3trimd024hdf(23986:35978,:),1);
a4trimd024hdfpn = a4trimd024hdf./mean(a4trimd024hdf(23986:35978,:),1);
a5trimd024hdfpn = a5trimd024hdf./mean(a5trimd024hdf(23986:35978,:),1);
a6trimd024hdfpn = a6trimd024hdf./mean(a6trimd024hdf(23986:35978,:),1);
a7trimd024hdfpn = a7trimd024hdf./mean(a7trimd024hdf(23986:35978,:),1);
a8trimd024hdfpn = a8trimd024hdf./mean(a8trimd024hdf(23986:35978,:),1);
a9trimd024hdfpn = a9trimd024hdf./mean(a9trimd024hdf(23986:35978,:),1);
%a10trimd024hdfpn = a10trimd024hdf./mean(a10trimd024hdf(23986:35978,:),1);
a11trimd024hdfpn = a11trimd024hdf./mean(a11trimd024hdf(23986:35978,:),1);
a12trimd024hdfpn = a12trimd024hdf./mean(a12trimd024hdf(23986:35978,:),1);

%% d0 traces in various combinations
% for subtracting linear pre-trend (just do this on d1 and beyond lights
% off) - consider using just the first 1.5 hours of the window, or expand
% the window to be more than -3 hours from start
all470d0rawttr = [a1trimd024hdf(:,[2 6 10]) a2trimd024hdf(:,[2 6 10]) a3trimd024hdf(:,[2 6 10])...
    a4trimd024hdf(:,[2 6 10]) a5trimd024hdf(:,[2 6 10]) a6trimd024hdf(:,[2 6 10]) a7trimd024hdf(:,[2 6 10])...
    a8trimd024hdf(:,[2 6 10]) a9trimd024hdf(:,[2 6 10])];
all470d0rawwt = [a11trimd024hdf(:,[2 6 10]) a12trimd024hdf(:,[2 6 10])];
excl470d0rawttr = [a1trimd024hdf(:,[2 6 10]) a2trimd024hdf(:,[2 6 10]) a3trimd024hdf(:,[2 6 10])...
    a4trimd024hdf(:,[2 6 10]) a5trimd024hdf(:,[6 10]) a6trimd024hdf(:,[2 6 10]) a7trimd024hdf(:,[2 6 10])...
    a8trimd024hdf(:,[2 6 10]) a9trimd024hdf(:,[2 6 10])];
excl470d0rawwt = [a11trimd024hdf(:,[2 6 10]) a12trimd024hdf(:,[6 10])];

all470d0rawttravg = mean(all470d0rawttr,2);
all470d0rawwtavg = mean(all470d0rawwt,2);
excl470d0rawttravg = mean(excl470d0rawttr,2);
excl470d0rawwtavg = mean(excl470d0rawwt,2);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
f470d0ttr = fit((1:length(excl470d0rawttravg))',excl470d0rawttravg,g);
f470d0wt = fit((1:length(excl470d0rawwtavg))',excl470d0rawwtavg,g);

%%
all415d0rawttr = [a1trimd024hdf(:,[3 7 11]) a2trimd024hdf(:,[3 7 11]) a3trimd024hdf(:,[3 7 11])...
    a4trimd024hdf(:,[3 7 11]) a5trimd024hdf(:,[3 7 11]) a6trimd024hdf(:,[3 7 11]) a7trimd024hdf(:,[3 7 11])...
    a8trimd024hdf(:,[3 7 11]) a9trimd024hdf(:,[3 7 11])];
all415d0rawwt = [a11trimd024hdf(:,[3 7 11]) a12trimd024hdf(:,[3 7 11])];
excl415d0rawttr = [a1trimd024hdf(:,[3 7 11]) a2trimd024hdf(:,[3 7 11]) a3trimd024hdf(:,[3 7 11])...
    a4trimd024hdf(:,[3 7 11]) a5trimd024hdf(:,[7 11]) a6trimd024hdf(:,[3 7 11]) a7trimd024hdf(:,[3 7 11])...
    a8trimd024hdf(:,[3 7 11]) a9trimd024hdf(:,[3 7 11])];
excl415d0rawwt = [a11trimd024hdf(:,[3 7 11]) a12trimd024hdf(:,[7 11])];

all415d0rawttravg = mean(all415d0rawttr,2);
all415d0rawwtavg = mean(all415d0rawwt,2);
excl415d0rawttravg = mean(excl415d0rawttr,2);
excl415d0rawwtavg = mean(excl415d0rawwt,2);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
f415d0ttr = fit((1:length(excl415d0rawttravg))',excl415d0rawttravg,g);
f415d0wt = fit((1:length(excl415d0rawwtavg))',excl415d0rawwtavg,g);

%%
all470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[4 8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12])];
all470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[4 8 12])];
excl470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12])];
excl470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[8 12])];

all470d0bleachttrpn = [a1trimd024hdfpn(:,[4 8 12]) a2trimd024hdfpn(:,[4 8 12]) a3trimd024hdfpn(:,[4 8 12])...
    a4trimd024hdfpn(:,[4 8 12]) a5trimd024hdfpn(:,[4 8 12]) a6trimd024hdfpn(:,[4 8 12]) a7trimd024hdfpn(:,[4 8 12])...
    a8trimd024hdfpn(:,[4 8 12]) a9trimd024hdfpn(:,[4 8 12])];
all470d0bleachwtpn = [a11trimd024hdfpn(:,[2 8 12]) a12trimd024hdfpn(:,[4 8 12])];
excl470d0bleachttrpn = [a1trimd024hdfpn(:,[4 8 12]) a2trimd024hdfpn(:,[4 8 12]) a3trimd024hdfpn(:,[4 8 12])...
    a4trimd024hdfpn(:,[4 8 12]) a5trimd024hdfpn(:,[8 12]) a6trimd024hdfpn(:,[4 8 12]) a7trimd024hdfpn(:,[4 8 12])...
    a8trimd024hdfpn(:,[4 8 12]) a9trimd024hdfpn(:,[4 8 12])];
excl470d0bleachwtpn = [a11trimd024hdfpn(:,[2 8 12]) a12trimd024hdfpn(:,[8 12])];


all470d0bleachttravg = mean(all470d0bleachttr,2);
all470d0bleachwtavg = mean(all470d0bleachwt,2);
excl470d0bleachttravg = mean(excl470d0bleachttr,2);
excl470d0bleachwtavg = mean(excl470d0bleachwt,2);

%% another way to get peak normalized traces
all470d0bleachttrpn2 = all470d0bleachttr./mean(all470d0bleachttr(23986:35978,:),1);
all470d0bleachwtpn2 = all470d0bleachwt./mean(all470d0bleachwt(23986:35978,:),1);
excl470d0bleachttrpn2 = excl470d0bleachttr./mean(excl470d0bleachttr(23986:35978,:),1);
excl470d0bleachwtpn2 = excl470d0bleachwt./mean(excl470d0bleachwt(23986:35978,:),1);

[~,idxttr] = sort(-mean(excl470d0bleachttr(23986:35978,:),1),2);
excl470d0bleachttrpn2s = excl470d0bleachttrpn2(:,idxttr);
excl470d0bleachttrpn2sd = excl470d0bleachttrpn2s(1:49:47972,:);

[~,idxwt] = sort(-mean(excl470d0bleachwt(23986:35978,:),1),2);
excl470d0bleachwtpn2s = excl470d0bleachwtpn2(:,idxwt);
excl470d0bleachwtpn2sd = excl470d0bleachwtpn2s(1:49:47972,:);

%% dark-dark heatmaps
% peak-normalize
a8trimd124hdfpn = a8trimd124hdf./abs(mean(a8trimd124hdf(23986:35978,:),1));
a9trimd124hdfpn = a9trimd124hdf./abs(mean(a9trimd124hdf(23986:35978,:),1));
a10trimd124hdfpn = a10trimd124hdf./abs(mean(a10trimd124hdf(23986:35978,:),1));
a8trimd224hdfpn = a8trimd224hdf./abs(mean(a8trimd224hdf(23986:35978,:),1));
a9trimd224hdfpn = a9trimd224hdf./abs(mean(a9trimd224hdf(23986:35978,:),1));
a10trimd224hdfpn = a10trimd224hdf./abs(mean(a10trimd224hdf(23986:35978,:),1));
a8trimd324hdfpn = a8trimd224hdf./abs(mean(a8trimd324hdf(23986:35978,:),1));
a9trimd324hdfpn = a9trimd224hdf./abs(mean(a9trimd324hdf(23986:35978,:),1));
a10trimd324hdfpn = a10trimd224hdf./abs(mean(a10trimd324hdf(23986:35978,:),1));

% combine
all470d1ddbleachttrpn = [a8trimd124hdfpn(:,[4 8 12]) a9trimd124hdfpn(:,[4 8 12]) a10trimd124hdfpn(:,[4 8 12])];
all470d2ddbleachttrpn = [a8trimd224hdfpn(:,[4 8 12]) a9trimd224hdfpn(:,[4 8 12]) a10trimd224hdfpn(:,[4 8 12])];
all470d3ddbleachttrpn = [a8trimd324hdfpn(:,[4 8 12]) a9trimd324hdfpn(:,[4 8 12]) a10trimd324hdfpn(:,[4 8 12])];

% subsample
all470d1ddbleachttrpnd = all470d1ddbleachttrpn(1:49:47972,:);
all470d2ddbleachttrpnd = all470d2ddbleachttrpn(1:49:47972,:);
all470d3ddbleachttrpnd = all470d3ddbleachttrpn(1:49:47972,:);

%%
all415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[5 9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13])];
all415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[5 9 13])];
excl415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13])];
excl415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[9 13])];

all415d0bleachttravg = mean(all415d0bleachttr,2);
all415d0bleachwtavg = mean(all415d0bleachwt,2);
excl415d0bleachttravg = mean(excl415d0bleachttr,2);
excl415d0bleachwtavg = mean(excl415d0bleachwt,2);

%% day 1 lights off (10x TTR-mNeon vs 2x WT)
idx11 = find(a1(:,1)>6.12E7)+143912;
idx12 = find(a2(:,1)>6.12E7)+143912;
idx13 = find(a3(:,1)>6.12E7)+143912;
idx14 = find(a4(:,1)>6.12E7)+143912;
idx15 = find(a5(:,1)>6.12E7)+143912;
idx16 = find(a6(:,1)>6.12E7)+143912;
idx17 = find(a7(:,1)>6.12E7)+143912;
% recall that runs 8-10 are dark-dark, so first two days are true lights off and
% the next two days are dark-dark "anticipated lights off"
idx18 = find(a8(:,1)>6.12E7)+143912;
idx19 = find(a9(:,1)>6.12E7)+143912;
idx110 = find(a10(:,1)>6.12E7); % this is the one where first lights off is actually second!!
idx111 = find(a11(35052:end,1)>3.6E6)+35051+143912;
idx112 = find(a12(36982:end,1)>3.6E6)+36981+143912;

% looks like it's actualy 143912 readings per 24 hours, not 144000
% runs 1-3 only had one lights off!
%a1trimd124h = a1(idx11(1):(idx11(1)+143911),:);
%a2trimd124h = a2(idx12(1):(idx12(1)+143911),:);
%a3trimd124h = a3(idx13(1):(idx13(1)+143911),:);
a4trimd124h = a4(idx14(1):(idx14(1)+143911),:);
a5trimd124h = a5(idx15(1):(idx15(1)+143911),:);
a6trimd124h = a6(idx16(1):(idx16(1)+143911),:);
a7trimd124h = a7(idx17(1):(idx17(1)+143911),:);
a8trimd124h = a8(idx18(1):(idx18(1)+143911),:);
a9trimd124h = a9(idx19(1):(idx19(1)+143911),:);
a10trimd124h = a10(idx110(1):(idx110(1)+143911),:);
a11trimd124h = a11(idx111(1):(idx111(1)+143911),:);
a12trimd124h = a12(idx112(1):(idx112(1)+143911),:);

%% NEW 6/27/2022 calcuate dF/F realtive to first hour mean
% F is average of the RAW signal (not the fit - I think)
% idea is to pull out column X-mean(column X)/mean(column X-2) for X = 4,8,12
% this seems to really blunt the lights off bump
for i = 1:3
    %a1trimd124hdff(:,i) = (a1trimd124h(:,i*4)-mean(a1trimd124h(1:5996,i*4),1))/mean(a1trimd124h(1:5996,(i*4-2)),1);
    %a2trimd124hdff(:,i) = (a2trimd124h(:,i*4)-mean(a2trimd124h(1:5996,i*4),1))/mean(a2trimd124h(1:5996,(i*4-2)),1);
    %a3trimd124hdff(:,i) = (a3trimd124h(:,i*4)-mean(a3trimd124h(1:5996,i*4),1))/mean(a3trimd124h(1:5996,(i*4-2)),1);
    a4trimd124hdff(:,i) = (a4trimd124h(:,i*4)-mean(a4trimd124h(1:5996,i*4),1))/mean(a4trimd124h(1:5996,(i*4-2)),1);
    a5trimd124hdff(:,i) = (a5trimd124h(:,i*4)-mean(a5trimd124h(1:5996,i*4),1))/mean(a5trimd124h(1:5996,(i*4-2)),1);
    a6trimd124hdff(:,i) = (a6trimd124h(:,i*4)-mean(a6trimd124h(1:5996,i*4),1))/mean(a6trimd124h(1:5996,(i*4-2)),1);
    a7trimd124hdff(:,i) = (a7trimd124h(:,i*4)-mean(a7trimd124h(1:5996,i*4),1))/mean(a7trimd124h(1:5996,(i*4-2)),1);
    a8trimd124hdff(:,i) = (a8trimd124h(:,i*4)-mean(a8trimd124h(1:5996,i*4),1))/mean(a8trimd124h(1:5996,(i*4-2)),1);
    a9trimd124hdff(:,i) = (a9trimd124h(:,i*4)-mean(a9trimd124h(1:5996,i*4),1))/mean(a9trimd124h(1:5996,(i*4-2)),1);
    a10trimd124hdff(:,i) = (a10trimd124h(:,i*4)-mean(a10trimd124h(1:5996,i*4),1))/mean(a10trimd124h(1:5996,(i*4-2)),1);
    a11trimd124hdff(:,i) = (a11trimd124h(:,i*4)-mean(a11trimd124h(1:5996,i*4),1))/mean(a11trimd124h(1:5996,(i*4-2)),1);
    a12trimd124hdff(:,i) = (a12trimd124h(:,i*4)-mean(a12trimd124h(1:5996,i*4),1))/mean(a12trimd124h(1:5996,(i*4-2)),1);
end

%% NEW 6/27/2022 mean dF/F traces (bleach corrected)

all470d1bleachttrdff = [a4trimd124hdff a5trimd124hdff a6trimd124hdff a7trimd124hdff...
    a8trimd124hdff a9trimd124hdff a10trimd124hdff];
all470d1bleachwtdff = [a11trimd124hdff a12trimd124hdff];
excl470d1bleachttrdff = [a4trimd124hdff a5trimd124hdff(:,2) a6trimd124hdff a7trimd124hdff...
    a8trimd124hdff a9trimd124hdff a10trimd124hdff];
excl470d1bleachwtdff = [a11trimd124hdff a12trimd124hdff(:,[2 3])];

%% calculate df relative to first one hour mean
% lights off at 17990 (idx(1) + 17988)!
% lights on at 89946 (idx(1) + 89945)!
%a1trimd124hdf = a1trimd124h-mean(a1trimd124h(1:5996,:),1);
%a2trimd124hdf = a2trimd124h-mean(a2trimd124h(1:5996,:),1);
%a3trimd124hdf = a3trimd124h-mean(a3trimd124h(1:5996,:),1);
a4trimd124hdf = a4trimd124h-mean(a4trimd124h(1:5996,:),1);
a5trimd124hdf = a5trimd124h-mean(a5trimd124h(1:5996,:),1);
a6trimd124hdf = a6trimd124h-mean(a6trimd124h(1:5996,:),1);
a7trimd124hdf = a7trimd124h-mean(a7trimd124h(1:5996,:),1);
a8trimd124hdf = a8trimd124h-mean(a8trimd124h(1:5996,:),1);
a9trimd124hdf = a9trimd124h-mean(a9trimd124h(1:5996,:),1);
a10trimd124hdf = a10trimd124h-mean(a10trimd124h(1:5996,:),1);
a11trimd124hdf = a11trimd124h-mean(a11trimd124h(1:5996,:),1);
a12trimd124hdf = a12trimd124h-mean(a12trimd124h(1:5996,:),1);

%% d1 traces in various combinations
% for subtracting linear pre-trend (just do this on d1 and beyond lights
% off) - consider using just the first 1.5 hours of the window, or expand
% the window to be more than -3 hours from start
all470d1rawttr = [a4trimd124hdf(:,[2 6 10]) a5trimd124hdf(:,[2 6 10]) a6trimd124hdf(:,[2 6 10]) a7trimd124hdf(:,[2 6 10])...
    a8trimd124hdf(:,[2 6 10]) a9trimd124hdf(:,[2 6 10]) a10trimd124hdf(:,[2 6 10])];
all470d1rawwt = [a11trimd124hdf(:,[2 6 10]) a12trimd124hdf(:,[2 6 10])];
excl470d1rawttr = [a4trimd124hdf(:,[2 6 10]) a5trimd124hdf(:,6) a6trimd124hdf(:,[2 6 10]) a7trimd124hdf(:,[2 6 10])...
    a8trimd124hdf(:,[2 6 10]) a9trimd124hdf(:,[2 6 10]) a10trimd124hdf(:,[2 6 10])];
excl470d1rawwt = [a11trimd124hdf(:,[2 6 10]) a12trimd124hdf(:,[6 10])];

all470d1rawttravg = mean(all470d1rawttr,2);
all470d1rawwtavg = mean(all470d1rawwt,2);
excl470d1rawttravg = mean(excl470d1rawttr,2);
excl470d1rawwtavg = mean(excl470d1rawwt,2);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
f470d1ttr = fit((1:length(excl470d1rawttravg))',excl470d1rawttravg,g);
f470d1wt = fit((1:length(excl470d1rawwtavg))',excl470d1rawwtavg,g);

%% NOT UPDATED FOR d1
all415d0rawttr = [a1trimd024hdf(:,[3 7 11]) a2trimd024hdf(:,[3 7 11]) a3trimd024hdf(:,[3 7 11])...
    a4trimd024hdf(:,[3 7 11]) a5trimd024hdf(:,[3 7 11]) a6trimd024hdf(:,[3 7 11]) a7trimd024hdf(:,[3 7 11])...
    a8trimd024hdf(:,[3 7 11]) a9trimd024hdf(:,[3 7 11]) a10trimd024hdf(:,[3 7 11])];
all415d0rawwt = [a11trimd024hdf(:,[3 7 11]) a12trimd024hdf(:,[3 7 11])];
excl415d0rawttr = [a1trimd024hdf(:,[3 7 11]) a2trimd024hdf(:,[3 7 11]) a3trimd024hdf(:,[3 7 11])...
    a4trimd024hdf(:,[3 7 11]) a5trimd024hdf(:,[7 11]) a6trimd024hdf(:,[3 7 11]) a7trimd024hdf(:,[3 7 11])...
    a8trimd024hdf(:,[3 7 11]) a9trimd024hdf(:,[3 7 11]) a10trimd024hdf(:,[3 7 11])];
excl415d0rawwt = [a11trimd024hdf(:,[3 7 11]) a12trimd024hdf(:,[7 11])];

all415d0rawttravg = mean(all415d0rawttr,2);
all415d0rawwtavg = mean(all415d0rawwt,2);
excl415d0rawttravg = mean(excl415d0rawttr,2);
excl415d0rawwtavg = mean(excl415d0rawwt,2);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-10000000,0,0,0],...
               'Upper',[10000000,10000000,1000000000,1],...
               'StartPoint',[1 1 1 1]);
g = fittype('a+(b*exp(-c*(x^d)/1000000))','options',fo);
f415d0ttr = fit((1:length(excl415d0rawttravg))',excl415d0rawttravg,g);
f415d0wt = fit((1:length(excl415d0rawwtavg))',excl415d0rawwtavg,g);

%%
all470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[4 8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12]) a10trimd024hdf(:,[4 8 12])];
all470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[4 8 12])];
excl470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12]) a10trimd024hdf(:,[4 8 12])];
excl470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[8 12])];

all470d0bleachttravg = mean(all470d0bleachttr,2);
all470d0bleachwtavg = mean(all470d0bleachwt,2);
excl470d0bleachttravg = mean(excl470d0bleachttr,2);
excl470d0bleachwtavg = mean(excl470d0bleachwt,2);

%%
all415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[5 9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13]) a10trimd024hdf(:,[5 9 13])];
all415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[5 9 13])];
excl415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13]) a10trimd024hdf(:,[5 9 13])];
excl415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[9 13])];

all415d0bleachttravg = mean(all415d0bleachttr,2);
all415d0bleachwtavg = mean(all415d0bleachwt,2);
excl415d0bleachttravg = mean(excl415d0bleachttr,2);
excl415d0bleachwtavg = mean(excl415d0bleachwt,2);

%% day 2 lights off (10x TTR-mNeon vs 2x WT) DARK-DARK

% recall that runs 8-10 are dark-dark, so first two days are true lights off and
% the next two days are dark-dark "anticipated lights off"
idx28 = find(a8(:,1)>6.12E7)+143912+143912;
idx29 = find(a9(:,1)>6.12E7)+143912+143912;
idx210 = find(a10(:,1)>6.12E7)+143912; % this is the one where first lights off is actually second!!

% looks like it's actualy 143912 readings per 24 hours, not 144000
% runs 1-3 only had one lights off!

a8trimd224h = a8(idx28(1):(idx28(1)+143911),:);
a9trimd224h = a9(idx29(1):(idx29(1)+143911),:);
a10trimd224h = a10(idx210(1):(idx210(1)+143911),:);

%% NEW 6/27/2022 calcuate dF/F realtive to first hour mean
% F is average of the RAW signal (not the fit - I think)
% idea is to pull out column X-mean(column X)/mean(column X-2) for X = 4,8,12
% this seems to really blunt the lights off bump
for i = 1:3
    a8trimd224hdff(:,i) = (a8trimd224h(:,i*4)-mean(a8trimd224h(1:5996,i*4),1))/mean(a8trimd224h(1:5996,(i*4-2)),1);
    a9trimd224hdff(:,i) = (a9trimd224h(:,i*4)-mean(a9trimd224h(1:5996,i*4),1))/mean(a9trimd224h(1:5996,(i*4-2)),1);
    a10trimd224hdff(:,i) = (a10trimd224h(:,i*4)-mean(a10trimd224h(1:5996,i*4),1))/mean(a10trimd224h(1:5996,(i*4-2)),1);
end

%% NEW 6/27/2022 mean dF/F traces (bleach corrected)

all470d1bleachttrdff = [a8trimd224hdff a9trimd224hdff a10trimd224hdff];
excl470d1bleachttrdff = [a8trimd224hdff a9trimd224hdff a10trimd224hdff];
%% calculate df relative to first two hour mean
% lights off at 17990 (idx(1) + 17988)!
% lights on at 89946 (idx(1) + 89945)!
a8trimd224hdf = a8trimd224h-mean(a8trimd224h(1:11993,:),1);
a9trimd224hdf = a9trimd224h-mean(a9trimd224h(1:11993,:),1);
a10trimd224hdf = a10trimd224h-mean(a10trimd224h(1:11993,:),1);

%% d2 DARK-DARK traces in various combinations
% for subtracting linear pre-trend (just do this on d1 and beyond lights
% off) - consider using just the first 1.5 hours of the window, or expand
% the window to be more than -3 hours from start
all470d2ddrawttr = [a8trimd224hdf(:,[2 6 10]) a9trimd224hdf(:,[2 6 10]) a10trimd224hdf(:,[2 6 10])];
excl470d2ddrawttr = [a8trimd224hdf(:,[2 6 10]) a9trimd224hdf(:,[2 6 10]) a10trimd224hdf(:,[2 6 10])];

all470d2ddrawttravg = mean(all470d2ddrawttr,2);
excl470d2ddrawttravg = mean(excl470d2ddrawttr,2);


%% 
all415d2ddrawttr = [a8trimd224hdf(:,[3 7 11]) a9trimd224hdf(:,[3 7 11]) a10trimd224hdf(:,[3 7 11])];
excl415d2ddrawttr = [a8trimd224hdf(:,[3 7 11]) a9trimd224hdf(:,[3 7 11]) a10trimd224hdf(:,[3 7 11])];

all415d2ddrawttravg = mean(all415d2ddrawttr,2);
excl415d2ddrawttravg = mean(excl415d2ddrawttr,2);

%% NOT UPDATED FOR d2 DARK-DARK
all470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[4 8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12]) a10trimd024hdf(:,[4 8 12])];
all470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[4 8 12])];
excl470d0bleachttr = [a1trimd024hdf(:,[4 8 12]) a2trimd024hdf(:,[4 8 12]) a3trimd024hdf(:,[4 8 12])...
    a4trimd024hdf(:,[4 8 12]) a5trimd024hdf(:,[8 12]) a6trimd024hdf(:,[4 8 12]) a7trimd024hdf(:,[4 8 12])...
    a8trimd024hdf(:,[4 8 12]) a9trimd024hdf(:,[4 8 12]) a10trimd024hdf(:,[4 8 12])];
excl470d0bleachwt = [a11trimd024hdf(:,[2 8 12]) a12trimd024hdf(:,[8 12])];

all470d0bleachttravg = mean(all470d0bleachttr,2);
all470d0bleachwtavg = mean(all470d0bleachwt,2);
excl470d0bleachttravg = mean(excl470d0bleachttr,2);
excl470d0bleachwtavg = mean(excl470d0bleachwt,2);

%%
all415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[5 9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13]) a10trimd024hdf(:,[5 9 13])];
all415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[5 9 13])];
excl415d0bleachttr = [a1trimd024hdf(:,[5 9 13]) a2trimd024hdf(:,[5 9 13]) a3trimd024hdf(:,[5 9 13])...
    a4trimd024hdf(:,[5 9 13]) a5trimd024hdf(:,[9 13]) a6trimd024hdf(:,[5 9 13]) a7trimd024hdf(:,[5 9 13])...
    a8trimd024hdf(:,[5 9 13]) a9trimd024hdf(:,[5 9 13]) a10trimd024hdf(:,[5 9 13])];
excl415d0bleachwt = [a11trimd024hdf(:,[3 9 13]) a12trimd024hdf(:,[9 13])];

all415d0bleachttravg = mean(all415d0bleachttr,2);
all415d0bleachwtavg = mean(all415d0bleachwt,2);
excl415d0bleachttravg = mean(excl415d0bleachttr,2);
excl415d0bleachwtavg = mean(excl415d0bleachwt,2);

%% day 3 lights off (10x TTR-mNeon vs 2x WT) DARK-DARK

% recall that runs 8-10 are dark-dark, so first two days are true lights off and
% the next two days are dark-dark "anticipated lights off"
% here just leave in 12 hours post!!
idx38 = find(a8(:,1)>6.12E7)+143912+143912+143912;
idx39 = find(a9(:,1)>6.12E7)+143912+143912+143912;
idx310 = find(a10(:,1)>6.12E7)+143912+143912; % this is the one where first lights off is actually second!!

% looks like it's actualy 143912 readings per 24 hours, not 144000
% runs 1-3 only had one lights off!

a8trimd324h = a8(idx38(1):(idx38(1)+71956),:);
a9trimd324h = a9(idx39(1):(idx39(1)+71956),:);
a10trimd324h = a10(idx310(1):(idx310(1)+71956),:);

%% calculate df relative to first two hour mean
% lights off at 17990 (idx(1) + 17988)!
% lights on at 89946 (idx(1) + 89945)!
a8trimd324hdf = a8trimd324h-mean(a8trimd324h(1:11993,:),1);
a9trimd324hdf = a9trimd324h-mean(a9trimd324h(1:11993,:),1);
a10trimd324hdf = a10trimd324h-mean(a10trimd324h(1:11993,:),1);

%% d3 DARK-DARK traces in various combinations
% for subtracting linear pre-trend (just do this on d1 and beyond lights
% off) - consider using just the first 1.5 hours of the window, or expand
% the window to be more than -3 hours from start
all470d3ddrawttr = [a8trimd324hdf(:,[2 6 10]) a9trimd324hdf(:,[2 6 10]) a10trimd324hdf(:,[2 6 10])];
excl470d3ddrawttr = [a8trimd324hdf(:,[2 6 10]) a9trimd324hdf(:,[2 6 10]) a10trimd324hdf(:,[2 6 10])];

all470d3ddrawttravg = mean(all470d3ddrawttr,2);
excl470d3ddrawttravg = mean(excl470d3ddrawttr,2);

%%
all415d3ddrawttr = [a8trimd324hdf(:,[3 7 11]) a9trimd324hdf(:,[3 7 11]) a10trimd324hdf(:,[3 7 11])];
excl415d3ddrawttr = [a8trimd324hdf(:,[3 7 11]) a9trimd324hdf(:,[3 7 11]) a10trimd324hdf(:,[3 7 11])];

all415d3ddrawttravg = mean(all415d3ddrawttr,2);
excl415d3ddrawttravg = mean(excl415d3ddrawttr,2);



%% day 0 lights off z-scoring
idx01 = find(a1(:,1)>6.12E7);
idx02 = find(a2(:,1)>6.12E7);
idx03 = find(a3(:,1)>6.12E7);
idx04 = find(a4(:,1)>6.12E7);
idx05 = find(a5(:,1)>6.12E7);
idx06 = find(a6(:,1)>6.12E7);
idx07 = find(a7(:,1)>6.12E7);

% z score calculation
a1trimd024h = a1(idx01(1):(idx01(1)+143999),:);
a2trimd024h = a2(idx02(1):(idx02(1)+143999),:);
a3trimd024h = a3(idx03(1):(idx03(1)+143999),:);
a4trimd024h = a4(idx04(1):(idx04(1)+143999),:);
a5trimd024h = a5(idx05(1):(idx05(1)+143999),:);
a6trimd024h = a6(idx06(1):(idx06(1)+143999),:);
a7trimd024h = a7(idx07(1):(idx07(1)+143999),:);

for i = 2:4
    a1trimd024h(:,i) = (a1trimd024h(:,i)-mean(a1trimd024h(:,i)))/std(a1trimd024h(:,i));
    a2trimd024h(:,i) = (a2trimd024h(:,i)-mean(a2trimd024h(:,i)))/std(a2trimd024h(:,i));
    a3trimd024h(:,i) = (a3trimd024h(:,i)-mean(a3trimd024h(:,i)))/std(a3trimd024h(:,i));
    a4trimd024h(:,i) = (a4trimd024h(:,i)-mean(a4trimd024h(:,i)))/std(a4trimd024h(:,i));
    a5trimd024h(:,i) = (a5trimd024h(:,i)-mean(a5trimd024h(:,i)))/std(a5trimd024h(:,i));
    a6trimd024h(:,i) = (a6trimd024h(:,i)-mean(a6trimd024h(:,i)))/std(a6trimd024h(:,i));
    a7trimd024h(:,i) = (a7trimd024h(:,i)-mean(a7trimd024h(:,i)))/std(a7trimd024h(:,i));
end

a1trimd0zpr = a1trimd024h(1:48000,:);
a2trimd0zpr = a2trimd024h(1:48000,:);
a3trimd0zpr = a3trimd024h(1:48000,:);
a4trimd0zpr = a4trimd024h(1:48000,:);
a5trimd0zpr = a5trimd024h(1:48000,:);
a6trimd0zpr = a6trimd024h(1:48000,:);
a7trimd0zpr = a7trimd024h(1:48000,:);

% dz calculation
a1trimd0z = a1trimd0zpr - mean(a1trimd0zpr(1:12000,:),1);
a2trimd0z = a2trimd0zpr - mean(a2trimd0zpr(1:12000,:),1);
a3trimd0z = a3trimd0zpr - mean(a3trimd0zpr(1:12000,:),1);
a4trimd0z = a4trimd0zpr - mean(a4trimd0zpr(1:12000,:),1);
a5trimd0z = a5trimd0zpr - mean(a5trimd0zpr(1:12000,:),1);
a6trimd0z = a6trimd0zpr - mean(a6trimd0zpr(1:12000,:),1);
a7trimd0z = a7trimd0zpr - mean(a7trimd0zpr(1:12000,:),1);

m1d0z = [a1trimd0z(:,2) a2trimd0z(:,2) a3trimd0z(:,2) a4trimd0z(:,2) a5trimd0z(:,2) a6trimd0z(:,2) a7trimd0z(:,2)];
m2d0z = [a1trimd0z(:,3) a2trimd0z(:,3) a3trimd0z(:,3) a4trimd0z(:,3) a5trimd0z(:,3) a6trimd0z(:,3) a7trimd0z(:,3)];
m3d0z = [a1trimd0z(:,4) a2trimd0z(:,4) a3trimd0z(:,4) a4trimd0z(:,4) a5trimd0z(:,4) a6trimd0z(:,4) a7trimd0z(:,4)];
malld0z = [m1d0z m2d0z m3d0z];
% with exclusions
m1d0zexcl = [a1trimd0z(:,2) a3trimd0z(:,2) a4trimd0z(:,2) a5trimd0z(:,2) a6trimd0z(:,2) a7trimd0z(:,2)];
m2d0zexcl = [a1trimd0z(:,3) a2trimd0z(:,3) a3trimd0z(:,3) a4trimd0z(:,3) a5trimd0z(:,3) a6trimd0z(:,3) a7trimd0z(:,3)];
m3d0zexcl = [a2trimd0z(:,4) a3trimd0z(:,4) a4trimd0z(:,4) a5trimd0z(:,4) a6trimd0z(:,4) a7trimd0z(:,4)];

%m1d0zexcl = [a1trimd0z(:,2) a3trimd0z(:,2) a4trimd0z(:,2) a5trimd0z(:,2) a6trimd0z(:,2) a7trimd0z(:,2)];
%m2d0zexcl = [a1trimd0z(:,3) a2trimd0z(:,3) a3trimd0z(:,3) a4trimd0z(:,3) a5trimd0z(:,3) a6trimd0z(:,3) a7trimd0z(:,3)];
%m3d0zexcl = [a2trimd0z(:,4) a3trimd0z(:,4) a4trimd0z(:,4) a5trimd0z(:,4) a6trimd0z(:,4) a7trimd0z(:,4)];
malld0zexcl = [m1d0zexcl m2d0zexcl m3d0zexcl];

%% day 0 lights on z-scoring
a1trimd0znpr = a1trimd024h(72001:120000,:);
a2trimd0znpr = a2trimd024h(72001:120000,:);
a3trimd0znpr = a3trimd024h(72001:120000,:);
a4trimd0znpr = a4trimd024h(72001:120000,:);
a5trimd0znpr = a5trimd024h(72001:120000,:);
a6trimd0znpr = a6trimd024h(72001:120000,:);
a7trimd0znpr = a7trimd024h(72001:120000,:);

% dz calculation
a1trimd0zn = a1trimd0znpr - mean(a1trimd0znpr(1:12000,:),1);
a2trimd0zn = a2trimd0znpr - mean(a2trimd0znpr(1:12000,:),1);
a3trimd0zn = a3trimd0znpr - mean(a3trimd0znpr(1:12000,:),1);
a4trimd0zn = a4trimd0znpr - mean(a4trimd0znpr(1:12000,:),1);
a5trimd0zn = a5trimd0znpr - mean(a5trimd0znpr(1:12000,:),1);
a6trimd0zn = a6trimd0znpr - mean(a6trimd0znpr(1:12000,:),1);
a7trimd0zn = a7trimd0znpr - mean(a7trimd0znpr(1:12000,:),1);

m1d0zn = [a1trimd0zn(:,2) a2trimd0zn(:,2) a3trimd0zn(:,2) a4trimd0zn(:,2) a5trimd0zn(:,2) a6trimd0zn(:,2) a7trimd0zn(:,2)];
m2d0zn = [a1trimd0zn(:,3) a2trimd0zn(:,3) a3trimd0zn(:,3) a4trimd0zn(:,3) a5trimd0zn(:,3) a6trimd0zn(:,3) a7trimd0zn(:,3)];
m3d0zn = [a1trimd0zn(:,4) a2trimd0zn(:,4) a3trimd0zn(:,4) a4trimd0zn(:,4) a5trimd0zn(:,4) a6trimd0zn(:,4) a7trimd0zn(:,4)];
malld0zn = [m1d0zn m2d0zn m3d0zn];

%% full 24 hour z scoring
% dz calculation
a1trimd024hz = a1trimd024h - mean(a1trimd024h(1:12000,:),1);
a2trimd024hz = a2trimd024h - mean(a2trimd024h(1:12000,:),1);
a3trimd024hz = a3trimd024h - mean(a3trimd024h(1:12000,:),1);
a4trimd024hz = a4trimd024h - mean(a4trimd024h(1:12000,:),1);
a5trimd024hz = a5trimd024h - mean(a5trimd024h(1:12000,:),1);
a6trimd024hz = a6trimd024h - mean(a6trimd024h(1:12000,:),1);
a7trimd024hz = a7trimd024h - mean(a7trimd024h(1:12000,:),1);

%% day 0 lights off
idx01 = find(a1(:,1)>6.12E7);
idx02 = find(a2(:,1)>6.12E7);
idx03 = find(a3(:,1)>6.12E7);
idx04 = find(a4(:,1)>6.12E7);
idx05 = find(a5(:,1)>6.12E7);
idx06 = find(a6(:,1)>6.12E7);
idx07 = find(a7(:,1)>6.12E7);

a1trimd0pr = a1(idx01(1):(idx01(1)+47999),:);
a2trimd0pr = a2(idx02(1):(idx02(1)+47999),:);
a3trimd0pr = a3(idx03(1):(idx03(1)+47999),:);
a4trimd0pr = a4(idx04(1):(idx04(1)+47999),:);
a5trimd0pr = a5(idx05(1):(idx05(1)+47999),:);
a6trimd0pr = a6(idx06(1):(idx06(1)+47999),:);
a7trimd0pr = a7(idx07(1):(idx07(1)+47999),:);

% df calculation
a1trimd0 = a1trimd0pr - mean(a1trimd0pr(1:12000,:),1);
a2trimd0 = a2trimd0pr - mean(a2trimd0pr(1:12000,:),1);
a3trimd0 = a3trimd0pr - mean(a3trimd0pr(1:12000,:),1);
a4trimd0 = a4trimd0pr - mean(a4trimd0pr(1:12000,:),1);
a5trimd0 = a5trimd0pr - mean(a5trimd0pr(1:12000,:),1);

m1d0 = [a1trimd0(:,2) a2trimd0(:,2) a3trimd0(:,2) a4trimd0(:,2) a5trimd0(:,2)];
m2d0 = [a1trimd0(:,3) a2trimd0(:,3) a3trimd0(:,3) a4trimd0(:,3) a5trimd0(:,3)];
m3d0 = [a1trimd0(:,4) a2trimd0(:,4) a3trimd0(:,4) a4trimd0(:,4) a5trimd0(:,4)];
malld0 = [m1d0 m2d0 m3d0];
% with exclusions
m1d0excl = [a1trimd0(:,2) a3trimd0(:,2) a4trimd0(:,2) a5trimd0(:,2)];
m2d0excl = [a1trimd0(:,3) a2trimd0(:,3) a3trimd0(:,3) a4trimd0(:,3) a5trimd0(:,3)];
m3d0excl = [a2trimd0(:,4) a3trimd0(:,4) a4trimd0(:,4) a5trimd0(:,4)];
malld0excl = [m1d0excl m2d0excl m3d0excl];
%% day 1 lights off
idx11 = find(a1((idx01(1)+42000):end,1)>6.12E7)+idx01(1)+41999;
idx12 = find(a2((idx02(1)+42000):end,1)>6.12E7)+idx02(1)+41999;
idx13 = find(a3((idx03(1)+42000):end,1)>6.12E7)+idx03(1)+41999;
idx14 = find(a4((idx04(1)+42000):end,1)>6.12E7)+idx04(1)+41999;
idx15 = find(a5((idx05(1)+42000):end,1)>6.12E7)+idx05(1)+41999;

%a1trimd1pr = a1(idx11(1):end,:);
%a2trimd1pr = a2(idx12(1):end,:);
a3trimd1pr = a3(idx13(1):(idx13(1)+47999),:);
a4trimd1pr = a4(idx14(1):(idx14(1)+47999),:);
a5trimd1pr = a5(idx15(1):(idx15(1)+47999),:);

% df calculation
%a1trimd1 = a1trimd1pr - mean(a1trimd1pr(1:12000,:),1);
%a2trimd1 = a2trimd1pr - mean(a2trimd1pr(1:12000,:),1);
a3trimd1 = a3trimd1pr - mean(a3trimd1pr(1:12000,:),1);
a4trimd1 = a4trimd1pr - mean(a4trimd1pr(1:12000,:),1);
a5trimd1 = a5trimd1pr - mean(a5trimd1pr(1:12000,:),1);

m1d1 = [a3trimd1(:,2) a4trimd1(:,2) a5trimd1(:,2)];
m2d1 = [a3trimd1(:,3) a4trimd1(:,3) a5trimd1(:,3)];
m3d1 = [a3trimd1(:,4) a4trimd1(:,4) a5trimd1(:,4)];
malld1 = [m1d1 m2d1 m3d1];
% with exclusions
m1d1excl = [a3trimd1(:,2) a4trimd1(:,2) a5trimd1(:,2)];
m2d1excl = [a3trimd1(:,3) a4trimd1(:,3) a5trimd1(:,3)];
m3d1excl = [a3trimd1(:,4) a4trimd1(:,4)];
malld1excl = [m1d1excl m2d1excl m3d1excl];
%% day 2 lights off
idx24 = find(a4((idx14(1)+42000):end,1)>6.12E7)+idx14(1)+41999;
idx25 = find(a5((idx15(1)+42000):end,1)>6.12E7)+idx15(1)+41999;

a4trimd2pr = a4(idx24(1):(idx24(1)+47999),:);
a5trimd2pr = a5(idx25(1):(idx25(1)+47999),:);

% df calculation
a4trimd2 = a4trimd2pr - mean(a4trimd2pr(1:12000,:),1);
a5trimd2 = a5trimd2pr - mean(a5trimd2pr(1:12000,:),1);

m1d2 = [a4trimd2(:,2) a5trimd2(:,2)];
m2d2 = [a4trimd2(:,3) a5trimd2(:,3)];
m3d2 = [a4trimd2(:,4) a5trimd2(:,4)];
malld2 = [m1d2 m2d2 m3d2];

m1d2excl = [a5trimd2(:,2)];
m2d2excl = [a4trimd2(:,3) a5trimd2(:,3)];
m3d2excl = [a4trimd2(:,4)];
malld2excl = [m1d2excl m2d2excl m3d2excl];

%%
%% day 0 lights on
idx01n = find(a1(75000:end,1)>1.80E7)+74999;
idx02n = find(a2(50000:end,1)>1.80E7)+49999;
idx03n = find(a3(62000:end,1)>1.80E7)+61999;
idx04n = find(a4(50000:end,1)>1.80E7)+49999;
idx05n = find(a5(80000:end,1)>1.80E7)+79999;

a1trimd0npr = a1(idx01n(1):(idx01n(1)+47999),:);
a2trimd0npr = a2(idx02n(1):(idx02n(1)+47999),:);
a3trimd0npr = a3(idx03n(1):(idx03n(1)+47999),:);
a4trimd0npr = a4(idx04n(1):(idx04n(1)+47999),:);
a5trimd0npr = a5(idx05n(1):(idx05n(1)+47999),:);

% df calculation
a1trimd0n = a1trimd0npr - mean(a1trimd0npr(1:12000,:),1);
a2trimd0n = a2trimd0npr - mean(a2trimd0npr(1:12000,:),1);
a3trimd0n = a3trimd0npr - mean(a3trimd0npr(1:12000,:),1);
a4trimd0n = a4trimd0npr - mean(a4trimd0npr(1:12000,:),1);
a5trimd0n = a5trimd0npr - mean(a5trimd0npr(1:12000,:),1);

m1d0n = [a1trimd0n(:,2) a2trimd0n(:,2) a3trimd0n(:,2) a4trimd0n(:,2) a5trimd0n(:,2)];
m2d0n = [a1trimd0n(:,3) a2trimd0n(:,3) a3trimd0n(:,3) a4trimd0n(:,3) a5trimd0n(:,3)];
m3d0n = [a1trimd0n(:,4) a2trimd0n(:,4) a3trimd0n(:,4) a4trimd0n(:,4) a5trimd0n(:,4)];
malld0n = [m1d0n m2d0n m3d0n];
% with exclusions
m1d0nexcl = [a1trimd0n(:,2) a2trimd0n(:,2) a3trimd0n(:,2) a4trimd0n(:,2)];
m2d0nexcl = [a2trimd0n(:,3) a3trimd0n(:,3) a4trimd0n(:,3) a5trimd0n(:,3)];
m3d0nexcl = [a1trimd0n(:,4) a2trimd0n(:,4) a3trimd0n(:,4) a4trimd0n(:,4)];
malld0nexcl = [m1d0nexcl m2d0nexcl m3d0nexcl];
%% day 1 lights on
%idx11n = find(a1((idx01n(1)+72000):end,1)>1.80E7)+idx01n(1)+71999;
%idx12n = find(a2((idx02n(1)+72000):end,1)>1.80E7)+idx02n(1)+71999;
%idx13n = find(a3((idx03n(1)+72000):end,1)>1.80E7)+idx03n(1)+71999;
idx14n = find(a4((idx04n(1)+115000):end,1)>1.80E7)+idx04n(1)+114999;
idx15n = find(a5((idx05n(1)+115000):end,1)>1.80E7)+idx05n(1)+114999;

%a1trimd1pr = a1(idx11(1):end,:);
%a2trimd1pr = a2(idx12(1):end,:);
%a3trimd1npr = a3(idx13n(1):(idx13n(1)+47999),:);
a4trimd1npr = a4(idx14n(1):(idx14n(1)+47999),:);
a5trimd1npr = a5(idx15n(1):(idx15n(1)+47999),:);

% df calculation
%a1trimd1 = a1trimd1pr - mean(a1trimd1pr(1:12000,:),1);
%a2trimd1 = a2trimd1pr - mean(a2trimd1pr(1:12000,:),1);
%a3trimd1n = a3trimd1npr - mean(a3trimd1npr(1:12000,:),1);
a4trimd1n = a4trimd1npr - mean(a4trimd1npr(1:12000,:),1);
a5trimd1n = a5trimd1npr - mean(a5trimd1npr(1:12000,:),1);

m1d1n = [a4trimd1n(:,2) a5trimd1n(:,2)];
m2d1n = [a4trimd1n(:,3) a5trimd1n(:,3)];
m3d1n = [a4trimd1n(:,4) a5trimd1n(:,4)];
malld1n = [m1d1n m2d1n m3d1n];
% with exclusions
m1d1nexcl = [a4trimd1n(:,2) a5trimd1n(:,2)];
m2d1nexcl = [a4trimd1n(:,3) a5trimd1n(:,3)];
m3d1nexcl = [a4trimd1n(:,4)];
malld1nexcl = [m1d1nexcl m2d1nexcl m3d1nexcl];

%%
%% delta off-on
% ceiling effect is correct to assess between lights on drop and next
% lights off bump - i.e. does a poor decrease at lights on predict a poor
% increase at next lights off?
% or - should we say that a poor lights off bump predicts a poor lights on
% decrease?

%% d1 off y - d0 on x

a3trimd1m = mean(a3trimd1(36001:48000,2:4),1)';
a4trimd1m = mean(a4trimd1(36001:48000,2:4),1)';
a5trimd1m = mean(a5trimd1(36001:48000,2:3),1)';

a3trimd0nm = mean(a3trimd0n(36001:48000,2:4),1)';
a4trimd0nm = mean(a4trimd0n(36001:48000,2:4),1)';
a5trimd0nm = mean(a5trimd0n(36001:48000,2:3),1)';

a3d10pairs = [a3trimd0nm a3trimd1m];
a4d10pairs = [a4trimd0nm a4trimd1m];
a5d10pairs = [a5trimd0nm a5trimd1m];
alld10pairs = [a3d10pairs; a4d10pairs; a5d10pairs];

%% d2 off y - d1 on x

a4trimd2m = mean(a4trimd2(36001:48000,2:4),1)';
a5trimd2m = mean(a5trimd2(36001:48000,2:3),1)';

a4trimd1nm = mean(a4trimd1n(36001:48000,2:4),1)';
a5trimd1nm = mean(a5trimd1n(36001:48000,2:3),1)';

a4d21pairs = [a4trimd1nm a4trimd2m];
a5d21pairs = [a5trimd1nm a5trimd2m];
alld21pairs = [a4d21pairs; a5d21pairs];

alldpairs = [alld10pairs; alld21pairs];
%% plots
% lights off mouse 1 all sessions
stdshade(m1d0excl')
hold on
xline(18000);
hold off
figure,stdshade(m1d1excl')
hold on
xline(18000);
hold off
figure,stdshade(m1d2excl')
hold on
xline(18000);
hold off

% lights off mouse 2 all sessions
stdshade(m2d0excl')
hold on
xline(18000);
hold off
figure,stdshade(m2d1excl')
hold on
xline(18000);
hold off
figure,stdshade(m2d2excl')
hold on
xline(18000);
hold off

% lights off mouse 3 all sessions
stdshade(m3d0excl')
hold on
xline(18000);
hold off
figure,stdshade(m3d1excl')
hold on
xline(18000);
hold off
figure,stdshade(m3d2excl')
hold on
xline(18000);
hold off

% lights off all mice all sessions
stdshade(malld0excl')
hold on
xline(18000);
hold off
figure,stdshade(malld1excl')
hold on
xline(18000);
hold off
figure,stdshade(malld2excl')
hold on
xline(18000);
hold off

figure,stdshade([malld0excl'; malld1excl'; malld2excl']);
hold on
xline(18000);
hold off

% lights on mouse 1 all sessions
stdshade(m1d0nexcl')
hold on
xline(18000);
hold off
figure,stdshade(m1d1nexcl')
hold on
xline(18000);
hold off

% lights on mouse 2 all sessions
stdshade(m2d0nexcl')
hold on
xline(18000);
hold off
figure,stdshade(m2d1nexcl')
hold on
xline(18000);
hold off

% lights on mouse 3 all sessions
stdshade(m3d0nexcl')
hold on
xline(18000);
hold off
figure,stdshade(m3d1nexcl')
hold on
xline(18000);
hold off

% lights on all mice all sessions
stdshade(malld0nexcl')
hold on
xline(18000);
hold off
figure,stdshade(malld1nexcl')
hold on
xline(18000);
hold off

figure,stdshade([malld0nexcl'; malld1nexcl']);
hold on
xline(18000);
hold off

% lights on drop predicts next lights off rise
%scatter(alldpairs(:,1),alldpairs(:,2),'black')
%hold on
f = fit(alldpairs(:,1),alldpairs(:,2),'poly1');
plot(f,alldpairs(:,1),alldpairs(:,2));

% example timecourses
% run4mouse3
ha = fill([18000 18000 90000 90000],[-3 4 4 -3],[0.8 0.8 0.8]);
hold on
hb = fill([162000 162000 234000 234000],[-3 4 4 -3],[0.8 0.8 0.8]);
hc = fill([306000 306000 378000 378000],[-3 4 4 -3],[0.8 0.8 0.8]);
plot(a4(idx04(1):end,4),'black');
hold off

% run5mouse2
ha = fill([18000 18000 90000 90000],[-5 35 35 -5],[0.8 0.8 0.8]);
hold on
hb = fill([162000 162000 234000 234000],[-5 35 35 -5],[0.8 0.8 0.8]);
hc = fill([306000 306000 378000 378000],[-5 35 35 -5],[0.8 0.8 0.8]);
plot(a5(idx05(1):end,3),'black');
hold off

%% 220418 T32 progress talk plots
ha = fill([-3 -3 0 0],[-3 1 1 -3],[0.8 0.8 0.8]);
hold on
stdshade([malld0nexcl'; malld1nexcl'],0.3,'b',vec);
xlabel('Time relative to lights on (hrs)');
ylabel('Fluorescence (arbitrary units)');
hold off

ha = fill([0 0 5 5],[-1 6 6 -1],[0.8 0.8 0.8]);
hold on
stdshade([malld0excl'; malld1excl'; malld2excl'],0.3,'b',vec);
xlabel('Time relative to lights off (hrs)');
ylabel('Fluorescence (arbitrary units)');
hold off

ha = fill([18000 18000 90000 90000],[-5 35 35 -5],[0.8 0.8 0.8]);
hold on
hb = fill([162000 162000 234000 234000],[-5 35 35 -5],[0.8 0.8 0.8]);
hc = fill([306000 306000 378000 378000],[-5 35 35 -5],[0.8 0.8 0.8]);
plot(a5(idx05(1):end,3),'black');
axis([0 404000 -5 35]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
hold off


valsx = a5(idx05(1):end,3);
trimx = ((26/403988)*(1:403988)+3)';
ha = fill([18000 18000 90000 90000],[-10 20 20 -10],[0.8 0.8 0.8]);
hold on
hb = fill([162000 162000 234000 234000],[-10 20 20 -10],[0.8 0.8 0.8]);
hc = fill([306000 306000 378000 378000],[-10 20 20 -10],[0.8 0.8 0.8]);
plot((valsx-trimx),'black');
axis([0 404000 -10 20]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
hold off