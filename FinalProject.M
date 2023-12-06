clear;
% import brace data sets
[step_stats_b1, angles_b1, moments_b1, grf_b1, startStance_b1, power_b1] = ProcessData("brace1.mat");
[step_stats_b2, angles_b2, moments_b2, grf_b2, startStance_b2, power_b2] = ProcessData("brace2.mat");
[step_stats_b3, angles_b3, moments_b3, grf_b3, startStance_b3, power_b3] = ProcessData("brace3.mat");

% import braceless data sets
[step_stats_w1, angles_w1, moments_w1, grf_w1, startStance_w1, power_w1] = ProcessData("walking1.mat");
[step_stats_w2, angles_w2, moments_w2, grf_w2, startStance_w2, power_w2] = ProcessData("walking2.mat");
[step_stats_w3, angles_w3, moments_w3, grf_w3, startStance_w3, power_w3] = ProcessData("walking3.mat");

% establish stance phases for all sets from visual inspection
bStnc = [[145,240];
        [187,282];
        [186,285]];
wStnc = [[255,341];
        [249,338];
        [281,372]];

figure
titles = {'Foot Angle','Shank Angle','Thigh Angle'};
for i = 1:3
    subplot(3,1,i)
    hold on
    % plot brace data sets in blue
    plot(linspace(0,100,length(angles_b1(bStnc(1,1):bStnc(1,2),i))),rad2deg(angles_b1(bStnc(1,1):bStnc(1,2),i)),"Color","#0a4a7a")
    plot(linspace(0,100,length(angles_b2(bStnc(2,1):bStnc(2,2),i))),rad2deg(angles_b2(bStnc(2,1):bStnc(2,2),i)),"Color","#3787c4")
    %plot(linspace(0,100,length(angles_b3(bStnc(3,1):bStnc(3,2),i))),rad2deg(angles_b3(bStnc(3,1):bStnc(3,2),i)),"Color","#8ccafa")
    % plot braceless data sets in red
    plot(linspace(0,100,length(angles_w1(wStnc(1,1):wStnc(1,2),i))),rad2deg(angles_w1(wStnc(1,1):wStnc(1,2),i)),"Color","#852908")
    %plot(linspace(0,100,length(angles_w2(wStnc(2,1):wStnc(2,2),i))),rad2deg(angles_w2(wStnc(2,1):wStnc(2,2),i)),"Color","#c7552c")
    plot(linspace(0,100,length(angles_w3(wStnc(3,1):wStnc(3,2),i))),rad2deg(angles_w3(wStnc(3,1):wStnc(3,2),i)),"Color","#faa98c")

    title(titles(i),'FontWeight','normal')
    if i == 2
        ylabel('Angle (Degrees)','FontWeight','bold')
    end
    hold off
    if i == 1
        legend('Brace1','Brace2','Walking1','Walking3','Orientation','horizontal','Location','northoutside')
    end
end
sgtitle('Global Member Angles During the Stance Phase','FontWeight','bold')
xlabel('Percentage of Stance Phase','FontWeight','bold')

figure
titles = {'Hip Power','Knee Power','Ankle Power'};
for i = 1:3
    subplot(3,1,i)
    hold on
    % plot brace data sets in blue
    plot(linspace(0,100,length(power_b1(bStnc(1,1):bStnc(1,2),i))),rad2deg(power_b1(bStnc(1,1):bStnc(1,2),i)),"Color","#0a4a7a")
    plot(linspace(0,100,length(power_b2(bStnc(2,1):bStnc(2,2),i))),rad2deg(power_b2(bStnc(2,1):bStnc(2,2),i)),"Color","#3787c4")
    %plot(linspace(0,100,length(power_b3(bStnc(3,1):bStnc(3,2),i))),rad2deg(power_b3(bStnc(3,1):bStnc(3,2),i)),"Color","#8ccafa")
    % plot braceless data sets in red
    plot(linspace(0,100,length(power_w1(wStnc(1,1):wStnc(1,2),i))),rad2deg(power_w1(wStnc(1,1):wStnc(1,2),i)),"Color","#852908")
    %plot(linspace(0,100,length(power_w2(wStnc(2,1):wStnc(2,2),i))),rad2deg(power_w2(wStnc(2,1):wStnc(2,2),i)),"Color","#c7552c")
    plot(linspace(0,100,length(power_w3(wStnc(3,1):wStnc(3,2),i))),rad2deg(power_w3(wStnc(3,1):wStnc(3,2),i)),"Color","#faa98c")

    title(titles(i),'FontWeight','normal')
    if i == 2
        ylabel('Power (watts)','FontWeight','bold')
    end
    hold off
    if i == 1
        legend('Brace1','Brace2','Walking1','Walking3','Orientation','horizontal','Location','northoutside')
    end
end
sgtitle('Joint Power During the Stance Phase','FontWeight','bold')
xlabel('Percentage of Stance Phase','FontWeight','bold')
