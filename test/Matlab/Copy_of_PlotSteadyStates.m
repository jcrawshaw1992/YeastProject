% script to calculate steady state ring of chain of springs with Pressure
% force P radially outwards.

clear
close all

N_D = [10,20 40 80];
N_Z = [6, 13, 26, 52, 104, 208];



for i = 1:3 %length(N_D)
    for j = 1:3
        %base_directory  = '/Users/ozzy/Desktop/Roley/mobernabeu/testoutput/CylinderValidation/';
        base_directory  = '/Users/ozzy/Desktop/CylinderValidation/'
        %base_directory = '/tmp/ozzy/testoutput/CylinderValidation/';
        directory  = [ base_directory num2str(N_D(i)) '_' num2str(N_Z(i+j-1))]
        file = [directory '/results_from_time_0/results.viznodes'];
        results = load(file);
        time(3*(i-1)+j,:) = results(:,1);
        
        for k=1:length(time)
            centre_x = mean(results(k,2:3:end));
            centre_y = mean(results(k,3:3:end));
            radius(3*(i-1)+j,k) = mean(sqrt((results(k,2:3:end)-centre_x).^2 + (results(k,3:3:end)-centre_y).^2));
        end
    end
end


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(time(1,:),radius(1,:),'k--',time(2,:),radius(2,:),'k',time(3,:),radius(3,:),'k-.','linewidth',2.0)
hold on;
plot(time(4,:),radius(4,:),'b--',time(5,:),radius(5,:),'b',time(6,:),radius(6,:),'b-.','linewidth',2.0)
plot(time(7,:),radius(7,:),'r--',time(8,:),radius(8,:),'r',time(9,:),radius(9,:),'r-.','linewidth',2.0)
%plot(time(10,:),radius(10,:),'g--',time(11,:),radius(11,:),'g',time(12,:),radius(12,:),'g:')

R_long = 0.001640934122444; 
R_norm = 0.001868409976191;
R_short = 0.002728920489522;
time_ends = [time(7,1),time(7,end)];
plot(time_ends, [R_long,R_long],'g-.','linewidth',2.0);
plot(time_ends, [R_norm,R_norm],'g','linewidth',2.0);
plot(time_ends, [R_short,R_short],'g--','linewidth',2.0);

axis([time_ends,1.0e-3,3.0e-3])



legend('10x6','10x13','10x26','20x13','20x26', '20x52', '40x26', '40x52', '40x104','Location', 'EastOutside')
print( '-dpng', 'RegularCylinderRadius');


% Now load iregular meshes
sampling_band_width = 0.5e-3;
sampling_band_centre = 0.0e-3;

results = load('/Users/ozzy/Desktop/CylinderValidation/581/results_from_time_0/results.viznodes');
time(13,:) = results(:,1);
for k=1:length(time)
    x = results(k,2:3:end);
    y = results(k,3:3:end);
    z = results(k,4:3:end);
    I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
    x(I) = [];
    y(I) = [];
    z(I) = [];
    radius(13,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
end
figure;
plot3(x-mean(x),y-mean(y),z,'r.');
hold on


results = load('/Users/ozzy/Desktop/CylinderValidation/2103/results_from_time_0/results.viznodes');
time(14,:) = results(:,1);
for k=1:length(time)
    x = results(k,2:3:end);
    y = results(k,3:3:end);
    z = results(k,4:3:end);
    I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
    x(I) = [];
    y(I) = [];
    z(I) = [];
    radius(14,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
end
plot3(x-mean(x),y-mean(y),z,'b.');

results = load('/Users/ozzy/Desktop/CylinderValidation/8334/results_from_time_0/results.viznodes');
time(15,:) = results(:,1);
for k=1:length(time)
    x = results(k,2:3:end);
    y = results(k,3:3:end);
    z = results(k,4:3:end);
    I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
    x(I) = [];
    y(I) = [];
    z(I) = [];
    radius(15,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
end
plot3(x-mean(x),y-mean(y),z,'g.');
%view(2);

results = load('/Users/ozzy/Desktop/CylinderValidation/33331/results_from_time_0/results.viznodes');
time(16,:) = results(:,1);
for k=1:length(time)
    x = results(k,2:3:end);
    y = results(k,3:3:end);
    z = results(k,4:3:end);
    I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
    x(I) = [];
    y(I) = [];
    z(I) = [];
    radius(16,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
end
plot3(x-mean(x),y-mean(y),z,'y.');
%view(2);

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
axes('Parent',fig,'FontSize',24);
plot(time(2,:),radius(2,:),'k','linewidth',2.0)
hold on;
plot(time(5,:),radius(5,:),'b','linewidth',2.0)
plot(time(8,:),radius(8,:),'r','linewidth',2.0)
%plot(time(10,:),radius(10,:),'g--',time(11,:),radius(11,:),'g',time(12,:),radius(12,:),'g:')

plot(time(13,:),radius(13,:),'k--','linewidth',2.0);
plot(time(14,:),radius(14,:),'b--','linewidth',2.0);
plot(time(15,:),radius(15,:),'r--','linewidth',2.0);
plot(time(16,:),radius(16,:),'y--','linewidth',2.0);

axis([time_ends,1.4e-3,2.2e-3])
legend('10x13','20x26', '40x52', '581','2103','8334','33331','Location', 'EastOutside')
print( '-dpng', 'IrregularCylinderRadius');