% script to calculate steady state ring of chain of springs with Pressure
% force P radially outwards.

% clear
 close all

sampling_band_width = 0.5; % 0.5e-3
sampling_band_centre = 0.0;


N_D = [20, 40, 80, 160];
N_Z = [15, 30, 60, 120, 240, 480];


for i = 1:4 %length(N_D)
    for j = 1:3
        %base_directory  = '/Users/ozzy/Desktop/Roley/mobernabeu/testoutput/CylinderValidation/';
        base_directory  = '/Users/ozzy/Desktop/CylinderValidation/Equilateral/'
        %base_directory = '/tmp/ozzy/testoutput/CylinderValidation/';
        directory  = [ base_directory num2str(N_D(i)) '_' num2str(N_Z(i+j-1))]
        file = [directory '/results_from_time_0/results.viznodes'];
        results = load(file);
        time(3*(i-1)+j,:) = results(:,1);

        for k=1:length(time)
            x = results(k,2:3:end)*1e3;
            y = results(k,3:3:end)*1e3;
            z = results(k,4:3:end)*1e3;
            I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
            x(I) = [];
            y(I) = [];
            z(I) = [];
            radius(3*(i-1)+j,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
        end
    end
end


R_squashed = 0.001679411112578*1e3;
R_equil = 0.001975257505321*1e3;
R_stretched = 0.003099915191342*1e3;
time_ends = [time(1,1),time(1,end)];



fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04

hold on;
plot(time(1,:),radius(1,:),'k')
plot(time(4,:),radius(4,:),'k')
plot(time(7,:),radius(7,:),'k')
plot(time(10,:),radius(10,:),'k')
plot(time_ends, [R_stretched,R_stretched],'k:','linewidth',2.0);


axis([time_ends,1.0,3.5])

% legend('n_y = 20, n_x = 15',...
%        'n_y = 40, n_x = 30',...
%        'n_y = 80, n_x = 60',...
%        'n_y = 160, n_x = 120',...
%        'Analytical',...
%        'Location', 'SouthEast')
       
xlabel('Time (hours)')
ylabel('Radius (mm)')
saveaspngandeps(-1,'StretchedRegularCylinderRadius',7, 7/5, 9);


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04

hold on;
plot(time(2,:),radius(2,:),'k')
plot(time(5,:),radius(5,:),'k')
plot(time(8,:),radius(8,:),'k')
plot(time(11,:),radius(11,:),'k')
plot(time_ends, [R_equil,R_equil],'k:','linewidth',2.0);


axis([time_ends,1.4,2.2])

% legend('n_y = 20 n_x = 30',...
%        'n_y = 40 n_x = 60',...
%        'n_y = 80 n_x = 120',...
%        'n_y = 160 n_x = 240',...
%        'Analytical',...
%        'Location', 'NorthEast')

xlabel('Time (hours)')
ylabel('Radius (mm)')
saveaspngandeps(-1,'EquilateralRegularCylinderRadius',7, 7/5, 9);


fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
hold on;
plot(time(3,:),radius(3,:),'k')
plot(time(6,:),radius(6,:),'k')
plot(time(9,:),radius(9,:),'k')
plot(time(12,:),radius(12,:),'k')
plot(time_ends, [R_squashed,R_squashed],'k:','linewidth',2.0);


axis([time_ends,1.4,1.8])

% legend('n_y = 20 n_x = 60',...
%        'n_y = 40 n_x = 120',...
%        'n_y = 80 n_x = 240',...
%        'n_y = 160 n_x = 480',...
%        'Analytical',...
%        'Location', 'NorthEast')
   
xlabel('Time (hours)')
ylabel('Radius (mm)')
saveaspngandeps(-1,'SquashedRegularCylinderRadius',7, 7/5, 9);



% Now load iregular meshes

% results = load('/Users/ozzy/Desktop/CylinderValidation/581/results_from_time_0/results.viznodes');
% time(13,:) = results(:,1);
% for k=1:length(time)
%     x = results(k,2:3:end)*1e3;
%     y = results(k,3:3:end)*1e3;
%     z = results(k,4:3:end)*1e3;
%     I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
%     x(I) = [];
%     y(I) = [];
%     z(I) = [];
%     radius(13,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
% end
% figure;
% plot3(x,y,z,'k.');
% hold on
% 
% 
% results = load('/Users/ozzy/Desktop/CylinderValidation/2103/results_from_time_0/results.viznodes');
% time(14,:) = results(:,1);
% for k=1:length(time)
%     x = results(k,2:3:end)*1e3;
%     y = results(k,3:3:end)*1e3;
%     z = results(k,4:3:end)*1e3;
%     I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
%     x(I) = [];
%     y(I) = [];
%     z(I) = [];
%     radius(14,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
% end
% plot3(x,y,z,'b.');
% 
% 
% results = load('/Users/ozzy/Desktop/CylinderValidation/8338/results_from_time_0/results.viznodes');
% time(15,:) = results(:,1);
% for k=1:length(time)
%     x = results(k,2:3:end)*1e3;
%     y = results(k,3:3:end)*1e3;
%     z = results(k,4:3:end)*1e3;
%     I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
%     x(I) = [];
%     y(I) = [];
%     z(I) = [];
%     radius(15,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
% end
% plot3(x,y,z,'r.');
% 
% results = load('/Users/ozzy/Desktop/CylinderValidation/33331/results_from_time_0/results.viznodes');
% time(16,:) = results(:,1);
% for k=1:length(time)
%     x = results(k,2:3:end)*1e3;
%     y = results(k,3:3:end)*1e3;
%     z = results(k,4:3:end)*1e3;
%     I = find(or(z>sampling_band_centre+sampling_band_width,z<sampling_band_centre-sampling_band_width));
%     x(I) = [];
%     y(I) = [];
%     z(I) = [];
%     radius(16,k) = mean(sqrt((x-mean(x)).^2 + (y-mean(y)).^2));
% end
% plot3(x,y,z,'cy.');
% view(2);

fig = figure('Renderer', 'painters'); % So it plots in ubuntu 8.04
%axes('Parent',fig,'FontSize',24);
plot(time(2,:),radius(2,:),'k','linewidth',2.0)
hold on;
plot(time(5,:),radius(5,:),'k','linewidth',2.0)
plot(time(8,:),radius(8,:),'k','linewidth',2.0)
plot(time(11,:),radius(11,:),'k','linewidth',2.0)


plot(time(13,:),radius(13,:),'k--','linewidth',2.0);
plot(time(14,:),radius(14,:),'k--','linewidth',2.0);
plot(time(15,:),radius(15,:),'k--','linewidth',2.0);
plot(time(16,:),radius(16,:),'k--','linewidth',2.0);

plot(time_ends, [R_equil,R_equil],'k:','linewidth',2.0);




axis([time_ends,1.4,2.2])
%legend('20x30','40x60', '80x120', '160x240','581','2103','8334','33331', 'Location', 'EastOutside')
%print( '-dpng', 'IrregularCylinderRadius');
xlabel('Time (hours)')
ylabel('Radius (mm)')
saveaspngandeps(-1,'IrregularCylinderRadius',7, 7/5, 9);