%Aleksandra Denisin (adenisin@stanford.edu), Pruitt Laboratory
%2/22/2016
%Code written for analysis of force-displacement curves obtained by AFM
%microindentation

function [E_kPa_H_2um,Rsq_H_2um, nameoutput,SlopeCanti_V_nm]=E_Hertz
close all
clear all
clear global xdata ydata F R
global R k xdata ydata F
%First we use a force-displacement curve on glass to calculate the Optical
%Lever Sensitivity. Make sure to save multiple trials of the glass
%force-displacment curve in a single Excel file. The first column is height
%of the stage, the second column is voltage output of the photodiode
%(cantilever deflection). This repeats with odd columns being stage
%position and even columns containing cantilever deflection.

filename = input('Type the filename containing the force-displacement curve on glass (for OLS). Use single quotations around filename:');
%user will type in something like 'StifferGelCalibrationGlass_10.xls'
[numeric_calibration, text, raw] = xlsread(filename);
R = input ('Type in the radius of the spherical indentor used (in meters):');
k = input('Type in the spring constant of your cantilever (in N/m):'); 
%Need to figure out the slope of the calibration curve to change voltage to
%distance We will measure the first inflection point of the calibration
%graph approaching the sample and take our measurement from there to the
%max
CantiK_nm = k/1E9; %N/nm change units to match displacement
Siz = size(numeric_calibration);
measurements = Siz(1,2)/2;

InflectionIndex = zeros(measurements,1);
MaxIndex = zeros(measurements,1);
MaxIndex_x = zeros(measurements,1);
Slope = zeros(measurements,1);
SlopeCanti = zeros(measurements,1);

%the data off the park system is in um on x axis. let's change to nm since
%the rest of the code is in nm...
%for i = 2:2:2*measurements;
%   numeric_calibration(:,i-1) = numeric_calibration(:,i-1)*1000;%conversion to nm
%end

for i = 2:2:2*measurements;
    [maxVal, idx] = max(numeric_calibration(:,i));
    x_idx = numeric_calibration(idx,i-1);
    MaxIndex(i/2,1) = maxVal;
    MaxIndex_x(i/2,1) = idx;
    %we analyse just the first points until the max is reached then find
    %the inflection point on the approach
    dy = zeros(size(numeric_calibration(idx-5:idx,i)));
    dy(1:end-1) = diff(numeric_calibration(idx-5:idx,i));
    dy(end) = dy(end-1);
    
    dy_x = numeric_calibration(idx-5:idx,i-1);
    %     figure1 = plot(numeric_calibration(:,1) ,numeric_calibration(:,2));
    %     hold on;
    %     plot(dy_x,dy);
    
    avg_dy = abs(mean(dy));
    for z = 20: length(dy);
        if (dy(z,1)>avg_dy);
            inflection = numeric_calibration(z,1);
            inflection_x = dy_x(z,1);
            InflectionIndex (i/2,1)= (min(find( numeric_calibration(:,i-1)== inflection_x)))+10;
            break
        end
    end
    
    p = polyfit(numeric_calibration(idx-5:idx,i-1),numeric_calibration(idx-5:idx,i),1);  % p returns 2 coefficients fitting r = a_1 * x + a_2
    Slope (i/2,1) = p(1); %V/nm
end
   
SlopeAvg = (abs(mean(Slope)));
SlopeCanti_V_nm= SlopeAvg; %V/nm
SlopeCanti_V_N = SlopeAvg/CantiK_nm; %V/nm*(nm/N) = V/N
filename = input('Type the filename containing the force-distance curve for your sample. Use single quotation marks around the filename:');
%user will type in something like 'StifferGel_10.xls' remember to include
%the '
[numeric_Data, text, raw] = xlsread(filename);

Siz = size(numeric_Data(100:end, :)); %throw out first 100 points to get rid of noise

%the data off the park system is in um on x axis. let's change to nm since
%the rest of the code is in nm...
%for i = 2:2:Siz(1,2);
%    numeric_Data(:,i-1) = numeric_Data(:,i-1)*1000;%conversion to nm
%end

%We will use Siz to determine how may datapoints we must go through to
%evaluate all of the AFM indentation trials
nameoutput = filename(1:end-6);
for i = 2:2:Siz(1,2);
    %to use F vs delta formulas later on, we need to flip our
    %data over the Y axis this means that f(x) changes to f(-x)
    nm_vs_Disp(:,i-1) = -numeric_Data(100:end,i-1); %cut the first 100 points to get rid of noise
    nm_vs_Disp(:,i) = numeric_Data(100:end,i)./(SlopeCanti_V_nm); %volt*(nm/V) = nm
end

figure(1)
xlabel('Sample Height (nm)');
ylabel('Cantilever Deflection (nm)');
title(sprintf(nameoutput));
grid on;
figureHandle = gcf;
%# make all text in the figure to size 14
set(findall(figureHandle,'type','text'),'fontSize',14)
hold on;
for j = 2:2:Siz(1,2);
    TurnAround = find(nm_vs_Disp(:,j)==(max(nm_vs_Disp(:,j))));
    figure1 = plot(nm_vs_Disp(1:TurnAround,j-1),nm_vs_Disp(1:TurnAround,j),nm_vs_Disp(TurnAround+1:end,j-1),nm_vs_Disp(TurnAround+1: end,j), 'Linewidth', 2);
end
N_vs_Disp = zeros(Siz(1,1),Siz(1,2));

%convert data to N vs m
for i = 2:2:Siz(1,2);
    N_vs_Disp(:, i-1) = nm_vs_Disp(:,i-1).*1E-9; %now 'sample height' is in meters!
    N_vs_Disp (:, i) = (nm_vs_Disp(:,i).*k)./1E9;%now cantilever deflection is changed to force experienced by cantilever (N)
end

figure(2)
xlabel('Sample Height (m)');
ylabel('Force (N)');
grid on;
figureHandle = gcf;
%# make all text in the figure to size 14
set(findall(figureHandle,'type','text'),'fontSize',14)
hold on;

for l = 2:2:Siz(1,2);
    figure2 = plot(N_vs_Disp(1:TurnAround,l-1), N_vs_Disp(1:TurnAround,l), N_vs_Disp(TurnAround+1:end,l-1), N_vs_Disp(TurnAround+1: end,l), 'Linewidth', 2);
    hold on
end

%make empty matrices to store the values needed for calculating
%Hertz and Oliver Pharr models
E_H = zeros(Siz(1,2)/2,1);
E_kPa_red_H = zeros(Siz(1,2)/2,1);
Rsq_H = zeros(Siz(1,2)/2,1);
Inf_index = zeros(Siz(1,2)/2,1);
EndPt = zeros(Siz(1,2)/2,1);
StartPt = zeros(Siz(1,2)/2,1);
d2 = zeros(Siz(1,2)/2,1);
z2 = zeros(Siz(1,2)/2,1);
d1 = zeros(Siz(1,2)/2,1);
z1 = zeros(Siz(1,2)/2,1);
cp = zeros(Siz(1,2)/2,1);
firstpt=zeros(Siz(1,2)/2,1);
Start = zeros(Siz(1,2)/2,1);
End = zeros(Siz(1,2)/2,1);
E_H_2um = zeros(Siz(1,2)/2,1);
E_kPa_H_2um = zeros(Siz(1,2)/2,1);
Rsq_H_2um = zeros(Siz(1,2)/2,1);

%Hertz fitting:
%we follow the method by the MacKay + Kumar to flatten out our graph
%0] rescale data such that z position (x) starts with 0 and increases. The
%deflection (y) is unchanged. Make sure data is nm vs nm
for m = 2:2:Siz(1,2);
    nm_vs_Disp_approach(:, m) = nm_vs_Disp(1:TurnAround, m);
    if nm_vs_Disp(1, m-1) <0
        nm_vs_Disp_approach(:, m-1) = nm_vs_Disp(1:TurnAround, m-1)-nm_vs_Disp_approach(1, m-1);
    end
    if nm_vs_Disp(1, m-1) <0
        nm_vs_Disp_approach(:, m-1) = nm_vs_Disp(1:TurnAround, m-1)+abs(nm_vs_Disp(1, m-1));
    end
    if nm_vs_Disp(1, m-1) ==0
        nm_vs_Disp_approach(:, m-1) = nm_vs_Disp(1:TurnAround, m-1);
    end
end

%plot z position on x axis and vertical deflection on y axis
figure(4)
xlabel('Sample Height z position(nm)');
ylabel('Deflection (nm)');
grid on;
hold on;
figureHandle = gcf;
%# make all text in the figure to size 20
set(findall(figureHandle,'type','text'),'fontSize',20)
for m = 2:2:Siz(1,2);
    figure4=plot(nm_vs_Disp_approach(:, m-1), nm_vs_Disp_approach(:, m));
    hold on
end

%1] define baseline of the curve by choosing two points on flat part at least
%2 um apart
BaselineEnd = 1000*input('Look at the graph and determine a point before the inflection. Enter in rounded microns:'); 

%from now on, do everything for 1 deflection curve at a time
for z = 2:2:Siz(1,2)
    %calculate slope betwen these points. multiply slope of z position data and subtract these values form vertical
    %deflection to straighten curve so baseline is perfectly horizontal
    val = BaselineEnd; %value to find
    tmp = abs(nm_vs_Disp_approach(:,z-1)-val);
    [idx idx] = min(tmp); %index of closest value
    Inf_index (z/2, 1)= idx;
    closest = nm_vs_Disp_approach(idx,z-1); %closest value
    EndPt(z/2, 1)=  Inf_index (z/2, 1);
    val = closest-500; %0.5 um away from indentation point identified by user
    tmp = abs(nm_vs_Disp_approach(:,z-1)-val);
    [idx idx] = min(tmp); %index of closest value
    StartPt (z/2, 1)= idx;
    %closest = nm_vs_Disp_approach(idx,z-1); %closest value
    BaselinePoints = [StartPt, EndPt];
    
    %let's  plot which points we are baselining along with the fit
    figure(5)
    xlabel('Sample Height z position(nm)');
    ylabel('Deflection (nm)');
    grid on;
    figureHandle = gcf;
    %# make all text in the figure to size 14
    set(findall(figureHandle,'type','text'),'fontSize',14)
    hold on;
    figure5 = plot(nm_vs_Disp_approach(BaselinePoints(z/2,1):BaselinePoints(z/2,2), z-1),nm_vs_Disp_approach(BaselinePoints(z/2,1):BaselinePoints(z/2,2), z));
    p_baseline = polyfit(nm_vs_Disp_approach(BaselinePoints(z/2,1):BaselinePoints(z/2,2),z-1),nm_vs_Disp_approach(BaselinePoints(z/2,1):BaselinePoints(z/2,2),z),1);
    % p returns 2 coefficients fitting r = a_1 * x + a_2
    S_baseline(z/2,1) = p_baseline(1); %N/nm
    fit_x = nm_vs_Disp_approach(BaselinePoints(z/2,1):BaselinePoints(z/2,2),z-1);
    fit_y = fit_x.*p_baseline(1,1)+p_baseline(1,2);
    figure5 = plot(fit_x, fit_y, 'r');
    
    %use StartPt as the 1st point in these series, get rid of others
    nm_vs_Disp_flat{:,z-1} = nm_vs_Disp_approach(StartPt(z/2,1):end, z-1); %no changes to the x axis which is still z position or displacement
    nm_vs_Disp_flat{:,z} = nm_vs_Disp_approach(StartPt(z/2,1):end, z)-(nm_vs_Disp_approach(StartPt(z/2,1):end, z-1).*S_baseline(z/2,1)); % we multiply the x axis by the slope and subtract it from the original y
    
    figure(6)
    xlabel('Sample Height z position(nm)');
    ylabel('Deflection (nm) (flat & baselined)');
    grid on;
    figureHandle = gcf;
    %# make all text in the figure to size 20
    set(findall(figureHandle,'type','text'),'fontSize',20)
    hold on;
    figure6 = plot(nm_vs_Disp_flat{:, z-1},nm_vs_Disp_flat{:, z});
    
    %calculate average deflection along baseline and subtract it from
    %vertical deflection data to normalize so that deflection starts from zero
    averageDef = mean(nm_vs_Disp_flat{:,z}(1:EndPt(z/2, 1)-StartPt(z/2, 1)));
    nm_vs_Disp_flat_avg{:, z-1} = nm_vs_Disp_flat{:, z-1};
    nm_vs_Disp_flat_avg {:, z} =  nm_vs_Disp_flat{:, z}-averageDef;
    figure6 = plot(nm_vs_Disp_flat_avg {:, z-1},nm_vs_Disp_flat_avg {:, z},'r');
    
    %2] calculate the contact point (z position of the probe at which it first
    %indents the material). MacKay & Kumar suggest choosing vertical
    %deflections 5–15 nm for the first point and 100 nm for the second point
    %(or the maximum deflection)
    %cp = (z2 - d1)-(z1-d1)(d2/d1)^n / 1-(d2/d1)^n 
    % where cp is the z position at the contact point, d 1 and z 1 refer
    % to the vertical deflection and z position of the first data point
    % to be fit, d 2 and z 2 refer to the vertical deflection and z position
    % of the second data point to be fit,
    
    %we will use d2 and z2 as max deflection
    n = 2/3; %for spherical tip
    val = 1; % 1 nm vertical deflection as suggested by MacKay and Kumar
    tmp = abs(nm_vs_Disp_flat_avg{:,z}-val);
    [idx idx] = min(tmp); %index of closest value
    firstpt (z/2, 1)= idx;
    closest = nm_vs_Disp_flat_avg{:,z}(idx); %closest value
    d1(z/2, 1)=  nm_vs_Disp_flat_avg{:, z}(firstpt(z/2, 1));
    z1(z/2, 1)= nm_vs_Disp_flat_avg{:,z-1} (firstpt(z/2, 1));
    
    if nm_vs_Disp_flat_avg{:,z}(end)>100
        val = 100;
        tmp = abs(nm_vs_Disp_flat_avg{:,z}-val);
        [idx idx] = min(tmp); %index of closest value
        secondpt (z/2, 1)= idx;
        closest = nm_vs_Disp_flat_avg{:,z}(idx); %closest value
        d2(z/2, 1)=  nm_vs_Disp_flat_avg{:, z}(secondpt(z/2, 1));
        z2(z/2, 1)= nm_vs_Disp_flat_avg{:,z-1} (secondpt(z/2, 1));
    end
    if nm_vs_Disp_flat_avg{:,z}(end)<100
        d2(z/2, 1)= nm_vs_Disp_flat_avg{:, z}(end);
        z2(z/2, 1)= nm_vs_Disp_flat_avg{:, z-1}(end);
    end
    
    %now find d1 and z1 from where 15 nm deflection was detected
    d1(z/2, 1)=  nm_vs_Disp_flat_avg{:, z}(firstpt(z/2, 1));
    z1(z/2, 1)= nm_vs_Disp_flat_avg{:,z-1} (firstpt(z/2, 1));
    cp(z/2, 1) = ((z2(z/2, 1)-d2(z/2, 1))-(z1(z/2, 1)-d1(z/2, 1))*(d2(z/2, 1)/d1(z/2, 1))^n)/(1-((d2(z/2, 1)/d1(z/2, 1))^n));
    
   %check contact point by overlaying on a graph
    figure(7)
    xlabel('Sample Height z position(nm)');
    ylabel('Deflection (nm) (flat & baselined)');
    grid on;
    figureHandle = gcf;
    %# make all text in the figure to size 20
    set(findall(figureHandle,'type','text'),'fontSize',20)
    hold on;
    figure7 = plot(cp,zeros(Siz(1,2)/2,1), 'ro');
    figure7 = plot(nm_vs_Disp_flat_avg{:,z-1},nm_vs_Disp_flat_avg{:, z});
    
    %For each data point after contact point, calculate indentation force (F =
    %(k*d) and inddentation depth (delta= z-d). change delta to meters and F to
    %N
    cp_index = zeros(Siz(1,2)/2,1);
    val = cp(z/2, 1);
    tmp = abs(nm_vs_Disp_flat_avg{:, z-1}-val);
    [idx idx] = min(tmp); %index of closest value
    cp_index (z/2, 1)= idx;
    closest = nm_vs_Disp_flat_avg{:,z-1}(idx); %closest value
    
    N_vs_delta_cp{:,z-1} = (nm_vs_Disp_flat_avg{:, z-1}(cp_index(z/2, 1):end)-nm_vs_Disp_flat_avg{:, z}(cp_index(z/2, 1):end))./1E9; % this is delta (z-d) in m
    N_vs_delta_cp{:,z} = (nm_vs_Disp_flat_avg{:, z}(cp_index(z/2, 1):end)).*k/1E9; % this is in N, N/m*nm*m/1E9nm = N
    
    figure(8)
    xlabel('delta (z-d) in m');
    ylabel('Force (N)');
    grid on;
    figureHandle = gcf;
    %# make all text in the figure to size 14
    set(findall(figureHandle,'type','text'),'fontSize',14)
    hold on;
    figure8 = plot(N_vs_delta_cp{:,z-1},N_vs_delta_cp{:,z});
    
    %need to normalize such that delta begins at 0
    N_vs_delta_cp_normalized{:,z-1} = N_vs_delta_cp{:,z-1}-N_vs_delta_cp{:,z-1}(1);
    N_vs_delta_cp_normalized{:,z}= N_vs_delta_cp{:,z};
    
    %3] fit Hertz model to normalized data
    %we have small indentations (d<R/10). So we can use:
    %F = (4/3)*(Er)*delta^(3/2)*sqrt(R), from Radmacher 2007
    %Fit 2 um of indentation
    val = N_vs_delta_cp_normalized{:,z-1}(1)+1E-6; %end at 1 um - 2um indentation
    tmp = abs(N_vs_delta_cp_normalized{:,z-1}-val);
    [idx idx] = min(tmp); %index of closest value
    End (z/2, 1)= idx;
    closest = N_vs_delta_cp_normalized{:,z-1}(idx); %closest value
    
    xdata = N_vs_delta_cp_normalized{:,z-1}(1:End(z/2,1),:);
    ydata = N_vs_delta_cp_normalized{:,z}(1:End(z/2,1),:);
    
    E_H_2um(z/2, 1) = fminbnd(@funHertz, 0, 1E6, optimset('TolX',1e-8));
    %Print the final fitting parameter
    error = ydata-F;
    sse = sum(error.^2);
    ssr = sum((F-mean(ydata)).^2);
    sst = sum((ydata-mean(ydata)).^2);
    Rsq_H_2um(z/2, 1) =1-(sse/sst); 
    
    %Plotting the fitted result with the measured data
    figure(9)
    xlabel('delta (z-d) in m');
    ylabel('Force (N)');
    grid on;
    figureHandle = gcf;
    %# make all text in the figure to size 20
    set(findall(figureHandle,'type','text'),'fontSize',20)
    hold on
    %figure6 = plot(N_vs_delta_cp{:,z-1}, N_vs_delta_cp{:,z}, 'g');
    figure9=plot(xdata,ydata, xdata,F,'--r', 'LineWidth', 3);
    E_kPa_red_H_2um (z/2, 1) = E_H_2um (z/2, 1)/1000; %kPa
end

%we actually found the Ereduced
% so we need to multiply by 1 - v^2 to find E
%E = Er*(1-v^2)
v = 0.48; %for polyacrylamide accourding to Boudou 2006
E_kPa_H_2um = E_kPa_red_H_2um.*(1-v^2);
end


