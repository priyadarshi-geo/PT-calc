%% Temperature determination using Amph-Grt and Amph-Plag geothermometer
% Last Update: 20.05.2016

clc; clear; close all; format compact; format short;
fprintf('_________________ AMPHIBOLE GEOTHERMOMETRY __________________\n\n');

%% Different Geothermometric expressions
    % 1 -- Graham and Powell 1984, JMG
    % 2 -- Powell 1985, JMG
    % 3 -- Perchuk et al. 1985, JMG
    % 4 -- Ravna 2000, Lithos
    
%% Reading dataset from file
fprintf (' Calculates temprature using GRT_AMPH and PLAG-AMPH assemblage\n\n');
fprintf (' Please form a ''.dat'' file for the dataset. Enter Fe, Mg, Ca, Mn for\n');
fprintf (' Garnet(12-O) then Fe, Mg total apfu followed by Fe, Mg in M1-M3 sites\n');
fprintf (' (site C acc. to Leake et al. 1997) for amphibole(23-O)and then Ca, Na\n');
fprintf (' for Plagioclase(8-O). Enter each dataset in a new line.\n\n');
fprintf (' [CAUTION: Enter the filenames without extensions]\n\n');      

filename1 = input(' Enter the file-name[e.g. Data_rim]: ','s');
filename = [filename1,'.dat'];
fid   =  fopen(filename,'r');
if fid == -1                                                                           
    fprintf('     >>ERROR: Problem in opening the file!\n');
else
    data  = fscanf(fid,'%f',[8 inf]);
    data  = data';
    [r,c] = size(data) ;
    dataset = r;
end
fclose(fid);

% storing Fe,Mg,Ca,Mn,Na values for Grt, Amph & Plag in vectors
Fe_G = data(:,1); Mg_G = data(:,2); Ca_G = data(:,3); Mn_G = data(:,4);
Fe_A = data(:,5); Mg_A = data(:,6); Fe_A2 = data(:,7); Mg_A2 = data(:,8); 
%Ca_Pl = data(:,9);Na_Pl = data(:,10);

%Pressure = input(' Enter Pressure(kbar): ');

%% Geothermometer formulations
%% Ravna(2000,Lithos)

T_Rav = zeros(dataset,1);
for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XMn_Grt = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_A2(i)/Mg_A2(i)); ln_Kd = log(Kd);
    
    T_Rav(i) = ((1504 + 1784*(XCa_Grt+XMn_Grt))/(ln_Kd + 0.720))-273;
end
% rounding off to nearest integer 
T_Rav = round(T_Rav);

%% Perchuk et al.(1985,JMG); Powell(1985,JMG); Graham and Powell(1984,JMG) 

T_Per = zeros(dataset,1);
T_Po = zeros(dataset,1);
T_GP = zeros(dataset,1);

for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_A(i)/Mg_A(i)); ln_Kd = log(Kd);
    
    T_Per(i) = (3330/(ln_Kd + 2.333))-273;
    T_Po(i) = ((2580 + 3340*XCa_Grt)/(ln_Kd + 2.20))-273;
    T_GP(i) = ((2880 + 3280*XCa_Grt)/(ln_Kd + 2.426))-273;
end
% rounding off to nearest integer 
T_Per = round(T_Per); T_Po = round(T_Po); T_GP = round(T_GP);

%% Printing certain results on screen

SlNo = [1:dataset]';

% Forming matrix of results per geothermometer per data
results = [SlNo,round(T_Rav),round(T_Per),round(T_Po),round(T_GP)];

% Forming matrix of max, min and average values per geothermometer
Rav_max = max(T_Rav); Rav_min = min(T_Rav); Rav_avg = mean(T_Rav);
Per_max = max(T_Per); Per_min = min(T_Per); Per_avg = mean(T_Per);
Po_max = max(T_Po); Po_min = min(T_Po); Po_avg = mean(T_Po);
GP_max = max(T_GP); GP_min = min(T_GP); GP_avg = mean(T_GP);


res_mat = zeros(4,3);
res_mat(1,:) = [Rav_max,Rav_min,Rav_avg];
res_mat(2,:) = [Per_max,Per_min,Per_avg];
res_mat(3,:) = [Po_max,Po_min,Po_avg];
res_mat(4,:) = [GP_max,GP_min,GP_avg];


mean_T = mean([Rav_avg,Per_avg,Po_avg,GP_avg]);

% printing on screen
fprintf ('\nSummary of Results:\n\n');

names = ['S.No.',' Rav-00',' Per-85',' Pow-85',' G&P-84'];

% printing the individual values and avg_value for each data
fprintf('%s \n',names);
for kk=1:dataset
    for ll=1:5
        if ll==1
            fprintf('%4.0d  ',results(kk,ll));
        else
            fprintf('%5.0d  ',results(kk,ll));
        end    
    end
    fprintf('\n');
end

% printing the max, min and avg value for all data per geothermometer
fprintf('\n');
for kk=1:3
    if kk==1
        fprintf(' MAX');
    elseif kk==2
        fprintf(' MIN');
    else
        fprintf(' AVG');
    end
    for ll=1:4
        fprintf('%7.0d',round(res_mat(ll,kk)));
    end    
    fprintf('\n');
end

fprintf('\nTotal dataset used for computation: %d \n', dataset);
fprintf ('Average temperature of all geothermometers: %0.2f °C\n',mean_T);
fprintf('------------------------------------------------------------------------\n');

%% Printing results to an output file
format short
f_name=[filename1,'_res_','.txt'];
fidww = fopen(f_name,'wt');
fprintf(fidww,'%s \n',names);
fclose(fidww);

fidww = fopen(f_name,'at');
for kk=1:dataset
    for ll=1:5
        if ll==1
            fprintf(fidww,'%4.0d  ',results(kk,ll));
        else
            fprintf(fidww,'%7.0d  ',results(kk,ll));
        end    
    end
    fprintf(fidww,'\n');
end

fprintf(fidww,'\n');
for kk=1:3
    if kk==1
        fprintf(fidww,' Max');
    elseif kk==2
        fprintf(fidww,' Min');
    else
        fprintf(fidww,' Avg');
    end

    for ll=1:4
        fprintf(fidww,'%8.0d',round(res_mat(ll,kk)));
    end
    fprintf(fidww,'\n');
end
fprintf(fidww,'\nTotal dataset used for computation: %d \n', dataset);
fprintf(fidww,'_____________________________________________________________________________________\n');
fclose(fidww);

%% Plotting the result
plot(1:dataset,T_Rav,'.-',1:dataset,T_Per,'s-',1:dataset,T_Po,'*-',...
    1:dataset,T_GP,'d-');
axis tight;
legend('Rav-00','Per-85','Pow-85','G&P-84');
xlabel('Dataset'); ylabel('Temperature [°C]');
title(['Temperature']);
box on; grid on;



    
    
    
    
    
    
    
    