%% Temperature determination using Grt-Cpx geothermometer
% Last Update: 21.06.2016

clc; clear; close all; format compact; format short;
fprintf('_________________GARNET-ORTHOPYROXENE GEOTHERMOMETRY__________________\n\n');

%% Different Geothermometric expressions
    % 1 -- Lee & Ganguly 1988, JP
    % 2 -- Harley 1984, CMP
    % 3 -- Bhattacharya & Sen 1991, JP
    % 4 -- Ganguly 1996, -- only plotting results
    
%% Reading dataset from file
fprintf (' Please form ''.dat'' files for the dataset. Enter Fe, Mg, Ca, Mn\n');
fprintf (' for Garnet(12-O) first & then for Orthopyroxene(6-O) in a single\n');
fprintf (' line. Enter each dataset of Grt-Opx in a new line. \n\n');
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
% storing Fe,Mg,Ca,Mn values for Grt and Cpx in vectors
Fe_G = data(:,1); Mg_G = data(:,2); Ca_G = data(:,3); Mn_G = data(:,4);
Fe_O = data(:,5); Mg_O = data(:,6); Ca_O = data(:,7); Mn_O = data(:,8);

Pressure = input(' Enter Pressure(kbar): ');
R_cal = 1.9872036; % in cal/K/mol

%% Geothermometer formulations
%% Bhattacharya et al. (1991,JP) formulation
T_Bh = zeros(dataset,1); P_Bh = Pressure*1e3; % in bars
for i=1:dataset
    XMg_Gt = Mg_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XFe_Gt = Fe_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XCa_Gt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XMn_Gt = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));

    XMg_Opx = Mg_O(i)/(Fe_O(i)+ Mg_O(i));
    
    A_term = -1220*XFe_Gt*XMg_Gt - 441*XCa_Gt*(XMg_Gt-XFe_Gt)...
                -136*XMg_Gt*XMg_Gt + 746*XFe_Gt*XFe_Gt;
    
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_O(i)/Mg_O(i)); ln_Kd = log(Kd);      
    T_Bh(i) = ((1611 + 0.021*P_Bh + 906*XCa_Gt + A_term + 477*(2*XMg_Opx-1))/...
            (ln_Kd+0.796))-273;
end

%% Lee & Ganguly(1988,JP) formulation

T_LG = zeros(dataset,1); P_LG = Pressure;
for i=1:dataset
    
    XCa_Gt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XMn_Gt = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    
    WCa_Gt = 3000; WMn_Gt = 3000; %in cal/mol
    non_ideal_term = (WCa_Gt*XCa_Gt+WMn_Gt*XMn_Gt)/R_cal;
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_O(i)/Mg_O(i)); ln_Kd = log(Kd);
       
    %T_LG_i(i) = ((1971 + 11.91*P_LG)/(ln_Kd+0.96))-273;
    T_LG(i) = ((1971 + 11.91*P_LG + non_ideal_term)/(ln_Kd+0.96))-273;
     
end
% rounding off to nearest integer 
T_LG = round(T_LG);

%% Harley(1984,CMP) formulation

T_H = zeros(dataset,1); P_H = Pressure; % in kbar
for i=1:dataset
    XCa_Gt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_O(i)/Mg_O(i)); ln_Kd = log(Kd);
    
    T_H(i) = ((3740 + 1400*XCa_Gt + 22.86*P_H)/(R_cal*ln_Kd+1.96))-273;
end
% rounding off to nearest integer 
T_H = round(T_H);

%% Results of Ganguly et al. 1996
ask = input(' Want to enter results of Ganguly et al. 96 ? [y-yes/n-No] : ','s');
if strcmpi (ask,'y') == 1
    filename2 = input(' Enter the file-name having results of Ganguly-96: ','s');
    filename2 = [filename2,'.dat'];
    fid   =  fopen(filename2,'r');
    if fid == -1
        fprintf('     >>ERROR: Problem in opening the file!\n');
    else
        data_Ganguly  = fscanf(fid,'%f',[10 inf]);
        data_Ganguly  = data_Ganguly';
    end
    T_Gang = data_Ganguly(:,2);
else
    T_Gang = zeros(dataset,1);
end

%% Printing certain results on screen

SlNo = [1:dataset]';

% Forming matrix of results per geothermometer per data
results = [SlNo,round(T_Bh),round(T_LG),round(T_H),round(T_Gang)];

% Forming matrix of max, min and average values per geothermometer
Bh_max = max(T_Bh); Bh_min = min(T_Bh); Bh_avg = mean(T_Bh);
LG_max = max(T_LG); LG_min = min(T_LG); LG_avg = mean(T_LG);
H_max = max(T_H); H_min = min(T_H); H_avg = mean(T_H);
Gan_max = max(T_Gang); Gan_min = min(T_Gang); Gan_avg = mean(T_Gang);

res_mat = zeros(4,3);
res_mat(1,:) = [Bh_max,Bh_min,Bh_avg];
res_mat(2,:) = [LG_max,LG_min,LG_avg];
res_mat(3,:) = [H_max,H_min,H_avg];
res_mat(4,:) = [Gan_max,Gan_min,Gan_avg];

mean_T = mean([Bh_avg,LG_avg,H_avg,Gan_avg]);

% printing on screen
fprintf ('\nSummary of Results:\n\n');

names = ['S.No.',' B&S-91',' L&G-88',' Har-84',' Gan-96'];

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
f_name=[filename1,'_res_',num2str(Pressure),'.txt'];
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
plot(1:dataset,T_Bh,'x-',1:dataset,T_LG,'p-',1:dataset,T_H,'.-',...
    1:dataset,T_Gang,'o-');
axis tight;
legend('B&S-91','L&G-88','Har-84','Gan-96');
xlabel('Dataset'); ylabel('Temperature [°C]');
title(['Temperature [',num2str(Pressure),' kbar]']);
box on; grid on;



    
    
    
    
    
    
    
    