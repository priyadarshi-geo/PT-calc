%% Pressure determination using Grt-Cpx-Plag-Qtz geobarometer
% Last update: 01.11.2017 - GCPS_P5_v2
% Previous code: 22.06.2016 - GCPS_P5_v1
% slightly changed in the coding algorithm
% Plagioclase activity is restricted to Newton 1983 model
clc; clear; close all; format compact; format short;

fprintf('_________________GRT-CPX-PlAG-QTZ GEOBAROMETRY__________________\n\n');

%% Different Geobarometric expressions
    % 1 -- Newton & Perkins 1982
    % 2 -- Eckert et al. 1991
    % 3 -- Eckert et al. 1991 & Bhattacharya et al. 1991
    % 4 -- Eckert et al. 1991 & Moecher et al. 1988
    % 5 -- log10(K) of Moecher et al. 1988 for CMAS sys (Di-reaction)
 
%% Reading dataset from file
fprintf (' Please form ''.dat'' files for the dataset. Enter Fe,Mg,Ca,Mn\n');
fprintf (' for Garnet(12-O) followed by Ca,Na for Plagioclase(8-O) and \n');
fprintf (' then XCa and XMg for Clinopyroxene(6-O) in single line.\n');
fprintf (' # Enter each dataset in a new line. \n\n');

fprintf (' [CAUTION: FOR CPX- Entries are XCa(M2-site) and XMg(M1-site)]\n');
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
Fe_Gt = data(:,1); Mg_Gt = data(:,2); Ca_Gt = data(:,3); Mn_Gt = data(:,4);
Ca_Pl = data(:,5);Na_Pl = data(:,6);
XCa_Cpx = data(:,7); XMg_Cpx = data(:,8);

Temp = input(' Enter Temperature [°C]: ');

%% Geobarometer formulations

P_NP       = zeros(dataset,1);
P_Eck      = zeros(dataset,1);
P_Eck_B    = zeros(dataset,1);
P_Eck_M    = zeros(dataset,1);
Di_K_M     = zeros(dataset,1);
T_Kel      = Temp+273;
R          = 8.314; % in J/mol/K

for i=1:dataset
    
    %% activities of GARNET end-members; after Newton & Perkins 1982 
    XFe_Grt = Fe_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i));
    XMg_Grt = Mg_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i));
    XCa_Grt = Ca_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i));
    
    W_CaMg  = 13807.2 - 6.276*T_Kel;
    RTln_gammaCa = W_CaMg*(XMg_Grt*XMg_Grt+XMg_Grt*XFe_Grt); 
    RTln_gammaMg = W_CaMg*(XCa_Grt*XCa_Grt+XCa_Grt*XFe_Grt);
    
    a_Grs = XCa_Grt*exp(RTln_gammaCa/(R*T_Kel));
    a_Prp = XMg_Grt*exp(RTln_gammaMg/(R*T_Kel));
    
    %% activitiy of ANORTHITE
    XAn = Ca_Pl(i)/(Ca_Pl(i)+ Na_Pl(i));
    XAb = Na_Pl(i)/(Ca_Pl(i)+ Na_Pl(i));
    
    % after Newton & Perkins, 1982
    % ln_gammaAn = (XAb*XAb)*((1019.076+(4751.663*XAn))/T_Kel);
    % after Newton 1983 AJS
    ln_gammaAn = (XAb*XAb)*((1032+(4727*XAn))/T_Kel);
    a_An = ((XAn*(1+XAn)*(1+XAn))/4)* exp(ln_gammaAn);
   
    %% activitiy of DIOPSIDE; after Wood & Banno 1973
    a_Di = XCa_Cpx(i) * XMg_Cpx(i);
     
    % Kd calculation for An + Di = 2/3 Grs + 1/3 Prp
    Kd = (a_Grs*a_Grs*a_Prp)/(a_An*a_Di); ln_Kd = log(Kd);
    
    %% Newton & Perkins (1982,Am Min) formulation
    P_NP(i) = 0.675 + 0.017179*T_Kel + 0.0035962*T_Kel*ln_Kd;
    
    %% Eckert et al. (1991,Am Min) formulation
    P_Eck(i) = 2.60 + 0.01718*T_Kel + 0.003596*T_Kel*ln_Kd;
    
    %% Bhattarcharya et al.1991,JP Garnet activity model 
    %  using Eckert et al. formulation
   
    % activities of GARNET end-members; after Bhattarcharya et al.1991
    XFe_Grt = Fe_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    XMg_Grt = Mg_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    XCa_Grt = Ca_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    
    ln_gamma_Grs = 0.906*XMg_Grt*XMg_Grt + XFe_Grt*XMg_Grt*(0.465-0.610*...
        (XFe_Grt-XMg_Grt));
    
    ln_gamma_Prp = 0.906*XCa_Grt*XCa_Grt + XFe_Grt*(0.746-1.22*XMg_Grt)+...
        XCa_Grt*XFe_Grt*(1.346-0.610*XMg_Grt);
    
    a_Grs = XCa_Grt*exp(ln_gamma_Grs);
    a_Prp = XMg_Grt*exp(ln_gamma_Prp);
    
    % Kd calculation for An + Di = 2/3 Grs + 1/3 Prp
    Kd = (a_Grs*a_Grs*a_Prp)/(a_An*a_Di); ln_Kd = log(Kd);
    
    % Eckert et al. 1991 formulation
    P_Eck_B(i) = 2.60 + 0.01718*T_Kel + 0.003596*T_Kel*ln_Kd;
    
    %% Moecher et al.1988,CMP log(10) K isopleths
    %  using Eckert et al. formulation
    %  Garnet activity model after Anovitz and Essene 1987
    
    % activities of GARNET end-members; after Anovitz & Essene 1987
    XFe_Grt = Fe_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    XMg_Grt = Mg_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    XCa_Grt = Ca_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    XMn_Grt = Mn_Gt(i)/(Fe_Gt(i)+ Mg_Gt(i)+ Ca_Gt(i) + Mn_Gt(i));
    
    RTln_gCa = (XMg_Grt*XMg_Grt)*(4047-1.5*T_Kel-6094*XCa_Grt)...
              +(XFe_Grt*XFe_Grt)*(150-1.5*T_Kel+7866*XCa_Grt)...
              +(XMg_Grt*XFe_Grt)*(3290-3*T_Kel+886*XCa_Grt+2300*(XMg_Grt-XFe_Grt))...
              + 4640 *(1-2*XCa_Grt)...
              +(XMn_Grt*XFe_Grt)*(2117-1.5*T_Kel+3933*XCa_Grt-1967*(1-2*XCa_Grt))...
              +(XMg_Grt*XMn_Grt)*(2524-1.5*T_Kel-3047*XCa_Grt+1524*(1-2*XCa_Grt))...
              + 2300*(XMg_Grt*XFe_Grt*XMn_Grt);
          
    RTln_gMg = (XCa_Grt*XCa_Grt)*(1000-1.5*T_Kel+6094*XMg_Grt)...
              +(XFe_Grt*XFe_Grt)*(2500-4600*XMg_Grt)+3000*(XMn_Grt*XMn_Grt)...
              +(XCa_Grt*XFe_Grt)*(1757+747*XMg_Grt-3933*(XCa_Grt-XFe_Grt))...
              + 4640 *(1-2*XMg_Grt)...
              +(XMn_Grt*XFe_Grt)*(4350-2300*XMg_Grt+1150*(1-2*XMg_Grt))...
              +(XMn_Grt*XCa_Grt)*(5524+3047*XMg_Grt-1524*(1-2*XMg_Grt))...
              + 3933*(XCa_Grt*XFe_Grt*XMn_Grt);

    RTln_gFe = (XMg_Grt*XMg_Grt)*(200+4600*XFe_Grt)...
              +(XCa_Grt*XCa_Grt)*(4083-1.5*T_Kel-7866*XFe_Grt)...
              +(XCa_Grt*XMg_Grt)*(943-1633*XFe_Grt-3047*(XMg_Grt-XCa_Grt))...
              - 4640 *(1-2*XFe_Grt)...
              +(XMn_Grt*XCa_Grt)*(2117-1.5*T_Kel-3933*XFe_Grt+1917*(1-2*XFe_Grt))...
              +(XMn_Grt*XMg_Grt)*(-1650+2300*XFe_Grt-1150*(1-2*XFe_Grt))...
              + 3048*(XCa_Grt*XMg_Grt*XMn_Grt);

    a_Grs   = XCa_Grt*exp(RTln_gCa/(R*T_Kel));
    a_Prp   = XMg_Grt*exp(RTln_gMg/(R*T_Kel));
    a_Alm   = XFe_Grt*exp(RTln_gFe/(R*T_Kel));  
          
    % Kd calculation for An + Di = 2/3 Grs + 1/3 Prp as per Eckert et al. 
    % 1991 formulation
    Kd    = (a_Grs*a_Grs*a_Prp)/(a_An*a_Di); ln_Kd = log(Kd);  
    P_Eck_M(i) = 2.60 + 0.01718*T_Kel + 0.003596*T_Kel*ln_Kd;
    
    % log10 (K) of Moecher et al. 1988
    a_Grs = (a_Grs)^3;  a_Prp = (a_Prp)^3;  a_Alm = (a_Alm)^3;  
    Kd_Di_Moe = (a_Grs*a_Grs*a_Prp)/(a_An*a_An*a_An*a_Di*a_Di*a_Di);
    Di_K_M(i) = log10(Kd_Di_Moe);
  
end

%% Printing certain results on screen

SlNo = [1:dataset]';
% Forming matrix of results per geobarometer per data
results = [SlNo,P_NP,P_Eck,P_Eck_B,P_Eck_M,Di_K_M];

% Forming matrix of max, min and average values per geobarometer
NP_max = max(P_NP); NP_min = min(P_NP); NP_avg = mean(P_NP);
Eck_max = max(P_Eck); Eck_min = min(P_Eck); Eck_avg = mean(P_Eck);
E_B_max = max(P_Eck_B); E_B_min = min(P_Eck_B); E_B_avg = mean(P_Eck_B);
E_M_max = max(P_Eck_M); E_M_min = min(P_Eck_M); E_M_avg = mean(P_Eck_M);

res_mat = zeros(4,3);
res_mat(1,:) = [NP_max,NP_min,NP_avg];
res_mat(2,:) = [Eck_max,Eck_min,Eck_avg];
res_mat(3,:) = [E_B_max,E_B_min,E_B_avg];
res_mat(4,:) = [E_M_max,E_M_min,E_M_avg];

mean_P = mean([NP_avg,Eck_avg,E_B_avg,E_M_avg]);

% printing on screen
fprintf ('\nSummary of Results:\n\n');
names = [' S.No.',' NP-82',' Eck-91',' Eck_Bh',' Eck_Mo',' Di_K_M'];

% printing the individual values for each data
fprintf('%s \n',names);
for kk=1:dataset
    for ll=1:6
        if ll == 1
            fprintf('%4.0d  ',results(kk,ll));
        else
            fprintf('%5.1f  ',results(kk,ll));
        end
    end
    fprintf('\n');
end

% printing the max, min and avg value for all data per geobarometer
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
        fprintf('%7.1f',(res_mat(ll,kk)));
    end
    fprintf('\n');
end

fprintf('\nTotal dataset used for computation: %d \n', dataset);
fprintf ('Average Pressure of all geobarometers: %0.2f kbars\n',mean_P);

%% Printing results to an output file
format short
f_name=[filename1,'_res_',num2str(Temp),'.txt'];
names = [' S.No.',' NP-82',' Eck-91',' Eck_Bh',' Eck_Mo',' Di_K_M'];
fidww = fopen(f_name,'wt');
fprintf(fidww,'%s \n',names);
fclose(fidww);

fidww = fopen(f_name,'at');
for kk=1:dataset
    for ll=1:6
        if ll == 1
            fprintf(fidww,'%4.0d   ',results(kk,ll));
        else
            fprintf(fidww,'%6.1f   ',results(kk,ll));
        end
    end
    fprintf(fidww,'\n');
end

fprintf(fidww,'\n');
for kk=1:3
    if kk==1
        fprintf(fidww,' MAX');
    elseif kk==2
        fprintf(fidww,' MIN');
    else
        fprintf(fidww,' AVG');
    end
    for ll=1:4  
        fprintf(fidww,'%7.1f   ',res_mat(ll,kk)); 
    end
    fprintf(fidww,'\n');
end
fprintf(fidww,'\nTotal dataset used for computation: %d \n', dataset);

fprintf(fidww,'_____________________________________________________________________________________\n');
fclose(fidww);

%% Plotting the result
plot(1:dataset,P_NP,'s-',1:dataset,P_Eck,'.-',1:dataset,P_Eck_B,'d-',...
    1:dataset,P_Eck_M,'x-');
axis tight;
legend('NP-82','Eck-91','Eck_B_h','Eck_M_o');
xlabel('Dataset'); ylabel('P [kbar]');
title(['Pressure [',num2str(Temp),' °C]']);
box on; grid on;

fprintf('------------------------------------------------------------------------\n');
fprintf ('Eck_Bh: Grt_activity-Bhattacharya et al. 1991\n');
fprintf ('Eck_Mo: Grt_activity-Moecher et al. 1988\n');
fprintf ('Plag_acitivty: Newton 1983 in all calculations \n');
fprintf ('Cpx_activity: Wood and Banno 1973 in all calculations \n');
fprintf ('Di_K_M: log10(K) of 2Grs+Prp=3An+3Di from Moecher et al. 1988 \n');
fprintf ('  \n\n');
    
    
    
    
    
    
    
    