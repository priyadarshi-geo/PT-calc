%% Temperature determination using Grt-Cpx geothermometer
% Last Update: 20.06.2016

clc; clear; close all; format compact; format short;
fprintf('_________________GARNET-CLINOPYROXENE GEOTHERMOMETRY__________________\n\n');

%% Different Geothermometric expressions
    % 1 -- Nakamura 2009, JMG
    % 2 -- Ravna 2000, JMG
    % 3 -- Ai 1994, CMP
    % 4 -- Pattinson & Newton 1989,CMP
    % 5 -- Krogh 1988, CMP
    % 6 -- Powell 1985, JMG
    % 7 -- Ganguly 1979, GCA
    % 8 -- Ellis & Green 1979, CMP
    % 9 -- Ganguly 1996, CMP -- only plotting results
    
%% Reading dataset from file
fprintf (' Please form ''.dat'' files for the dataset. Enter Fe, Mg, Ca, Mn\n');
fprintf (' for Garnet(12-O) first & then for Clinopyroxene(6-O) in a single\n');
fprintf (' line. Enter each dataset of Grt-Cpx in a new line. \n\n');
fprintf (' [CAUTION: Enter the filenames without extensions]\n\n');      

filename1 = input(' Enter the file-name[e.g. Data_rim]: ','s');
filename = [filename1,'.dat'];
fid   =  fopen(filename,'r');
if fid == -1                                                                           
    fprintf('     >>ERROR: Problem in opening the file!\n');
else
    data  = fscanf(fid,'%f',[10 inf]);
    data  = data';
    [r,c] = size(data) ;
    dataset = r;
end
fclose(fid);
% storing Fe,Mg,Ca,Mn values for Grt and Cpx in vectors
Fe_G = data(:,1); Mg_G = data(:,2); Ca_G = data(:,3); Mn_G = data(:,4);
Fe_C = data(:,5); Mg_C = data(:,6); Ca_C = data(:,7); Mn_C = data(:,8);
% storing Al & Fe3+ values of Cpx for Nakamura 2009 in vectors
Al_C = data(:,9); Fe3_C = data(:,10); 

Pressure = input(' Enter Pressure(kbar): ');

%% Geothermometer formulations
%% Nakamura(2009,JMG) formulation

if max(Fe3_C)== 0.000 %|| min(Fe3_C)==0.000
    T_Nak = zeros(dataset,1);
    fprintf ('\n >> Fe3-Cpx values are invalid/not given.\n\n');
else
    T_Nak = zeros(dataset,1);
    P_Nak = Pressure;
    for i=1:dataset
        Xalm = Fe_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
        Xprp = Mg_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
        Xgrs = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
        Xsps = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
        
        XMg_Cpx = Mg_C(i)/(Fe_C(i)+ Mg_C(i)+ Al_C(i)+ Fe3_C(i));
        XFe_Cpx = Fe_C(i)/(Fe_C(i)+ Mg_C(i)+ Al_C(i)+ Fe3_C(i));
        
        Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
        
        A = 0.5*Xgrs*(Xprp-Xalm-Xsps); B = 0.5*Xgrs*(Xprp-Xalm+Xsps);
        C = 0.5*(Xgrs+Xsps)*(Xprp-Xalm);
        
        num = 2784+14.52*P_Nak+(2601+1.44*P_Nak)*(2*Xgrs*Xprp-A)+...
            (1183+6.98*P_Nak)*(Xgrs*Xgrs-A)-105*(2*Xgrs*Xalm+B)+...
            (814.6+3.61*P_Nak)*(Xgrs*Xgrs+B)-(254.6+8.42*P_Nak)*...
            (2*Xprp*Xalm-Xalm*Xalm+C)-83.6*(Xprp*Xprp-2*Xprp*Xalm+C)+...
            1388*Xsps-462*(XMg_Cpx-XFe_Cpx);
        
        den = ln_Kd+1.431+0.695*(2*Xgrs*Xprp+Xgrs*Xgrs-2*A)+...
            0.203*(Xgrs*Xgrs-2*Xgrs*Xalm)+0.922*Xsps;
        
        T_Nak(i) = (num/den)-273;
    end
    % rounding off to nearest integer
    T_Nak = round(T_Nak);
end

%% Ravna(2000,JMG) formulation

% for verification use file: Ravna.dat; (Results are almost identical)
% [verification is done with Temperature obtained in Table 4 & Table 5 in
% Ravna 2000(JMG,18,211-219);Grt-Cpx apfu are given in respective tables]
% Results can aldo verify Ai,Krogh,Pow and EG formulaiotns using same data

T_Ravna = zeros(dataset,1);
P_Ravna = Pressure*0.1;
for i=1:dataset
    XMg_Grt = Mg_G(i)/(Fe_G(i)+ Mg_G(i));
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    XMn_Grt = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_Ravna(i) = ((1939.9+(3270*XCa_Grt)-(1396*XCa_Grt*XCa_Grt)+...
                  (3319*XMn_Grt)-(3535*XMn_Grt*XMn_Grt)+(1105*XMg_Grt)-...
                  (3561*XMg_Grt*XMg_Grt)+(2324*XMg_Grt*XMg_Grt*XMg_Grt)+...
                  (169.4*P_Ravna))/(ln_Kd+1.223))-273;
end
% rounding off to nearest integer 
T_Ravna = round(T_Ravna);
                    
%% Ai(1994,CMP) formulation

% for verification use file: Ai.dat (Results does not match exactly)
% [verification is done with Temperature for Sample XM46 & XM48 from Ai et 
% al.1994(Table 3-Diamondiferrous lherzolite Shee et al.1982,CMP,81,79-87)]

T_Ai = zeros(dataset,1);
P_Ai = Pressure;
for i=1:dataset
    XMg_Grt = 100*Mg_G(i)/(Fe_G(i)+ Mg_G(i));
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
   
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_Ai(i) = (((-1629*XCa_Grt*XCa_Grt)+(3648.55*XCa_Grt)-(6.59*XMg_Grt)...
               +1987.98+(17.66*P_Ai))/(ln_Kd+1.076))-273;
end
% rounding off to nearest integer 
T_Ai = round(T_Ai);

%% Pattinson & Newton(1989,CMP) formulation

% for verification use file: P&N.dat (Results does not match exactly)
% [verification is done with Temperature for Sample XM46 & XM48 from Ai et 
% al.1994(Table 3-Diamondiferrous lherzolite Shee et al.1982,CMP,81,79-87)]

T_PN = zeros(dataset,1);
P_PN = Pressure;
for i=1:dataset
    X = Mg_G(i)/(Fe_G(i)+ Mg_G(i));
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i));
    
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    if X>0.124999 && X<0.6000001
        if XCa_Grt <= 0.20 
            Xca1 = 0; Xca2 = 0; T2 = 0;
            a0=3.606; b0=-5.172; c0=2.317; d0=0.1742;
            a1=26370; b1=-32460; c1=11050; d1=1012;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
    
        elseif XCa_Grt >0.20 && XCa_Grt< 0.3
            Xca1 = 0.2;
            a0=3.606; b0=-5.172; c0=2.317; d0=0.1742;
            a1=26370; b1=-32460; c1=11050; d1=1012;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
            Xca2 = 0.3;
            a0=15.87; b0=-20.30; c0=7.468; d0=-0.1479;
            a1=43210; b1=-53230; c1=18120; d1=776;
            T2=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
        elseif XCa_Grt ==0.3 
            Xca1 = 0; Xca2 = 0; T2 = 0;
            a0=15.87; b0=-20.30; c0=7.468; d0=-0.1479;
            a1=43210; b1=-53230; c1=18120; d1=776;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
   
        elseif XCa_Grt >0.3 && XCa_Grt< 0.4
            Xca1 = 0.3;
            a0=15.87; b0=-20.30; c0=7.468; d0=-0.1479;
            a1=43210; b1=-53230; c1=18120; d1=776;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
            Xca2 = 0.4;
            a0=14.64; b0=-18.72; c0=6.940; d0=-0.2583;
            a1=44900; b1=-55250; c1=18820; d1=712;
            T2=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
        elseif XCa_Grt ==0.4 
            Xca1 = 0; Xca2 = 0; T2 = 0;
            a0=14.64; b0=-18.72; c0=6.940; d0=-0.2583;
            a1=44900; b1=-55250; c1=18820; d1=712;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
    
        elseif XCa_Grt >0.4 && XCa_Grt< 0.5
            Xca1 = 0.4;
            a0=14.64; b0=-18.72; c0=6.940; d0=-0.2583;
            a1=44900; b1=-55250; c1=18820; d1=712;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
            Xca2 = 0.5;
            a0=9.408; b0=-12.37; c0=4.775; d0=-0.2331;
            a1=38840; b1=-47880; c1=16300; d1=859;
            T2=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;
            
        elseif XCa >= 0.5
            Xca1 = 0; Xca2 = 0; T2 = 0;
            a0=9.408; b0=-12.37; c0=4.775; d0=-0.2331;
            a1=38840; b1=-47880; c1=16300; d1=859;
            T1=(((a1*X*X*X+b1*X*X+c1*X+d1)/(ln_Kd+a0*X*X*X+b0*X*X+c0*X+d0))+...
                5.5*(P_PN-15))-273;    
        end
        
        if T2 == 0
            T_PN(i) = T1;
        else
            T_PN(i) = T1+(T2-T1)*(XCa_Grt-Xca1)/(Xca2-Xca1);
        end
    else
        fprintf(' ERR: Mg#-Garnet falls out of experimental range\n');
    end

end
% rounding off to nearest integer 
T_PN = round(T_PN);

%% Krogh(1988,CMP) formulation

% for verification use file: Krogh.dat; Result: T = 675 °C at 18.2 kbar
% [verification is done with Temperature for Sample 1338-Ii from Krogh et 
% al.1990(JMG,8,289-309; Column 1 of Table 8); Grt & Cpx compositions for 
% the sample are given in Table-1 & Table-2 respectively]

T_Krogh = zeros(dataset,1);
P_Krogh = Pressure;
for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i)); 
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_Krogh(i) = (((-6173*XCa_Grt*XCa_Grt)+(6731*XCa_Grt)+1879+10*P_Krogh)/...
        (ln_Kd+1.393))-273;
end
% rounding off to nearest integer 
T_Krogh = round(T_Krogh);

%% Powell(1985,JMG) formulation

T_Po = zeros(dataset,1);
P_Po = Pressure;
for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    XFe_Grt = Fe_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    XMg_Grt = Mg_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    
    XFe_Cpx = Fe_C(i)/(Fe_C(i)+ Mg_C(i));
    XMg_Cpx = Mg_C(i)/(Fe_C(i)+ Mg_C(i));
    
    Kd = (XFe_Grt/XMg_Grt)/(XFe_Cpx/XMg_Cpx); ln_Kd = log(Kd);
    %Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_Po(i) = ((2790+10*P_Po+3140*XCa_Grt)/(ln_Kd+1.735))-273;
end
% rounding off to nearest integer 
T_Po = round(T_Po);

%% Ganguly(1979,GCA) formulation

T_G = zeros(dataset,1);
P_G = Pressure;
for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i)); 
    XMn_Grt = Mn_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)+ Mn_G(i)); 
    
    ln_psi = ((1586*XCa_Grt))+((1308*XMn_Grt));
    Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_G(i) = ((4801+(11.07*Pressure)+ln_psi)/(ln_Kd+2.930))-273;
end
% rounding off to nearest integer 
T_G = round(T_G);

%% Ellis & Green(1979,CMP) formulation

% for verification use file: EG.dat (Result differ very very silghtly)
% [verification is done with Temperature for Sample 9-5-2C,264-3,163-K,
% 248-6 from Ellis & Green 1976(Table 4); Grt & Cpx compositions are from 
% Evans et al. 1978(EPSL, first 2 data) and 1979(Am. Min, last 2 data)]

T_EG = zeros(dataset,1);
P_EG = Pressure;
for i=1:dataset
    XCa_Grt = Ca_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i)); 
    XFe_Grt = Fe_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    XMg_Grt = Mg_G(i)/(Fe_G(i)+ Mg_G(i)+ Ca_G(i));
    
    XFe_Cpx = Fe_C(i)/(Fe_C(i)+ Mg_C(i));
    XMg_Cpx = Mg_C(i)/(Fe_C(i)+ Mg_C(i));
    
    Kd = (XFe_Grt/XMg_Grt)/(XFe_Cpx/XMg_Cpx); ln_Kd = log(Kd);
    %Kd = (Fe_G(i)/Mg_G(i))/(Fe_C(i)/Mg_C(i)); ln_Kd = log(Kd);
    
    T_EG(i) = ((3104*XCa_Grt+3030+10.86*P_EG)/(ln_Kd+1.9034))-273;
end
% rounding off to nearest integer 
T_EG = round(T_EG);

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

%% Determining lnKd=A+B/T format for Mg-Fe exchange;
%Fe_Mg_Gt = Fe_G./Mg_G; Fe_Mg_Cpx = Fe_C./Mg_C;
%lnKd1 = log(Fe_Mg_Gt./Fe_Mg_Cpx); lnKd2 = log(1./(Fe_Mg_Gt./Fe_Mg_Cpx));

%% Printing certain results on screen

SlNo = [1:dataset]';

% Forming matrix of results per geothermometer per data
results = [SlNo,round(T_Nak),round(T_Ravna),round(T_Ai),round(T_PN),round(T_Krogh),...
    round(T_Po),round(T_G),round(T_EG),round(T_Gang)];

% Forming matrix of max, min and average values per geothermometer
Nak_max = max(T_Nak); Nak_min = min(T_Nak); Nak_avg = mean(T_Nak);
Ravna_max = max(T_Ravna); Ravna_min = min(T_Ravna); Ravna_avg = mean(T_Ravna);
Ai_max = max(T_Ai); Ai_min = min(T_Ai); Ai_avg = mean(T_Ai);
PN_max = max(T_PN); PN_min = min(T_PN); PN_avg = mean(T_PN);
Krogh_max = max(T_Krogh); Krogh_min = min(T_Krogh); Krogh_avg = mean(T_Krogh);
Po_max = max(T_Po); Po_min = min(T_Po); Po_avg = mean(T_Po);
G_max = max(T_G); G_min = min(T_G); G_avg = mean(T_G);
EG_max = max(T_EG); EG_min = min(T_EG); EG_avg = mean(T_EG);
Gan_max = max(T_Gang); Gan_min = min(T_Gang); Gan_avg = mean(T_Gang);

res_mat = zeros(9,3);
res_mat(1,:) = [Nak_max,Nak_min,Nak_avg];
res_mat(2,:) = [Ravna_max,Ravna_min,Ravna_avg];
res_mat(3,:) = [Ai_max,Ai_min,Ai_avg];
res_mat(4,:) = [PN_max,PN_min,PN_avg];
res_mat(5,:) = [Krogh_max,Krogh_min,Krogh_avg];
res_mat(6,:) = [Po_max,Po_min,Po_avg];
res_mat(7,:) = [G_max,G_min,G_avg];
res_mat(8,:) = [EG_max,EG_min,EG_avg];
res_mat(9,:) = [Gan_max,Gan_min,Gan_avg];

mean_T = mean([Nak_avg,Ravna_avg,Ai_avg,PN_avg,Krogh_avg,Po_avg,G_avg,EG_avg,Gan_avg]);

% printing on screen
fprintf ('\nSummary of Results:\n\n');

names = ['S.No.',' Nak-09',' Rav-00','  Ai-94','  PN-89','  Kro-88',' Pow-85',...
    ' Gan-79',' EG-79',' Gan-96'];

% printing the individual values and avg_value for each data
fprintf('%s \n',names);
for kk=1:dataset
    for ll=1:10
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
    for ll=1:9
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
    for ll=1:10
        if ll==1
            fprintf(fidww,'%4.0d  ',results(kk,ll));
        else
            fprintf(fidww,'%7.0d  ',results(kk,ll));
        end    
    end
    fprintf(fidww,'\n');
end

% fprintf(fidww,'\n');
% for kk=1:3
%     if kk==1
%         fprintf(fidww,' Max');
%     elseif kk==2
%         fprintf(fidww,' Min');
%     else
%         fprintf(fidww,' Avg');
%     end
% 
%     for ll=1:9
%         fprintf(fidww,'%8.0d',round(res_mat(ll,kk)));
%     end
%     fprintf(fidww,'\n');
% end
fprintf(fidww,'\nTotal dataset used for computation: %d \n', dataset);
fprintf(fidww,'_____________________________________________________________________________________\n');
fclose(fidww);

%% Plotting the result
plot(1:dataset,T_Nak,1:dataset,T_Ravna,1:dataset,T_Ai,1:dataset,T_PN,1:dataset,T_Krogh,...
    1:dataset,T_Po,1:dataset,T_G,1:dataset,T_EG,1:dataset,T_Gang);
axis tight;
legend('Nak-09','Rav-00','Ai-94','PN-89','Kro-88','Pow-85','G-79','E&G-79','Gan-96');
xlabel('Dataset'); ylabel('Temperature [°C]');
title(['Temperature [',num2str(Pressure),' kbar]']);
box on; grid on;



    
    
    
    
    
    
    
    