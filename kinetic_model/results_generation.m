clc;clear;
model = 'kinetic_model';
open(model);
load_system(model);
files_path=glob('../data/WithLabel/1*.xlsx');
liq_gas_ratio=[];
for file_num =1:size(files_path,1)
    [filepath,name,ext] = fileparts(files_path(file_num));
    name
    input_df = xlsread(char(files_path(file_num)),1);
    n_set = size(input_df,1);
    ini_time = 600; %s
    run_time = 43*n_set;
    rho_liq = 1004; % kg/m3
    x0=[500.6,1,1005.5,4380,1,1263];
    
    %Variable: Fliq, Fvap, Ca, Tliq_in, P_bot, P_top
    
    P_bot = (input_df(:,31)+1.013).*1E5; %Pa
    P_bot = timeseries(P_bot,(ini_time:43:ini_time+run_time-1));
    P_bot =addsample(P_bot, 'Data', P_bot.data(1), 'Time', 0);
    
    P_top = (input_df(:,37)+1.013).*1E5; %Pa
    P_top = timeseries(P_top,(ini_time:43:ini_time+run_time-1));
    P_top =addsample(P_top, 'Data', P_top.data(1), 'Time', 0);
    
    avgP = ((P_bot+P_top)/2); %Pa
    
    Tliq_in = 0.5*(input_df(:,65)+input_df(:,66))+273.15; %K
    Tliq_in = timeseries(Tliq_in,(ini_time:43:ini_time+run_time-1));
    Tliq_in =addsample(Tliq_in, 'Data', Tliq_in.data(1), 'Time', 0);
    
    Tvap_in = input_df(:,43)+273.15; %K
    Tvap_in = timeseries(Tvap_in,(ini_time:43:ini_time+run_time-1));
    Tvap_in =addsample(Tvap_in, 'Data', Tvap_in.data(1), 'Time', 0);
    
    Vm = 8.314*Tvap_in/avgP; % m3/mol
    
    FCO2_in = input_df(:,17); % m3/h
    FCO2_in = timeseries(FCO2_in,(ini_time:43:ini_time+run_time-1));
    FCO2_in =addsample(FCO2_in, 'Data', FCO2_in.data(1), 'Time', 0);
    
    PCO2_in = (input_df(:,36)+1.013).*1E5; %Pa
    PCO2_in = timeseries(PCO2_in,(ini_time:43:ini_time+run_time-1));
    PCO2_in =addsample(PCO2_in, 'Data', PCO2_in.data(1), 'Time', 0);
    
    TCO2_in = input_df(:,87)+273.15; %K
    TCO2_in = timeseries(TCO2_in,(ini_time:43:ini_time+run_time-1));
    TCO2_in =addsample(TCO2_in, 'Data', TCO2_in.data(1), 'Time', 0);
    
    FN2_in = input_df(:,19); % m3/h
    FN2_in = timeseries(FN2_in,(ini_time:43:ini_time+run_time-1));
    FN2_in =addsample(FN2_in, 'Data', FN2_in.data(1), 'Time', 0);
    
    PN2_in = (input_df(:,37)+1.013).*1E5; %Pa
    PN2_in = timeseries(PN2_in,(ini_time:43:ini_time+run_time-1));
    PN2_in =addsample(PN2_in, 'Data', PN2_in.data(1), 'Time', 0);
    
    TN2_in = input_df(:,76)+273.15; %K
    TN2_in = timeseries(TN2_in,(ini_time:43:ini_time+run_time-1));
    TN2_in =addsample(TN2_in, 'Data', TN2_in.data(1), 'Time', 0);
    
    Fliq_in = input_df(:,7)/rho_liq/3600;
    Fliq_in = timeseries(Fliq_in,(ini_time:43:ini_time+run_time-1));
    Fliq_in =addsample(Fliq_in, 'Data', Fliq_in.data(1), 'Time', 0);
    
    Fvap = FCO2_in*(PCO2_in/P_bot)*(Tvap_in/TCO2_in)+FN2_in*(PN2_in/P_bot)*(Tvap_in/TN2_in); %m3/h
    Fvap = Fvap/3600; %m3/s
    
    Origin_CA = input_df(:,3)/100;
    label = input_df(:,end);
    first_six=false;
    first_conc=0;
    for i = 1:n_set
        if label(i) ~= 6
            Origin_CA(i)=nan;
        elseif(first_six == false)
            first_conc=Origin_CA(i);
            first_six=true;
        end
    end
    prev_conc=first_conc;
    for i = 1:n_set
        if label(i) ~= 6
            Origin_CA(i)=prev_conc;
        else
            prev_conc=Origin_CA(i);
        end
    end
    Origin_CA = timeseries(Origin_CA,(ini_time:43:ini_time+run_time-1));
    Origin_CA =addsample(Origin_CA, 'Data', first_conc, 'Time', 0);
    C_Ain = Origin_CA/Vm;
    
    new_initial = ones(1,6);
    new_initial(1)=Fliq_in.Data(1);
    new_initial(2)=Fvap.Data(1);
    new_initial(3)= C_Ain.Data(1); %CA mol/m3
    new_initial(4) = 2465.619; %CB
    new_initial(5) = 0; %#CC
    new_initial(6) = Tliq_in.Data(1);
    initial_mat = new_initial;
    initial_mat(3)=0;
    input_str =['[' num2str(x0(1)) ';' num2str(x0(2)) ';' num2str(x0(3)) ';' num2str(x0(4)) ';' ...
            num2str(x0(5)) ';' num2str(x0(6)) ';]'];
    ini_cond_str=['[' num2str(initial_mat(1)) ';' num2str(initial_mat(2)) ';' num2str(initial_mat(3)) ';' ...
            num2str(initial_mat(4)) ';' num2str(initial_mat(5)) ';' num2str(initial_mat(6)) ';]'];
    
    set_param('kinetic_model/m file','parameters', [input_str ',' ini_cond_str]);
    set_param('kinetic_model/m file1','parameters', [input_str ',' ini_cond_str]);
    set_param('kinetic_model/m file2','parameters', [input_str ',' ini_cond_str]);
    set_param('kinetic_model/m file3','parameters', [input_str ',' ini_cond_str]);
    set_param('kinetic_model/m file4','parameters', [input_str ',' ini_cond_str]);
    set_param(model,'StartTime','0','StopTime',num2str(ini_time+run_time));

    liq_gas_ratio=[liq_gas_ratio, Fliq_in.Data(1)/Fvap.Data(1)];
    
    simOut = sim(model);
    CA_log = zeros(n_set,6);
    CA_log(:,1)=CA1(end-run_time+1:43:end);
    CA_log(:,2)=CA2(end-run_time+1:43:end);
    CA_log(:,3)=CA3(end-run_time+1:43:end);
    CA_log(:,4)=CA4(end-run_time+1:43:end);
    CA_log(:,5)=CA5(end-run_time+1:43:end);
    CA_log(:,6)=Origin_CA.Data(end-n_set+1:end);
    
    writematrix(CA_log,[name,'.csv']);
end