clc
clear
close all

%% -- Control Variables --


ERB = [0 0 2000];

global TM_table

PicoTower = [1700 0 300];

max_users = 2;

uovertT_1 = 1;

powerMacro = 40; %(dbm)
powerPico = 7; %(dbm)
powerNoise = -115; % (dbm)
T = 0.001; %(s)
n = 4;
d0 = 1; %(m)
alpha = 0.9;
beta = 1.1; 
B = 5000; %(Hz)

%% -- Pos Excecution Control --

shouldPlotMap = 0;
% 1 for printing map
% 0 for omitting the plot

%% -- Plot elements --
if(shouldPlotMap)
    % ------------------ BEGIN function setUsers -------------------------
    
    radius_usr1 = 2000;
    radius_usr2 = 300;
    
    % randRadius = [radius*sqrt(rand(max_users,1))];
    randRadius = [radius_usr1*sqrt(rand(1,1)); radius_usr2*sqrt(rand(1,1))];

    % generate random angle for TM
    angle = 2*pi*rand(max_users,1);

    % TM_table =>
    %               |----------------------X,Y---------------------|-----------PicoTower---------|
    TM_table    =   [randRadius.*cos(angle) randRadius.*sin(angle)]+[0 0;PicoTower(1) PicoTower(2)];
    % ------------------ END function setUsers -------------------------
    
    scatter(ERB(1),ERB(2),'*')
    hold on
    scatter(PicoTower(1),PicoTower(2),'*')

    for x = 1:1:max_users
        scatter(TM_table(x,1), TM_table(x,2))
    end
    circle(ERB(1),ERB(2),ERB(3));
    circle(PicoTower(1),PicoTower(2),PicoTower(3));
    hold off
end

%% -- MAIN  --
result1matrix = [];
result2matrix = [];
result3matrix = [];
result4matrix = [];
    
aux = [];
max_iterations = 100;

timepercentage_2 = 0.5;
timepercentage_4 = 0.5;
angle_4 = 70;


% Loop
for loop = 1:max_iterations
 
    %% Primeira questão item 1
    
    %Create new users randomly spread through the area
    
    % ------------------ BEGIN function setUsers -------------------------
    
    radius_usr1 = 2000;
    radius_usr2 = 300;
    
    % randRadius = [radius*sqrt(rand(max_users,1))];
    randRadius = [radius_usr1*sqrt(rand(1,1)); radius_usr2*sqrt(rand(1,1))];

    % generate random angle for TM
    angle = 2*pi*rand(max_users,1);

    % TM_table =>
    %               |----------------------X,Y---------------------|-----------PicoTower---------|
    TM_table    =   [randRadius.*cos(angle) randRadius.*sin(angle)]+[0 0;PicoTower(1) PicoTower(2)];
    % ------------------ END function setUsers -------------------------

    
    
    %SINR for the MACRO CELL
    SINR = powerMacro / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    
    capacity_aux = min((uovertT_1)*alpha*B*log2(1+SINR/beta),10^7);
    aux(1) = capacity_aux;
    
    %SINR for the PICO CELL
    SINR = powerPico / (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) + powerNoise);
    
    capacity_aux = min((uovertT_1)*alpha*B*log2(1+SINR/beta),10^7);
    aux(2) = capacity_aux;
    
    result1matrix = [result1matrix; aux];

    %% Primeira questão item 2

    aux = [];
    
    
    %SINR for the MACRO CELL

    SINR = powerMacro / (38.4 + 10*n*log10(hypot(TM_table(1,1),TM_table(1,2))/d0) + powerNoise);
    capacity_aux = min((timepercentage_2)*alpha*B*log2(1+SINR/beta),10^7);
    aux(1) = capacity_aux;
    
    
    %SINR for the PICO CELL
    SINR = powerPico / (38.4 + 10*n*log10(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))/d0) + powerNoise);
    capacity_aux = min((1 - timepercentage_2)*alpha*B*log2(1+SINR/beta),10^7);
    aux(2) = capacity_aux;
    
    result2matrix = [result2matrix; aux];

%% Primeira questão item 3


    aux = [];
    angle = 10;

    %SINR for the MACRO CELL
    
    SINR = (powerMacro + min(20,360/angle)) / (38.4 + 10*n*log10(hypot(TM_table(1,1),TM_table(1,2))/d0) + powerNoise);
    capacity_aux = min((1)*alpha*B*log2(1+SINR/beta),10^7);
    aux(1) = capacity_aux;
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (38.4 + 10*n*log10(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))/d0) + powerNoise);
    capacity_aux = min((1)*alpha*B*log2(1+SINR/beta),10^7);
    aux(2) = capacity_aux;
    
    result3matrix = [result3matrix; aux];

%% Primeira questão item 4

    aux = [];
    
    %SINR for the MACRO CELL
    SINR = (powerPico + min(20,360/angle_4)) / (38.4 + 10*n*log10(hypot(TM_table(1,1),TM_table(1,2))/d0) + powerNoise);
    capacity_aux = min((timepercentage_4)*alpha*B*log2(1+SINR/beta),10^7);
    aux(1) = capacity_aux ;
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle_4)) / (38.4 + 10*n*log10(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2))/d0) + powerNoise);
    capacity_aux = min((timepercentage_4)*alpha*B*log2(1+SINR/beta),10^7);
    aux(1) = aux(1) + capacity_aux;
    
    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle_4)) / (38.4 + 10*n*log10(hypot(TM_table(2,1),TM_table(2,2))/d0) + powerNoise);
    capacity_aux = min((1- timepercentage_4)*alpha*B*log2(1+SINR/beta),10^7);
    aux(2) = capacity_aux;
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle_4)) / (38.4 + 10*n*log10(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))/d0) + powerNoise);
    capacity_aux = min((1- timepercentage_4)*alpha*B*log2(1+SINR/beta),10^7);
    aux(2) = aux(2) + capacity_aux;
    
    result4matrix = [result4matrix; aux];
    
end


% cdfCompare(result1matrix, 1)
    date_str = datestr(now);
    
    filename =   strcat('fig/', num2str(1), '-', strrep(date_str(13:end),':','_'),  '.png');
    
    aux = figure;
    cdfplot(real(result1matrix(:,1)))
    hold on
    cdfplot(real(result1matrix(:,2)))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(1), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(1), ' macroTime: ', num2str(timepercentage_4) , ' angle: ', num2str(min(20,(360/angle_4)) ) ) )
    else
        title(strcat('Cenário ', int2str(1)))
    end
    
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
    
    saveas(aux,filename,'png')
% 
% cdfCompare(result2matrix, 2)
    date_str = datestr(now);
    
    filename =   strcat('fig/', num2str(2), '-', strrep(date_str(13:end),':','_'),  '.png');
    
    aux = figure;
    cdfplot(real(result2matrix(:,1)))
    hold on
    cdfplot(real(result2matrix(:,2)))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(2), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(2), ' macroTime: ', num2str(timepercentage_4) , ' angle: ', num2str(min(20,(360/angle_4)) ) ) )
    else
        title(strcat('Cenário ', int2str(2)))
    end
    
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
    
    saveas(aux,filename,'png')
% 
% cdfCompare(result3matrix, 3)
    date_str = datestr(now);
    
    filename =   strcat('fig/', num2str(3), '-', strrep(date_str(13:end),':','_'),  '.png');
    
    aux = figure;
    cdfplot(real(result3matrix(:,1)))
    hold on
    cdfplot(real(result3matrix(:,2)))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(3), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(3), ' macroTime: ', num2str(timepercentage_4) , ' angle: ', num2str(min(20,(360/angle_4)) ) ) )
    else
        title(strcat('Cenário ', int2str(3)))
    end
    
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
    
    saveas(aux,filename,'png')
% 
% cdfCompare(result4matrix, 4)
    date_str = datestr(now);
    
    filename =   strcat('fig/', num2str(4), '-', strrep(date_str(13:end),':','_'),  '.png');
    
    aux = figure;
    cdfplot(real(result4matrix(:,1)))
    hold on
    cdfplot(real(result4matrix(:,2)))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(4), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(4), ' macroTime: ', num2str(timepercentage_4) , ' angle: ', num2str(min(20,(360/angle_4)) ) ) )
    else
        title(strcat('Cenário ', int2str(4)))
    end
    
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
    
    saveas(aux,filename,'png')


%% -- Function declarations --

% -- Set Users  --


% Desconsiderar cdigo abaixo

data = [];
% ------------------ BEGIN Normalize data -------------------------
result = (data - min(data)) / ( max(data) - min(data) );
% ------------------ END  Normalize data -------------------------

% ------------------ BEGIN function setUsers -------------------------
% ------------------ END function setUsers -------------------------