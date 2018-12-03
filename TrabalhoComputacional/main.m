clc
clear
close all

set(groot,'defaultLineLineWidth',2)

%% -- Control Variables --

interestingPlot = [];
reference = [2000 0];


global TM_table

global MacroTower
global PicoTower

MacroTower = [0 0 2000];
PicoTower = [1700 0 300];

global max_users

max_users = 2;

global powerMacro
global powerPico 
global T 
global n 
global d0 
global powerNoise 
global alpha 
global beta 
global B 

powerMacro = 40-30; %(db)
powerPico = 7-30; %(db)
powerNoise = -115-30; % (db)
T = 0.001; %(s)
n = 4;
d0 = 1; %(m)
alpha = 0.9;
beta = 1.1; 
B = 5*10^6; %(Hz)

max_iterations = 10000;

global cenario3_interf
cenario3_interf = [];

global timepercentage_2;
global timepercentage_4;
global angle_3;

timepercentage_2 = 1;
timepercentage_4 = 1;
angle_3 = 18;


%% -- Pos Excecution Control --

shouldPlotMap = 1;
% 1 for printing map
% 0 for omitting the plot

%% -- Plot elements --
if(shouldPlotMap)
    setUsers()
    plotScenario()
end

%% -- MAIN  --

for loopout = 1:18

    result1matrix = [];
    result2matrix = [];
    result3matrix = [];
    result4matrix = [];
    
    cenario3_interf = [cenario3_interf; 0];

    aux = [];

    timepercentage_2 = timepercentage_2 - 0.05;
    timepercentage_4 = timepercentage_4 - 0.05;

    angle_3 = angle_3 + loopout;
    disp(angle_3)


    % Loop
    for loop = 1:max_iterations

        %Create new users randomly spread through the area
        setUsers()

    %% Primeira questão item 1
        %

        %SINR for the MACRO CELL
        potencia_recebida1 = powerMacro - pathLoss(hypot(TM_table(1,1),TM_table(1,2)));
        potencia_interferente1 = powerPico - pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2)));
        SINR1 = potencia_recebida1 - lin2db(db2lin(potencia_interferente1) + db2lin(powerNoise));

        aux(1) = capacity(1, SINR1);


        %SINR for the PICO CELL
        potencia_recebida2 = powerPico - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2)));
        potencia_interferente2 = powerMacro - pathLoss(hypot(TM_table(2,1),TM_table(2,2)));
        SINR2 = potencia_recebida2 - lin2db(db2lin(potencia_interferente2) + db2lin(powerNoise));
        
        aux(2) = capacity(1, SINR2);

        result1matrix = [result1matrix; aux];

     %% Primeira questão item 2
        % SINR é a mesma SNR já que há multiplexação no tempo
        
        aux = [];

        %SINR for the MACRO CELL
        SINR = powerMacro - pathLoss(hypot(TM_table(1,1),TM_table(1,2))) - powerNoise;
        aux(1) = capacity(timepercentage_2, SINR);

        %SINR for the PICO CELL
        SINR = powerPico - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) - powerNoise;
        aux(2) = capacity(1 - timepercentage_2, SINR);

        result2matrix = [result2matrix; aux];

    %% Primeira questão item 3
        aux = [];
        
        angle_user_macro_from_macro = atan2d(TM_table(1,2), TM_table(1,1));
        angle_user_pico_from_macro = atan2d(TM_table(2,2), TM_table(2,1));
        
        angle_user_macro_from_pico = atan2d(TM_table(1,2)-PicoTower(2), TM_table(1,1)-PicoTower(1));
        angle_user_pico_from_pico = atan2d(TM_table(2,2)-PicoTower(2), TM_table(2,1)-PicoTower(1));
        
        if(abs(angle_user_macro_from_macro - angle_user_pico_from_macro) <= angle_3/2)
            interference_macro2pico = true;
        else
            interference_macro2pico = false;
        end
            
        if(abs(angle_user_macro_from_pico - angle_user_pico_from_pico) <= angle_3/2)
            interference_pico2macro = true;
        else
            interference_pico2macro = false;
        end
        
        if(interference_macro2pico || interference_pico2macro)
           cenario3_interf(loopout) = cenario3_interf(loopout) + 1;
        end
            
        %SINR for the MACRO CELL
        if(interference_pico2macro)
            potencia_recebida1 = (powerMacro + min(20,360/angle_3)) -  pathLoss( hypot(TM_table(1,1),TM_table(1,2)));
            potencia_interferente1 = (powerPico + min(20,360/angle_3)) - pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2)));
            SINR1 = potencia_recebida1 - lin2db(db2lin(potencia_interferente1) + db2lin(powerNoise));
        else
            SINR1 = (powerMacro + min(20,360/angle_3)) -  pathLoss( hypot(TM_table(1,1),TM_table(1,2))) - powerNoise ;
        end
        aux(1) = capacity(1, SINR1);

        %SINR for the PICO CELL
        if(interference_macro2pico)
            potencia_recebida2 = (powerPico + min(20,360/angle_3)) - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2)));
            potencia_interferente2 = (powerMacro + min(20,360/angle_3)) -  pathLoss(hypot(TM_table(2,1),TM_table(2,2)));
            SINR2 = potencia_recebida2 - lin2db(db2lin(potencia_interferente2) + db2lin(powerNoise));
        else
            SINR2 = (powerPico + min(20,360/angle_3)) - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) - powerNoise;
        end
        aux(2) = capacity(1, SINR2);

        result3matrix = [result3matrix; aux];

    %% Primeira questão item 4

        aux = [];

        %SINR for the user MACRO CELL
        SINR = lin2db(...
                      db2lin((powerMacro) - pathLoss(hypot(TM_table(1,1),TM_table(1,2)))) ...
                    + db2lin((powerPico) - pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2))))...
                ) ...
                - powerNoise;
            
        aux(1) = capacity(timepercentage_4, SINR);

        %SINR for the user PICO CELL
        SINR = lin2db(...
                      db2lin((powerMacro) - pathLoss(hypot(TM_table(2,1),TM_table(2,2)))) ...
                    + db2lin((powerPico) - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))))...
                ) ...
                - powerNoise;
            
        aux(2) = capacity(1 - timepercentage_4, SINR);

        result4matrix = [result4matrix; aux];

    end

cdfCompare(result1matrix, 1)

cdfCompare(result2matrix, 2)

interestingPlot(loopout,:) = mean(result2matrix,1);

cdfCompare(result3matrix, 3)

cdfCompare(result4matrix, 4)

close all
end

% %plot InterestingPlot (scenario 3)
% plot(0.95:-0.05:0.1,interestingPlot(:,1),'-x');
% hold on; 
% plot(0.95:-0.05:0.1,interestingPlot(:,2),'-x');
% hold off;

%% -- Function declarations --

% -- Set Users  --
function setUsers()
    
    global max_users
    global PicoTower
    radius_Macro = 2000;
    radius_Pico = 300;
    
    % randRadius = [radius*sqrt(rand(max_users,1))];
    randRadius = [radius_Macro*sqrt(rand(1,1)); radius_Pico*sqrt(rand(1,1))];

    % generate random angle for TM
    angle = 2*pi*rand(max_users,1);

    % TM_table =>
    %               |----------------------X,Y---------------------|-----------PicoTower---------|
    global TM_table
    TM_table    =   [randRadius.*cos(angle) randRadius.*sin(angle)]+[0 0;PicoTower(1) PicoTower(2)];
    
end

function plotScenario()
    global MacroTower
    global PicoTower
    global max_users
    global TM_table
    global angle_3
    markerSize = 100;
    
    %Plot antennas
    scatter(MacroTower(1),MacroTower(2),markerSize,'d')
    hold on
    scatter(PicoTower(1),PicoTower(2),markerSize ,'d', 'filled')
    
    % Plot Users
    for x = 1:1:max_users
        scatter(TM_table(x,1), TM_table(x,2), markerSize, 'filled')
    end
    
    % Plot antennas range
    circle(MacroTower(1),MacroTower(2),MacroTower(3));
    circle(PicoTower(1),PicoTower(2),PicoTower(3));
    
%     title('Without Beamforming')
    xlabel('Distance (m)')
    ylabel('Distance (m)')
    grid on
    hold on
    
    
    %plot beamforming angle 
    angle1 = atan2d(TM_table(1,2), TM_table(1,1))+angle_3/2;
    xplot = [0 cosd(angle1)*MacroTower(3)];
    yplot = [0 sind(angle1)*MacroTower(3)];
    plot(xplot, yplot, 'LineWidth',1)
    angle2 = atan2d(TM_table(1,2), TM_table(1,1))-angle_3/2;
    xplot = [0 cosd(angle2)*MacroTower(3)];
    yplot = [0 sind(angle2)*MacroTower(3)];
    plot(xplot, yplot, 'LineWidth',1)
    
    angle3 = atan2d(TM_table(2,2)-PicoTower(2), TM_table(2,1)-PicoTower(1))+angle_3/2;
    xplot = [PicoTower(1) PicoTower(1)+cosd(angle3)*MacroTower(3)];
    yplot = [PicoTower(2) sind(angle3)*MacroTower(3)];
    plot(xplot, yplot, 'LineWidth',1)
    angle4 = atan2d(TM_table(2,2)-PicoTower(2), TM_table(2,1)-PicoTower(1))-angle_3/2;
    xplot = [PicoTower(1) PicoTower(1)+cosd(angle4)*MacroTower(3)];
    yplot = [PicoTower(2) sind(angle4)*MacroTower(3)];
    plot(xplot, yplot, 'LineWidth',1)
    axis([-2500 2500 -2500 2500])
    
    legend('Antenna Macro','Antenna Pico','User Macro','User Pico')
    
    grid on
    hold off
end

function C = capacity(uovertT, SINR) %(bits/s)
    
    global B
    global alpha
    global beta
    C = min((uovertT)*alpha*B*log2(1+db2lin(SINR)/beta),10^7);

end

function pathloss = pathLoss(d)

    global n
    global d0
    pathloss = 38.4 + 10*n*log10(d/d0);

end

function cdfCompare(datamatrix, extra_param)
    
    global timepercentage_2;
    global timepercentage_4;
    global angle_3;

    
    date_str = datestr(now);
    
    filename =   strcat('fig/', num2str(extra_param), '-', strrep(date_str(13:end),':','_'),  '.png');
    
    aux = figure;
    cdfplot(datamatrix(:,1))
    hold on
    cdfplot(datamatrix(:,2))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(extra_param), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif (extra_param == 3)
        title(strcat('Cenário ', int2str(extra_param), ' gain: ', num2str(min(20,(360/angle_3)) ), ' angle', num2str(angle_3) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(extra_param), ' macroTime: ', num2str(timepercentage_4) ) )
    else
        title(strcat('Cenário ', int2str(extra_param)))
%         set(gca, 'XScale', 'log')
    end
    
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
    
    saveas(aux,filename,'png')
    
end

function result = normalizeData(data)
    result = (data - min(data)) / ( max(data) - min(data) );
end

function result = dbm2lin(input)
    result = 10.^((input-30)/10);
end

function result = lin2dbm(input)
    result = 10*log10(input) + 30;
end

function result = db2lin(input)
    result = 10.^((input)/10);
end

function result = lin2db(input)
    result = 10*log10(input);
end