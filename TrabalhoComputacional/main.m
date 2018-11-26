clc
clear
close all

%% -- Control Variables --

interestingPlot = [];

ERB = [0 0 2000];

global TM_table

global PicoTower

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

powerMacro = 40 - 30; %(db)
powerPico = 7 - 30; %(db)
powerNoise = -115 - 30; % (db)
T = 0.001; %(s)
n = 4;
d0 = 1; %(m)
alpha = 0.9;
beta = 1.1; 
B = 5000; %(Hz)

%% -- Pos Excecution Control --

shouldPlotMap = 1;
% 1 for printing map
% 0 for omitting the plot

%% -- Plot elements --
if(shouldPlotMap)
    setUsers()
    
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

max_iterations = 1000;

global timepercentage_2;
global timepercentage_4;
global angle_3;

timepercentage_2 = 1;
timepercentage_4 = 1;
angle_3 = 18;

for loopout = 1:1

    result1matrix = [];
    result2matrix = [];
    result3matrix = [];
    result4matrix = [];

    aux = [];

    timepercentage_2 = timepercentage_2 - 0.05;
    timepercentage_4 = timepercentage_4 - 0.05;

    angle_3 = angle_3 + loopout;


    % Loop
    for loop = 1:max_iterations

        %Create new users randomly spread through the area
%         setUsers()

        %% Primeira questão item 1
        %


        %SINR for the MACRO CELL
        potencia_recebida = powerMacro - pathLoss(hypot(TM_table(1,1),TM_table(1,2)));
        potencia_interferente = powerPico - pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2)));
        SINR = potencia_recebida - lin2db(db2lin(potencia_interferente) + db2lin(powerNoise));

        aux(1) = capacity(1, db2lin(SINR));

        %SINR for the PICO CELL
        potencia_recebida = powerMacro - pathLoss(hypot(TM_table(2,1),TM_table(2,2)));
        potencia_interferente = powerPico - pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2)));
        SINR = potencia_recebida - lin2db(db2lin(potencia_interferente) + db2lin(powerNoise));
        
        aux(2) = capacity(1, db2lin(SINR));

        result1matrix = [result1matrix; aux];

        %% Primeira questão item 2
        % SINR é 
        aux = [];

        %SINR for the MACRO CELL
        SINR = powerMacro - (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) - powerNoise);
        aux(1) = capacity(timepercentage_2, db2lin(SINR));

        %SINR for the PICO CELL
        SINR = powerPico - (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) - powerNoise);
        aux(2) = capacity(1-timepercentage_2, db2lin(SINR));

        result2matrix = [result2matrix; aux];

    %% Primeira questão item 3


        aux = [];
        angle = 10;

        %SINR for the MACRO CELL
        SINR = (powerMacro + min(20,360/angle)) - (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) - powerNoise);
        aux(1) = capacity(1, db2lin(SINR));

        %SINR for the PICO CELL
        SINR = (powerPico + min(20,360/angle)) - (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) - powerNoise);
        aux(2) = capacity(1, db2lin(SINR));

        result3matrix = [result3matrix; aux];

    %% Primeira questão item 4

        aux = [];

        %SINR for the MACRO CELL
        SINR = (powerPico) - (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) - powerNoise);
        aux(1) = capacity(timepercentage_4, db2lin(SINR));

        %SINR for the PICO CELL
        SINR = (powerPico) - (pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2))) - powerNoise);
        aux(1) = aux(1) + capacity(timepercentage_4, db2lin(SINR));

        %SINR for the MACRO CELL
        SINR = (powerMacro) - (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) - powerNoise);
        aux(2) = capacity(1 - timepercentage_4, db2lin(SINR));

        %SINR for the PICO CELL
        SINR = (powerPico) - (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) - powerNoise);
        aux(2) = aux(2) + capacity(1 - timepercentage_4, db2lin(SINR));

        result4matrix = [result4matrix; aux];

    end

cdfCompare(result1matrix, 1)

% cdfCompare(result2matrix, 2)

interestingPlot(loopout,:) = mean(result2matrix,1);

% cdfCompare(result3matrix, 3)

% cdfCompare(result4matrix, 4)

% close all
end

%% -- Function declarations --

% -- Set Users  --
function setUsers()
    
    global max_users
    global PicoTower
    radius_usr1 = 2000;
    radius_usr2 = 300;
    
    % randRadius = [radius*sqrt(rand(max_users,1))];
    randRadius = [radius_usr1*sqrt(rand(1,1)); radius_usr2*sqrt(rand(1,1))];

    % generate random angle for TM
    angle = 2*pi*rand(max_users,1);

    % TM_table =>
    %               |----------------------X,Y---------------------|-----------PicoTower---------|
    global TM_table
    TM_table    =   [randRadius.*cos(angle) randRadius.*sin(angle)]+[0 0;PicoTower(1) PicoTower(2)];
    
end

function C = capacity(uovertT, SINR) %(bits/s)
    
    global B
    global alpha
    global beta
    C = min((uovertT)*alpha*B*log2(1+dbm2lin(SINR)/beta),10^7);

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
%     cdfplot(datamatrix(:,2))
    
    if (extra_param == 2 )
        title(strcat('Cenário ', int2str(extra_param), ' macroTime: ', num2str(timepercentage_2 ) ) )
    elseif (extra_param == 3)
        title(strcat('Cenário ', int2str(extra_param), ' angle: ', num2str(min(20,(360/angle_3)) ) ) )
    elseif(extra_param == 4)
        title(strcat('Cenário ', int2str(extra_param), ' macroTime: ', num2str(timepercentage_4) ) )
    else
        title(strcat('Cenário ', int2str(extra_param)))
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

function result = db2lin(input)
    result = 10.^((input)/10);
end

function result = lin2dbm(input)
    result = 10*log10(input) + 30;
end

function result = lin2db(input)
    result = 10*log10(input);
end