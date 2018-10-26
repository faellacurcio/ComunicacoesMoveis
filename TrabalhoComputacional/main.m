clc
clear
close all

%% -- Control Variables --


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
result1matrix = [];
result2matrix = [];
result3matrix = [];
result4matrix = [];
    
aux = [];
max_iterations = 100;


% Loop
for loop = 1:max_iterations
 
    %% Primeira questão item 1
    
    %Create new users randomly spread through the area
    setUsers()
    
    
    %SINR for the MACRO CELL
    SINR = powerMacro / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    
    aux(1) = capacity(1, SINR);
    
    %SINR for the PICO CELL
    SINR = powerPico / (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) + powerNoise);
    
    aux(2) = capacity(1, SINR);
    
    result1matrix = [result1matrix; aux];

    %% Primeira questão item 2

    aux = [];
    timepercentage = 0.5;
    
    %SINR for the MACRO CELL
    SINR = powerMacro / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = powerPico / (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) + powerNoise);
    aux(2) = capacity(1-timepercentage, SINR);
    
    result2matrix = [result2matrix; aux];

%% Primeira questão item 3


    aux = [];
    angle = 10;

    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(1, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) + powerNoise);
    aux(2) = capacity(1, SINR);
    
    result3matrix = [result3matrix; aux];

%% Primeira questão item 4

    aux = [];
    timepercentage = 0.5;
    
    %SINR for the MACRO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1)-PicoTower(1),TM_table(1,2)-PicoTower(2))) + powerNoise);
    aux(1) = aux(1) + capacity(timepercentage, SINR);
    
    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1)-PicoTower(1),TM_table(2,2)-PicoTower(2))) + powerNoise);
    aux(2) = aux(2) + capacity(timepercentage, SINR);
    
    result4matrix = [result4matrix; aux];
    
end

cdfCompare(result1matrix)
cdfCompare(result2matrix)
cdfCompare(result3matrix)
cdfCompare(result4matrix)

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
    global T
    global beta
    C = min((uovertT/T)*alpha*B*log2(1+SINR/beta),10^7);
        
end

function pathloss = pathLoss(d)

    global n
    global d0
    pathloss = 38.4 + 10*n*log10(d/d0);

end

function cdfCompare(datamatrix)
    figure
    cdfplot((datamatrix(:,1)))
    hold on
    cdfplot((datamatrix(:,2)))
    legend('User 1 (macro)','User 2 (pico)','Location','best')
    hold off
end

function result = normalizeData(data)
    result = (data - min(data)) / ( max(data) - min(data) );
end
