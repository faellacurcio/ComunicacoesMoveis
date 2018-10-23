clc
clear
close all

%% -- Control Variables --


ERB = [0 0 2000];

global TM_table

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
    scatter(ERB(1),ERB(2),'*')
    hold on
    scatter(PicoCell(1),PicoCell(2),'*')

    for x = 1:1:max_users
        scatter(TM_table(x,1), TM_table(x,2))
    end
    circle(ERB(1),ERB(2),ERB(3));
    circle(PicoCell(1),PicoCell(2),PicoCell(3));
    hold off
end

%% -- MAIN  --
result1matrix = [];
aux = [];
max_iterations = 100;

%% Primeira quest�o item 1

for loop = 1:max_iterations
    
    %Create new users randomly spread through the area
    setUsers()
    
    %SINR for the MACRO CELL
    SINR = powerMacro / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(1, SINR);
    
    %SINR for the PICO CELL
    SINR = powerPico / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = capacity(1, SINR);
    
    result1matrix = [result1matrix; aux];
end

%% Primeira quest�o item 2

result2matrix = [];
aux = [];
timepercentage = 0.5;

for loop = 1:max_iterations
    
    %Create new users randomly spread through the area
    setUsers()
    
    %SINR for the MACRO CELL
    SINR = powerMacro / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = powerPico / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = capacity(1-timepercentage, SINR);
    
    result2matrix = [result2matrix; aux];
end

%% Primeira quest�o item 3

result3matrix = [];
aux = [];
angle = 10;

for loop = 1:max_iterations
    
    %Create new users randomly spread through the area
    setUsers()
    
    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(1, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = capacity(1, SINR);
    
    result3matrix = [result3matrix; aux];
end

%% Primeira quest�o item 4

result4matrix = [];
aux = [];
timepercentage = 0.5;

for loop = 1:max_iterations
    
    %Create new users randomly spread through the area
    setUsers()
    
    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(1,1),TM_table(1,2))) + powerNoise);
    aux(1) = aux(1) + capacity(timepercentage, SINR);
    
    %SINR for the MACRO CELL
    SINR = (powerMacro + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = capacity(timepercentage, SINR);
    
    %SINR for the PICO CELL
    SINR = (powerPico + min(20,360/angle)) / (pathLoss(hypot(TM_table(2,1),TM_table(2,2))) + powerNoise);
    aux(2) = aux(2) + capacity(timepercentage, SINR);
    
    result4matrix = [result4matrix; aux];
    
end
%% -- Function declarations --

% -- Set Users  --
function setUsers()
    
    global max_users
    
    radius_usr1 = 2000;
    radius_usr2 = 300;
    
    PicoCell = [1700 0 300];
    
    % randRadius = [radius*sqrt(rand(max_users,1))];
    randRadius = [radius_usr1*sqrt(rand(1,1)); radius_usr2*sqrt(rand(1,1))];

    % generate random angle for TM
    angle = 2*pi*rand(max_users,1);

    % TM_table =>
    %               |----------------------X,Y---------------------|-----------Picocell---------|
    global TM_table
    TM_table    =   [randRadius.*cos(angle) randRadius.*sin(angle)]+[0 0;PicoCell(1) PicoCell(2)];
end

function C = capacity(uovetT, SINR) %(bits/s)
    global B
    global alpha
    global T
    global beta
    C = min((uovetT/T)*alpha*B*log2(1+SINR/beta),10^7);
end

function pathloss = pathLoss(d)
    global n
    global d0
    pathloss = 38.4 + 10*n*log10(d/d0);
end


