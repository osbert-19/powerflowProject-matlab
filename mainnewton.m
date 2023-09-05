% Set base MVA, accuracy, and maximum number of iterations and am usinf 
basemva = 100;  % Base MVA

% Load flow convergence criteria
accuracy = 0.001;  

% Maximum number of iterations
maxiter = 10;  

% Form the bus admittance matrix
lfybus  


%function is used to start a timer that measures the elapsed time =>"tic"
tic
lfnewton   % Load flow solution by Newton-Raphson method

toc %"toc" function is used to read the elapsed time from the timer started by the "tic"

busout   % Prints the power flow solution on the screen

lineflow   % Computes and displays the line flow and losses


