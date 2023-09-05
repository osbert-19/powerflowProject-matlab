
% NAME: BARIGYE OSBERT  STUDENT NUMBER :2021171055  POWER ENGINEERING FINAL


% Clear all variables and command window
clear global
clc

% Set the bus system number to 30 and store it in "ch" variable 
ch = 30;

% Load data for bus system 30 from the data30.m file 
data30

% Perform load flow analysis using Newton-Raphson method
mainnewton

