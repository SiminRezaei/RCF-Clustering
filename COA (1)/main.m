
%% 
% Coati Optimization Algorithm: A New Bio-Inspired Metaheuristic Algorithm for Solving Optimization Problems
% Knowledge-Based Systems
% Mohammad Dehghani, Zeinab Montazeri and Pavel Trojovský1
% Department of Mathematics, Faculty of Science, University of Hradec Králové, 50003 Hradec Králové, Czech Republic
% " Optimizer"
%%
clc
clear
close all
%%
%%
Fun_name='F5'; % number of test functions: 'F1' to 'F23'
SearchAgents=30;                      % number of Coati (population members) 
Max_iterations=1000;                  % maximum number of iteration
[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name); % Object function information
[Best_score,Best_pos,COA_curve]=COA(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);  % Calculating the solution of the given problem using COA 
%%
display(['The best solution obtained by COA for ' [num2str(Fun_name)],'  is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by COA  for ' [num2str(Fun_name)],'  is : ', num2str(Best_score)]);
        