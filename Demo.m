close all;
clc;clear;
starminmax = [1 0;
%               5 -5;
%               5 -5;
%               5 -5;
%               5 -5;
%               5 -5;
%               5 -5;
%               5 -5;         
%               5 -5;
%               5 -5;
              1 0;
              1 0;              
              1 0;
              1 0;              
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;              
              1 0;
              1 0;              
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;
              1 0;
%               1 0;
%               1 0;              
%               1 0;
%               1 0;              
%               1 0;
%               1 0;
%               1 0;
%               1 0;
%               1 0;
% % % % % % % %               1 0;              
% % % % % % % % % %              
];
% % % 
fitnessfcn = @(x)DLTZ7(x);
% starminmax = [10 -10];
% fitnessfcn = @(x)DLTZ7(x);%目标函数
nvars = 20;%变量个数
nobj =3;%目标个数
narc = 50;%Pareto解个数
bh_option = struct('maxgen',600,'sizestar',200,'els_max',0.3,'els_min',0.1);
%maxgen 迭代次数
%sizestar 种群个数
%els_max 精英学习率最大值
%els_min 精英学习率最小值
[X,fval] = AMOBH(fitnessfcn,nvars,nobj,starminmax,narc,bh_option);
