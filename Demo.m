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
% fitnessfcn = @(x)DLTZ7(x);%Ŀ�꺯��
nvars = 20;%��������
nobj =3;%Ŀ�����
narc = 50;%Pareto�����
bh_option = struct('maxgen',600,'sizestar',200,'els_max',0.3,'els_min',0.1);
%maxgen ��������
%sizestar ��Ⱥ����
%els_max ��Ӣѧϰ�����ֵ
%els_min ��Ӣѧϰ����Сֵ
[X,fval] = AMOBH(fitnessfcn,nvars,nobj,starminmax,narc,bh_option);
