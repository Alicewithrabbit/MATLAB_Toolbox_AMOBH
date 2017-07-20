%%
% author: Alicewithrabbit (Wu Chong)
% email: imroxaswc@gmail.com
% description: It's a multi-objective optimization toolbox using black hole
% algorithm.
%
% fitnessfcn is the objective function
% nvars is the number of variables
% nobj is the number of objectives
% starminmax is the bound of every variables
% narc is the size of global archive
% bh_option is the setting of bh_algorithm
% Syntax:
%         [X,fval] = AMOBH(fitnessfcn,nvars,nobj,starminmax,narc)
%         [X,fval] = AMOBH(fitnessfcn,nvars,nobj,starminmax,narc,bh_option)
% In most cases, we suggest to set bh_option =
% struct('maxgen',600,'sizestar',300,'els_max',0.3,'els_min',0.1) to get
% good performance and time complexity. If the problem is complex we
% suggest struct('maxgen',600,'sizestar',300,'els_max',0.3,'els_min',0.1)
% version: 1.0.0
function  [X,fval] = AMOBH(fitnessfcn,nvars,nobj,starminmax,narc,bh_option)

global h;
global se;
global dse;
global APS;
global APF;

if nargin == 5
    bh_option = struct('maxgen',600,'sizestar',200,'els_max',0.3,'els_min',0.1);
end



if nobj == 1
    disp('Error:The number of objectives must larger than 1 ')
    return;
end
K = narc;
KT = 2*(nobj+1);
N = nvars;
M = nobj;
ELS_MAX = bh_option.els_max;%Maximum of learning rate;
ELS_MIN = bh_option.els_min;%Minimum of learning rate;
ELS = ELS_MAX;
P = 0.3;
step_x = (ELS_MAX - ELS_MIN)/bh_option.maxgen;%step size of learning;
delta_s = 2/(N*K)*log2(2);
status = 'convergence';
gArchive = [];



%% 初始化边界
if N == 1
    starmax(1) = starminmax(1,1);
    starmin(1) = starminmax(1,2);
else
    for i = 1:N
        starmax(i) = starminmax(i,1);
        starmin(i) = starminmax(i,2);
    end
end

solution_vector = [];
gArchive = [];
%%初始化种群
for i=1:bh_option.sizestar
    
    if N ~=1
        % ramdon inital stars
        for j = 1:N
            star(i,j) = (starmax(j) - starmin(j)) * rand + starmin(j);
        end
    else
        star(i,1) = (starmax(1) - starmin(1)) * rand + starmin(1);
    end
    % calculate the fitness
    solution_vector(i,:) = fitnessfcn(star(i,:));
    
end
p = [];
% 初始化Pareto解
for i = 1:M
    [p(i),q] = find(solution_vector(:,i) == min(solution_vector(:,i)));
end
p_new = unique(p);

m = length(p_new);

for i = 1:m
    gArchive(i,1:M) = solution_vector(p_new(i),:);
    gArchive(i,M+1:M+N) = star(p_new(i),:);
end

% gArchive = [solution_vector,star];


%%
% map the fitness into PCCS

n = size(gArchive,1);

if ~isfield(bh_option,'gui')
figure(1)
set(gcf, 'position', [300 0 300 900]);
j = 0:1:K;
i = 0:M;
[x,y] = meshgrid(i,j);
plot(x,y,'k')
hold on
plot(x',y','k');

title('Pareto Front')
set(gca,'xticklabel',{''});
set(gca,'yticklabel',{''});
xlabel('Objective')
ylabel('Label number')

axis([0 M 0 K])
end

for i = 1:n
    for j = 1:M
        L(i,j) = ceil(n*(gArchive(i,j) - min(gArchive(1:n,j)))/(max(gArchive(1:n,j)) -...
            min(gArchive(1:n,j))));
        
        
        if gArchive(i,j) == min(gArchive(1:n,j))
            L(i,j) = 1;
        end
        
        h = plot(-100,-100,'mp', 'MarkerSize',5);
        
    end
end



Entropy  = zeros(1,bh_option.maxgen);

Entropy_temp = Entropy(1);

% calculate the entropy
for t = 1:M
    for k = 1:n
        if(~isempty(find(L(:,t) == k, 1)))
            Entropy(1) = Entropy(1) - length(find(L(:,t) == k))/n*M*log2(length(find(L(:,t) == k))/(n*M));
        else
            Entropy(1) = Entropy(1) - 0;
        end
    end
end

delta_Entropy = Entropy - Entropy_temp;



sc = [];
gbest_all = [];
n1 = n;
if n == 1
    gbest_all = gArchive;
else
    for k = 1:n
        sc(k) = score(L,n,k);
    end
    if all(sc == 0)
        gbest_all = gArchive;
    else
        q = find(sc == max(sc));
        sel = randperm(length(q),1);
        if length(q) == 1
            gbest_all = gArchive(q(1),:);
        else
            
            for j = 1:length(q)
                gbest_all =  [gbest_all;gArchive(q(j),:)];
            end
        end
    end
    
    
end




%% 迭代寻优
for i = 1:bh_option.maxgen
    
    
    
    %%
    
    
    for j = 1:bh_option.sizestar
        gbest = gbest_all(randperm(size(gbest_all,1),1),:);
        
        if rand < ELS
            k = ceil(N * rand);
            
            gbest(:,k+M) = gbest(:,k+M) + (max(gArchive(:,k+M)) - min(gArchive(:,k+M)))*randn;
            gbest(:,1:M) = fitnessfcn(gbest(:,M+1:N+M));
            
        end
        
        if  N == 1
            star(j,1) = star(j,1) + rand*(gbest(1,1 + M) - star(j,1));
            if star(j,1) > starmax(1)
                star(j,1) = starmax(1);
            end
            if star(j,1) < starmin(1)
                star(j,1) = starmin(1);
            end
        else
%             for k = 1:N
%                 star(j,k) = star(j,k) + rand*(gbest(1,k + M) - star(j,k));
%             end
            star(j,:) = (gbest(1,M+1:N+M) - star(j,:)).*rand(1,N) + star(j,:);
            
            for k = 1:N
                if star(j,k) > starmax(k)
                    star(j,k) = starmax(k);
                end
                if star(j,k) < starmin(k)
                    star(j,k) = starmin(k);
                end
            end
            
        end
        % Constellation updated
        
        % calculate the fitness
        solution_vector(j,:) = fitnessfcn(star(j,:));
        
        
        
        
        % adaptive mutation
        if rand < P
            k = ceil(N * rand);
            
            star(j,k) = starmin(k) + (starmax(k) - starmin(k))*rand;
            solution_vector(j,:) = fitnessfcn(star(j,:));
            
        end
        
        %         if all(solution_vector(j,:) <= gbest(1, 1:M))&&any(solution_vector(j,:) <gbest(1, 1:M))
        %             I = find(ismember(gbest_all, gbest,'rows'),1);
        %             gbest(1,M + 1:N + M) = star(j,:);
        %             gbest(1, 1:M) = solution_vector(j,:);
        % %             if isempty(I)
        % %                 gbest_all = [gbest_all;gbest];
        % %             else
        % %                 gbest_all(I,:) = gbest;
        % %             end
        %         end
        
        
        %判断是否新解进入外部文档
        count = 0;
        for t = 1:n
            %新解优于旧解，进入并替代
            if all(solution_vector(j,:) <= gArchive(t,1:M))&&any(solution_vector(j,:) < gArchive(t,1:M)) && count == 0 %占优旧解
                gArchive(t,:) = [solution_vector(j,:),star(j,:)];
                count = 1;
            elseif all(solution_vector(j,:) <= gArchive(t,1:M))&&any(solution_vector(j,:) < gArchive(t,1:M))&& count ~= 0
                
                gArchive(t,:) = NaN;
                
            end
        end
        I = all(isnan(gArchive), 2);
        gArchive(I,:) =[];
        
        n = size(gArchive,1);
       
        if count == 0  %Non-dominant
            if n < K
                 %互为非占优，并且文档未满，直接进入文档
                if all(sum(repmat(solution_vector(j,:),size(gArchive,1),1) >= gArchive(:,1:M),2) < M)
                    gArchive(n + 1,:) = [solution_vector(j,:),star(j,:)];
                end
            else
                if all(sum(repmat(solution_vector(j,:),size(gArchive,1),1) >= gArchive(:,1:M),2) < M)
                    %互为非占优，文档已满，判断密度
                    flag = 0;
                    den = zeros(1,n + 1);
                    dis = zeros(1,n + 1);
                    L_temp = [];
                    gArchive_temp = [];
                    gArchive_temp = [gArchive(:,1:M);solution_vector(j,:)];
                    for t = 1:M
                        L_temp(:,t) = ceil((n + 1)*(gArchive_temp(:,t) - min(gArchive_temp(:,t)))/(max(gArchive_temp(:,t)) -...
                            min(gArchive_temp(:,t))));
                        
                        
                    end
                    l_label = find(L_temp==0);
                    L_temp(l_label) = 1;
                    %计算密度
                    for t = 1:n+1
                        
                        den(t) = sum(sum(repmat(L_temp(t,:),size(L_temp,1),1) == L_temp));
                        
                        
                    end
                    
                    
                    
                    c = find(den == max(den));
                    %若密度小于最大密度，进入文档
                    select = c(randperm(length(c),1));
                    if den(end) < den(select)
                        gArchive(select,:) = [solution_vector(j,:),star(j,:)];
                        
                    end
                    
                    
                    
                    
                end
                
            end
        end
        n = size(gArchive,1);
        
        
        
        
        
        
    end
    
    n = size(gArchive,1);
    x1 = [];y1 = [];
    L = [];
    if n ~= 1
        
        for j = 1:M
            L(:,j) = ceil(n*(gArchive(:,j) - min(gArchive(:,j)))/(max(gArchive(:,j)) -...
                min(gArchive(:,j))));
            
            
            
            x1 = [x1 ; j - 0.7 + rand(size(L,1),1)*0.5];
            y1 = [y1 ; L(:,j) - 0.7 + 0.5* rand(size(L,1),1)];
            
            
        end
        l_label1 = find(L==0);
        L(l_label1) = 1;
        
    else
        L = [];
        
        for j = 1:M
            L(1,j) = ceil(n*(gArchive(1,j) - min(gArchive(1:n,j)))/(max(gArchive(1:n,j)) -...
                min(gArchive(1:n,j))));
            
            
            if gArchive(1,j) == min(gArchive(1:n,j))
                L(1,j) = 1;
            end
            
            x1 = [x1 ; j - 0.7 + rand(length(L),1)*0.5];
            y1 = [y1 ; L(:,j) - 0.7 + 0.5* rand(length(L),1)];
            
            
            
        end
        
        
        
        
    end
    
    set(h,'Xdata',x1,'Ydata',y1);
    drawnow
    % calculate the entropy
    
    for t = 1:M
        for k = 1:n
            if(~isempty(find(L(:,t) == k, 1)))
                Entropy(i) = Entropy(i) - length(find(L(:,t) == k))/n*M*log2(length(find(L(:,t) == k))/(n*M));
            else
                Entropy(i) = Entropy(i) - 0;
            end
        end
    end
    
    delta_Entropy(i) = Entropy(i) - Entropy_temp;
    
    % check the status
    delta_c = 2/n*log2(2);
    Entropy_temp = Entropy(i);
    if abs(delta_Entropy(i)) > delta_c || abs(n - n1) >0
        status = 'convergence';
    end
    
    if abs(delta_Entropy(i)) < delta_c && n == K && n1 == K
        if abs(delta_Entropy(i)) > delta_s
            status = 'diversity';
        end
    end
    
    if abs(delta_Entropy(i)) < delta_s && n == n1
        status = 'stagnation';
    end
    %
    
    den1 = [];
    if n~=1
        for t = 1:n
            
            den1(t) = sum(sum(repmat(L(t,:),size(L,1),1) == L));
            
            
        end
    else
        
        den1 = 3;
        
    end
    
    %状态收敛时
    if(strcmp(status,'convergence'))
        
        
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for t = 1:KT/2-1
                I = find(den1 == den2(t),1);
                gbest_all =  [gbest_all;gArchive(I,:)];
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2+1
                I = find(sc == scc(k),1);
                gbest_all = [gbest_all;gArchive(I,:)];
            end
            
        end
        
        
        
    end
    %状态多样时
    if(strcmp(status,'diversity'))
        ELS = ELS - step_x*abs(delta_Entropy(i));
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for t = 1:KT/2+1
                I = find(den1 == den2(t),1);
                gbest_all =  [gbest_all;gArchive(I,:)];
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2-1
                I = find(sc == scc(k),1);
                gbest_all = [gbest_all;gArchive(I,:)];
            end
            
        end
        

        
    end
    %状态停滞时
    if(strcmp(status,'stagnation'))
        ELS = ELS + 2*step_x*abs(1 + delta_Entropy(i));
        sc = [];
        scc = [];
        den2 = [];
        gbest_all = [];
        n1 = n;
        for k = 1:n
            sc(k) = score(L,n,k);
        end
        if length(den1) < KT
            
            gbest_all = gArchive;
            
        else
            den2 = sort(den1);
            for t = 1:KT/2
                I = find(den1 == den2(t),1);
                gbest_all =  [gbest_all;gArchive(I,:)];
            end
            scc = sort(sc,'descend');
            for k = 1:KT/2
                I = find(sc == scc(k),1);
                gbest_all = [gbest_all;gArchive(I,:)];
            end
            
        end

        
    end
    
    
    v = size(gbest_all,1);
    R = [];
    %计算黑洞半径
    if v == 1
        R = abs(gbest_all(v,1:M)./(sum(solution_vector)));
    else
        for j = 1:v
            R = [R;abs(gbest_all(j,1:M)./sum(solution_vector))];
        end
    end
    % 判断是否进入黑洞边界
    for j = 1:bh_option.sizestar
        %         S = min(abs(solution_vector(j,:) - gbest(1,1:M)));
        if v == 1;
           
            if any(sum(abs(repmat(solution_vector(j,:),size(gbest_all,1),1) - gbest_all(:,1:M)) > R,2) ==0) %进入吸收边界:吸收并重新生成
                if N == 1
                    star(j,1) = (starmax(1) - starmin(1)) * rand + starmin(1);
                else
                    
                    star(j,:) = (starmax - starmin).*rand(1,N) + starmin;
                end
                
            end
        else
            for  t = 1:v
                
                if any(sum(abs(repmat(solution_vector(j,:),size(gbest_all,1),1) - gbest_all(:,1:M)) > R,2) ==0) %进入吸收边界:吸收并重新生成
                    if N == 1
                        star(j,1) = (starmax(1) - starmin(1)) * rand + starmin(1);
                    else
                        
                        star(j,:) = (starmax - starmin).*rand(1,N) + starmin;
                        
                        
                    end
                    
                end
            end
            
        end
    end
    
end



set(h,'Xdata',x1,'Ydata',y1);
drawnow



fval = gArchive(:,1:M);
Ln = M+N;
if Ln ~= M + 1
    X = gArchive(:,M+1:Ln);
else
    X = gArchive(:,M+1);
end

APF = fval;
APS = X;

se = Entropy;
dse = delta_Entropy;

if ~isfield(bh_option,'gui')
figure(2)
set(gcf, 'position', [100 0 650 500]);
plot(1:bh_option.maxgen,Entropy,'-','LineWidth',0.5)
hold on
plot(1:bh_option.maxgen,delta_Entropy,'-r','LineWidth',0.5,'MarkerSize',1.5);
xlabel('t')
ylabel('H')
title('Entropy')
legend('Entropy','Delta entropy')
switch(M)
    case 2
        ploting2
    case 3
        ploting3
    otherwise
        
end
end