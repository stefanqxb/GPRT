function [ hyp2, fitgen, timewin, counter_time, counter_pre ] = PSO_hypa( hypt, func, test_targ ,varargin )


pso_option = struct('c1',5,'c2',4,'maxgen', 20,'sizepop', 7, ...
                    'k',1.5,'wV',1,'wP',1, ...
                     'popcov1max',10^(-2),'popcov1min',0,'popcov2max',10^(-2),'popcov2min', 0 ,'popmeanmax',10^(-2),'popmeanmin', -10);

% c1:��ʼΪ1.5,pso�����ֲ���������
% c2:��ʼΪ1.7,pso����ȫ����������
% maxgen:��ʼΪ200,����������
% sizepop:��ʼΪ20,��Ⱥ�������
% k:��ʼΪ0.6(k belongs to [0.1,1.0]),���ʺ�x�Ĺ�ϵ(V = kX)
% wV:��ʼΪ1(wV best belongs to [0.8,1.2]),���ʸ��¹�ʽ���ٶ�ǰ��ĵ���ϵ��
% wP:��ʼΪ1,��Ⱥ���¹�ʽ���ٶ�ǰ��ĵ���ϵ��
% v:��ʼΪ3,SVM Cross Validation����
% popcmax:��ʼΪ100,SVM ����c�ı仯�����ֵ.
% popcmin:��ʼΪ0.1,SVM ����c�ı仯����Сֵ.
% popgmax:��ʼΪ1000,SVM ����g�ı仯�����ֵ.
% popgmin:��ʼΪ0.01,SVM ����c�ı仯����Сֵ.

Vc1max = pso_option.k*pso_option.popcov1max;
Vc1min = -Vc1max ;
Vc2max = pso_option.k*pso_option.popcov2max;
Vc2min = -Vc2max ;
Vmeanmax = pso_option.k*pso_option.popmeanmax;
Vmeanmin = -Vmeanmax ;



Z = hypt;
eps = 2^(-4);

time = zeros(pso_option.sizepop,1);
counter_time = zeros(pso_option.sizepop,pso_option.maxgen);
counter_pre = zeros(pso_option.sizepop,pso_option.maxgen);

%% ������ʼ���Ӻ��ٶ�
for i=1:pso_option.sizepop
   
    % ���������Ⱥ���ٶ�
    hyp(i,1) = (pso_option.popcov1max-pso_option.popcov1min)*rand+pso_option.popcov1min;  
    hyp(i,2) = (pso_option.popcov2max-pso_option.popcov2min)*rand+pso_option.popcov2min;
    hyp(i,3) = (pso_option.popmeanmax-pso_option.popmeanmin)*rand+pso_option.popmeanmin;
    V(i,1) = Vc1max*rands(1,1);  
    V(i,2) = Vc2max*rands(1,1);
    V(i,3) = Vmeanmax*rands(1,1);
   
    % �����ʼ��Ӧ��    %%%% GPML model
    [f0,~] = feval(func, rewrap(Z,hyp(i,:)'), varargin{:});
    gamma = corrcoef(f0,test_targ);
    fitness(i) = gamma(1,2);
    diff_mat = abs(f0 - test_targ);
    max_t = max(test_targ);
    min_t = min(test_targ);
    step = 100;
    thist = hist(diff_mat,step);
    time(i) = time_95_diff(thist,max_t,min_t,step);
    time(i) = time(i) / max_t;
end



[global_fitness,bestindex]=max(fitness); % ȫ�ּ�ֵ
local_fitness=fitness;   % ���弫ֵ��ʼ��

global_time = time(bestindex);
local_time = time;

global_x=hyp(bestindex,:);   % ȫ�ּ�ֵ��
local_x=hyp;    % ���弫ֵ���ʼ��


% ÿһ����Ⱥ��ƽ����Ӧ��
fitgen = zeros(pso_option.maxgen,1); 
timewin = zeros(pso_option.maxgen,1);

%% ����Ѱ��
for i=1:pso_option.maxgen
    i   
    for j=1:pso_option.sizepop
       
        %�ٶȸ���        
        V(j,:) = pso_option.wV*V(j,:) + pso_option.c1*rand*(local_x(j,:) - hyp(j,:)) + pso_option.c2*rand*(global_x - hyp(j,:)); 
        if V(j,1) > Vc1max
            V(j,1) = Vc1max;
        end
        if V(j,1) < Vc1min
            V(j,1) = Vc1min;
        end
        if V(j,2) > Vc2max
            V(j,2) = Vc2max;
        end
        if V(j,2) < Vc2min
            V(j,2) = Vc2min;
        end
        if V(j,3) > Vmeanmax
            V(j,3) = Vmeanmax;
        end
        if V(j,3) < Vmeanmin
            V(j,3) = Vmeanmin;
        end
        
        %��Ⱥ����
        hyp(j,:)=hyp(j,:) + pso_option.wP*V(j,:);
        if hyp(j,1) > pso_option.popcov1max
            hyp(j,1) = pso_option.popcov1max;
        end
        if hyp(j,1) < pso_option.popcov1min
            hyp(j,1) = pso_option.popcov1min;
        end
        if hyp(j,2) > pso_option.popcov2max
            hyp(j,2) = pso_option.popcov2max;
        end
        if hyp(j,2) < pso_option.popcov2min
            hyp(j,2) = pso_option.popcov2min;
        end
        if hyp(j,3) > pso_option.popmeanmax
            hyp(j,3) = pso_option.popmeanmax;
        end
        if hyp(j,3) < pso_option.popmeanmin
            hyp(j,3) = pso_option.popmeanmin;
        end
        
        % ����Ӧ���ӱ���
       % if rand>0.5
       %     k=ceil(2*rand);
       %     if k == 1
       %         hyp(j,k) = (20-1)*rand+1;
       %     end
       %     if k == 2
       %         hyp(j,k) = (pso_option.popcov2max-pso_option.popcov2min)*rand + pso_option.popcov2min;
       %     end            
       % end
       
        %��Ӧ��ֵ
        [f0,~] = feval(func, rewrap(Z,hyp(j,:)'), varargin{:}); 
        gamma = corrcoef(f0,test_targ);
        fitness(j) = gamma(1,2);
        counter_pre(j,i) = fitness(j);
        diff_mat = abs(f0 - test_targ);
        max_t = max(test_targ);
        min_t = min(test_targ);
        step = 100;
        thist = hist(diff_mat,step);
        time(j) = time_95_diff(thist,max_t,min_t,step);
        time(j) = time(j) / max_t;
        counter_time(j,i) = time(j); 
        
      
        %�������Ÿ���
        if fitness(j) > local_fitness(j)
            local_x(j,:) = hyp(j,:);
            local_fitness(j) = fitness(j);
            local_time(j) = time(j);
        end

        if abs(fitness(j) - local_fitness(j)) <= eps && time(j,1) <= local_time(j,1)
            local_x(j,:) = hyp(j,:);
            local_fitness(j) = fitness(j);
            local_time(j) = time(j);
        end        
        
        %Ⱥ�����Ÿ���
        if fitness(j) > global_fitness
            global_x = hyp(j,:);
            global_fitness = fitness(j);
            global_time = time(j);
        end

        if abs(fitness(j) - global_fitness) <= eps  && time(j,1) <= global_time(1)
            global_x = hyp(j,:);
            global_fitness = fitness(j);
            global_time = time(j);
        end
        
    end
  
    fitgen(i) = global_fitness;
    timewin(i) = global_time;
end


hypbest = global_x;


hyp2 = rewrap(Z,hypbest');

end

function v = unwrap(s)
% Extract the numerical values from "s" into the column vector "v". The
% variable "s" can be of any type, including struct and cell array.
% Non-numerical elements are ignored. See also the reverse rewrap.m. 
v = [];   
if isnumeric(s)
  v = s(:);                        % numeric values are recast to column vector
elseif isstruct(s)
  v = unwrap(struct2cell(orderfields(s))); % alphabetize, conv to cell, recurse
elseif iscell(s)
  for i = 1:numel(s)             % cell array elements are handled sequentially
    v = [v; unwrap(s{i})];
  end
end % other types are ignored

end

function [s,v] = rewrap(s, v)
% Map the numerical elements in the vector "v" onto the variables "s" which can
% be of any type. The number of numerical elements must match; on exit "v"
% should be empty. Non-numerical entries are just copied. See also unwrap.m.
if isnumeric(s)
  if numel(v) < numel(s)
    error('The vector for conversion contains too few elements')
  end
  s = reshape(v(1:numel(s)), size(s));            % numeric values are reshaped
  v = v(numel(s)+1:end);                        % remaining arguments passed on
elseif isstruct(s) 
  [s p] = orderfields(s); p(p) = 1:numel(p);      % alphabetize, store ordering
  [t v] = rewrap(struct2cell(s), v);                 % convert to cell, recurse
  s = orderfields(cell2struct(t,fieldnames(s),1),p);  % conv to struct, reorder
elseif iscell(s)
  for i = 1:numel(s)             % cell array elements are handled sequentially 
    [s{i} v] = rewrap(s{i}, v);
  end
end
end
