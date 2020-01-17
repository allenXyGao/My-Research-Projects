clc
clear
%主程序（模拟退火）
%利用两个函数的模拟退火logistic程序    这两个函数分别是hdiff.m(水位梯度函数)，和logitsalt.m（正确率函数）
%分三段流量，分别建模，
possibility=0.5;   %概率
threshold=0.5;%阈值
T=200;%历史序列长度
flow1=3500;     %  70%流量
flow2=2500;     %  90%流量
data=xlsread('2011to2013.xls');
%这个数据是补充了多潮汐信号的数据 从第五列开始
data1=data(8786-T:13129,:);   %训练数据  (3692:12465)    (1:8760)  (1000:11713)  
%8586:13129-----保证训练2012/9/1-2013/2/28  11514:end ----2013年
%8586:end---2012 9 1到2013年底

% data=xlsread('newtrain.xlsx');
% data1=data(600:end,:);      %训练数据  (3692:12465)    (1:8760)   1000:11713
data1(:,1)=data1(:,1)+data1(:,2);
data1(:,2)=data1(:,3);
data1(:,3)=data1(:,4);
% data1(:,4)=[];
data1(:,4)=data1(:,5);
data1(:,5)=data1(:,6);
data1(:,6)=data1(:,7);
data1(:,7)=data1(:,8);
data1(:,8)=data1(:,9);
data1(:,9)=data1(:,10);
data1(:,10)=data1(:,11);
data1(:,11)=[];

q1=data1(:,1); %逐时流量
X1=q1(T+1:end);

%计算过去一周的平均流量，再分类 ： 2500以下  ；2500~3500  ；3500以上
for i=1:length(X1)
    average(i)=sum(q1(T-168+i:i+T-1))/168;  %过去一周流量的平均
end

x1=find(average<=flow2);          %2500以下
x2=find(average<=flow1&average>flow2);      %2500~3500
x3=find(average>flow1);        %3500以上

dd=data1(T+1:end,:);

%第三类
% V=[0.5,30,0.1,20,0.1,20,0.1,20,0.1,20,0.1,20,0.1,20,0.1,20];
% [ y3,a300,a301,a310,a311 ]= logitsalt( data1,x3,V);


lb=[0,1, 0,1, 0,1, 0,1,0, 1, 0,1, 0,1, 0,1]; % 参数取值下界
ub=[1,60,1,40,1,40,1,40,1,40,1,40,1,40,1,40]; % 参数取值上界
% 冷却表参数
MarkovLength=200; % 马可夫链长度
DecayScale=0.9 ; % 衰减参数
StepFactor=0.2; % Metropolis步长因子
sp=0.2;  %搜索步长因子
Temperature0=90; % 初始温度
Temperatureend=0.1; % 最终温度
Boltzmann_con=1; % Boltzmann常数
% AcceptPoints=0.0; % Metropolis过程中总接受点
% 随机初始化参数
range=ub-lb;
Par_ini=rand(size(lb)).*range+lb
Par_cur=Par_ini; % 用Par_cur表示当前解
Par_best_cur=Par_cur; % 用Par_best_cur表示当前最优解
Par_best=rand(size(lb)).*range+lb; % 用Par_best表示冷却中的最好解
% 每迭代一次退火(降温)一次，直到满足迭代条件为止
t=Temperature0;
itr_num=0; % 记录迭代次数
k=1;
while t>Temperatureend
    itr_num=itr_num+1;
    AcceptPoints=0.0; %每个Metropolis过程中接受点
    t=DecayScale*t; % 温度更新（降温)
    for i=1:MarkovLength
        % 在此当前参数点附近随机选下一点
        p=0;
        while p==0
            %Par_new=Par_cur+StepFactor.*range.*(rand(size(lb))-0.5);
            u=rand(size(lb));
            eta=sign(u-0.5).*t.*((1+1./t).^abs(2*u-1)-1);
            delta=sp*eta;
            Par_new=Par_cur+delta;
            % 防止越界
            if sum(Par_new>ub)+sum(Par_new<lb)==0
                p=1;
            end
        end
        % 检验当前解是否为全局最优解  logitsalt( data1,x3,V)  
        if (logitsalt(data1,x3,Par_new)<logitsalt(data1,x3,Par_best))
            % 保留上一个最优解
            Par_best_cur=Par_best;
            % 此为新的最优解
            Par_best=Par_new;
        end
        % Metropolis过程
        if (logitsalt(data1,x3,Par_new)-logitsalt(data1,x3,Par_cur)>0)
            % 接受新解
            Par_cur=Par_new;
            AcceptPoints=AcceptPoints+1;
        else
            changer=-1*(logitsalt(data1,x3,Par_new)-logitsalt(data1,x3,Par_cur))/(Boltzmann_con*t);
            p1=exp(changer);
            if p1>rand
                Par_cur=Par_new;
                AcceptPoints=AcceptPoints+1;
            end
        end
    end
    accp=AcceptPoints/MarkovLength;
    sp=I(accp)*sp;
    k=k+1;
end
%% 结果显示
disp(['最小值点:',num2str(Par_best)]);
Objval_best= logitsalt(data1,x3,Par_best);
disp(['最小值为:',num2str(Objval_best)]);
