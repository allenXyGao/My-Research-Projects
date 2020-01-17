%台风来临时的logistic regression
clc
clear
%分三段流量，分别建模，
possibility=0.5;   %概率
threshold=0.5;%阈值
c1alpha=0.2;
T=200;%历史序列长度
data=xlsread('1.xlsx');
data1=data((300:1000),:);   %训练数据  (3692:12465)    (1:8760)
data1(:,1)=data1(:,1)+data1(:,2);
data1(:,2)=data1(:,3);
data1(:,3)=data1(:,4);
data1(:,4)=[];
dd=data1(T+1:end,:);
% figure(1)
% plot(dd(:,1))
% figure(2)
% plot(dd(:,2))
x1=1:length(dd);
c1x1=log(dd(x1,1));  %t时刻
%c1x1=c1x1*mg1;
c1x2=dd(x1,2);    
%c1x2=c1x2*mg2;
c1y=dd(x1,3);
for i=1:length(c1y)
    if c1y(i)>threshold    %不能取水
        c1y(i)=0;
    elseif c1y(i)<=threshold   %能取水
        c1y(i)=1; 
    end
end

for i=1:T
    midc1x3(i,:)=data1(T+x1-i,1)';
end

%权重  
c1N1=0.3;    %第一套         
c1M1=40;    %第一套
c1j1=1;    
t=1:T;
c1w=(c1N1^c1M1)/gamma(c1M1).*(c1j1.*t).^(c1M1-1).*exp(-c1N1*c1j1.*t);  %gamma函数
c1W=sum(c1w);
c1weight1=c1w/c1W;
c1weight1=c1weight1';
for j=1:length(c1x1)
         c1x3(j)=midc1x3(:,j)'*c1weight1;
end
c1x3=log(c1x3');  %历史流量因素


%下面是历史水位梯度因素
for i=1:T
    midc1x4(i,:)=data1(T+x1-i,2)';   %过去水位  按照t-1,t-2,t-3,......,t-T的顺序排列
end
[c1x4m,c1x4n]=size(midc1x4);

for i=1:c1x4n
    for j=2:c1x4m-1
        if (midc1x4(j,i)>midc1x4(j-1,i)&midc1x4(j,i)>midc1x4(j+1,i))||(midc1x4(j,i)<midc1x4(j-1,i)&midc1x4(j,i)<midc1x4(j+1,i))  %判断峰与谷
            c1feng(j,i)=midc1x4(j,i);
        else
            c1feng(j,i)=0;
        end
    end
end

 %确定峰和谷的具体位置
c1place=zeros(c1x4m,c1x4n);
for i=1:c1x4n
    c1placeN(i)=length(find(c1feng(:,i)~=0));  
     c1place(1:c1placeN(i),i)=find(c1feng(:,i)~=0);
end
c1dimension=min(c1placeN);
c1place=c1place(1:c1dimension,:);


c1N2=0.11;    %第二套           % 4个系数分别为 0.3 40  0.1  1.7 正确率 0.9053     0.11  1.9      
c1M2=1.9;    %第二套
c1j2=1;
t2=1:c1dimension;
c1w2=(c1N2^c1M2)/gamma(c1M2).*(c1j2.*t2).^(c1M2-1).*exp(-c1N2*c1j2.*t2);
c1W2=sum(c1w2);
c1weight2=c1w2/c1W2;
c1weight2=c1weight2';

for i=1:c1x4n
    for j=1:size(c1place,1)-1
        c1x4MID(j,i)=(c1feng(c1place(j,i),i)-c1feng(c1place(j+1,i),i))/(c1place(j+1,i)-c1place(j,i));  %（峰-谷）/时间差  or  （谷-峰）/时间差
    end
end

for i=1:c1x4n
    c1add(i)=(c1x2(i)-c1feng(c1place(1,i),i))/c1place(1,i);     %这是（当前水位-最近的峰（or谷））/时间差
end
c1x4MID=[c1add;c1x4MID]; %合并

for i=1:c1x4n
    c1x4(i)=c1x4MID(:,i)'*c1weight2;
end
c1x4=c1x4';

%这是简单的梯度版本，单纯地h（t+1）-h(t),直接忽视这段，保留下仅仅是为了能有些启发。
% c1weight2=zeros(T,1);
% for i=1:length(c1weight2)
%     c1weight2(i)=c1alpha*(1-c1alpha)^(i-1);
% end
% c1weight2=c1weight2/sum(c1weight2);
% for k=1:length(c1x1)    % 水位梯度的加权平均
%     A(:,k)=diff(midc1x4(:,k));   %梯度
%      c1x4(k)=A(:,k)'*c1weight2(1:end-1);
% end
% c1x4=c1x4';
%c1x4=c1x4*mg4;


c1X=[c1x3,c1x4];
c1theta = glmfit(c1X, [c1y ones(length(c1y),1)], 'binomial', 'link', 'logit');
 c1p=1./(1+exp(-[ones(length(c1y),1),c1X]*c1theta));
 
 c1P=c1p;  %为了不破坏原结果c1p，建立c1P可用来作后续分析
 
 for i=1:length(c1P)      %此处是依据概率来分2值，如果概率>=0.5,则视为1；否则视为0
    if c1P(i)>=possibility
       c1P(i)=1;
    elseif c1P(i)<possibility
       c1P(i)=0;
    end
end
c1pp=c1y;
c1delta=c1P-c1pp;
[c1M,c1N]=find(c1delta);
c1successfulrate=(length(c1P)-length(c1M))/length(c1P)

a00=length(find(c1y==0&c1P==0)); 
a01=length(find(c1y==0&c1P==1));
a10=length(find(c1y==1&c1P==0));
a11=length(find(c1y==1&c1P==1));
%a01是实际不可以取水但预测错成可以取水，a10是实际可以取水的，但预测错了不可以取水，a11实际和预测都可以取水
%detail是综合
detail=[a00,a01;
  a10,a11  ]

% % %找问题
% faultc1=[c1x2(c1M),c1x3(c1M),c1x4(c1M),c1y(c1M),c1p(c1M)];

 
figure(1)
subplot(2,1,1)
plot(1:T,c1weight1,'b-','LineWidth',1.5)
legend('the weight of flow')
title('below 2500')
subplot(2,1,2)
plot(1:c1dimension,c1weight2,'r-','LineWidth',1.5)
legend('the weight of sea level')
title('below 2500')

