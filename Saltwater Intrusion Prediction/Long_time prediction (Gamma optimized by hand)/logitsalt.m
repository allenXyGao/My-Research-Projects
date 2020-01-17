function [ y,a00,a01,a10,a11 ] = logitsalt( data1,x,V)
% y是成功率，
% a00是预测和实际都不可以取水，a01是实际不可以取水但预测错成可以取水，a10是实际可以取水的，但预测错了不可以取水，a11实际和预测都可以取水
% data1是数据集， T是历史持续时间，x是按照过去一周流量的分类   
% V=[c3N1,c3M1,c3d1N,c3d1M,c3d2N,c3d2M,c3d3N,c3d3M,c3d4N,c3d4M,c3d5N,c3d5M,c3d6N,c3d6M,c3d7N,c3d7M];
% V的范围如下：
% V(1): 0-1; 
% V(2): 1-50
% V(3): 0-1
% V(4)：1-40
% V(5),V(7),V(9),...V(15)同V(3)    V(6),V(8),...V(16)同V(4)
possibility=0.5;   %概率
threshold=0.5;%阈值
T=200; 
dd=data1(T+1:end,:);
c3x1=dd(x,1);
c3y=dd(x,3);

%将盐度值二值化（大于0.5时不能取水，小于或等于则能取水）
for i=1:length(c3y)
    if c3y(i)>threshold    %不能取水
        c3y(i)=0;
    elseif c3y(i)<=threshold   %能取水
        c3y(i)=1; 
    end
end

% 历史流量因素
for i=1:T
    midc3x3(i,:)=data1(T+x-i,1)';
end
j1=1;
t=1:T;
c3w=(V(1).^V(2))./gamma(V(2)).*(j1.*t).^(V(2)-1).*exp(-V(1)*j1.*t);  %gamma函数权重公式
c3W=sum(c3w);
c3weight1=c3w/c3W;  %归一化
c3weight1=c3weight1';
for j=1:length(c3x1)
    c3x3(j)=midc3x3(:,j)'*c3weight1;   %原始数据乘上等数据长度的权重
end
c3x3=c3x3';  %最终计算而得的历史流量因素（这个变量的名称是c3x3，不太顺眼是因为这个程序是基于无退火算法程序改来的，沿用了之前的变量名称，别的不顺眼的名称也是这个原因）

%下面是水位因素，我们把单一的水位数据经过拆分得到多组潮汐信号，需要哪几组潮汐信号是根据你的需求，下面的程序包含了7组（但主要因素可能只有其中三至四组）
%多潮汐信号第一列
n1=4;  % 潮汐信号的第一列位于原始数据的第四列，下同
% c3d1N=0.5;
% c3d1M=8;
c3d1N=V(3);
c3d1M=V(4);
c3x4d1=hdiff(data1,T,x,n1,c3d1N,c3d1M);
%多潮汐信号第二列
n2=5;
% c3d2N=0.5;
% c3d2M=5;
c3d2N=V(5);
c3d2M=V(6);
c3x4d2=hdiff(data1,T,x,n2,c3d2N ,c3d2M);
%多潮汐信号第三列
 n3=6;
% c3d3N=0.5;
% c3d3M=5;
c3d3N=V(7);
c3d3M=V(8);
c3x4d3=hdiff(data1,T,x,n3,c3d3N,c3d3M);
%多潮汐信号第四列
 n4=7;
% c3d4N=0.5;
% c3d4M=5;
c3d4N=V(9);
c3d4M=V(10);
c3x4d4=hdiff(data1,T,x,n4,c3d4N ,c3d4M );
%多潮汐信号第五列
n5=8;
c3d5N=V(11);
c3d5M=V(12);
c3x4d5=hdiff(data1,T,x,n5,c3d5N ,c3d5M );
%多潮汐信号第六列
n6=9;
c3d6N=V(13);
c3d6M=V(14);
c3x4d6=hdiff(data1,T,x,n6,c3d6N ,c3d6M );
%多潮汐信号第七列
n7=10;
c3d7N=V(15);
c3d7M=V(16);
c3x4d7=hdiff(data1,T,x,n7,c3d7N ,c3d7M);

c3X=[c3x3,c3x4d1,c3x4d2,c3x4d3,c3x4d4,c3x4d5,c3x4d6,c3x4d7];

c3theta = glmfit(c3X, [c3y ones(length(c3y),1)], 'binomial', 'link', 'logit');
 c3p=1./(1+exp(-[ones(length(c3y),1),c3X]*c3theta));
 
 c3P=c3p;
 
 for i=1:length(c3P)      %此处是依据概率来分2值，如果概率>=0.5,则视为1；否则视为0
    if c3P(i)>=possibility
       c3P(i)=1;
    elseif c3P(i)<possibility
       c3P(i)=0;
    end
end
c3pp=c3y;
c3delta=c3P-c3pp;
[c3M,c3N]=find(c3delta);
c3successfulrate=(length(c3P)-length(c3M))/length(c3P);
y=c3successfulrate;

a00=length(find(c3y==0&c3P==0));  %见开头注释
a01=length(find(c3y==0&c3P==1));
a10=length(find(c3y==1&c3P==0));
a11=length(find(c3y==1&c3P==1));


end

