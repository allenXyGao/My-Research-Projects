function [ y ] = hdiff(data1,T,x,n,N,M)
%水位梯度函数
%输出项得到的是一个向量y（调和分析后多潮汐信号之一）   
%输入项 data1是数据集 T是历史持续时间长度 x是分类 n是数据集中列数  N是gamma函数参数  M是gamma函数参数
dd=data1(T+1:end,:);
c3x2=dd(x,n);
for i=1:T
    midc3x4(i,:)=data1(T+x-i,n)';   %过去流量  按照t-1,t-2,t-3,......,t-T的顺序排列
end
[c3x4m,c3x4n]=size(midc3x4);

for i=1:c3x4n
    for j=2:c3x4m-1
        if (midc3x4(j,i)>midc3x4(j-1,i)&midc3x4(j,i)>midc3x4(j+1,i))||(midc3x4(j,i)<midc3x4(j-1,i)&midc3x4(j,i)<midc3x4(j+1,i))  %判断峰和谷
            c3feng(j,i)=midc3x4(j,i);  %如果是则保留原数据
        else
            c3feng(j,i)=0;    %如果不是则设置为0， 这样方便后面计算处理
        end
    end
end

c3place=zeros(c3x4m,c3x4n);
for i=1:c3x4n
    c3placeN(i)=length(find(c3feng(:,i)~=0));
     c3place(1:c3placeN(i),i)=find(c3feng(:,i)~=0);
end
c3dimension=min(c3placeN);  %取min是为了统一数据矩阵的维数
c3place=c3place(1:c3dimension,:);


c3N2=N;         
c3M2=M;         
c3j2=1;
t2=1:c3dimension;
c3w2=(c3N2.^c3M2)./gamma(c3M2).*(c3j2.*t2).^(c3M2-1).*exp(-c3N2*c3j2.*t2);   %gamma函数
c3W2=sum(c3w2);
c3weight2=c3w2/c3W2; %归一化
c3weight2=c3weight2';

for i=1:c3x4n
    for j=1:size(c3place,1)-1
        c3x4MID(j,i)=(c3feng(c3place(j,i),i)-c3feng(c3place(j+1,i),i))/(c3place(j+1,i)-c3place(j,i));   %（峰-谷）/时间差  or  （谷-峰）/时间差
    end
end

for i=1:c3x4n
    c3add(i)=(c3x2(i)-c3feng(c3place(1,i),i))/c3place(1,i);  %这是（当前水位-最近的峰（or谷））/时间差
end
c3x4MID=[c3add;c3x4MID]; %合并

for i=1:c3x4n
    c3x4(i)=c3x4MID(:,i)'*c3weight2;  %乘上权重
end
c3x4=c3x4';
y=c3x4;
end

