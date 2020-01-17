clc
clear
close all
%分三段流量，分别建模，
possibility=0.5;   %概率
threshold=0.5;%阈值
T=200;%历史序列长度
flow1=3500;     %  70%流量
flow2=2500;     %  90%流量
data=xlsread('2011to2013.xls');
data1=data(8786-T:13129,:);   %训练数据  (3692:12465)    (1:8760)  (1000:11713)   %c2time 代表未来的20个小时
%8586:13129-----保证训练2012/9/1-2013/2/28  11514:end ----2013年
%8586:end---2012 9 1到2013年底

% data=xlsread('newtrain.xlsx');
% data1=data(600:end,:);      %训练数据  (3692:12465)    (1:8760)   1000:11713
data1(:,1)=data1(:,1)+data1(:,2);
data1(:,2)=data1(:,3);
data1(:,3)=data1(:,4);
data1(:,4)=[];


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

if isempty(x1)    
    c1successfulrate=0;
elseif ~isempty(x1)
% %第一类 2500以下
c1x1=log(dd(x1,1));
c1x2=dd(x1,2);
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
c1N1=0.1;    %第一套         
c1M1=30;    %第一套
c1j1=1;    
t=1:T;
c1w=(c1N1^c1M1)/gamma(c1M1).*(c1j1.*t).^(c1M1-1).*exp(-c1N1*c1j1.*t);
c1W=sum(c1w);
c1weight1=c1w/c1W;
c1weight1=c1weight1';

for j=1:length(c1x1)
         c1x3(j)=midc1x3(:,j)'*c1weight1;
end
%c1x3=log(c1x3');
c1x3=c1x3';


for i=1:T
    midc1x4(i,:)=data1(T+x1+1-i,2)';   %过去水位  按照t,t-1,t-2,t-3,......,t-T的顺序排列
end
[c1x4m,c1x4n]=size(midc1x4);

for i=1:c1x4n
    for j=2:c1x4m-1
        if (midc1x4(j,i)>midc1x4(j-1,i)&midc1x4(j,i)>midc1x4(j+1,i))||(midc1x4(j,i)<midc1x4(j-1,i)&midc1x4(j,i)<midc1x4(j+1,i))
            c1feng(j,i)=midc1x4(j,i);
        else
            c1feng(j,i)=0;
        end
    end
end

c1place=zeros(c1x4m,c1x4n);
for i=1:c1x4n
    c1placeN(i)=length(find(c1feng(:,i)~=0));
     c1place(1:c1placeN(i),i)=find(c1feng(:,i)~=0);
end
c1dimension=min(c1placeN);
c1place=c1place(1:c1dimension,:);

c1N2=0.11;    %第二套             
c1M2=1.8;    %第二套            
c1j2=1;
t2=1:c1dimension;
c1w2=(c1N2^c1M2)/gamma(c1M2).*(c1j2.*t2).^(c1M2-1).*exp(-c1N2*c1j2.*t2);
c1W2=sum(c1w2);
c1weight2=c1w2/c1W2;
c1weight2=c1weight2';

for i=1:c1x4n
    for j=1:size(c1place,1)-1
        c1x4MID(j,i)=(c1feng(c1place(j,i),i)-c1feng(c1place(j+1,i),i))/(c1place(j+1,i)-c1place(j,i));
    end
end

for i=1:c1x4n
      c1add(i)=(c1x2(i)-c1feng(c1place(1,i),i))/(c1place(1,i)-1);
end
c1x4MID=[c1add;c1x4MID];

for i=1:c1x4n
    c1x4(i)=c1x4MID(:,i)'*c1weight2;
end
c1x4=c1x4';

c1X=[c1x3,c1x4];

c1theta = glmfit(c1X, [c1y ones(length(c1y),1)], 'binomial', 'link', 'logit');
 c1p=1./(1+exp(-[ones(length(c1y),1),c1X]*c1theta));
 
 c1P=c1p;
 
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
c1successfulrate=(length(c1P)-length(c1M))/length(c1P);
end


%第二类 3500以下+2500以上
c2x1=log(dd(x2,1));
c2x2=dd(x2,2);
c2y=dd(x2,3);
for i=1:length(c2y)
    if c2y(i)>threshold    %不能取水
        c2y(i)=0;
    elseif c2y(i)<=threshold   %能取水
        c2y(i)=1; 
    end
end

for i=1:T
    midc2x3(i,:)=data1(T+x2-i,1)';
end
c2N1=0.294;    %第一套
c2M1=28;      %第一套
j2=1;
t=1:T;
c2w=(c2N1^c2M1)/gamma(c2M1).*(j2.*t).^(c2M1-1).*exp(-c2N1*j2.*t);
c2W=sum(c2w);
c2weight1=c2w/c2W;
c2weight1=c2weight1';

for j=1:length(c2x1)
         c2x3(j)=midc2x3(:,j)'*c2weight1;
end
c2x3=c2x3';

for i=1:T
    midc2x4(i,:)=data1(T+x2+1-i,2)';   %过去流量  按照t-1,t-2,t-3,......,t-T的顺序排列
end
[c2x4m,c2x4n]=size(midc2x4);

for i=1:c2x4n
    for j=2:c2x4m-1
        if (midc2x4(j,i)>midc2x4(j-1,i)&midc2x4(j,i)>midc2x4(j+1,i))||(midc2x4(j,i)<midc2x4(j-1,i)&midc2x4(j,i)<midc2x4(j+1,i))
            c2feng(j,i)=midc2x4(j,i);
        else
            c2feng(j,i)=0;
        end
    end
end

c2place=zeros(c2x4m,c2x4n);
for i=1:c2x4n
    c2placeN(i)=length(find(c2feng(:,i)~=0));
     c2place(1:c2placeN(i),i)=find(c2feng(:,i)~=0);
end
c2dimension=min(c2placeN);
c2place=c2place(1:c2dimension,:);

c2N2=0.54;      %第二套     
c2M2=6.78;     %第二套
c2j2=1;
t2=1:c2dimension;
c2w2=(c2N2^c2M2)/gamma(c2M2).*(c2j2.*t2).^(c2M2-1).*exp(-c2N2*c2j2.*t2);
c2W2=sum(c2w2);
c2weight2=c2w2/c2W2;
c2weight2=c2weight2';

for i=1:c2x4n
    for j=1:size(c2place,1)-1
        c2x4MID(j,i)=(c2feng(c2place(j,i),i)-c2feng(c2place(j+1,i),i))/(c2place(j+1,i)-c2place(j,i));
    end
end

for i=1:c2x4n
      c2add(i)=(c2x2(i)-c2feng(c2place(1,i),i))/(c2place(1,i)-1);
end
c2x4MID=[c2add;c2x4MID];

for i=1:c2x4n
    c2x4(i)=c2x4MID(:,i)'*c2weight2;
end
c2x4=c2x4';

c2X=[c2x3,c2x4];

c2theta = glmfit(c2X, [c2y ones(length(c2y),1)], 'binomial', 'link', 'logit');
 c2p=1./(1+exp(-[ones(length(c2y),1),c2X]*c2theta));
 
 c2P=c2p;
 
 for i=1:length(c2P)      %此处是依据概率来分2值，如果概率>=0.5,则视为1；否则视为0
    if c2P(i)>=possibility
       c2P(i)=1;
    elseif c2P(i)<possibility
       c2P(i)=0;
    end
end
c2pp=c2y;
c2delta=c2P-c2pp;
[c2M,c2N]=find(c2delta);
c2successfulrate=(length(c2P)-length(c2M))/length(c2P);

% 
%第三类  3500以上
c3x1=log(dd(x3,1));
c3x2=dd(x3,2);
c3y=dd(x3,3);
for i=1:length(c3y)
    if c3y(i)>threshold    %不能取水
        c3y(i)=0;
    elseif c3y(i)<=threshold   %能取水
        c3y(i)=1; 
    end
end

for i=1:T
    midc3x3(i,:)=data1(T+x3-i,1)';
end

%  0.45 15  0.29 2.5 0.9475   (2012/9/1-2013/2/28)
c3N1=0.45;     %第一套         
c3M1=14.5;    %第二套
j1=1;
t=1:T;
c3w=(c3N1^c3M1)/gamma(c3M1).*(j1.*t).^(c3M1-1).*exp(-c3N1*j1.*t);
c3W=sum(c3w);
c3weight1=c3w/c3W;
c3weight1=c3weight1';


for j=1:length(c3x1)
         c3x3(j)=midc3x3(:,j)'*c3weight1;
end
c3x3=c3x3';

for i=1:T
    midc3x4(i,:)=data1(T+x3+1-i,2)';   %历史水位  按照t-1,t-2,t-3,......,t-T的顺序排列
end
[c3x4m,c3x4n]=size(midc3x4);

for i=1:c3x4n
    for j=2:c3x4m-1
        if (midc3x4(j,i)>midc3x4(j-1,i)&midc3x4(j,i)>midc3x4(j+1,i))||(midc3x4(j,i)<midc3x4(j-1,i)&midc3x4(j,i)<midc3x4(j+1,i))
            c3feng(j,i)=midc3x4(j,i);
        else
            c3feng(j,i)=0;
        end
    end
end

c3place=zeros(c3x4m,c3x4n);
for i=1:c3x4n
    c3placeN(i)=length(find(c3feng(:,i)~=0));
     c3place(1:c3placeN(i),i)=find(c3feng(:,i)~=0);
end
c3dimension=min(c3placeN);
c3place=c3place(1:c3dimension,:);

for i=1:c3x4n
    for j=1:size(c3place,1)-1
        c3x4MID(j,i)=(c3feng(c3place(j,i),i)-c3feng(c3place(j+1,i),i))/(c3place(j+1,i)-c3place(j,i));
    end
end

for i=1:c3x4n
    c3add(i)=(c3x2(i)-c3feng(c3place(1,i),i))/(c3place(1,i)-1);
end
c3x4MID=[c3add;c3x4MID];

c3N2=0.304;         %第一套  0.29  2.5
c3M2=2.45;         %第二套
c3j2=1;
t2=1:c3dimension;
c3w2=(c3N2^c3M2)/gamma(c3M2).*(c3j2.*t2).^(c3M2-1).*exp(-c3N2*c3j2.*t2);
c3W2=sum(c3w2);
c3weight2=c3w2/c3W2;
c3weight2=c3weight2';

for i=1:c3x4n
    c3x4(i)=c3x4MID(:,i)'*c3weight2;
end
c3x4=c3x4';

c3X=[c3x3,c3x4];

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


% 
% %找问题
%faultc1=[c1x2(c1M),c1x3(c1M),c1x4(c1M),c1y(c1M),c1p(c1M)];
faultc2=[c2x2(c2M),c2x3(c2M),c2x4(c2M),c2y(c2M),c2p(c2M)];
faultc3=[c3x2(c3M),c3x3(c3M),c3x4(c3M),c3y(c3M),c3p(c3M)];
% figure(1)
% plot(exp(c1x3),c1p,'+',exp(c1x3),c1y,'.','markersize',1)
% legend('概率','实际')
% figure(2)
% plot(exp(c2x3),c2p,'+',exp(c2x3),c2y,'.','markersize',1)
% legend('概率','实际')
% figure(3)
% plot(exp(c3x3),c3p,'+',exp(c3x3),c3y,'.','markersize',1)
% legend('概率','实际')

% c3pick=[c3X,c3y];
% xlswrite('c3pick.xls', c3pick)


%画10张图找错误
% c3leveln=sort(c3x3(c3M));
% c3leveln=c3leveln(end-10);
% c3Mmax10=find(c3x3(c3M)>c3leveln);
% c3fault10=c3M(c3Mmax10);     %3500以上错误的序号中历史流量最大的10个
% c3faultflow=c3x3(c3fault10);
% 
% for i=1:4*24+1  %错误的数据的水位
%     c3faultlevel(:,i)=c3x2(c3fault10+1-i);
% end
% 
% for i=1:length(c3fault10)
%     c3faultlevel(i,:)=flipud(c3faultlevel(i,:)');
% end
% 
% c3right=zeros(length(c3Mmax10),20);
% for i=1:length(c3Mmax10)
%     c3length(i)=length(find(abs(c3x3-c3faultflow(i))<=0.001));
%     c3right(i,1:c3length(i))=find(abs(c3x3-c3faultflow(i))<=0.001);
%     for j=1:length(c3right(i,:))
%         if c3right(i,j)==c3fault10(i)
%             c3right(i,j)=0;
%         end
%     end
% end
% c3rightpick=c3right(:,2);  %找到非错误数据中的项
% c3rightflow=c3x3(c3rightpick);
% 
% for i=1:4*24+1  %   正确的水位
%     c3rightlevel(:,i)=c3x2(c3rightpick+1-i);    %从过去到现在
% end
% 
% for i=1:length(c3fault10)
%     c3rightlevel(i,:)=flipud(c3rightlevel(i,:)');
% end
% 
% for j=1:length(c3fault10)
%     figure(j) 
%     str=['fault=',num2str(c3faultflow(j)),'right=',num2str(c3rightflow(j))];    %画10张图
%      plot(1:4*24+1,c3faultlevel(j,:),1:4*24+1,c3rightlevel(j,:))
%     title(str);
%    legend('错误','正确')
% end

totalsuccess=length(c2y)/length(dd)*c2successfulrate+length(c3y)/length(dd)*c3successfulrate;
%  totalsuccess=length(c2y)/(length(dd)-20)*c2successfulrate+length(c3y)/(length(dd)-20)*c3successfulrate;
 
% figure(1)
% subplot(2,1,1)
% plot(1:T,c1weight1,'b-','LineWidth',1.5)
% legend('the weight of flow')
% title('below 2500')
% subplot(2,1,2)
% plot(1:c1dimension,c1weight2,'r-','LineWidth',1.5)
% legend('the weight of sea level')
% title('below 2500')

figure(2)
subplot(2,1,1)
plot(1:T,c2weight1,'b-','LineWidth',1.5)
legend('the weight of flow')
title('between 2500 and 3500')
subplot(2,1,2)
 plot(1:c2dimension,c2weight2,'r-','LineWidth',1.5)
%plot(1:T,c2weight2,'r-','LineWidth',1.5)
legend('the weight of sea level')
title('between 2500 and 3500')

figure(3)
subplot(2,1,1)
plot(1:T,c3weight1,'b-','LineWidth',1.5)
legend('the weight of flow')
title('above 3500')
subplot(2,1,2)
plot(1:c3dimension,c3weight2,'r-','LineWidth',1.5)
% plot(1:T,c3weight2,'r-','LineWidth',1.5)
legend('the weight of sea level')
title('above 3500')

figure(4)
plot(c2p)
hold on
plot(c2y,'r')


figure(5)
plot(c3p)
hold on
plot(c3y,'r')
% 
% 
an=[c1successfulrate,c2successfulrate,c3successfulrate,totalsuccess]


 for i=1:length(c2y)
     if c2y(i)==0&c2P(i)==0
         c2aaaa(i)=0;
     elseif c2y(i)==0&c2P(i)==1
         c2aaaa(i)=0.1;
     elseif c2y(i)==1&c2P(i)==0
         c2aaaa(i)=10;
     elseif c2y(i)==1&c2P(i)==1
         c2aaaa(i)=11;
     end
 end

     c2num00=length(find(c2aaaa==0));
     c2num01=length(find(c2aaaa==0.1));
     c2num10=length(find(c2aaaa==10));
     c2num11=length(find(c2aaaa==11));
     
      [c2num00  c2num10
       c2num01   c2num11 ]


 for i=1:length(c3y)
     if c3y(i)==0&c3P(i)==0
         c3aaaa(i)=0;
     elseif c3y(i)==0&c3P(i)==1
         c3aaaa(i)=0.1;
     elseif c3y(i)==1&c3P(i)==0
         c3aaaa(i)=10;
     elseif c3y(i)==1&c3P(i)==1
         c3aaaa(i)=11;
     end
 end

     c3num00=length(find(c3aaaa==0));
     c3num01=length(find(c3aaaa==0.1));
     c3num10=length(find(c3aaaa==10));
     c3num11=length(find(c3aaaa==11));
     
     [c3num00  c3num10
       c3num01   c3num11 ]
