function [S_mc] = MC(S0,r,sigma)
%MC 模特卡洛模拟
  
dt=1/240;                   % 一天的年单位时间
expTerm=r*dt;                 % 漂移项dt
stddev=sigma*sqrt(dt);       % 波动项o:dz(t) 
nDays1 = 23;                    
nTrials = 100;             % 要模拟的总次数
S = ones(1,nTrials)*S0;

for j=1:nTrials 
    n = randn(1,nDays1);           %生成nDays个标准正态分布随机数 
    for i=1:nDays1 
        dS = S(i,j)*(expTerm+stddev*n(i));   % 模拟计算股票价格的增量
        S(i+1,j) = S(i,j)+dS;                      %计算股票价格
    end
end
   S_mc = S(4,:);                 % 计算每天模拟的股票价格的均值，作为价格的估值 

end

