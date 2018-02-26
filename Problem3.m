load('data.mat'); %contain raw .txt data 

Train = table2array(trainingdata);
Test = table2array(testdata);

[n,m] = size(Train);
tag = Train(:,m);
X1 = Train(tag == 1,1:m-1);
X2 = Train(tag == 2,1:m-1);

C1 = 0.6;
C2 = 0.4;

mu1 = mean(X1)
mu2 = mean(X2)

%independent Si
S1 = cov(X1);
S2 = cov(X2);
err_independent = 0;
for i = 1 : length(Test)
    x = Test(i,1:m-1);
    g1 = -1/2*log(det(S1))-1/2*(x-mu1)*S1^-1*(x-mu1)'+log(C1);
    g2 = -1/2*log(det(S2))-1/2*(x-mu2)*S2^-1*(x-mu2)'+log(C2);
    if g1 - g2 >0 && Test(i,m)==2
        err_independent = err_independent +1;
    end
    if g1 - g2 <0 && Test(i,m)==1
        err_independent = err_independent +1;
    end    
end

%share Si, S1 = S2;
S = cov([X1;X2]);
err_share = 0;
for i = 1 : length(Test)
    x = Test(i,1:m-1);
    g1 = -1/2*(x-mu1)*S^-1*(x-mu1)'+log(C1);
    g2 = -1/2*(x-mu2)*S^-1*(x-mu2)'+log(C2);
    if g1 - g2 >0 && Test(i,m)==2
        err_share = err_share +1;
    end
    if g1 - g2 <0 && Test(i,m)==1
        err_share = err_share +1;
    end
end

%S1 and S2 are diagnal
S1_dia = diag(diag(cov(X1)));
S2_dia = diag(diag(cov(X2)));
err_diagnal = 0;
for i = 1 : length(Test)
    x = Test(i,1:m-1);
    g1 = -1/2*log(det(S1_dia))-1/2*(x-mu1)*S1_dia^-1*(x-mu1)'+log(C1);
    g2 = -1/2*log(det(S2_dia))-1/2*(x-mu2)*S2_dia^-1*(x-mu2)'+log(C2);
    if g1 - g2 >0 && Test(i,m)==2
        err_diagnal = err_diagnal +1;
    end
    if g1 - g2 <0 && Test(i,m)==1
        err_diagnal = err_diagnal +1;
    end    
end

%S1 and S2 are diagnal, S1 = alpha1 I, S2 = alpha2 I
err_dia2 = 0;
S1_dia2 = sum(diag(cov(X1)))/(m-1) * eye(m-1);
S2_dia2 = sum(diag(cov(X2)))/(m-1) * eye(m-1);
alpha1 = sum(diag(cov(X1)))/(m-1)
alpha2 = sum(diag(cov(X2)))/(m-1)

for i = 1 : length(Test)
    x = Test(i,1:m-1);
    g1 = -1/2*log(det(S1_dia2))-1/2*(x-mu1)*S1_dia2^-1*(x-mu1)'+log(C1);
    g2 = -1/2*log(det(S2_dia2))-1/2*(x-mu2)*S2_dia2^-1*(x-mu2)'+log(C2);
    if g1 - g2 >0 && Test(i,m)==2
        err_dia2 = err_dia2 +1;
    end
    if g1 - g2 <0 && Test(i,m)==1
        err_dia2 = err_dia2 +1;
    end
end
