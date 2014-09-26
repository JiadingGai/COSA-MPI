close all;
clear all;

MU1 = [1 2];
SIGMA1 = [2 0.25; 0.25 .5];
%SIGMA1 = [2 0; 0 2];
MU2 = [-3 -5];
SIGMA2 = [1 0.95; 0.95 1];
MU3 = [6 -5];
SIGMA3 = [6 0.95; 0.95 1];

X = [mvnrnd(MU1,SIGMA1,1500);mvnrnd(MU2,SIGMA2,1500);mvnrnd(MU3,SIGMA3,1500);];

N = 4500;
K = 3;
D = 2;
% M-step
m = rand(D,K);
s = zeros(D,D,K);
for k=1:K
    s(:,:,k) = eye(D,D);
end
p = ones(K,1) / N;
% E-step
P = zeros(K,N);

for i = 1:600
    % E-step
    for k = 1:K
        mu = m(:,k)';
        Sigma = squeeze( s(:,:,k) );
        P(k,:) = mvnpdf(X,mu,Sigma);
    end
    pmat = repmat(p,1,N);      
    P = P .* pmat;
   
    denorm = repmat(sum(P,1),K,1);
    P = P ./ (denorm);
       
    % M-step
    m = P * X;
    m = m ./ repmat(sum(P,2),1,D);
    m = m';
   
    for k = 1:K
        mvector = repmat( m(:,k)',N,1 );
        tempP = repmat( P(k,:),D,1 );
        temp1 = tempP .* (X-mvector)';
        temp = temp1 * (X-mvector);
        temp2 = temp / (sum(P(k,:)));
        s(:,:,k) = temp2;
    end
    p = sum(P,2) / N;
end

figure;
scatter(X(:,1),X(:,2),10,'.');
hold on;
plot(m(1,1),m(2,1),'r*');
hold on;
plot(m(1,2),m(2,2),'g*');
hold on;
plot(m(1,3),m(2,3),'b*');


options = statset('Display','final');
obj = gmdistribution.fit(X,3,'Options',options);
h = ezcontour(@(x,y)pdf(obj,[x y]),[-12 12],[-12 12]);


idx = cluster(obj,X);
cluster1 = X(idx == 1,:);
cluster2 = X(idx == 2,:);
cluster3 = X(idx == 3,:);

figure;
h1 = scatter(cluster1(:,1),cluster1(:,2),10,'r.');
hold on;
h2 = scatter(cluster2(:,1),cluster2(:,2),10,'g.');
hold on;
h3 = scatter(cluster3(:,1),cluster3(:,2),10,'b.');


