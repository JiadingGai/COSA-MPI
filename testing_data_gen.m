n = 2000;
n0 = 10;
MU1 = zeros(n,1);
COVAR1 = eye(n,n);

MU2 = zeros(n,1);
MU2(1:n0) = 1.5;
COVAR2 = eye(n,n);
COVAR2(1:n0,1:n0) = 0.2*eye(n0,n0);

x1 = gsamp(MU1,COVAR1,85);
x2 = gsamp(MU2,COVAR2,15);

x = [x2;x1];

mean = sum(x,1) / 100;
the_std = std(x);

% save matrix x to e1.txt in a row-wise manner
% ie., the second element read from e1.txt is
% the entry on the first row of the second col-
% umn of x.
dlmwrite('e1.txt', x, 'delimiter', '\t', 'precision',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 3;
n0 = 2;
MU1 = zeros(n,1);
COVAR1 = eye(n,n);

MU2 = zeros(n,1);
MU2(1:n0) = 1.5;
COVAR2 = eye(n,n);
COVAR2(1:n0,1:n0) = 0.04*eye(n0,n0);

x1 = gsamp(MU1,COVAR1,7);
x2 = gsamp(MU2,COVAR2,3);

x = [x2;x1];

mean = sum(x,1) / 10;
the_std = std(x);

x = (x-repmat(mean,10,1)) ./ repmat(the_std,10,1);
% save matrix x to e1.txt in a row-wise manner
% ie., the second element read from e1.txt is
% the entry on the first row of the second col-
% umn of x.
dlmwrite('e2.txt', x, 'delimiter', '\t', 'precision',6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 4000;
n0 = 500;
MU1 = zeros(n,1);
COVAR1 = eye(n,n);

MU2 = zeros(n,1);
MU2(1:n0) = 1.5;
COVAR2 = eye(n,n);
COVAR2(1:n0,1:n0) = 0.2*eye(n0,n0);

x1 = gsamp(MU1,COVAR1,85+50);
x2 = gsamp(MU1,COVAR2,50+50);
x3 = gsamp(MU2,COVAR2,15+50);

x = [x3;x2;x1];

mean = sum(x,1) / 300;
the_std = std(x);

% save matrix x to e1.txt in a row-wise manner
% ie., the second element read from e1.txt is
% the entry on the first row of the second col-
% umn of x.
dlmwrite('e01.txt', x, 'delimiter', '\t', 'precision',6);