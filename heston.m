close all;
clear all;

% nsims - number of simultations
% s_init - initial asset price
% v_init - initial variance
% c - [mu kappa sigma theta]
% N - number of grid points
% dt - time step

nsims = 1000;
s_init = 100;
v_init = 0.03;
tbounds = [0 1];
c = [0.01 0.1 0.1 0.02];
N = 1000;
dt = (tbounds(2)-tbounds(1))/(N-1);

% vector of times
tvec = linspace(tbounds(1), tbounds(2), N);

% pre-allocate output
svec = zeros(nsims,N);

% set initial asset price
svec(:,1) = s_init;

% pre-allocate output
vvec = zeros(nsims,N);

% set initial variance
vvec(:,1) = v_init;

% correlation coefficient
rho = -0.5;

% pre-allocate correlated weiner increments
dW_correlated = zeros(2,N-1,nsims);

% correlation matrix
corr = [1 , rho ; rho , 1];

% Cholesky decomposition
lower_triangular_matrix = chol(corr,"lower");

% uncorrelated Weiner increments
dW_uncorrelated = sqrt(dt)*randn(2,N-1,nsims);

% iteratively compute correlated increments
for i = 1:nsims
    correlated_increments = lower_triangular_matrix*dW_uncorrelated(:,:,i);
    dW_correlated(:,:,i) = correlated_increments;
end

% stochastic Runge-Kutta method
for i = 1:nsims
    for j=2:numel(tvec)
        t = tbounds(1)+(j-1).*dt;
        s = svec(i,j-1);
        v = vvec(i,j-1);
        mu = c(1);
        kappa = c(2);
        theta = c(4);
        sig = c(3);
        dW1 = dW_correlated(1,:,i);
        dW2 = dW_correlated(2,:,i);
        k1s = dt.*mu.*s + (dW1(j-1)-(plusminusone().*sqrt(dt))).*sqrt(v).*s;
        k1v = dt.*kappa.*(theta-v) + (dW2(j-1)-(plusminusone().*sqrt(dt))).*sig.*sqrt(v);
        k2s = dt.*mu.*(s+k1s) + (dW1(j-1)+(plusminusone().*sqrt(dt))).*sqrt(v+k1v).*(s+k1s);
        k2v = dt.*kappa.*(theta-(v+k1v)) + (dW2(j-1)+(plusminusone().*sqrt(dt))).*sig.*sqrt(v+k1v);
        svec(i,j) = s + 0.5.*(k1s+k2s);
        vvec(i,j) = v + 0.5.*(k1v+k2v);
        if vvec(i,j) <= 0
            vvec(i,j) = 0;
        end
        if svec(i,j) <= 0
            svec(i,j) = 0;
        end
    end
end

figure;
hold on;

% plot all simulations
for i = 1:nsims
    plot(tvec, svec(i,:));
end

xlabel('Time');
ylabel('Asset Price');
title('Asset Price Paths');
grid on;

figure;
hold on;

% plot all simulations
for i = 1:nsims
    plot(tvec, vvec(i,:));
end

xlabel('Time');
ylabel('Variance');
title('Variance Paths');
grid on;
