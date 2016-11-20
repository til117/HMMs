% Artificial Intelligence
% HW1: HMM

clear all, close all, clc  %#ok<*DUALC,*CLALL>

% User Input

    % HMM
    HMM = 4;

%% HMM 1
if HMM == 1;

% Assume that A = NxN and B = NxM and P = 1xN
A = [4 4 0.2 0.5 0.3 0.0 0.1 0.4 0.4 0.1 0.2 0.0 0.4 0.4 0.2 0.3 0.0 0.5];
B = [4 3 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0 0.2 0.6 0.2];
P = [1 4 0.0 0.0 0.0 1.0];

N = A(2);
M = B(2);
Out1 = zeros(1,N);
Out = zeros(1,M);

% Column multiplication
for e = 1:N
    for i = 1:N
        Out1(e) = Out1(e) + P(i+2) * A(2+(i-1)*N+e);
    end
end

for e = 1:M
    for i = 1:N
        Out(e) = Out(e) + Out1(i) * B(2+(i-1)*M+e);
    end
end
Out = [[1 M] Out] %#ok<*NOPTS>

end

%% HMM 2
if HMM == 2

A1 = [0.0 0.8 0.1 0.1; 0.1 0.0 0.8 0.1; 0.1 0.1 0.0 0.8; 0.8 0.1 0.1 0.0];
B1 = [0.9 0.1 0.0 0.0; 0.0 0.9 0.1 0.0; 0.0 0.0 0.9 0.1; 0.1 0.0 0.0 0.9];

% Assume that A = NxN and B = NxM and P = 1xN
A = [4 4 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.8 0.1 0.1 0.0];
B = [4 4 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.1 0.0 0.0 0.9];
P = [1 4 1.0 0.0 0.0 0.0];
O = [8 0 1 2 3 0 1 2 3];

N = A(2);
M = B(2);
T = O(1);

alpha1 = zeros(N,1);
first_obs = O(2) + 1;
for i = 1:N
    alpha1(i) = P(i+2) * B(2+(i-1)*M+first_obs);
end

prev_alpha = alpha1;
alpha = zeros(N,1);
for t = 2:T
    % Row multiplication
    obs = O(t+1) + 1;
    alpha = zeros(N,1);
    for i = 1:N
        for j = 1:N
            alpha(i) = alpha(i) + prev_alpha(j)*A(2+i+(j-1)*N);
        end
        alpha(i) = alpha(i) * B(2+(i-1)*M+obs);
    end
    prev_alpha = alpha;
    alpha
end
%alpha

Prob = 0;
for i = 1:N
    Prob = Prob + alpha(i);
end
Prob

end

%% HMM 3
if HMM == 3

A = [4 4 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.8 0.1 0.1 0.0];  
B = [4 4 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.1 0.0 0.0 0.9];  
P = [1 4 1.0 0.0 0.0 0.0];  
O = [4 1 1 2 2];  

% A = [4 4 0.6 0.1 0.1 0.2 0.0 0.3 0.2 0.5 0.8 0.1 0.0 0.1 0.2 0.0 0.1 0.7]; 
% B = [4 4 0.6 0.2 0.1 0.1 0.1 0.4 0.1 0.4 0.0 0.0 0.7 0.3 0.0 0.0 0.1 0.9]; 
% P = [1 4 0.5 0.0 0.0 0.5]; 
% O = [4 2 0 3 1];

N = A(2);
M = B(2);
T = O(1);

possible_paths = 0;
argmax_state_index = zeros(1,T);
delta = zeros(N*T,1);
first_obs = O(2) + 1;
for i = 1:N
    delta(i) = P(i+2) * B(2+(i-1)*M+first_obs);
end

state_index = zeros(N*T,1);
for t = 2:T
    obs = O(t+1) + 1;
    for i = 1:N
        for j = 1:N
            tmp = delta((t-2)*N+j)*A(2+i+(j-1)*N) * B(2+(i-1)*M+obs);
            if tmp == delta((t-1)*N+i)
                if tmp ~= 0
                    possible_paths = possible_paths + 1;
                end
            end
            if tmp > delta((t-1)*N+i)
                delta((t-1)*N+i) = tmp;
                state_index((t-1)*N+i) = j;
            end
        end
    end
    

end

out = zeros(T,1);
for t = 1:T
    comp = (t-1)*N+1;
    out(t) = 1;
    for i = 1:N-1
        if delta((t-1)*N+1+i) > delta(comp)
            out(t) = i+1;
            comp = (t-1)*N+i+1
        end
    end
end
    
delta
state_index
out
possible_paths = 2^possible_paths
end


%% HMM 4
if HMM == 4
tic
% A = [4 4 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.1 0.1 0.1 0.0 0.8 0.8 0.1 0.1 0.0];
% B = [4 4 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.9 0.1 0.1 0.0 0.0 0.9];
% P = [1 4 1.0 0.0 0.0 0.0];
% O = [8 0 1 2 3 0 1 2 3];

A = [4 4 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4]; 
B = [4 4 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4 0.2 0.2 0.2 0.2 0.4]; 
P = [1 4 0.241896 0.266086 0.249153 0.242864]; 
O = [1000 0 1 2 3 3 0 0 1 1 1 2 2 2 3 0 0 0 1 1 1 2 3 3 0 0 0 1 1 1 2 3 3 0 1 2 3 0 1 1 1 2 3 3 0 1 2 2 3 0 0 0 1 1 2 2 3 0 1 1 2 3 0 1 2 2 2 2 3 0 0 1 2 3 0 1 1 2 3 3 3 0 0 1 1 1 1 2 2 3 3 3 0 1 2 3 3 3 3 0 1 1 2 2 3 0 0 0 0 0 0 0 0 0 1 1 1 1 1 2 2 2 3 3 3 3 0 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 3 3 3 0 1 2 3 0 1 1 1 2 3 0 1 1 2 2 2 2 2 3 0 1 1 1 2 2 2 2 3 0 0 0 0 0 1 1 1 1 2 2 3 3 0 1 2 3 3 0 0 0 0 0 0 1 1 2 2 3 0 0 1 1 1 1 1 1 2 3 3 0 0 1 1 1 2 3 0 0 1 2 3 0 1 1 2 3 3 0 0 0 1 2 3 3 3 0 1 1 1 1 2 3 3 3 3 3 3 0 1 2 2 2 2 2 2 3 0 1 1 1 2 2 3 3 3 3 0 1 2 3 0 0 0 1 1 2 2 3 0 0 0 0 0 0 0 1 2 2 2 3 3 3 3 0 0 1 2 2 2 3 3 3 0 0 1 2 2 3 0 0 0 0 1 1 1 2 3 3 3 3 3 3 3 3 0 1 2 3 0 0 1 2 3 3 3 0 0 0 0 0 1 1 1 1 2 3 0 0 0 1 2 2 3 3 0 0 0 1 1 1 1 1 2 3 3 3 3 0 1 1 1 2 2 3 0 1 2 3 3 3 3 0 0 0 0 1 2 3 3 0 1 2 2 3 3 0 0 1 1 2 3 3 0 1 2 2 3 3 3 0 0 1 1 2 3 3 3 3 0 0 1 1 2 3 3 0 1 2 3 0 1 1 2 2 3 0 1 2 3 3 0 1 1 1 2 2 2 3 3 0 0 1 1 1 1 1 2 3 3 3 0 1 1 2 2 2 2 3 3 0 0 1 2 3 0 1 1 2 2 2 2 3 0 0 1 2 2 3 0 0 0 0 0 1 1 1 2 3 0 0 1 2 3 3 0 0 0 1 2 2 2 3 3 0 0 0 1 2 2 2 2 2 3 0 1 1 2 3 0 0 1 1 1 2 2 3 0 0 0 0 1 1 1 2 2 3 0 1 1 1 2 2 2 3 3 0 0 1 2 2 3 3 3 0 1 1 2 3 0 0 0 0 0 1 2 2 2 3 3 3 0 0 0 1 2 3 0 1 1 2 3 3 3 0 1 2 2 2 3 0 0 1 1 1 1 2 3 3 0 0 0 0 1 2 3 3 3 0 0 0 1 1 2 3 0 1 1 1 1 2 2 2 2 2 2 3 0 0 0 0 1 2 2 2 2 3 0 1 2 2 3 0 1 2 3 0 1 2 3 0 0 0 1 1 2 2 3 3 0 1 1 1 1 2 2 3 3 0 1 1 1 2 2 2 3 3 3 0 1 1 2 3 3 0 1 2 3 0 0 0 0 1 2 3 0 0 0 0 0 0 1 2 2 3 3 0 0 1 2 3 0 1 2 2 3 0 0 0 1 1 2 2 2 2 2 3 3 3 3 3 0 1 2 2 3 3 3 3 3 0 0 1 1 2 2 3 0 0 1 2 2 3 3 3 0 0 0 1 2 2 2 2 3 3 0 1 2 3 0 0 1 1 1 2 2 3 0 0 1 1 2 2 2 3 3 0 0 1 1 1 1 1 2 3 3 3 0 1 2 2 2 2 3 3 3 3 3 3 0 0 0 0 0 0 1 2 3 0 0 1 1 1 2 3 0 0 1 1 2 2 2 2 3 3 3 0 1 1 2 2 2 3 3 0 0 0 0 0 0 1 2 2 3 3 0 0 0 0 0 0 1 2 3 3 3 0 1 1 1 2 2 2 2 2 3 3 3 0 1 2 2 2 3 3 3 3 0 0 0 0 1 2 3 3 3 3 3 3 0 0 1 1 1 1 2 3 0 1 2 3 0 1 1 2 3 3 3 0 0 0 0 1 1 2 3 3 3 3 0 0 1 1 1 2 2 2 2 2 2 3 3 0 0 0 1 2 3 0 0 1 1 2 2 3 3 3 3 3 0 0 1 2 2 2 2 3 0 0 1 1 1 1 1 2 3 3 0 0 1 1 1 2 3 3 3 0 0];

N = A(2);
M = B(2);
T = O(1);

iter = 0;
maxiter = 100;
logprob = 0;
oldlogprob = - inf;

% figure
% axis([0 30 0 10])

while iter < maxiter && logprob > oldlogprob
    %Initialize Parameters
    alpha = zeros(N*T,1);
    alphaT_sum = 0;
    beta = zeros(N*T,1);
    first_obs = O(2) + 1;
    gamma = zeros(N*T,1);
    digamma = zeros(N*N*T,1);
    
    %Initialize first alpha
    c0 = 0;
    for i = 1:N
        alpha(i) = P(i+2) * B(2+(i-1)*M+first_obs);
        c0 = c0 + alpha(i);
    end
    c0 = 1 / c0;
    for i = 1:N
        alpha(i) = alpha(i) * c0;
    end
    
    % Calculate Alpha
    c = zeros(T,1);
    for t = 2:T
        c(t) = 0;
        obs = O(t+1) + 1;
        for i = 1:N
            for j = 1:N
                alpha((t-1)*N+i) = alpha((t-1)*N+i) + alpha((t-2)*N+j)*A((j-1)*N+i+2);
            end
            alpha((t-1)*N+i) = alpha((t-1)*N+i) * B(2+(i-1)*M+obs);
            c(t) = c(t) + alpha((t-1)*N+i);
        end
        if c(t) == 0; disp('Warning'); end
        c(t) = 1 / c(t);
        for i = 1:N
            alpha((t-1)*N+i) = alpha((t-1)*N+i) * c(t);
        end
    end
    
    for i = 1:N
        alphaT_sum = alphaT_sum + alpha(T*N-N+i);
    end
    
    % Initialize beta
    for i = 1:N
        beta((T-1)*N+i) = c(T);
    end
    
    % Calculate Beta
    for tt = 1:T-1
        t = T - tt;
        obs_plus = O(t+2) + 1;
        for i = 1:N
            for j = 1:N
                beta((t-1)*N+i) = beta((t-1)*N+i) + A((i-1)*N+j+2)*beta(t*N+j)*B(2+(j-1)*M+obs_plus);
            end
            beta((t-1)*N+i) = beta((t-1)*N+i) * c(t);
        end
    end
    
    % Calculate Gamma
    for t = 1:T
        for i = 1:N
            gamma((t-1)*N+i) = alpha((t-1)*N+i) * beta((t-1)*N+i) / alphaT_sum;
        end
    end
    
    % Calculate DiGamma
    for t = 1:T-1
        obs_plus = O(t+2) + 1;
        for i = 1:N
            for j = 1:N
                digamma((t-1)*N*N+(i-1)*N+j) = alpha((t-1)*N+i) * A((i-1)*N+j+2)*beta(t*N+j)*B(2+(j-1)*M+obs_plus) / alphaT_sum;
            end
        end
    end
    
    % Control that Gammas are Correct
    control_gamma = zeros(N*T,1);
    for t = 1:T-1
        for i = 1:N
            for j = 1:N
                control_gamma((t-1)*N+i) = control_gamma((t-1)*N+i) + digamma((t-1)*N*N+(i-1)*N+j);
            end
        end
    end
    
    % Calculate P
    for i = 1:N
        P(i) = gamma(i);
    end
    % Calculate A and B
    sum_gamma = zeros(N,1);
    sum_gamma_b = zeros(N,1);
    sum_digamma = zeros(N*N,1);
    for i = 1:N
        for t = 1:T-1
            sum_gamma(i) = sum_gamma(i) + gamma((t-1)*N+i);
        end
    end
    for i = 1:N
        for j = 1:N
            for t = 1:T-1
                sum_digamma((i-1)*N+j) = sum_digamma((i-1)*N+j) + digamma((t-1)*N*N+(i-1)*N+j);
            end
        end
    end
    for i = 1:N
        for j = 1:N
            A((i-1)*N+j+2) = sum_digamma((i-1)*N+j) / sum_gamma(i);
        end
    end

    for i = 1:N
        for k = 1:M
            for t = 1:T-1
                obs = O(t+1)+1;
                if obs == k
                    sigmund = 1;
                else
                    sigmund = 0;
                end
                sum_gamma_b(i) = sum_gamma_b(i) + gamma((t-1)*N+i) * sigmund;
            end
        end
    end

    for i = 1:N
        for k = 1:M
            B((i-1)*M+k+2) = sum_gamma_b(i) / sum_gamma(i);
        end
    end
    
    % LogProb
%     oldlogprob = logprob;
%     logprob = 0;
%     for t = 1:T
%         logprob = logprob + log(c(t));
%     end
%     logprob = -logprob;
    
    
    iter = iter + 1
end




toc
end












