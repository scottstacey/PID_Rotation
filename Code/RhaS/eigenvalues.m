%% This script mainstains all parameters as in the RhaS model but iterates over
 % lots of values of alpha_1 and alpha_2, finding the imaginary and real
 % parts of the slowest eigenvalue, to determine how the buffering
 % parameters affect the response of the controller. I export the data from
 % this into R where I am more comfortable with creating nice graphics.
 % Note this is for the RNA model with no degradation of Z_1 or Z_2
clear; clc;
syms lambda

%% Fixed parameters from the antithetic model:
U        = 250000; 
k_1      = 0.1*60*60;
k_2      = 0.06*60*60;
k_3      = 1.5*60*60;
theta    = 0.025*60*60;       
K_U      = 178000;          
K_X      = 2600;             
delta    = 0.00039*60*60;   


%% n for the while loops:
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:10000
    for jj = 1:5
       data(ii,jj)=0;
    end
end

%% Nested while loops which iterate over all possible values and combinations 
 % of alpha_1 and alpha_2 between 0.00001 and 10.00001
while n < 100
    alpha_1 = ((0.00001 + (n * 0.05)) * 60 * 60);
    n_1 = 0;
    while n_1 < 100
        alpha_2 = ((0.00001 + (n_1 * 0.05)) * 60 * 60);
        S = lambda^4 + ((alpha_2 + 2*delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*lambda^3) + (delta*(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U)))) + (delta + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))))*lambda^2 + delta*(((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*theta + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))))*lambda + (2*K_X*alpha_1*k_2*k_3*theta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X +  ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))^2);
        S = vpa(S, 6);
        S = solve(S, lambda);
        data(((n*100) + (n_1+1)), 1) = alpha_1/60/60;
        data(((n*100) + (n_1+1)), 2) = alpha_2/60/60;
        data(((n*100) + (n_1+1)), 3) = abs(imag(S(4)));
        data(((n*100) + (n_1+1)), 4) = real(S(4));
        data(((n*100) + (n_1+1)), 5) = real(S(3));
        n_1 = n_1 + 1;
    end
    n = n + 1;
end



delta    = 0.00039*60*60*10;

%% n for the while loops:
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:10000
    for jj = 1:5
       data1(ii,jj)=0;
    end
end

%% Nested while loops which iterate over all possible values and combinations 
 % of alpha_1 and alpha_2 between 0.00001 and 10.00001
while n < 100
    alpha_1 = ((0.00001 + (n * 0.05)) * 60 * 60);
    n_1 = 0;
    while n_1 < 100
        alpha_2 = ((0.00001 + (n_1 * 0.05)) * 60 * 60);
        S = lambda^4 + ((alpha_2 + 2*delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*lambda^3) + (delta*(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U)))) + (delta + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))))*lambda^2 + delta*(((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*theta + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))))*lambda + (2*K_X*alpha_1*k_2*k_3*theta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X +  ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))^2);
        S = vpa(S, 6);
        S = solve(S, lambda);
        data1(((n*100) + (n_1+1)), 1) = alpha_1/60/60;
        data1(((n*100) + (n_1+1)), 2) = alpha_2/60/60;
        data1(((n*100) + (n_1+1)), 3) = abs(imag(S(4)));
        data1(((n*100) + (n_1+1)), 4) = real(S(4));
        data1(((n*100) + (n_1+1)), 5) = real(S(3));
        n_1 = n_1 + 1;
    end
    n = n + 1;
end



delta    = 0.00039*60*60*0.1;

%% n for the while loops:
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:10000
    for jj = 1:5
       data2(ii,jj)=0;
    end
end

%% Nested while loops which iterate over all possible values and combinations 
 % of alpha_1 and alpha_2 between 0.00001 and 10.00001
while n < 100
    alpha_1 = ((0.00001 + (n * 0.05)) * 60 * 60);
    n_1 = 0;
    while n_1 < 100
        alpha_2 = ((0.00001 + (n_1 * 0.05)) * 60 * 60);
        S = lambda^4 + ((alpha_2 + 2*delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*lambda^3) + (delta*(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U)))) + (delta + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))))*lambda^2 + delta*(((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*theta + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))))*lambda + (2*K_X*alpha_1*k_2*k_3*theta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X +  ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))^2);
        S = vpa(S, 6);
        S = solve(S, lambda);
        data2(((n*100) + (n_1+1)), 1) = alpha_1/60/60;
        data2(((n*100) + (n_1+1)), 2) = alpha_2/60/60;
        data2(((n*100) + (n_1+1)), 3) = abs(imag(S(4)));
        data2(((n*100) + (n_1+1)), 4) = real(S(4));
        data2(((n*100) + (n_1+1)), 5) = real(S(3));
        n_1 = n_1 + 1;
    end
    n = n + 1;
end


delta    = 0.00039*60*60*0.5;

%% n for the while loops:
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:10000
    for jj = 1:5
       data3(ii,jj)=0;
    end
end

%% Nested while loops which iterate over all possible values and combinations 
 % of alpha_1 and alpha_2 between 0.00001 and 10.00001
while n < 100
    alpha_1 = ((0.00001 + (n * 0.05)) * 60 * 60);
    n_1 = 0;
    while n_1 < 100
        alpha_2 = ((0.00001 + (n_1 * 0.05)) * 60 * 60);
        S = lambda^4 + ((alpha_2 + 2*delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*lambda^3) + (delta*(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U)))) + (delta + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))))*lambda^2 + delta*(((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*theta + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))))*lambda + (2*K_X*alpha_1*k_2*k_3*theta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X +  ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))^2);
        S = vpa(S, 6);
        S = solve(S, lambda);
        data3(((n*100) + (n_1+1)), 1) = alpha_1/60/60;
        data3(((n*100) + (n_1+1)), 2) = alpha_2/60/60;
        data3(((n*100) + (n_1+1)), 3) = abs(imag(S(4)));
        data3(((n*100) + (n_1+1)), 4) = real(S(4));
        data3(((n*100) + (n_1+1)), 5) = real(S(3));
        n_1 = n_1 + 1;
    end
    n = n + 1;
end

delta    = 0.00039*60*60*2;

%% n for the while loops:
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:10000
    for jj = 1:5
       data4(ii,jj)=0;
    end
end

%% Nested while loops which iterate over all possible values and combinations 
 % of alpha_1 and alpha_2 between 0.00001 and 10.00001
while n < 100
    alpha_1 = ((0.00001 + (n * 0.05)) * 60 * 60);
    n_1 = 0;
    while n_1 < 100
        alpha_2 = ((0.00001 + (n_1 * 0.05)) * 60 * 60);
        S = lambda^4 + ((alpha_2 + 2*delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*lambda^3) + (delta*(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U)))) + (delta + theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))))*lambda^2 + delta*(((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*theta + theta*((k_1*U)/(theta*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))*(alpha_2 + delta + 4*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))))*lambda + (2*K_X*alpha_1*k_2*k_3*theta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))*((2*alpha_1*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1)))^2 + delta*(sqrt(((alpha_2 + delta)* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(alpha_1))) -2*alpha_2* ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X +  ((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))^2);
        S = vpa(S, 6);
        S = solve(S, lambda);
        data4(((n*100) + (n_1+1)), 1) = alpha_1/60/60;
        data4(((n*100) + (n_1+1)), 2) = alpha_2/60/60;
        data4(((n*100) + (n_1+1)), 3) = abs(imag(S(4)));
        data4(((n*100) + (n_1+1)), 4) = real(S(4));
        data4(((n*100) + (n_1+1)), 5) = real(S(3));
        n_1 = n_1 + 1;
    end
    n = n + 1;
end

