syms lambda

k_1   = 0.1*60*60;
k_2   = 0.06*60*60;
k_3   = 1.5*60*60;
K_U   = 178000;
K_X   = 2600;
delta = 0.00039*60*60;
theta = 0.025*60*60;
n=0;


%% Creating an empty matrix for storing the data
for ii = 1:2500
    for jj = 1:4
       data_U(ii,jj)=0;
       data_k_1(ii,jj)=0;
       data_k_2(ii,jj)=0;
       data_k_3(ii,jj)=0;
       data_K_X(ii,jj)=0;
       data_delta(ii,jj)=0;
       data_theta(ii,jj)=0;
       data_K_U(ii,jj)=0;
    end
end
% %% U loop
% while n < 2500
%     U = ((0.0001 + (n * 1000)) * 60 * 60);
%     S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
%     theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
%     + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
%     + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
%     + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
%     + k_3*U - k_1*U)))^2);
%     S = vpa(S, 6);
%     A = solve(S, lambda);
%     data_U((n+1),1) = U;
%     data_U((n+1),2) = A(1);
%     data_U((n+1),3) = A(2);
%     data_U((n+1),4) = A(3);
%     n = n + 1;
% end
% n = 0;
% U = 250000;

%% k_1 loop
while n < 2500
    k_1 = ((0.0001 + (n * 0.01)) * 60 * 60);
    S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
    theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
    + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
    + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
    + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
    + k_3*U - k_1*U)))^2);
    S = vpa(S, 6);
    A = solve(S, lambda);
    data_k_1((n+1),1) = k_1;
    data_k_1((n+1),2) = A(1);
    data_k_1((n+1),3) = A(2);
    data_k_1((n+1),4) = A(3);
    n = n + 1;
end
n = 0;
k_1 = 0.1*60*60;

%% k_2 loop
while n < 2500
    k_2 = ((0.0001 + (n * 0.01)) * 60 * 60);
    S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
    theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
    + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
    + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
    + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
    + k_3*U - k_1*U)))^2);
    S = vpa(S, 6);
    A = solve(S, lambda);
    data_k_2((n+1),1) = k_2;
    data_k_2((n+1),2) = A(1);
    data_k_2((n+1),3) = A(2);
    data_k_2((n+1),4) = A(3);
    n = n + 1;
end
n = 0;
k_2 = 0.06*60*60;

%% k_3 loop
while n < 2500
    k_3 = ((0.0001 + (n * 0.01)) * 60 * 60);
    S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
    theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
    + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
    + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
    + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
    + k_3*U - k_1*U)))^2);
    S = vpa(S, 6);
    A = solve(S, lambda);
    data_k_3((n+1),1) = k_3;
    data_k_3((n+1),2) = A(1);
    data_k_3((n+1),3) = A(2);
    data_k_3((n+1),4) = A(3);
    n = n + 1;
end
n = 0;
k_3 = 1.5*60*60;


% % K_U loop
% while n < 2500
%     K_U = ((0.0001 + (n * 1000)));
%     S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
%     theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
%     + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
%     + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
%     + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
%     + k_3*U - k_1*U)))^2);
%     S = vpa(S, 6);
%     A = solve(S, lambda);
%     data_K_U((n+1),1) = K_U;
%     data_K_U((n+1),2) = A(1);
%     data_K_U((n+1),3) = A(2);
%     data_K_U((n+1),4) = A(3);
%     n = n + 1;
% end
% n = 0;
% K_U = 250000;

% %% K_X loop
% while n < 2500
%     K_X = ((0.0001 + (n * 1000)));
%     S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
%     theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
%     + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
%     + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
%     + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
%     + k_3*U - k_1*U)))^2);
%     S = vpa(S, 6);
%     A = solve(S, lambda);
%     data_K_X((n+1),1) = K_X;
%     data_K_X((n+1),2) = A(1);
%     data_K_X((n+1),3) = A(2);
%     data_K_X((n+1),4) = A(3);
%     n = n + 1;
% end
% n = 0;
% K_X = 2600;

%% delta loop
while n < 2500
    delta = ((0.0001 + (n * 0.0001)) * 60 * 60);
    S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
    theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
    + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
    + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
    + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
    + k_3*U - k_1*U)))^2);
    S = vpa(S, 6);
    A = solve(S, lambda);
    data_delta((n+1),1) = delta;
    data_delta((n+1),2) = A(1);
    data_delta((n+1),3) = A(2);
    data_delta((n+1),4) = A(3);
    n = n + 1;
end
n = 0;
delta = 0.00039*60*60;

%% theta loop
while n < 2500
    theta = ((0.0001 + (n * 0.0001)) * 60 * 60);
    S = lambda^3 + lambda^2*(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) + ...
    theta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))) + delta) ...
    + lambda*(theta*delta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)) ...
    + theta*delta*((k_1 * U)/(theta*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2))*(K_U + U))))...
    + (theta*k_2*k_3*K_X*((delta*((k_1*U*K_X)/(k_3*K_U + k_3*U - k_1*U)))/(k_2)))/((K_X + ((k_1*U*K_X)/(k_3*K_U ...
    + k_3*U - k_1*U)))^2);
    S = vpa(S, 6);
    A = solve(S, lambda);
    data_theta((n+1),1) = theta;
    data_theta((n+1),2) = A(1);
    data_theta((n+1),3) = A(2);
    data_theta((n+1),4) = A(3);
    n = n + 1
end
n = 0;
theta = 0.025*60*60;

%% Creating a plot of imaginary parts of eigenvalues vs parameter values 
% figure(1);
% set(gca, 'fontsize', 12);
% plot(data_U(:,1), abs(imag(data_U(:,2))), 'r', data_U(:,1), abs(imag(data_U(:,3))), 'g', data_U(:,1), abs(imag(data_U(:,4))), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('U concentration (nM)');
% ylabel('Imaginary part of eigenvalues');
% title('Antithetic Controller U Variation');
% 
% figure(2);
% set(gca, 'fontsize', 12);
% plot(data_U(:,1), real(data_U(:,2)), 'r', data_U(:,1), real(data_U(:,3)), 'g', data_U(:,1), real(data_U(:,4)), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('U concentration (nM)');
% ylabel('Real part of eigenvalue');
% title('Antithetic Controller U Variation');
% 
% figure(3);
% set(gca, 'fontsize', 12);
% plot(data_K_U(:,1), abs(imag(data_K_U(:,2))), 'r', data_K_U(:,1), abs(imag(data_K_U(:,3))), 'g', data_K_U(:,1), abs(imag(data_K_U(:,4))), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('K_U');
% ylabel('Imaginary part of eigenvalues');
% title('Antithetic Controller K_U Variation');
% 
% figure(4);
% set(gca, 'fontsize', 12);
% plot(data_K_U(:,1), real(data_K_U(:,2)), 'r', data_K_U(:,1), real(data_K_U(:,3)), 'g', data_K_U(:,1), real(data_K_U(:,4)), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('K_U');
% ylabel('Real part of eigenvalue');
% title('Antithetic Controller K_U Variation');
% 
% figure(5);
% set(gca, 'fontsize', 12);
% plot(data_K_U(:,1), abs(imag(data_K_X(:,2))), 'r', data_K_X(:,1), abs(imag(data_K_X(:,3))), 'g', data_K_X(:,1), abs(imag(data_K_X(:,4))), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('K_X');
% ylabel('Imaginary part of eigenvalues');
% title('Antithetic Controller K_X Variation');
% 
% figure(6);
% set(gca, 'fontsize', 12);
% plot(data_K_X(:,1), real(data_K_X(:,2)), 'r', data_K_X(:,1), real(data_K_X(:,3)), 'g', data_K_X(:,1), real(data_K_X(:,4)), 'b', 'LineWidth', 3);
% legend('lambda_1', 'lambda_2', 'lambda_3');
% xlabel('K_X');
% ylabel('Real part of eigenvalue');
% title('Antithetic Controller K_X Variation');

figure(7);
set(gca, 'fontsize', 12);
plot(data_k_1(:,1), abs(imag(data_k_1(:,2))), 'r', data_k_1(:,1), abs(imag(data_k_1(:,3))), 'g', data_k_1(:,1), abs(imag(data_k_1(:,4))), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_1');
ylabel('Imaginary part of eigenvalues');
title('Antithetic Controller k_1 Variation');

figure(8);
set(gca, 'fontsize', 12);
plot(data_k_1(:,1), real(data_k_1(:,2)), 'r', data_k_1(:,1), real(data_k_1(:,3)), 'g', data_k_1(:,1), real(data_k_1(:,4)), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_1');
ylabel('Real part of eigenvalue');
title('Antithetic Controller k_1 Variation');


figure(9);
set(gca, 'fontsize', 12);
plot(data_k_2(:,1), abs(imag(data_k_2(:,2))), 'r', data_k_2(:,1), abs(imag(data_k_2(:,3))), 'g', data_k_2(:,1), abs(imag(data_k_2(:,4))), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_2');
ylabel('Imaginary part of eigenvalues');
title('Antithetic Controller k_2 Variation');

figure(10);
set(gca, 'fontsize', 12);
plot(data_k_2(:,1), real(data_k_2(:,2)), 'r', data_k_2(:,1), real(data_k_2(:,3)), 'g', data_k_2(:,1), real(data_k_2(:,4)), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_2');
ylabel('Real part of eigenvalue');
title('Antithetic Controller k_2 Variation');

figure(11);
set(gca, 'fontsize', 12);
plot(data_k_3(:,1), abs(imag(data_k_3(:,2))), 'r', data_k_3(:,1), abs(imag(data_k_3(:,3))), 'g', data_k_3(:,1), abs(imag(data_k_3(:,4))), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_3');
ylabel('Imaginary part of eigenvalues');
title('Antithetic Controller k_3 Variation');

figure(12);
set(gca, 'fontsize', 12);
plot(data_k_3(:,1), real(data_k_3(:,2)), 'r', data_k_3(:,1), real(data_k_3(:,3)), 'g', data_k_3(:,1), real(data_k_3(:,4)), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('k_3');
ylabel('Real part of eigenvalue');
title('Antithetic Controller k_3 Variation');

figure(13);
set(gca, 'fontsize', 12);
plot(data_theta(:,1), abs(imag(data_theta(:,2))), 'r', data_theta(:,1), abs(imag(data_theta(:,3))), 'g', data_theta(:,1), abs(imag(data_theta(:,4))), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('theta');
ylabel('Imaginary part of eigenvalues');
title('Antithetic Controller theta Variation');

figure(14);
set(gca, 'fontsize', 12);
plot(data_theta(:,1), real(data_theta(:,2)), 'r', data_theta(:,1), real(data_theta(:,3)), 'g', data_theta(:,1), real(data_theta(:,4)), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('theta');
ylabel('Real part of eigenvalue');
title('Antithetic Controller theta Variation');

figure(15);
set(gca, 'fontsize', 12);
plot(data_delta(1:500,1), abs(imag(data_delta(1:500,2))), 'r', data_delta(1:500,1), abs(imag(data_delta(1:500,3))), 'g', data_delta(1:500,1), abs(imag(data_delta(1:500,4))), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('delta');
ylabel('Imaginary part of eigenvalues');
title('Antithetic Controller delta Variation');

figure(16);
set(gca, 'fontsize', 12);
plot(data_delta(:,1), real(data_delta(:,2)), 'r', data_delta(:,1), real(data_delta(:,3)), 'g', data_delta(:,1), real(data_delta(:,4)), 'b', 'LineWidth', 3);
legend('lambda_1', 'lambda_2', 'lambda_3');
xlabel('delta');
ylabel('Real part of eigenvalue');
title('Antithetic Controller delta Variation');


 
