syms lambda Z_2 theta Z_1 k_3 K_X X k_2 k_4 alpha_1 alpha_2 delta_1 delta_2 R R_star X_star

delta_3 = ((delta_1+delta_2)/2);

M = [(-theta*Z_2-lambda), (-theta*Z_1), (0), (0), (0), (0); ...
    (-theta*Z_2), (-theta*Z_1 - lambda), (0), ((k_3 * K_X)/((K_X + X)^2)), (0), (0);...
    (k_2), (0), (-delta_1-4*alpha_1*R-alpha_1*R_star-lambda), (2*alpha_2), (-alpha_1*R), (alpha_2);...
    (0), (0), (2*alpha_1*R), (-alpha_2-delta_1-lambda), (0), (0);
    (0), (0), (-alpha_1*R_star), (0), (-delta_2-alpha_1*R-lambda), (alpha_2);
    (0), (0), (alpha_1*R_star), (0), (alpha_1*R), (-alpha_2-delta_3-lambda)];



A = det(M)