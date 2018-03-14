%% Homework 10 Question 1
clc
clear all
close all
%% P(z) is a Half Band, Zero Phase Filter
syms a b c d z
g = expand((1+z)^4);
h = expand((1+(z^-1))^4);
i = expand((a*z^3)+(b*z^2)+(c*z)+d+(c*z^-1)+(b*z^-2)+(a*z^-3));
coeffs_z = collect((h*g*i), z);
display(coeffs_z)
% Determining the values of coefficients a b c and d 
% P0 = 1 , Even Powers of z have coefficients zero.
% Solving the four equations
eqns =[16*a + 56*b + 112*c + 70*d - 1 == 0, 56*a + 28*b + 8*c + d == 0, 56*a + 71*b + 64*c + 28*d == 0, 8*a + b == 0];
vars = [a b c d];
[sola, solb, solc, sold] = solve(eqns, vars);

% Substituting the values of a, b, c, d to obtain H0, H1, G0, G1
r1 = expand((sola*z^3)+(solb*z^2)+(solc*z)+sold+(solc*z^-1)+(solb*z^-2)+(sola*z^-3));
coeffs_z1 = collect((h*g*r1*(z^-7)), z); % Factorizing out z^k
factor_z1 = factor(coeffs_z1);
display(factor_z1)

% H0(z) and G0(z) are 8-tap symmetric linear phase lowpass filters 
coeffs_h0 = collect(((z^-7)*((z+1)^7)*(1/32)), z);
coeffs_g0 = collect(((5*z^6 - 40*z^5 + 131*z^4 - 208*z^3 + 131*z^2 - 40*z + 5)*(z+1)*(z^-7)*(-1/64)), z);
% H1(z) =  G0(-z)
% G1(z) = -H0(-z)
coeffs_h1 = collect(((5 - (35*(-z)^-1) + (91*z^-2) - (77*(-z)^-3) - (77*z^-4) + (91*(-z)^-5) - (35*z^-6) + (5*(-z)^-7))*(-1/64)), z);
coeffs_g1 = collect(((1 + (7*(-z)^-1) + (21*z^-2) + (35*(-z)^-3) + (35*z^-4) + (21*(-z)^-5) + (7*z^-6) + ((-z)^-7))*(1/32)*(-1)), z);
% Perfect Reconstruction Filters:
% A(z) = (1/2){(H0(-z)*G0(z))+(H1(-z)*G1(z))}
coeffs_h11 = collect(((5 - (35*(z)^-1) + (91*z^-2) - (77*(z)^-3) - (77*z^-4) + (91*(z)^-5) - (35*z^-6) + (5*(z)^-7))*(-1/64)),z);
display(coeffs_h11)
coeffs_h01 = collect(((1 + (7*(-z)^-1) + (21*z^-2) + (35*(-z)^-3) + (35*z^-4) + (21*(-z)^-5) + (7*z^-6) + ((-z)^-7))*(1/32)),z);
display(coeffs_h01)

% Coefficients of A(z)= 0
coeffs_az = collect (((coeffs_h01*coeffs_g0)+(coeffs_g1*coeffs_h11)),z);
display(coeffs_az)

% Coefficients of T(z)= C(z^-k)
coeffs_tz = collect(((1/2)*((coeffs_h0*coeffs_h11)-(coeffs_h1*coeffs_h01))),z);
display(coeffs_tz)

%% Magnitude plots for Analysis Filter
mag_h1 = [-5  -35  -91  -77  77 91  35 5]./64;
mag_h0 = [1 7  21 35  35 21  7 1]./32;

figure(1)
[h1,w1]=freqz(mag_h0);
hmag_h0 = abs(h1);
tmp = angle(h1);
tol = 0.85*pi;
hang1 = unwrap(tmp, tol);
plot(w1,hmag_h0,'--')
hold on;

[h2,w2]=freqz(mag_h1,1);
hmag_h1 = abs(h2);
tmp = angle(h2);
tol = 0.85*pi;
hang2 = unwrap(tmp, tol);
plot(w2,hmag_h1)
legend('|H0(w)|','|H1(w)|')
grid on;
xlabel('Frequency (radians)','fontsize',12);
ylabel('Magnitude |H(w)|','fontsize',12);
title('Magnitude Plot for Analysis Filter');

figure(3)
plot(w1,hang1,'--')
hold on;
plot(w2,hang2)
legend('H0(w)','H1(w)')
grid on;
title('Phase Plot for Analysis Filter')
xlabel('Frequency (radians)')
ylabel('Phase (radians)')


%% Magnitude plots for Synthesis Filter
mag_g1 = [-1 7 -21 35 -35 21 -7 1]./32;
mag_g0 = [-5 35 -91 77  77 -91 35 -5]./64;

figure(2)
[g1,w1]=freqz(mag_g0);
hmag_g0 = abs(g1);
tmp = angle(g1);
tol = 0.85*pi;
hang3 = unwrap(tmp, tol);
plot(w1,hmag_g0,'--')
hold on;

[g2,w2]=freqz(mag_g1,1);
hmag_g1 = abs(g2);
tmp = angle(g2);
tol = 0.85*pi;
hang4 = unwrap(tmp, tol);
plot(w2,hmag_g1) 
legend('|G0(w)|','|G1(w)|')
grid on;
xlabel('Frequency (radians)','fontsize',12);
ylabel('Magnitude |G(w)|','fontsize',12);
title('Magnitude Plot for Synthesis Filter');

figure(4)
plot(w1,hang3,'--')
hold on
plot(w2,hang4)
legend('G0(w)','G1(w)')
grid on;
title('Phase Plot for Synthesis Filter')
xlabel('Frequency (radians)')
ylabel('Phase (radians)')

