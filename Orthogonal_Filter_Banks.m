clc
clear all;
close all;

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
R = [sola solb solc sold solc solb sola];
Q = vpa(roots(R));

A = [1 1];
R1 = [1 -0.32887591778603086663314159546921];
R2 = [1 -(0.28409629819182161933311014347123+0.24322822591037988274771916445887i)];
R3 = [1 -(0.28409629819182161933311014347123-0.24322822591037988274771916445887i)];
R = conv(R1,conv(R2,R3));
h0 = conv(A,conv(A,conv(A,conv(A,R))));
H0 = h0;
G0 = flip(H0);
H1 = zeros(1,length(G0));
G1 = zeros(1,length(H0));
factor1 = 1/32;
factor2 = 1/64;
for i = 1:1:length(H0)
    if (mod(i,2) == 1)
        H1(i) = G0(i);
        G1(i) = H0(i);
    else
        H1(i) = -G0(i);
        G1(i) = -H0(i);
    end
end
H0 = factor1.*H0;
G0 = factor2.*G0;
H1 = (factor2).*H1;
G1 = (-factor1).*G1;

figure(1)
[h0,w0] = freqz(H0,1);
scale = sqrt(2)/max(abs(h0));
H0 = scale*H0;
[h0,w0] = freqz(H0,1);
tmp = angle(h0);
tol = 0.85*pi;
hang1 = unwrap(tmp, tol);
plot(w0,abs(h0))
hold on
[h1,w1] = freqz(H1,1);
scale = sqrt(2)/max(abs(h1));
H1 = scale*H1;
[h1,w1] = freqz(H1,1);
tmp = angle(h1);
tol = 0.85*pi;
hang2 = unwrap(tmp, tol);
plot(w1,abs(h1))
legend('|H0(w)|','|H1(w)|')
grid on;
xlabel('Frequency (radians)','fontsize',12);
ylabel('Magnitude |H(w)|','fontsize',12);
title('Magnitude Plot for Analysis Filter');

figure(2)
plot(w0,hang1)
hold on
plot(w1,hang2)
legend('H0(w)','H1(w)')
grid on;
title('Phase Plot for Analysis Filter')
xlabel('Frequency (radians)')
ylabel('Phase (radians)')

figure(3)
[g0,w0] = freqz(G0,1);
scale = sqrt(2)/max(abs(g0));
G0 = scale*G0;
[g0,w0] = freqz(G0,1);
tmp = angle(g0);
tol = 0.85*pi;
hang3 = unwrap(tmp, tol);
plot(w0,abs(g0))
hold on
[g1,w1] = freqz(G1,1);
scale = sqrt(2)/max(abs(g1));
G1 = scale*G1;
[g1,w1] = freqz(G1,1);
tmp = angle(g1);
tol = 0.85*pi;
hang4 = unwrap(tmp, tol);
plot(w1,abs(g1))
legend('|G0(w)|','|G1(w)|')
grid on;
xlabel('Frequency (radians)','fontsize',12);
ylabel('Magnitude |G(w)|','fontsize',12);
title('Magnitude Plot for Synthesis Filter');

figure(4)
plot(w0,hang3)
hold on
plot(w1,hang4)
legend('G0(w)','G1(w)')
grid on;
title('Phase Plot for Synthesis Filter')
xlabel('Frequency (radians)')
ylabel('Phase (radians)')

T1 = conv(H0,G0);
T2 = conv(H1,G1);
T = 1/2 * (T1 + T2);
Tz = poly2sym(T,z);
Tz = sym2poly(Tz);
A1 = conv(G0,-G1);
A2 = conv(G1,G0);
A = 1/2 * (A1 + A2);
Az = poly2sym(A,z);
Az = sym2poly(Az);
