% Homework10 Question 4
clc 
clear all
close all
h0 = [0.0106   -0.0329   -0.0308    0.1870    0.0280   -0.6309   -0.7148   -0.2304];
g0 = [-0.2304   -0.7148   -0.6309    0.0280    0.1870   -0.0308   -0.0329    0.0106];
h1 = [-0.2304    0.7148   -0.6309   -0.0280    0.1870    0.0308   -0.0329   -0.0106];
g1 = [-0.0106   -0.0329    0.0308    0.1870   -0.0280   -0.6309    0.7148   -0.2304];
x = sampdata();
[v000,v001,v01,v10,v11] = part1(x,h0,h1);
y = part2(v000,v001,v01,v10,v11,g0,g1);

figure(1)
subplot(2,1,1)
plot(1:1:length(x),x)
subplot(2,1,2)
plot(1:1:length(y),y)
grid on;

function [v0,v1] = analyse(h0,h1,x)
v01 = conv(x,h0);
v0 = downsample(v01,2);
v11 = conv(x,h1);
v1 = downsample(v11,2);
end

function[y]= synth(x1,x2,g0,g1)
y1 = upsample(x1,2);
y2 = conv(y1,g0);
y3 = upsample(x2,2);
y4 = conv(y3,g1);
y = y4 + y2;
end

function[v000,v001,v01,v10,v11] = part1(x,h0,h1)
 [a1,b1] = analyse(h0,h1,x);
 [a2,v01] = analyse(h0,h1,a1);
 [v000,v001] = analyse(h0,h1,a2);
 [v10,v11] = analyse(h0,h1,b1);
end

function[y] = part2(v000,v001,v01,v10,v11,g0,g1)
  c1 = synth(v000,v001,g0,g1);
  d1 = synth(c1,[zeros(1,(length(c1)-length(v01))) v01],g0,g1);
  d2 = synth(v10,v11,g0,g1);
  y = synth(d1,[zeros(1,(length(d1)-length(d2))) d2],g0,g1);
end