% y = sampdata;
%
%   This function returns a sample data sequence to be used as a test input
%       for digital signal processing problems.
%   The input data is filtered with a low pass filter to make it band
%       limited (cutoff frequency is set at 0.95*pi).

function y = sampdata

ntaps = 65;
f = [0.0 0.9 0.95 1.0];
mag = [ 1.0 1.0 0.7071 0.0];
b = fir2(ntaps, f, mag);
n1 = length(b);
len1 = 256 - n1 + 1;
%
%  Generate the input sample sequence.
%
data = 5*[zeros(1, 24) ones(1, 48) zeros(1, 48) -1*ones(1, 48) zeros(1,23)];
y = conv(b, data);
return
