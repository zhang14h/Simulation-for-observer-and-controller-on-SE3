clc;
clear;
close all;

fs = 10^(3);

t = 0:1/fs:2;
%     n(i) = randi([1,8]);
for i = 1:length(t)
%     x(i) = randi([1,2])*(-7+13*rand(1));
    x(i) =1;
end;
% x = vco(sin(t));
% x(1000) = 34;
% x(1500) = 34;

stft(x,fs,'window',kaiser(100));

view (-40,40);
colormap jet
