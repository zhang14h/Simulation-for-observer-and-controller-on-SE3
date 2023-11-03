clc
clear 
close all


load ('C:\Users\black\OneDrive\Desktop\filefrommac\localo+localc.csv');
load ('C:\Users\black\OneDrive\Desktop\filefrommac\localo+semic.csv');
load ('C:\Users\black\OneDrive\Desktop\filefrommac\semio+localc.csv');
load ('C:\Users\black\OneDrive\Desktop\filefrommac\semio+semic.csv');



y = [2 2.5 2.5 1.5 4.3 3.4 2.7 1.9;
     0 1.2 -1.1 -1.0 -0.4 0.8 -2.1 -1.6;
     0.2 0.35 0.25 0.2 0.4 0.6 1 1.2;
     1 1 1 1 1 1 1 1];
 
T = eye(4); 
