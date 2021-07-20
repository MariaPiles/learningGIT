
clear;clc;

% True system come from an ARMA process, not an AR!
%A_true=[0.5];
%B_true=[1];
%A_true=[0.5665 0.32 0.4839];
%B_true=[1 0.5432 -0.7431 -0.98542];
A_true=[1 0.5 0.3];
B_true=[1 0.5];

% Data
N = 200;
x=randn(N,1);
Pn=0.01; 
signal = filter(B_true,[1 A_true],x) + sqrt(Pn)*randn(size(x));

% Identify the system with the gamma-filter
mu = 0.1;  % [0,1] where mu=0 is the FIR=AR filter and mu=0 is IIR
x = signal(1:N-1); y = signal(2:N);
[W,mubest,M] = gamma_ls(x,y,x,y,1)
