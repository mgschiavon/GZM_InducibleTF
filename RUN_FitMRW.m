%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% RUN: Find parameters that best fit the data %

% Created by Mariana Gómez-Schiavon
% August 2019

%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m
%   See also FN_FitError.m
%   See also FN_FitMRW.m

clear

%% Inputs & conditions
% Data
ExID = 'GEMc';	% Experiment (TF) to consider
    load('DATA_synTF.mat','xd');
    X = mean(xd.(ExID).X,1);
    H = xd.(ExID).H;
    D = xd.(ExID).Y;
    clear xd
S = [1:10];     % Random number seed(s) (1 per run)
I = 10000;      % Iterations per fitting run
% Kinetic parameters:
    p.mY = 0.997;
    p.kb = 0.6;
    p.nO = 1.6;
    p.KO = 0.99;
    p.KX = 15;
    p.aX = 0.00065;
    p.Im = 0.6;
    p.gY = 0.01;
% Parameters to fit:
    i = 1;
    f(i).par = 'mY';
    f(i).cov = 1e-5;
    f(i).lim = [0,1];
    i = i + 1;
    f(i).par = 'kB';
    f(i).cov = 1e-4;
    f(i).lim = [0,10];
    i = i + 1;
    f(i).par = 'nO';
    f(i).cov = 1e-4;
    f(i).lim = [0,10];
    i = i + 1;
    f(i).par = 'KO';
    f(i).cov = 1e-3;
    f(i).lim = [0,100];
    i = i + 1;
    f(i).par = 'KX';
    f(i).cov = 1e-3;
    f(i).lim = [0,100];
    i = i + 1;
    f(i).par = 'aX';
    f(i).cov = 1e-5;
    f(i).lim = [0,1];
    i = i + 1;
    f(i).par = 'nO';
    f(i).cov = 1e-4;
    f(i).lim = [0,10];
    clear i

%% Run fitting:
for s = S
    cat(2,'Running seed #',num2str(s))
    FN_FitMRW(X,H,p,D,s,f,I,ExID)
end

%% END