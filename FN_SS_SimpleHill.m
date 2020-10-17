%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Finding steady state solution     %
%           for the simple Hill model         %

% Created by Mariana Gómez-Schiavon
% September 2019

% FN_SS_SimpleHill : Finding steady state solution for the simple Hill
%                     model given a set of biophysical paramers and
%                     inducer (hormone) concentration.
%
%   Xs = FN_SS_SimpleHill(X,H,p)
%   X : Vector (array) of TF concentration ([nM])
%   H : Vector of inducer (hormone) concentration ([nM])
%   p : Structure with the kinetic parameters
%    .a : Related to basal synthesis of the naked promoter ([0,1])
%    .n : TF regulation nonlinearity (Hill coefficient)
%    .K : Related to TF-promoter dissociation rate ([nM])
%    .m : Maximum synthesis rate
%    .g : Degradation/dilution rate of the output ([1/min])
%
%   OUTPUT Ye : Matrix [length(H) x length(X)]
%
%   See also FN_SS_Mechanistic.m
%   See also FN_SS_HillxBasal.m
%   See also FN_FitError.m

function Ye = FN_SS_SimpleHill(X,H,p)
    OHFn = @(x,n,k,a) [a+((1-a)*(x.^n)./((x.^n)+(k^n)))];
    Ye = zeros(length(H),length(X));
    for h = 1:length(H)
        for i = 1:length(X)
            Ye(h,i) = p.m * OHFn(H(h)*X(i),p.n,p.K,p.a)/p.g;
        end
    end
    clear h i
end