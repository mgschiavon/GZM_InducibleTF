%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Finding steady state solution     %
%           for the allosteric model         %

% Created by Mariana Gómez-Schiavon
% July 2020

% FN_SS_Allosteric : Finding steady state solution for the allosteric
%                     model given a set of biophysical paramers and
%                     inducer (hormone) concentration.
%
%   Xs = FN_SS_Allosteric(X,H,p)
%   X : Vector (array) of TF concentration
%   H : Vector of inducer (hormone) concentration
%   p : Structure with the kinetic parameters
%    .KX : Related to Hormone-TF (X:H) dissociation rate ([nM])
%    .b  : Related to basal activity of free (non-active) TF ([0,1])
%    .a : Related to basal synthesis of the naked promoter ([0,1])
%    .n : TF regulation nonlinearity (Hill coefficient)
%    .K : Related to TF-promoter dissociation rate ([nM])
%    .m : Maximum synthesis rate
%    .g : Degradation/dilution rate of the output ([1/min])
%
%   OUTPUT Ye : Matrix [length(H) x length(X)]
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_SS_Mechanistic.m
%   See also FN_FitError.m

function Ye = FN_SS_Allosteric(X,H,p)
    OHFn = @(x,n,k,a) [a+((1-a)*(x.^n)./((x.^n)+(k^n)))];
    Ye = zeros(length(H),length(X));
    for h = 1:length(H)
        for i = 1:length(X)
            XT = X(i);
            Xa = roots([1,-(H(h)+XT+p.KX),H(h)*XT]);
            Xa = Xa([Xa<XT]);
            if(length(Xa)~=1)
                'error -- multiple solutions'
                Xa = NaN;
            end
            Ye(h,i) = p.m * OHFn(Xa+(p.b*(XT-Xa)),p.n,p.K,p.a)/p.g;
        end
    end
    clear h i
end