%%   [GZM] Inducible Transcription Factors   %%
% ------------------------------------------- %
% FUNCTION: Finding steady state solution     %
%           for the mechanistic model         %

% Created by Mariana Gómez-Schiavon
% August 2019

% FN_SS_Mechanistic : Finding steady state solution for the mechanistic
%                     model given a set of biophysical paramers and
%                     inducer (hormone) concentration.
%
%   Xs = FN_SS_Mechanistic(X,H,p)
%   X : Vector (array) of TF concentration
%   H : Vector of inducer (hormone) concentration
%   p : Structure with the kinetic parameters
%    .mY : Related to basal synthesis of the naked promoter ([0,1])
%    .kb : Related to efficiency rate of the transcription factor
%    .nO : Promoter occupancy nonlinearity (Hill coefficient)
%    .KO : Related to TF-promoter dissociation rate ([nM])
%    .KX : Related to Hormone-TF (X:H) dissociation rate ([nM])
%    .aX : Related to basal activity of free (non-active) TF ([0,1])
%    .Im : Maximum synthesis rate given the gene & translocation rate
%    .gY : Degradation/dilution rate of the output ([1/min])
%
%   OUTPUT Ye : Matrix [length(H) x length(X)]
%
%   See also FN_SS_SimpleHill.m
%   See also FN_SS_HillxBasal.m
%   See also FN_FitError.m

function Ye = FN_SS_Mechanistic(X,H,p)
    OHFn = @(x,n,k,a) [a+((1-a)*(x.^n)./((x.^n)+(k^n)))];
    OImx = @(x,I,b,m) [I*(1-(m*exp(-b*x/I)))];
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
            Xo = OHFn(Xa+(p.aX*(XT-Xa)),p.nO,p.KO,0);
            Ye(h,i) = OImx(Xo,p.Im*p.nM,p.kb,p.mY)/p.gY;
        end
    end
    clear h i
end