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
    X = mean(xd.(ExID).X,1)*10;
    H = xd.(ExID).H;
    D = xd.(ExID).Y;
    clear xd
S = [1:1000];     % Random number seed(s) (1 per run)
I = 20000;      % Iterations per fitting run
printAll = 0;   % Flag for printing full random walk
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
    f(i).lim = [1,10];
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
    f(i).cov = 1e-6;
    f(i).lim = [0,0.2];
    i = i + 1;
    f(i).par = 'Im';
    f(i).cov = 1e-4;
    f(i).lim = [0,2];
    clear i

%% Run fitting:
bestP = zeros(length(S),length(f));
minE  = zeros(length(S),1);
for s = S
    cat(2,'Running seed #',num2str(s))
    [bP,mE] = FN_FitMRW(X,H,p,D,s,f,I,ExID,printAll)
    bestP(s,:) = bP;
    minE(s) = mE;
    if(mod(s,10)==0)
        save('TEMP_MRWs.mat');
    end
end
clear s bP mE
save(cat(2,'MRW_',ExID,'.mat'));

%% Figures
if(printAll)
    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [1 0 18 10];
    fig.Position = fig.PaperPosition;
    C = colormap('parula');
    for s = S
        load(cat(2,'MRW_',ExID,'_s',num2str(s),'.mat'),'mrw');
        for i = 1:size(mrw.P,2)
            subplot(2,4,i)
            hold on;
            plot(mrw.P(:,i),'LineWidth',2,'Color',C(s*6,:))
                xlabel('Iterations')
                ylabel(f(i).par)
                xlim([0,I])
                set(gca,'YScale','log')
                box on
        end
        subplot(2,4,8)
        hold on;
        plot(mrw.e,'LineWidth',2,'Color',C(s*6,:))
            xlabel('Iterations')
            ylabel('Error')
            xlim([0,I])
            set(gca,'YScale','log')
            box on
    end
    clear s i a b
    print(gcf,cat(2,'MRW_',ExID,'_Runs.png'),'-dpng','-r300')

    fig = figure();
    fig.Units = 'inches';
    fig.PaperPosition = [0 0 18 8];
    fig.Position = fig.PaperPosition;
    C = colormap('jet');
    Hi = log2(H); Hi(length(H)) = Hi(length(Hi)-1) - 1;
    for s = S
        for i = 1:length(f)
            p.(f(i).par) = bestP(s,i);
        end
        Ye = FN_SS_Mechanistic(X,H,p);
        subplot(2,5,s)
        hold on
        for i = 1:length(X)
            plot(Hi,Ye(:,i),'Color',C(i*10,:),'LineWidth',2)
            plot(Hi,D(:,i),'Color',C(i*10,:),'LineStyle','none','Marker','o')
        end
                xlabel('Hormone')
                ylabel('Output')
                title(cat(2,'Error = ',num2str(FN_FitError(X,H,p,D))))
                xlim([min(Hi)-0.5 max(Hi)+0.5])
                set(gca,'YScale','log','YGrid','on',...
                    'XTick',Hi([length(H):-3:1]),'XTickLabel',H([length(H):-3:1]),...
                    'XTickLabelRotation',45)
                box on
    end
    clear Hi s i 
    print(gcf,cat(2,'MRW_',ExID,'_BestFits.png'),'-dpng','-r300')
end

%% END