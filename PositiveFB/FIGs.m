clear;
sim.mm = 'FSMv1';       % Label for motif file
sim.ex = 'Ex02';        % Label for parameters file
sim.an = 'ExSSs';       % Chose analysis type (Options: ExSSs)

%%
x = importdata(cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'.txt'),'\t',1);
HT = x.data(:,1);
XT0 = x.data(:,2);
Xi = x.data(:,3);
Xa = x.data(:,4);
H  = x.data(:,5);
clear x

%%
figure();
hold on;
    I = [XT0==min(XT0)];
    plot(H(I),Xi(I)+Xa(I),'LineWidth',2,...
        'DisplayName',cat(2,'X_{T0} = ',num2str(min(XT0),2)))
    I = [XT0==max(XT0)];
    plot(H(I),Xi(I)+Xa(I),'LineWidth',2,'LineStyle','--',...
        'DisplayName',cat(2,'X_{T0} = ',num2str(max(XT0),2)))
        legend('Location','northwest')
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        xlabel('Hormone')
        ylabel('Total X')
        xlim(10.^[-3 3])
%         ylim(10.^[-2 1])
        title(cat(2,'ZPM | ',sim.mm))
        box on
        print(gcf,cat(2,'OUT_',sim.an,'_',sim.mm,'_',sim.ex,'.png'),'-dpng','-r300')
    