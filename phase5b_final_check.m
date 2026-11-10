%% PHASE 5b : gene-level check of the OptKnock design + envelope figure
clc; clear; close all;
load(fullfile('phase1_out','iCN1361_loaded.mat'),'model');
changeCobraSolver('glpk','LP');
out = 'phase5_out'; if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'phase5b_log.txt')); diary on;
BIO='R_Biomass'; SUC='EX_SUC_e'; iB=findRxnIDs(model,BIO); iS=findRxnIDs(model,SUC);
mF = model;
mA = changeRxnBounds(model,'EX_BETA_D_FRUCTOSE_e',0,'b');
mA = changeRxnBounds(mA,'EX_CARBON_DIOXIDE_e',-10,'l');
mA = changeRxnBounds(mA,'EX_HYDROGEN_MOLECULE_e',-50,'l');
conds = {mF,mA}; names = {'Fructose','CO2_H2'};

D = { {},                                        'Wild type';
      {'H16_A2631','H16_A0670','H16_B1630'},     'SDH + ACK(both)  [3 genes]';
      {'H16_A2631','H16_B1871','H16_B1631'},     'SDH + PTA(both)  [3 genes]';
      {'H16_A2631','H16_A0547'},                 'SDH + SucCD      [2 genes]' };
for k = 1:2
    m0 = conds{k};
    fprintf('\n===== %s =====\n%-28s %8s %9s %9s %9s %9s\n', names{k}, ...
        'Design','Growth','SUCmin','SUCmax','MOMA_mu','MOMA_SUC');
    f = figure('Color','w','Visible','off'); hold on; C = lines(size(D,1));
    for d = 1:size(D,1)
        mk = m0; if ~isempty(D{d,1}), mk = deleteModelGenes(m0, D{d,1}); end
        mu = optimizeCbModel(mk).f;
        mm = changeRxnBounds(mk,BIO,0.999*mu,'l');
        a = optimizeCbModel(changeObjective(mm,SUC),'max');
        b = optimizeCbModel(changeObjective(mm,SUC),'min');
        mg = NaN; ms = NaN;
        if ~isempty(D{d,1})
            try, [sd,~,~,st] = linearMOMA(m0,mk,'max',false,false);
                 if st==1, mg = sd.x(iB); ms = sd.x(iS); end, catch, end
        end
        fprintf('%-28s %8.4f %9.4f %9.4f %9.4f %9.4f\n', D{d,2}, mu, b.f, a.f, mg, ms);
        gg = linspace(0,mu,12); lo = 0*gg; hi = lo;
        for j = 1:numel(gg)
            q = changeRxnBounds(mk,BIO,gg(j),'b');
            x1 = optimizeCbModel(changeObjective(q,SUC),'max'); hi(j)=x1.f;
            x2 = optimizeCbModel(changeObjective(q,SUC),'min'); lo(j)=x2.f;
        end
        plot(gg,hi,'-','Color',C(d,:),'LineWidth',2,'DisplayName',strtrim(D{d,2}));
        plot(gg,lo,'--','Color',C(d,:),'LineWidth',2,'HandleVisibility','off');
    end
    xlabel('Growth rate (h^{-1})'); ylabel('Succinate (mmol gDW^{-1} h^{-1})');
    title(['Envelopes, solid = max, dashed = min: ' strrep(names{k},'_','+')]);
    legend('Location','northeast'); grid on; box on; set(gca,'FontSize',12);
    try
        exportgraphics(f, fullfile(out,['Fig_final_envelope_' names{k} '.png']),'Resolution',200);
    catch
        print(f, fullfile(out,['Fig_final_envelope_' names{k}]), '-dpng', '-r200');
    end
    close(f);
end
diary off; fprintf('\nDONE.\n');
