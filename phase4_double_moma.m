%% =====================================================================
%  PHASE 4 : (a) corrected single-KO summary
%            (b) double-gene knockouts (all pairs of 16 targets) - FBA + MOMA
%            (c) best multi-gene designs + production envelopes
%  Run from C:\CnecatorGSM\Model.  ~10-30 min.
%  =====================================================================
clc; clear; close all;
load(fullfile('phase1_out','iCN1361_loaded.mat'),'model');
changeCobraSolver('glpk','LP');
out = 'phase4_out'; if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'phase4_log.txt')); diary on;

BIO='R_Biomass'; SUC='EX_SUC_e'; FRU='EX_BETA_D_FRUCTOSE_e';
mF = model;
mA = changeRxnBounds(model, FRU, 0, 'b');
mA = changeRxnBounds(mA, 'EX_CARBON_DIOXIDE_e', -10, 'l');
mA = changeRxnBounds(mA, 'EX_HYDROGEN_MOLECULE_e', -50, 'l');
conds = {mF, mA}; names = {'Fructose','CO2_H2'};
iB = findRxnIDs(model,BIO); iS = findRxnIDs(model,SUC);

%% ---- (a) corrected single-KO summary from Phase 3 files ------------
fprintf('===== SINGLE-KO SUMMARY (from phase3) =====\n');
for k = 1:2
    T = readtable(fullfile('phase3_out',['succinate_screen_' names{k} '.csv']));
    nl = T(~T.Lethal & ~isnan(T.SUC_guaranteed),:);
    fprintf('%-9s lethal %d | non-lethal with effect %d | no effect %d | guaranteed SUC>0: %d\n', ...
        names{k}, sum(T.Lethal), height(nl), height(T)-sum(T.Lethal)-height(nl), sum(nl.SUC_guaranteed>1e-6));
end

%% ---- (b) target panel and all pairs --------------------------------
tg = {'H16_A2631','SDH';      'H16_A0547','SucCD';    'H16_A2227','ICL-1';
      'H16_A2211','ICL-2';    'H16_A2217','MS';       'H16_A1437','PhaC1';
      'H16_A2003','PhaC2';    'H16_A1439','PhaB1';    'H16_B1871','PTA-1';
      'H16_B1631','PTA-2';    'H16_A0670','ACK-1';    'H16_B1630','ACK-2';
      'H16_A2325','OGDH';     'H16_A3711','PEPCK';    'H16_A1251','PYC';
      'H16_A2921','PPC'};
pairs = nchoosek(1:size(tg,1),2);
fprintf('\nTargets %d -> %d gene pairs\n', size(tg,1), size(pairs,1));

for k = 1:2
    m0 = conds{k}; wt = optimizeCbModel(m0,'max','one'); muWT = wt.f;
    R = cell(size(pairs,1),8);
    for p = 1:size(pairs,1)
        g = tg(pairs(p,:),1)'; lab = strjoin(tg(pairs(p,:),2)','+');
        mk = deleteModelGenes(m0, g);
        s = optimizeCbModel(mk); muKO = 0; if s.stat==1, muKO = s.f; end
        sMax = NaN; sMin = NaN; momaG = NaN; momaS = NaN;
        if muKO > 0.01*muWT
            a = optimizeCbModel(changeObjective(changeRxnBounds(mk,BIO,0.5*muWT,'l'),SUC),'max');
            b = optimizeCbModel(changeObjective(changeRxnBounds(mk,BIO,0.999*muKO,'l'),SUC),'min');
            if a.stat==1, sMax=a.f; end; if b.stat==1, sMin=b.f; end
            try
                [sd,~,~,st] = linearMOMA(m0, mk, 'max', false, false);
                if st==1, momaG = sd.x(iB); momaS = sd.x(iS); end
            catch, end
        end
        R(p,:) = {g{1}, g{2}, lab, muKO, sMax, sMin, momaG, momaS};
    end
    T = cell2table(R,'VariableNames',{'Gene1','Gene2','Label','FBA_growth','SUC_max_50pct', ...
        'SUC_guaranteed','MOMA_growth','MOMA_SUC'});
    T = sortrows(T, {'SUC_guaranteed','MOMA_SUC'}, {'descend','descend'}, 'MissingPlacement','last');
    writetable(T, fullfile(out,['double_KO_' names{k} '.csv']));
    fprintf('\n===== DOUBLE KO: %s (WT growth %.4f) - top 15 =====\n', names{k}, muWT);
    disp(T(1:15,3:end));
    fprintf('Pairs lethal: %d | pairs with guaranteed SUC > 0: %d | pairs with MOMA SUC > 0: %d\n', ...
        sum(T.FBA_growth < 0.01*muWT), sum(T.SUC_guaranteed > 1e-6), sum(T.MOMA_SUC > 1e-6));
end

%% ---- (c) best multi-gene designs -----------------------------------
D = { {},                                                      'Wild type';
      {'H16_A2631','H16_A0547'},                               'SDH+SucCD';
      {'H16_A2631','H16_A0547','H16_B1871','H16_B1631'},       'SDH+SucCD+PTA(both)';
      {'H16_A2631','H16_A0547','H16_A0670','H16_B1630'},       'SDH+SucCD+ACK(both)';
      {'H16_A2631','H16_A0547','H16_A1437','H16_A2003'},       'SDH+SucCD+PhaC(both)' };
fprintf('\n===== MULTI-GENE DESIGNS =====\n');
for k = 1:2
    m0 = conds{k}; muWT = optimizeCbModel(m0).f;
    fprintf('\n-- %s --\n%-24s %8s %10s %10s %10s %10s\n', names{k}, 'Design','FBA_mu','SUCmax50','SUCguar','MOMA_mu','MOMA_SUC');
    figure('Color','w'); hold on; C = lines(size(D,1));
    for d = 1:size(D,1)
        mk = m0; if ~isempty(D{d,1}), mk = deleteModelGenes(m0, D{d,1}); end
        s = optimizeCbModel(mk); muKO = 0; if s.stat==1, muKO = s.f; end
        a = optimizeCbModel(changeObjective(changeRxnBounds(mk,BIO,min(0.5*muWT,muKO),'l'),SUC),'max');
        b = optimizeCbModel(changeObjective(changeRxnBounds(mk,BIO,0.999*muKO,'l'),SUC),'min');
        mg = NaN; ms = NaN;
        if ~isempty(D{d,1})
            try, [sd,~,~,st] = linearMOMA(m0, mk, 'max', false, false);
                 if st==1, mg = sd.x(iB); ms = sd.x(iS); end, catch, end
        end
        fprintf('%-24s %8.4f %10.4f %10.4f %10.4f %10.4f\n', D{d,2}, muKO, a.f, b.f, mg, ms);
        % production envelope
        gg = linspace(0, muKO, 12); lo = zeros(size(gg)); hi = lo;
        for j = 1:numel(gg)
            mm = changeRxnBounds(mk, BIO, gg(j), 'b');
            x1 = optimizeCbModel(changeObjective(mm,SUC),'max'); hi(j) = x1.f;
            x2 = optimizeCbModel(changeObjective(mm,SUC),'min'); lo(j) = x2.f;
        end
        plot(gg,hi,'-','Color',C(d,:),'LineWidth',2,'DisplayName',[D{d,2} ' max']);
        plot(gg,lo,'--','Color',C(d,:),'LineWidth',2,'HandleVisibility','off');
    end
    xlabel('Growth rate (h^{-1})','Color','k'); ylabel('Succinate (mmol gDW^{-1} h^{-1})','Color','k');
    title(['Production envelopes (dashed = minimum): ' strrep(names{k},'_','+')],'Color','k');
    lg = legend('Location','northeast'); set(lg,'Color','w','TextColor','k');
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',12); grid on; box on;
    saveas(gcf, fullfile(out,['Fig_envelopes_designs_' names{k} '.png']));
end
diary off; fprintf('\nDONE. Upload all files from "%s" to Claude.\n', out);
