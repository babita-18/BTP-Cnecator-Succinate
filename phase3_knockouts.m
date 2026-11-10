%% =====================================================================
%  PHASE 3 : Single-gene knockout screen (FBA)
%   (a) lethal / essential genes   (b) succinate screen for every gene
%  Run from C:\CnecatorGSM\Model.  Takes ~5-20 minutes.
%  =====================================================================
clc; clear; close all;
load(fullfile('phase1_out','iCN1361_loaded.mat'),'model');
changeCobraSolver('glpk','LP');
out = 'phase3_out'; if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'phase3_log.txt')); diary on;

BIO='R_Biomass'; FRU='EX_BETA_D_FRUCTOSE_e'; CO2='EX_CARBON_DIOXIDE_e';
H2='EX_HYDROGEN_MOLECULE_e'; SUC='EX_SUC_e';

mF = model;
mA = changeRxnBounds(model, FRU, 0, 'b');
mA = changeRxnBounds(mA, CO2, -10, 'l');
mA = changeRxnBounds(mA, H2,  -50, 'l');
conds = {mF, mA}; names = {'Fructose','CO2_H2'};
genes = model.genes(~strcmp(model.genes,'Spontaneous'));
fprintf('Genes screened: %d\n', numel(genes));

%% ---- 0. linprog check + loop-free (parsimonious) WT fluxes ----------
for k = 1:2
    s = optimizeCbModel(conds{k},'max','one');   % 'one' = pFBA, removes loops
    fprintf('%-9s pFBA growth %.6f | ICDH-NADP %.3f | ICDH-NAD %.3f\n', names{k}, s.f, ...
        s.x(findRxnIDs(model,'R_ISOCITDEH_RXN')), s.x(findRxnIDs(model,'R_ISOCITRATE_DEHYDROGENASE_NAD_RXN')));
    vWT{k} = s.x;
end
if exist('linprog','file') == 2
    opt = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
    for k = 1:2
        m = conds{k};
        [~, fval, flag] = linprog(-m.c, [], [], m.S, zeros(size(m.S,1),1), m.lb, m.ub, opt);
        fprintf('linprog %-9s growth = %.6f h^-1 (exitflag %d)\n', names{k}, -fval, flag);
    end
else
    fprintf('linprog still missing -> install Optimization Toolbox\n');
end
save(fullfile(out,'WT_pFBA.mat'),'vWT');

%% ---- (a) Essential genes ------------------------------------------
for k = 1:2
    [grRatio, grKO, grWT] = singleGeneDeletion(conds{k}, 'FBA', genes, false);
    grRatio(isnan(grRatio)) = 0;
    lethal = grRatio < 0.01;
    sick   = grRatio >= 0.01 & grRatio < 0.9;
    fprintf('\n===== ESSENTIALITY: %s =====\n', names{k});
    fprintf('WT growth %.6f | lethal %d | reduced growth (<90%%) %d | no effect %d\n', ...
        grWT, sum(lethal), sum(sick), sum(grRatio>=0.9));
    T = table(genes, grKO, grRatio, lethal, 'VariableNames',{'Gene','KO_growth','Ratio','Lethal'});
    writetable(T, fullfile(out, ['essentiality_' names{k} '.csv']));
    E{k} = genes(lethal);
end
both = intersect(E{1},E{2}); onlyA = setdiff(E{2},E{1}); onlyF = setdiff(E{1},E{2});
fprintf('\nLethal in both: %d | only autotrophic: %d | only fructose: %d\n', numel(both), numel(onlyA), numel(onlyF));
fprintf('Paper benchmark (fructose): 156 predicted essential genes (Pearcy et al. 2022)\n');
fprintf('Only-autotrophic lethal genes:\n'); disp(onlyA');
writetable(table(onlyA,'VariableNames',{'Gene'}), fullfile(out,'lethal_only_autotrophic.csv'));

%% ---- (b) Succinate screen: every non-lethal gene -------------------
for k = 1:2
    m0 = conds{k}; muWT = optimizeCbModel(m0).f;
    res = nan(numel(genes),4);
    for i = 1:numel(genes)
        [mk, hasEffect] = deleteModelGenes(m0, genes{i});
        if ~hasEffect, res(i,:) = [muWT 0 NaN 0]; continue; end
        s = optimizeCbModel(mk); muKO = s.f;
        if s.stat ~= 1 || muKO < 0.01*muWT, res(i,:) = [0 NaN NaN 1]; continue; end
        mf  = changeRxnBounds(mk, BIO, 0.5*muWT, 'l');           % 50% growth floor
        a = optimizeCbModel(changeObjective(mf,SUC),'max');
        mc  = changeRxnBounds(mk, BIO, 0.999*muKO, 'l');          % at the KO's own max growth
        b = optimizeCbModel(changeObjective(mc,SUC),'min');       % guaranteed succinate
        sMax = NaN; if a.stat==1, sMax = a.f; end
        sMin = NaN; if b.stat==1, sMin = b.f; end
        res(i,:) = [muKO sMax sMin 0];
        if mod(i,200)==0, fprintf('  %s: %d/%d genes done\n', names{k}, i, numel(genes)); end
    end
    T = table(genes, res(:,1), res(:,2), res(:,3), logical(res(:,4)), ...
        'VariableNames',{'Gene','KO_growth','SUC_max_at50pct','SUC_guaranteed','Lethal'});
    writetable(T, fullfile(out, ['succinate_screen_' names{k} '.csv']));
    ok = T(~T.Lethal & ~isnan(T.SUC_max_at50pct),:);
    ok = sortrows(ok, {'SUC_guaranteed','SUC_max_at50pct'}, {'descend','descend'});
    fprintf('\n===== TOP 15 SINGLE KNOCKOUTS FOR SUCCINATE: %s =====\n', names{k}); disp(ok(1:min(15,end),:));
    fprintf('Genes giving guaranteed succinate > 0: %d\n', sum(ok.SUC_guaranteed > 1e-6));
    top{k} = ok;
end

%% ---- (c) Literature targets checked one by one ---------------------
lit = {'H16_A1437','phaC1'; 'H16_A1438','phaA'; 'H16_A1439','phaB1'; 'H16_A2003','phaC2';
       'H16_A2629','sdhC'; 'H16_A2630','sdhD'; 'H16_A2631','sdhA'; 'H16_A2632','sdhB';
       'H16_A0547','sucC'; 'H16_A0548','sucD'; 'H16_A2227','aceA1'; 'H16_A2211','aceA2'};
fprintf('\n===== LITERATURE TARGETS (single KO) =====\n');
for k = 1:2
    T = readtable(fullfile(out, ['succinate_screen_' names{k} '.csv']));
    fprintf('\n-- %s --\n', names{k});
    for j = 1:size(lit,1)
        r = find(strcmp(T.Gene, lit{j,1}));
        if isempty(r), fprintf('%-10s %-6s not in model\n', lit{j,1}, lit{j,2}); continue; end
        fprintf('%-10s %-6s growth %.4f | SUCmax@50%% %.4f | SUCguaranteed %.4f\n', ...
            lit{j,1}, lit{j,2}, T.KO_growth(r), T.SUC_max_at50pct(r), T.SUC_guaranteed(r));
    end
end

%% ---- (d) Figures ---------------------------------------------------
for k = 1:2
    ok = top{k}; n = min(10,height(ok));
    figure('Color','w'); barh(flipud(ok.SUC_max_at50pct(1:n)),'FaceColor',[0 0.45 0.74]);
    set(gca,'YTick',1:n,'YTickLabel',flipud(ok.Gene(1:n)),'Color','w','XColor','k','YColor','k','FontSize',12);
    xlabel('Max succinate at 50% growth (mmol gDW^{-1} h^{-1})','Color','k');
    title(['Top 10 single-gene knockouts: ' strrep(names{k},'_','+')],'Color','k'); grid on;
    saveas(gcf, fullfile(out, ['Fig_top10_' names{k} '.png']));
end
diary off; fprintf('\nDONE. Upload all files from "%s" to Claude.\n', out);
