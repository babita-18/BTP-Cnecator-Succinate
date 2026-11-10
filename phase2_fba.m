%% =====================================================================
%  PHASE 2 : Wild-type FBA (fructose + CO2/H2) and maximum succinate
%  Run from the same folder as phase1 (C:\CnecatorGSM\Model)
%  =====================================================================
clc; clear; close all;
load(fullfile('phase1_out','iCN1361_loaded.mat'),'model');
changeCobraSolver('glpk','LP');
out = 'phase2_out'; if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'phase2_log.txt')); diary on;

BIO='R_Biomass'; FRU='EX_BETA_D_FRUCTOSE_e'; CO2='EX_CARBON_DIOXIDE_e';
H2 ='EX_HYDROGEN_MOLECULE_e'; O2='EX_OXYGEN_MOLECULE_e'; SUC='EX_SUC_e'; PHB='EX_PHB_e';

%% ---- 1. Two growth conditions --------------------------------------
mF = model;                                   % heterotrophic: fructose 2.6
mA = model;                                   % autotrophic: CO2 + H2 only
mA = changeRxnBounds(mA, FRU, 0,   'b');
mA = changeRxnBounds(mA, CO2, -10, 'l');
mA = changeRxnBounds(mA, H2,  -50, 'l');
conds = {mF, mA}; names = {'Fructose','CO2+H2'};

%% ---- 2. FBA with COBRA (optimizeCbModel) ---------------------------
fprintf('\n===== WILD-TYPE FBA (COBRA) =====\n');
mu = zeros(1,2); V = zeros(numel(model.rxns),2);
for k = 1:2
    s = optimizeCbModel(conds{k},'max');
    mu(k) = s.f; V(:,k) = s.x;
    fprintf('%-9s growth = %.6f h^-1  (status %d)\n', names{k}, s.f, s.stat);
end

%% ---- 3. Same FBA written directly with linprog ---------------------
fprintf('\n===== WILD-TYPE FBA (MATLAB linprog) =====\n');
if exist('linprog','file') == 2
    opt = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
    for k = 1:2
        m  = conds{k};
        f  = -m.c;                            % linprog minimises, so use -c
        Aeq = m.S;  beq = zeros(size(m.S,1),1);
        [v, fval, flag] = linprog(f, [], [], Aeq, beq, m.lb, m.ub, opt);
        fprintf('%-9s growth = %.6f h^-1  (exitflag %d)  diff vs COBRA = %.2e\n', ...
            names{k}, -fval, flag, abs(-fval - mu(k)));
    end
else
    fprintf('linprog NOT found -> install "Optimization Toolbox" from Home > Add-Ons\n');
end

%% ---- 4. Key fluxes table --------------------------------------------
keyR = {FRU,CO2,H2,O2,BIO,'R_ATPM','R_HYDROGENASE_SH_NAD','R_MBH_HYDROGENASE', ...
 'R_RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN','R_PHOSPHORIBULOKINASE_RXN','R_KDPGALDOL_RXN', ...
 'R_RXN0_1134','R_CITSYN_RXN','R_ISOCITDEH_RXN','R_ISOCITRATE_DEHYDROGENASE_NAD_RXN', ...
 'R_2OXOGLUTDECARB_RXN','R_SUCCCOASYN_RXN','R_SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN', ...
 'R_FUMHYDR_RXN','R_MALATE_DEH_RXN','R_ISOCIT_CLEAV_RXN','R_MALSYN_RXN','R_PEPCARBOX_RXN', ...
 'R_PYRUVATE_CARBOXYLASE_RXN','R_41132_RXN','R_ACETYL_COA_ACETYLTRANSFER_RXN','R_RXN_5901', ...
 'R_RXN1_42','R_31175_RXN',PHB,SUC};
idx = findRxnIDs(model, keyR);
T = table(keyR', V(idx,1), V(idx,2), 'VariableNames',{'Reaction','Fructose','CO2_H2'});
disp(T); writetable(T, fullfile(out,'WT_key_fluxes.csv'));

%% ---- 5. Max succinate (and PHB) at 50% of WT growth -----------------
fprintf('\n===== MAX THEORETICAL PRODUCT AT 50%% GROWTH =====\n');
for k = 1:2
    m = changeRxnBounds(conds{k}, BIO, 0.5*mu(k), 'l');
    s1 = optimizeCbModel(changeObjective(m, SUC),'max');
    s2 = optimizeCbModel(changeObjective(m, PHB),'max');
    fprintf('%-9s growth floor %.4f | max SUC = %.4f | max PHB = %.4f mmol/gDW/h\n', ...
        names{k}, 0.5*mu(k), s1.f, s2.f);
end

%% ---- 6. Production envelope (succinate vs growth) -------------------
figure('Color','w'); hold on; cols = {[0 0.45 0.74],[0.85 0.33 0.1]};
for k = 1:2
    g = linspace(0, mu(k), 15); smax = zeros(size(g)); smin = smax;
    for j = 1:numel(g)
        m = changeRxnBounds(conds{k}, BIO, g(j), 'b');
        a = optimizeCbModel(changeObjective(m,SUC),'max'); smax(j) = a.f;
        b = optimizeCbModel(changeObjective(m,SUC),'min'); smin(j) = b.f;
    end
    plot(g, smax, '-o','Color',cols{k},'LineWidth',2,'DisplayName',[names{k} ' max']);
    plot(g, smin, '--','Color',cols{k},'LineWidth',1.5,'DisplayName',[names{k} ' min']);
end
xlabel('Growth rate (h^{-1})'); ylabel('Succinate secretion (mmol gDW^{-1} h^{-1})');
title('Wild-type production envelope: succinate vs growth'); legend('Location','northeast');
set(gca,'Color','w','XColor','k','YColor','k','FontSize',13); grid on; box on;
saveas(gcf, fullfile(out,'Fig_production_envelope.png'));

%% ---- 7. S matrix plot again, white background -----------------------
figure('Color','w','Position',[100 100 850 750]);
spy(model.S); h = get(gca,'Children'); set(h,'Color',[0 0.3 0.7],'MarkerSize',2);
set(gca,'Color','w','XColor','k','YColor','k','FontSize',13);
xlabel('Reactions (n = 1292)'); ylabel('Metabolites (m = 1265)');
title('Stoichiometric matrix S of iCN1361 (nz = 5732)','Color','k');
saveas(gcf, fullfile(out,'Fig_S_matrix_spy.png'));

%% ---- 8. Small S sample: PHB vs succinate branch ---------------------
R = {'R_ACETYL_COA_ACETYLTRANSFER_RXN','R_RXN_5901','R_RXN1_42','R_31175_RXN',PHB, ...
     'R_ISOCIT_CLEAV_RXN','R_SUCCCOASYN_RXN','R_SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN',SUC};
M = {'M_ACETYL_COA[c]','M_CO_A[c]','M_ACETOACETYL_COA[c]','M_CPD_650[c]','M_Poly_Hydroxybutyrate[c]', ...
     'M_THREO_DS_ISO_CITRATE[c]','M_GLYOX[c]','M_SUC[c]','M_SUC_COA[c]','M_FUM[c]'};
sub = full(model.S(findMetIDs(model,M), findRxnIDs(model,R)));
ST = array2table(sub,'RowNames',M,'VariableNames',matlab.lang.makeValidName(R));
fprintf('\n===== S SUB-MATRIX (PHB vs succinate) =====\n'); disp(ST);
writetable(ST, fullfile(out,'S_sample_PHB_succinate.csv'),'WriteRowNames',true);
BT = ST; BT{:,:} = double(sub~=0);
fprintf('\n===== BINARY (0/1) VERSION =====\n'); disp(BT);
writetable(BT, fullfile(out,'B_sample_PHB_succinate.csv'),'WriteRowNames',true);

save(fullfile(out,'phase2_results.mat'),'mu','V','keyR');
diary off; fprintf('\nDONE. Upload all files from "%s" to Claude.\n', out);
