%% =====================================================================
%  PHASE 5 : OptKnock (bilevel MILP)  outer: max succinate | inner: max growth
%  Run from C:\CnecatorGSM\Model.  Each run can take a few minutes.
%  If one run takes > 20 min, press Ctrl+C and send what was printed.
%  =====================================================================
clc; clear; close all;
load(fullfile('phase1_out','iCN1361_loaded.mat'),'model');
changeCobraSolver('glpk','LP'); changeCobraSolver('glpk','MILP');
out = 'phase5_out'; if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'phase5_log.txt')); diary on;

BIO='R_Biomass'; SUC='EX_SUC_e'; FRU='EX_BETA_D_FRUCTOSE_e';
mF = model;
mA = changeRxnBounds(model, FRU, 0, 'b');
mA = changeRxnBounds(mA, 'EX_CARBON_DIOXIDE_e', -10, 'l');
mA = changeRxnBounds(mA, 'EX_HYDROGEN_MOLECULE_e', -50, 'l');
conds = {mF, mA}; names = {'Fructose','CO2_H2'};

% Candidate reactions OptKnock may delete (central carbon metabolism)
cand = {'R_SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN','R_SUCCCOASYN_RXN','R_ISOCIT_CLEAV_RXN', ...
 'R_MALSYN_RXN','R_FUMHYDR_RXN','R_MALATE_DEH_RXN','R_2OXOGLUTDECARB_RXN','R_ISOCITDEH_RXN', ...
 'R_ISOCITRATE_DEHYDROGENASE_NAD_RXN','R_PEPCARBOX_RXN','R_PYRUVATE_CARBOXYLASE_RXN', ...
 'R_41132_RXN','R_RXN0_1134','R_ACETYL_COA_ACETYLTRANSFER_RXN','R_RXN_5901','R_RXN1_42', ...
 'R_PHOSACETYLTRANS_RXN','R_ACETATEKIN_RXN','R_MALIC_NADP_RXN','R_GABATRANSAM_RXN', ...
 'R_RXN0_5293_NAD','R_METHYLISOCITRATE_LYASE_RXN','R_NADPtranshydrogenase','R_PEPDEPHOS_RXN', ...
 'R_1TRANSKETO_RXN','R_2TRANSKETO_RXN','R_KDPGALDOL_RXN','R_PGLUCONDEHYDRAT_RXN'};
cand = cand(ismember(cand, model.rxns));
fprintf('Candidate reactions: %d\n', numel(cand));

allRes = {};
for k = 1:2
    m0 = conds{k}; muWT = optimizeCbModel(m0).f;
    constrOpt.rxnList = {BIO}; constrOpt.values = 0.3*muWT; constrOpt.sense = 'G';
    for K = 1:3
        prev = {};
        for alt = 1:3                                   % up to 3 different solutions
            opt.targetRxn = SUC; opt.numDel = K; opt.numDelSense = 'L'; opt.vMax = 1000;
            fprintf('\n[%s] OptKnock K<=%d, solution %d ... ', names{k}, K, alt); tic;
            try
                sol = OptKnock(m0, cand, opt, constrOpt, prev, false);
            catch ME
                fprintf('error: %s\n', ME.message); break;
            end
            fprintf('%.1f s\n', toc);
            if ~isfield(sol,'rxnList') || isempty(sol.rxnList), fprintf('  no new solution\n'); break; end
            ko = sol.rxnList(:)';
            % validate with plain FBA: growth, and min/max succinate at that growth
            mk = changeRxnBounds(m0, ko, 0, 'b');
            s = optimizeCbModel(mk); mu = s.f;
            mm = changeRxnBounds(mk, BIO, 0.999*mu, 'l');
            a = optimizeCbModel(changeObjective(mm,SUC),'max');
            b = optimizeCbModel(changeObjective(mm,SUC),'min');
            fprintf('  KO: %s\n  growth %.4f | SUC at max growth: min %.4f  max %.4f\n', ...
                strjoin(strrep(ko,'R_',''),' + '), mu, b.f, a.f);
            allRes(end+1,:) = {names{k}, K, alt, strjoin(ko,' + '), mu, b.f, a.f}; %#ok<SAGROW>
            prev{end+1} = ko; %#ok<SAGROW>
        end
    end
end
T = cell2table(allRes,'VariableNames',{'Condition','MaxKO','Alt','Knockouts','Growth','SUC_min','SUC_max'});
writetable(T, fullfile(out,'optknock_solutions.csv')); disp(T);

% Envelope: wild type vs best OptKnock design (per condition)
for k = 1:2
    Tk = T(strcmp(T.Condition,names{k}),:); if isempty(Tk), continue; end
    [~,i] = max(Tk.SUC_min); ko = strsplit(Tk.Knockouts{i},' + ');
    figure('Color','w'); hold on;
    for d = 1:2
        mk = conds{k}; if d==2, mk = changeRxnBounds(mk, ko, 0, 'b'); end
        mu = optimizeCbModel(mk).f; gg = linspace(0,mu,12); lo = 0*gg; hi = lo;
        for j = 1:numel(gg)
            mm = changeRxnBounds(mk,BIO,gg(j),'b');
            x1 = optimizeCbModel(changeObjective(mm,SUC),'max'); hi(j)=x1.f;
            x2 = optimizeCbModel(changeObjective(mm,SUC),'min'); lo(j)=x2.f;
        end
        c = [0 0.45 0.74; 0.85 0.33 0.1]; lab = {'Wild type','OptKnock design'};
        plot(gg,hi,'-','Color',c(d,:),'LineWidth',2,'DisplayName',[lab{d} ' (max)']);
        plot(gg,lo,'--','Color',c(d,:),'LineWidth',2,'DisplayName',[lab{d} ' (min)']);
    end
    xlabel('Growth rate (h^{-1})','Color','k'); ylabel('Succinate (mmol gDW^{-1} h^{-1})','Color','k');
    title(['OptKnock design vs wild type: ' strrep(names{k},'_','+')],'Color','k');
    lg = legend('Location','northeast'); set(lg,'Color','w','TextColor','k');
    set(gca,'Color','w','XColor','k','YColor','k','FontSize',12); grid on; box on;
    saveas(gcf, fullfile(out,['Fig_optknock_envelope_' names{k} '.png']));
end
diary off; fprintf('\nDONE. Upload all files from "%s" to Claude.\n', out);
