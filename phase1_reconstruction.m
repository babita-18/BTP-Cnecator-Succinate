%% =====================================================================
%  PHASE 1 : Network reconstruction + Stoichiometric matrix (Section 4.1)
%  Organism : Cupriavidus necator H16      Model : iCN1361
%  Run this file once. It writes everything into the folder "phase1_out".
%  Then upload the files from phase1_out back to Claude.
%  =====================================================================
clc; clear;

%% ---- 0. Settings: EDIT THESE TWO LINES ONLY -------------------------
modelFile = 'iCN1361.xml';   % <-- exact name of the .xml in CnecatorGSM/Model
solver    = 'glpk';          % 'glpk' comes with COBRA. Use 'gurobi' if installed
%  ---------------------------------------------------------------------

outDir = 'phase1_out';
if ~exist(outDir,'dir'), mkdir(outDir); end
diary(fullfile(outDir,'phase1_log.txt')); diary on;

%% ---- 1. Start COBRA and load the model ------------------------------
initCobraToolbox(false);
changeCobraSolver(solver,'LP');
model = readCbModel(modelFile);
fprintf('\n===== MODEL LOADED: %s =====\n', modelFile);

% ---- 1b. Repair missing fields (this SBML comes from ScrumPy) --------
fprintf('Fields present in model: %s\n', strjoin(fieldnames(model)', ', '));
nR = numel(model.rxns);
if ~isfield(model,'rxnNames'), model.rxnNames = model.rxns; end
if ~isfield(model,'metNames'), model.metNames = model.mets; end
if ~isfield(model,'genes'),    model.genes = {}; end
if ~isfield(model,'grRules')
    if isfield(model,'rules') && ~isempty(model.genes)
        try, model = creategrRulesField(model); catch, end
    end
end
if ~isfield(model,'grRules'), model.grRules = repmat({''}, nR, 1); end
if ~isfield(model,'rxnGeneMat') && ~isempty(model.genes)
    try, model = buildRxnGeneMat(model); catch, end
end
if ~isfield(model,'rxnGeneMat'), model.rxnGeneMat = sparse(nR, max(1,numel(model.genes))); end
fprintf('Genes found in SBML file      : %d\n', numel(model.genes));
fprintf('Reactions with a gene rule    : %d\n', sum(~cellfun(@isempty, model.grRules)));

%% ---- 2. Whole-network statistics ------------------------------------
S      = model.S;
[m, n] = size(S);
nGenes = numel(model.genes);
isEx   = (sum(S~=0,1) == 1)';                 % exchange = touches 1 metabolite
isRev  = model.lb < 0;
hasGPR = ~cellfun(@isempty, model.grRules);

fprintf('\n===== RECONSTRUCTION STATISTICS =====\n');
fprintf('Metabolites (m)          : %d\n', m);
fprintf('Reactions (n)            : %d\n', n);
fprintf('Genes                    : %d\n', nGenes);
fprintf('Exchange/boundary rxns   : %d\n', sum(isEx));
fprintf('Internal reactions       : %d\n', sum(~isEx));
fprintf('Reversible reactions     : %d\n', sum(isRev));
fprintf('Irreversible reactions   : %d\n', sum(~isRev));
fprintf('Reactions with GPR rule  : %d\n', sum(hasGPR));
fprintf('Reactions without gene   : %d\n', sum(~hasGPR));
if isfield(model,'comps'), fprintf('Compartments             : %s\n', strjoin(model.comps',', ')); end

%% ---- 3. Stoichiometric matrix S and binary matrix B -----------------
B   = double(S ~= 0);                         % 1 = metabolite takes part, 0 = not
rS  = rank(full(S));
fprintf('\n===== STOICHIOMETRIC MATRIX =====\n');
fprintf('Size of S (m x n)        : %d x %d\n', m, n);
fprintf('Total entries            : %d\n', m*n);
fprintf('Non-zero entries         : %d\n', nnz(S));
fprintf('Zero entries             : %d\n', m*n - nnz(S));
fprintf('Sparsity (%% zeros)       : %.3f %%\n', 100*(1 - nnz(S)/(m*n)));
fprintf('Negative entries (used)  : %d\n', nnz(S<0));
fprintf('Positive entries (made)  : %d\n', nnz(S>0));
fprintf('Avg metabolites/reaction : %.2f\n', nnz(S)/n);
fprintf('Rank of S                : %d\n', rS);
fprintf('Degrees of freedom (n-r) : %d\n', n - rS);
fprintf('Distinct coefficient values in S: %s\n', mat2str(unique(full(S(S~=0)))'));

figure('Color','w','Position',[100 100 900 750]);
spy(S,'k',1);
xlabel(sprintf('Reactions (n = %d)',n)); ylabel(sprintf('Metabolites (m = %d)',m));
title('Stoichiometric matrix S of iCN1361 (C. necator H16)');
set(gca,'FontSize',14);
saveas(gcf, fullfile(outDir,'Fig_S_matrix_spy.png'));

%% ---- 4. Hub metabolites (most connected) ----------------------------
deg = sum(B,2);
[ds, di] = sort(deg,'descend');
fprintf('\n===== TOP 20 HUB METABOLITES =====\n');
for k = 1:20
    fprintf('%-35s %-40s %d rxns\n', model.mets{di(k)}, model.metNames{di(k)}, ds(k));
end

%% ---- 5. Export the full reaction network ----------------------------
formulas = printRxnFormula(model, model.rxns, false);
if isfield(model,'subSystems')
    ss = model.subSystems;
    if iscell(ss) && ~isempty(ss) && iscell(ss{1}), ss = cellfun(@(x) strjoin(x,';'), ss, 'UniformOutput', false); end
else
    ss = repmat({''}, n, 1);
end
T = table(model.rxns, model.rxnNames, formulas, model.lb, model.ub, model.grRules, ss, isEx, ...
    'VariableNames', {'ID','Name','Formula','LB','UB','GPR','Subsystem','IsExchange'});
writetable(T, fullfile(outDir,'all_reactions.csv'));
writetable(table(model.mets, model.metNames, deg, 'VariableNames',{'ID','Name','Degree'}), ...
    fullfile(outDir,'all_metabolites.csv'));
writetable(table(model.genes,'VariableNames',{'Gene'}), fullfile(outDir,'all_genes.csv'));

%% ---- 6. Exchange reactions and default medium -----------------------
fprintf('\n===== EXCHANGE REACTIONS (default bounds) =====\n');
exIdx = find(isEx);
for k = exIdx'
    fprintf('%-40s LB=%8.2f UB=%8.2f  %s\n', model.rxns{k}, model.lb(k), model.ub(k), formulas{k});
end

%% ---- 7. Find the key pathway reactions (for the reaction table) -----
norm = @(c) regexprep(upper(c), '[^A-Z0-9]', '');
rN = norm(model.rxns); nN = norm(model.rxnNames);
keys = { ...
 'CBB: RuBisCO',                 {'RIBULOSEBISPHOSPHATECARBOX','RUBISCO'}
 'CBB: Phosphoribulokinase',     {'PHOSPHORIBULOKINASE'}
 'H2: Hydrogenase',              {'HYDROGENASE','HYDROG'}
 'PHB: beta-ketothiolase phaA',  {'ACETYLCOAACETYLTRANSFER','ACETOACETYLCOATHIOLASE'}
 'PHB: reductase phaB',          {'ACETOACETYLCOAREDUCT'}
 'PHB: synthase phaC',           {'PHB','HYDROXYBUTYRATE','POLYHYDROXY','PHA'}
 'PHB: depolymerase phaZ',       {'DEPOLYMER'}
 'Anaplerotic: PEP carboxylase', {'PEPCARBOX'}
 'Anaplerotic: PEP carboxykinase',{'PEPCARBOXYKIN'}
 'Anaplerotic: Pyruvate carboxyl',{'PYRUVATECARBOXYLASE'}
 'Node: Pyruvate dehydrogenase', {'PYRUVDEH','PYRUVATEDEHYDROGENASE'}
 'TCA: Citrate synthase',        {'CITSYN'}
 'TCA: Aconitase',               {'ACONITATEDEHYDR','ACONITATEHYDR'}
 'TCA: Isocitrate DH',           {'ISOCITDEH'}
 'TCA: 2-oxoglutarate DH',       {'2OXOGLUTARATEDEH','OXOGLUTARATEDEH'}
 'TCA: Succinyl-CoA synthetase', {'SUCCCOASYN'}
 'TCA: Succinate DH',            {'SUCCINATEDEHYDROGENASE','SUCCDEH','SUCCINATEDEH'}
 'TCA: Fumarase',                {'FUMHYDR'}
 'TCA: Malate DH',               {'MALATEDEH'}
 'Glyoxylate: Isocitrate lyase', {'ISOCITCLEAV'}
 'Glyoxylate: Malate synthase',  {'MALSYN'}
 'Succinate transport/exchange', {'SUC','SUCC'}
 'Biomass',                      {'BIOMASS','GROWTH'}
 'Maintenance ATPM',             {'ATPM'}
 'CO2 / H2 / O2 / fructose',     {'CARBONDIOXIDE','CO2','HYDROGENMOLECULE','H2TX','OXYGENMOLECULE','O2','FRUCTOSE','FRU'} };

fid = fopen(fullfile(outDir,'key_pathway_reactions.txt'),'w');
fprintf('\n===== KEY PATHWAY REACTIONS =====\n');
for k = 1:size(keys,1)
    hit = false(n,1);
    for q = keys{k,2}, hit = hit | contains(rN, q{1}) | contains(nN, q{1}); end
    idx = find(hit);
    fprintf(fid, '\n### %s  (%d hits)\n', keys{k,1}, numel(idx));
    fprintf('\n### %s  (%d hits)\n', keys{k,1}, numel(idx));
    for j = idx(1:min(end,15))'
        line = sprintf('%-40s | %s | GPR: %s | %s', model.rxns{j}, model.rxnNames{j}, model.grRules{j}, formulas{j});
        fprintf(fid, '%s\n', line); fprintf('%s\n', line);
    end
end
fclose(fid);

%% ---- 8. phaCAB operon genes and their reactions ---------------------
fprintf('\n===== phaCAB OPERON (H16_A1437 / A1438 / A1439) =====\n');
for g = {'A1437','A1438','A1439'}
    gi = find(contains(model.genes, g{1}));
    for j = gi'
        rx = find(model.rxnGeneMat(:,j));
        fprintf('Gene %s controls: %s\n', model.genes{j}, strjoin(model.rxns(rx)', ', '));
    end
    if isempty(gi), fprintf('Gene %s NOT found in model\n', g{1}); end
end

%% ---- 9. Small 0/1 and S sample around succinate (for Table 4.2) -----
sucMet = find(contains(norm(model.mets),'SUC') & ~contains(norm(model.mets),'SUCCCOA'));
if ~isempty(sucMet)
    rxSuc = find(any(S(sucMet,:)~=0,1));
    metSet = find(any(S(:,rxSuc)~=0,2));
    sub  = full(S(metSet, rxSuc(1:min(end,6))));
    keep = any(sub~=0,2); metSet = metSet(keep); sub = sub(keep,:);
    subT = array2table(sub, 'RowNames', model.mets(metSet), ...
        'VariableNames', matlab.lang.makeValidName(model.rxns(rxSuc(1:min(end,6)))));
    writetable(subT, fullfile(outDir,'S_submatrix_succinate.csv'),'WriteRowNames',true);
    writetable(array2table(double(sub~=0),'RowNames',model.mets(metSet), ...
        'VariableNames', subT.Properties.VariableNames), ...
        fullfile(outDir,'B_binary_submatrix_succinate.csv'),'WriteRowNames',true);
    fprintf('\nSuccinate sub-matrix saved (%d mets x %d rxns)\n', size(sub,1), size(sub,2));
end

%% ---- 10. Quick sanity FBA with default bounds -----------------------
sol = optimizeCbModel(model,'max');
fprintf('\n===== DEFAULT FBA =====\n');
fprintf('Objective reaction : %s\n', strjoin(model.rxns(model.c~=0)', ', '));
fprintf('Status             : %d (1 = optimal)\n', sol.stat);
fprintf('Objective value    : %.6f\n', sol.f);

save(fullfile(outDir,'iCN1361_loaded.mat'),'model');
diary off;
fprintf('\nDONE. Upload all files from folder "%s" to Claude.\n', outDir);
