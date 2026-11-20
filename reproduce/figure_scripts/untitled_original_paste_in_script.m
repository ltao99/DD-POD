% ==== Paste your full text block here ====
txt = [
"Frequency (Hz):              0.0010"
"Processes:                   40"
"Mode:                    TE"
"Mesh size (ngx, ngy, ngz):    100    100     60"
"Subdomains (nsy, nsz):          1     20"
"Iterations:                    91"
"Read Model (ONCE) (s):            0.1698"
"Nodes Time (ONCE) (s):            0.0008"
"Basis Time (ONCE) (s):         2778.7386"
"Source Time (ONCE) (s):           0.0192"
"Assembly Time (ONCE) (s)          4.04"
"Compute Coef Time (ONCE) (s):     0.0019"
"Analysis and Factor (s):          0.0111"
"Solve Time (s):                   0.0118"
"Communication Time (s):           0.4968"
"RHS Time (s):                     1.5007"
"x Projection Time (s):            1.078"
"Mapping Time (s):                 0.1214"
"Calc Error Time (s):              0.0315"
"Gather Time (s):                  0.0025"
"Output Time (s):                  0.0097"
"Overhead Time (s):                0.0000"
"Total Time (s):                   4.4583"
"Frequency (Hz):              0.0010"
"Processes:                   18"
"Mode:                    TM"
"Mesh size (ngx, ngy, ngz):    100    100    60"
"Subdomains (nsy, nsz):          1     1"
"Iterations:                    92"
"Read Model (ONCE) (s):            0.0000"
"Nodes Time (ONCE) (s):            0.0000"
"Basis Time (ONCE) (s):            0.0000"
"Source Time (ONCE) (s):           0.0000"
"Assembly Time (ONCE) (s)          0.0000"
"Compute Coef Time (ONCE) (s):     0.0000"
"Analysis and Factor (s):          0.0000"
"Solve Time (s):                   0.0000"
"Communication Time (s):           0.0000"
"RHS Time (s):                     0.0000"
"x Projection Time (s):            0.0000"
"Mapping Time (s):                 0.0000"
"Calc Error Time (s):              0.0000"
"Gather Time (s):                  0.0000"
"Output Time (s):                  0.0000"
"Overhead Time (s):                0.0000"
"Total Time (s):                  7.20"
];

% ---------- Normalize input ----------
if ischar(txt), txt = string(strsplit(txt, newline)); end
if isstring(txt) && isscalar(txt), txt = splitlines(txt); end
txt = txt(~(strlength(strtrim(txt))==0));
txt = cellstr(txt);

starts   = @(s,pfx) strncmpi(strtrim(s), pfx, length(pfx));
lastnum  = @(s) str2double(get_last_tokens(s,1));

% ---------- Steps to include ----------
keepKeys = { ...
    'Assembly Time (ONCE)', ...
    'RHS Time', ...
    'Analysis and Factor', ...
    'Solve Time', ...
    'x Projection Time', ...
    'Mapping Time', ...
    'Communication Time', ...
    'Total Time' ...
};

prettyNames = { ...
    'Assembly', ...
    'RHS Assembly', ...
    'Analysis/Factorization', ...
    'Solve', ...
    'Get Full Order Solution', ...
    'Mapping', ...
    'Communication', ...
    'Total Runtime' ...
};

% ---------- Split by Mode ----------
isMode = cellfun(@(s) starts(s,'Mode:'), txt);
modeIdx = find(isMode);
blockStarts = modeIdx(:);
blockEnds   = [modeIdx(2:end)-1; numel(txt)];

modeMaps = struct();
itersTE = NaN; itersTM = NaN;
lastTotal = NaN;

for b = 1:numel(blockStarts)
    lines = txt(blockStarts(b):blockEnds(b));
    parts = strsplit(strtrim(lines{1}));
    modeName = upper(string(parts{end}));

    % Iterations
    itLine = lines(cellfun(@(s) starts(s,'Iterations:'), lines));
    if ~isempty(itLine)
        its = lastnum(itLine{1});
        if modeName=="TE", itersTE = its; else, itersTM = its; end
    end

    % Step values (skip Overhead Time)
    M = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(lines)
        L = strtrim(lines{i});
        if contains(L,'Overhead Time'), continue; end
        if contains(L,'(s)') || starts(L,'Total Time')
            val = lastnum(L);
            key = L;
            k = strfind(key,'(s)');
            if ~isempty(k), key = key(1:k(1)-1); end
            key = strtrim(key);
            if endsWith(key,':'), key = key(1:end-1); end
            if contains(L,'(ONCE)') && ~contains(key,'(ONCE)')
                key = [key ' (ONCE)'];
            end
            if any(strcmp(key, keepKeys)), M(key) = val; end
            if starts(L,'Total Time'), lastTotal = val; end
        end
    end
    modeMaps.(char(modeName)) = M;
end

% ---------- Build values (1 bar per step) ----------
vals = zeros(numel(keepKeys),1);
for si = 1:numel(keepKeys)
    key = keepKeys{si};
    if strcmp(key,'Total Time')
        vals(si) = lastTotal;                       % only the last total (TM)
    elseif contains(key,'(ONCE)')
        if isfield(modeMaps,'TE') && isKey(modeMaps.TE, key)
            vals(si) = modeMaps.TE(key);
        elseif isfield(modeMaps,'TM') && isKey(modeMaps.TM, key)
            vals(si) = modeMaps.TM(key);
        else
            vals(si) = 0;
        end
    else
        v = 0;
        if isfield(modeMaps,'TE') && isKey(modeMaps.TE, key), v = v + modeMaps.TE(key); end
        if isfield(modeMaps,'TM') && isKey(modeMaps.TM, key), v = v + modeMaps.TM(key); end
        vals(si) = v;
    end
end

% ---------- Plot (identical canvas and margins, bigger fonts) ----------
FIG_W = 1000; FIG_H = 500; BAR_W = 0.70;

% Font sizes
TITLE_F = 40; AXIS_LABEL_F = 26; TICK_F = 22; VALUE_F = 20; LEGEND_F = 22;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure('Color','w','Units','pixels','Position',[100 100 FIG_W FIG_H]);

b = bar(vals, 'EdgeColor','none', 'BarWidth', BAR_W); hold on
b.FaceColor = [0.28 0.49 0.82];

% Grey for Total runtime
c = repmat(b.FaceColor, numel(vals), 1);
totIdx = find(strcmp(keepKeys,'Total Time'));
if ~isempty(totIdx), c(totIdx,:) = [0.35 0.35 0.35]; end
b.CData = c;

xticks(1:numel(prettyNames));
xticklabels(prettyNames);
xtickangle(25);

ylabel('\textbf{Time} (s)','FontSize',AXIS_LABEL_F);
%h = xlabel('\textbf{Solver step}','FontSize',AXIS_LABEL_F);

title(sprintf('\\textbf{DDG--POD($\\epsilon = 10^{-2}$) -- Iterations: TE = %g, TM = %g}', ...
               itersTE, itersTM), ...
      'FontSize', TITLE_F, 'Interpreter', 'latex');

set(gca,'FontSize',TICK_F);
ylim([0 230.1]);   % per your request
grid on; box on;

%pos = get(h, 'Position');   % current [x y z] position
%pos(2) = pos(2) - 5;      % move it upward (increase y-value)
%set(h, 'Position', pos);    % update position

% --- FORCE IDENTICAL AXES WIDTH/MARGINS ---
ax = gca;
ax.Units = 'normalized';
ax.ActivePositionProperty = 'Position';
ax.Position   = [0.10 0.20 0.80 0.70];   % [left bottom width height]
ax.LooseInset = [0 0 0 0];

% dashed separator before "Total Runtime"
if ~isempty(totIdx)
    xline(totIdx-0.5,'k:','LineWidth',1);
end

% Enlarge legend if one exists
lgd = findobj(gcf,'Type','Legend');
if ~isempty(lgd), set(lgd,'FontSize',LEGEND_F); end

% ---- Add values on top of bars (2 decimals) ----
for i = 1:numel(vals)
    text(i, vals(i) + 0.02, sprintf('%.2f', vals(i)), ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontSize',VALUE_F, 'Interpreter','latex');
end

% ---------- Helper ----------
function toks = get_last_tokens(s, k)
    parts = strsplit(strtrim(s));
    toks = string(parts(end-k+1:end));
end





%% ----------- CODE 2 (pad with fake bars to match 9 categories) -----------
% ===== Paste your full text block here =====
txt = [
"Frequency (Hz):                          0.0010"
"Processes:                               18"
"Mode:                                TE"
"Mesh size (ngx, ngy, ngz):               100     100     60"
"Assembly Time (s):                       16.16"
"Analysis and Factorization Time (s):      0.15"
"Solve Time (s):                          0.01"
"Output Time (s):                         0.2254"
"Total Time (s):                          17.78"
"Frequency (Hz):                          0.0010"
"Processes:                               18"
"Mode:                                TM"
"Mesh size (ngx, ngy, ngz):               100     100     60"
"Assembly Time (s):                       0.0000"
"Analysis and Factorization Time (s):     0.0000"
"Solve Time (s):                          0.0000"
"Output Time (s):                         0.0000"
"Total Time (s):                          17.78"
];


% ---------- Normalize input ----------
if ischar(txt), txt = string(strsplit(txt, newline)); end
if isstring(txt) && isscalar(txt), txt = splitlines(txt); end
txt = string(txt);
txt = txt(strlength(strtrim(txt))>0);
txt = cellstr(txt);

starts = @(s,pfx) strncmpi(strtrim(s), pfx, length(pfx));
lastnum_inline = @(s) str2double(regexp(strtrim(s), '\S+$', 'match', 'once'));

% ---------- Keys we actually read from the log ----------
keySpecs = { ...
    'Assembly Time',                    'Assembly'; ...
    'Analysis and Factorization Time',  'Analysis/Factorization'; ...
    'Solve Time',                       'Solve'; ...
    'Total Time',                       'Total Runtime' ...
};

% ---------- Parse TE/TM blocks ----------
isMode = cellfun(@(s) starts(s,'Mode:'), txt);
modeIdx = find(isMode);
blockStarts = modeIdx(:);
blockEnds   = [modeIdx(2:end)-1; numel(txt)];

modeMaps = struct(); lastTotal = NaN;
for b = 1:numel(blockStarts)
    lines = txt(blockStarts(b):blockEnds(b));
    mdparts  = strsplit(strtrim(lines{1}));
    modeName = upper(string(mdparts{end}));
    M = containers.Map('KeyType','char','ValueType','double');
    for i = 1:size(keySpecs,1)
        patt = keySpecs{i,1};
        hit = cellfun(@(s) starts(s,[patt ' (s):']), lines);
        if any(hit)
            v = lastnum_inline(lines{find(hit,1,'last')});
            M(patt) = v;
            if strcmp(patt,'Total Time'), lastTotal = v; end
        end
    end
    modeMaps.(char(modeName)) = M;
end

% ---------- Sum TE+TM for real metrics ----------
realVals = containers.Map('KeyType','char','ValueType','double');
realVals('Assembly')               = get_from_map(modeMaps,'TE','Assembly Time') ...
                                   + get_from_map(modeMaps,'TM','Assembly Time');
realVals('Analysis/Factorization') = get_from_map(modeMaps,'TE','Analysis and Factorization Time') ...
                                   + get_from_map(modeMaps,'TM','Analysis and Factorization Time');
realVals('Solve')                  = get_from_map(modeMaps,'TE','Solve Time') ...
                                   + get_from_map(modeMaps,'TM','Solve Time');
realVals('Total Runtime')          = lastTotal;

% ---------- Target order (match CODE 1: 9 categories) ----------
targetLabels = { ...
    'Assembly', ...
    'RHS Assembly', ...
    'Analysis/Factorization', ...
    'Solve', ...
    'Get Full Order Solution', ...
    'Mapping', ...
    'Communication', ...
    'Total Runtime' ...
};

% ---------- Build padded values (fake bars = 0; blank tick labels) ----------
vals = zeros(numel(targetLabels),1);
dispLabels = targetLabels;
for i = 1:numel(targetLabels)
    lbl = targetLabels{i};
    if isKey(realVals, lbl)
        vals(i) = realVals(lbl);
    else
        vals(i) = 0;        % fake bar
        dispLabels{i} = ''; % hide label of fake category
    end
end

% ---------- Plot (identical canvas and margins, bigger fonts) ----------
FIG_W = 1000; FIG_H = 500; BAR_W = 0.70;

% Font sizes
TITLE_F = 40; AXIS_LABEL_F = 26; TICK_F = 22; VALUE_F = 20; LEGEND_F = 22;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure('Color','w','Units','pixels','Position',[100 100 FIG_W FIG_H]);

b = bar(vals,'EdgeColor','none','BarWidth',BAR_W); hold on
b.FaceColor = [0.28 0.49 0.82];

% Color the Total Runtime bar gray (index 9)
totIdx = 9;
C = repmat(b.FaceColor, numel(vals), 1);
C(totIdx,:) = [0.35 0.35 0.35];
b.CData = C;

% dashed separator before "Total Runtime"
xline(totIdx-1.5,'k:','LineWidth',1);

xticks(1:numel(dispLabels));
xticklabels(dispLabels);
xtickangle(25);

ylabel('\textbf{Time} (s)','FontSize',AXIS_LABEL_F);
%h = xlabel('\textbf{Solver step}','FontSize',AXIS_LABEL_F);
title('\textbf{G--POD}','FontSize',TITLE_F,'Interpreter','latex');

%pos = get(h, 'Position');   % current [x y z] position
%pos(2) = pos(2) - 0.05 ;      % move it upward (increase y-value)
%set(h, 'Position', pos);    % update position

set(gca,'FontSize',TICK_F);
ylim([0 230]);   % per your request
grid on; box on;

% --- FORCE IDENTICAL AXES WIDTH/MARGINS ---
ax = gca;
ax.Units = 'normalized';
ax.ActivePositionProperty = 'Position';
ax.Position   = [0.10 0.20 0.80 0.70];   % same as Code 1
ax.LooseInset = [0 0 0 0];

% Enlarge legend if one exists
lgd = findobj(gcf,'Type','Legend');
if ~isempty(lgd), set(lgd,'FontSize',LEGEND_F); end

% Value labels only for non-fake bars
for i = 1:numel(vals)
    if vals(i) > 0
        text(i, vals(i)+0.02, sprintf('%.2f', vals(i)), ...
            'HorizontalAlignment','center','VerticalAlignment','bottom', ...
            'FontSize',VALUE_F,'Interpreter','latex');
    end
end

% ---------- Local helper (keep helpers at end of script) ----------
function v = get_from_map(m, fld, key)
    v = 0;
    if isfield(m, fld)
        mp = m.(fld);
        if isKey(mp, key)
            v = mp(key);
        end
    end
end


