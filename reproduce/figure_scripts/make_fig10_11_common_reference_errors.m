function make_common_reference_error_plots()
% Draft error figures with the global solver as the common reference.
% Only archived numerical responses are read. No forward solver is invoked.
folder = fileparts(mfilename('fullpath'));
dataFolder = fullfile(fileparts(folder),'sources','100x100x60(1airlayer)');
modes = [15 25 50 75 100 125 150 175 200];
frequencies = [0.001 10];
tokens = {'0.10000E-02','0.10000E+02'};
globalNames = {'rho_phase_1.00000E-03Hz_100x100x60_40.txt', ...
    'rho_phase_1.00000E+01Hz_100x100_60.txt'};
gpodFrequencyText = {' 0.0010','10.0000'};
suffixes = {'001','10'};
% Match Figure 7's border after scaling each EPS to its manuscript width.
% Reference EPS: 1512 pt wide, 0.5 pt axes, displayed at 8.3 cm.
axesLineWidths = 0.5*(506/1512)*(8.3./[8.3 8.1]);
qoiNames = {'rho_xy','phi_xy','rho_yx','phi_yx'};
yLabels = {'$\mathbf{RelL_{2}(\rho_{xy})\ [\%]}$', ...
    '$\mathbf{RelL_{2}(\phi_{xy})\ [\%]}$', ...
    '$\mathbf{RelL_{2}(\rho_{yx})\ [\%]}$', ...
    '$\mathbf{RelL_{2}(\phi_{yx})\ [\%]}$'};
frequencyTitles = {'$f=0.001\,\mathrm{Hz}$', ...
    '$f=10\,\mathrm{Hz}$'};
errorsG = zeros(9,4,2);
errorsDD = zeros(9,4,2);
resultRows = cell(72,5);
row = 0;

for fi = 1:2
    G = readmatrix(fullfile(dataFolder,globalNames{fi}),'FileType','text');
    assert(isequal(size(G),[10000 11]));
    assert(all(abs(G(:,3)-1/frequencies(fi)) < 1e-8));
    for pi = 1:2
        if pi == 1
            pol = 'xy';
            columns = [6 7];
        else
            pol = 'yx';
            columns = [8 9];
        end
        reference = G(:,columns);
        for ni = 1:numel(modes)
            gpodPattern = sprintf(['r%s_y_fase_freq=%sHz*' ...
                '_np2_NS%d_tol1.0e-02*'],pol,gpodFrequencyText{fi},modes(ni));
            ddpodPattern = sprintf(['GUESSr%s_*%sHz*' ...
                '_np40_NS%d_tol1.0e-02*'],pol,tokens{fi},modes(ni));
            gpod = read_one(dataFolder,gpodPattern);
            ddpod = read_one(dataFolder,ddpodPattern);
            assert(isequal(size(gpod),[10000 4]));
            assert(isequal(size(ddpod),[10000 4]));
            % The saved global coordinates have fewer printed digits.
            assert(max(abs(gpod(:,1:2)-G(:,1:2)),[],'all') < 1e-3);
            assert(max(abs(ddpod(:,1:2)-G(:,1:2)),[],'all') < 1e-3);
            [gr,gp] = errors_against_G(gpod(:,3:4),reference);
            [dr,dp] = errors_against_G(ddpod(:,3:4),reference);
            qi = (pi-1)*2 + (1:2);
            errorsG(ni,qi,fi) = [gr gp];
            errorsDD(ni,qi,fi) = [dr dp];
            for j = 1:2
                row = row+1;
                resultRows(row,:) = {frequencies(fi),modes(ni), ...
                    qoiNames{qi(j)},errorsG(ni,qi(j),fi),errorsDD(ni,qi(j),fi)};
            end
        end
    end

    % Standalone panels preserve the original MATLAB error-plot style.
    for qi = 1:4
        fig = figure('Visible','off','Color','w','Units','pixels', ...
            'Position',[100 100 600 500],'Renderer','painters');
        ax = axes(fig);
        draw_panel(ax,modes,errorsG(:,qi,fi),errorsDD(:,qi,fi),yLabels{qi},axesLineWidths(fi));
        style_common_reference_error_panel(fig,ax,qi,fi);
        base = fullfile(folder,sprintf('%s_error_%s_vs_G',qoiNames{qi},suffixes{fi}));
        export_draft(fig,base,150);
        close(fig);
    end

    % Composite previews are drafts only and do not replace manuscript assets.
    fig = figure('Visible','off','Color','w','Units','pixels', ...
        'Position',[100 100 1240 790],'Renderer','painters');
    layout = tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');
    for qi = 1:4
        ax = nexttile(layout);
        draw_panel(ax,modes,errorsG(:,qi,fi),errorsDD(:,qi,fi),yLabels{qi},axesLineWidths(fi));
        letters = {'a)','b)','c)','d)'};
        text(ax,-0.14,1.005,letters{qi},'Units','normalized', ...
            'Interpreter','none','FontName','Helvetica','FontSize',16, ...
            'VerticalAlignment','middle','Clipping','off');
    end
    sgtitle(layout,frequencyTitles{fi},'Interpreter','latex','FontSize',20);
    base = fullfile(folder,['errors_vs_G_' suffixes{fi}]);
    export_draft(fig,base,130);
    close(fig);
end

assert(row == 72);
results = cell2table(resultRows,'VariableNames', ...
    {'frequency_Hz','nPOD','QoI','GPOD_vs_G_percent','DDPOD_vs_G_percent'});
writetable(results,fullfile(folder,'common_reference_errors_matlab.csv'));
save(fullfile(folder,'common_reference_errors_matlab.mat'), ...
    'modes','frequencies','errorsG','errorsDD','results');
fprintf('Exported eight MATLAB panels and two composite previews.\n');
fprintf('Both methods use G as the reference at all nine mode counts.\n');
end

function data = read_one(folder,pattern)
files = dir(fullfile(folder,pattern));
assert(numel(files) == 1,'Expected exactly one saved response file.');
data = readmatrix(fullfile(files(1).folder,files(1).name),'FileType','text');
assert(all(isfinite(data),'all'));
end

function [rhoError,phiError] = errors_against_G(response,reference)
% Same phase transformation and normalisation as the recovered test_xy.m.
phi = mod(response(:,2)+90,360)-180;
phiG = mod(reference(:,2)+90,360)-180;
rhoError = 100*norm(response(:,1)-reference(:,1),2)/norm(reference(:,1),2);
phiError = 100*norm(phi-phiG,2)/(norm(mod(abs(phiG),90),2)+eps);
end

function draw_panel(ax,modes,gpod,ddpod,yLabel,axesLineWidth)
hold(ax,'on');
hG = plot(ax,modes,gpod,'-ro','LineWidth',2, ...
    'MarkerFaceColor','r','DisplayName','G-POD');
hDD = plot(ax,modes,ddpod,'-bs','LineWidth',2, ...
    'MarkerFaceColor','b','DisplayName','DD-POD');
set(ax,'TickLabelInterpreter','latex','FontSize',12,'FontWeight','bold', ...
    'Box','on','LineWidth',axesLineWidth);
% A small extension prevents clipping the newly computed 15.10% point.
xlim(ax,[5 210]);
ylim(ax,[0 16]);
xticks(ax,modes);
yticks(ax,[0 5 10 15]);
ylabel(ax,yLabel,'Interpreter','latex','FontSize',18,'FontWeight','bold');
xlabel(ax,'$\mathbf{\#\ Modes}$ ($n_{\mathrm{POD}}$)', ...
    'Interpreter','latex','FontSize',16,'FontWeight','bold');
grid(ax,'on');
legend(ax,[hG hDD],'Interpreter','latex','Location','northeast', ...
    'FontSize',14,'FontWeight','normal');
hold(ax,'off');
end

function export_draft(fig,base,resolution)
drawnow;
savefig(fig,[base '.fig']);
set(fig,'InvertHardcopy','off');
print(fig,[base '.eps'],'-depsc','-painters','-loose');
print(fig,[base '.png'],'-dpng',['-r' num2str(resolution)]);
end
