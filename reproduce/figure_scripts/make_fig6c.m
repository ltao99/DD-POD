% Fig. 6c (G-POD, 80x80x72, 0.001 Hz) rebuilt from the simulation log
%   old_code/MT_G_ROM/outputs/svdtiming_nprocs2_ngx80_ngy80_nsy1_nsz1_TM_freq  0.0010.txt
% bar order: Assembly | RHS Assembly | Analysis/Factorization | Solve | Get Full Order Solution | Mapping | Communication | Total
R='/Users/luis/Desktop/Code/Results/80x80x72 2/'; OUT='/Users/luis/Desktop/2nd_revision/figures_rescaled/';
v=[5.9669 0 0.0066 0.0001 0 0 0 7.3810];     % same bars as the original panel are left at zero (G-POD has no iterations)
h=openfig([R 'g_pod.fig'],'invisible'); h.Position=[1 1 1512 894]; ax=findobj(h,'Type','axes'); b=findobj(ax,'Type','bar');
fprintf('old bars=%s\n',mat2str(round(b.YData,4))); b.YData=v; top=ax.YLim(2);
for x=findobj(ax,'Type','text')'
  s=str2double(strjoin(string(x.String),'')); if isnan(s), continue; end
  xi=round(x.Position(1)); x.String=sprintf('%.2f',v(xi)); x.Position(2)=v(xi)+0.012*top; x.VerticalAlignment='bottom';
end
fprintf('new bars=%s  ylim=[0 %g] (left panel: axis unchanged)\n',mat2str(v),top);
h.PaperPositionMode='auto'; print(h,'-depsc','-vector','-loose',[OUT 'fig6c_G-POD.eps']); savefig(h,[OUT 'fig6c_G-POD.fig']); exportgraphics(h,[OUT 'fig6c_G-POD.png'],'Resolution',60); close(h);
