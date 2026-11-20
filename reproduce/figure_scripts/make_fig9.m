% Fig. 9 rebuilt at f = 10 Hz (all four panels) from the cluster timing files (TM mode = "XY").
% bar order: Assembly | RHS Assembly | Analysis/Factorization | Solve | Get Full Order Solution | Mapping | Communication | Total Runtime
R='/Users/luis/Desktop/Code/Results/100x100x60(1airlayer)/'; OUT='/Users/luis/Desktop/2nd_revision/figures_rescaled/';
% G      : old_code/MT_G_factorized/outputs/timing_nprocs40_ngx100_ngy100_nsy1_nsz1_TE_freq10.txt  (TM block)
% DD     : old_code/DD_ROM_FREQ/outputs/timingsnprocs40_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt
% G-POD  : old_code/MT_G_ROM/outputs/svdtiming_nprocs200_ngx100_ngy100_nsy1_nsz1_TM_freq 10.0000.txt  (unchanged from original panel)
% DD-POD : old_code/DD_ROM_FREQ/outputs/svdtimingsGUESSnprocs200_ngx100_ngy100_nsy1_nsz20_TM_freq 10.0000.txt
P={ 'g',          'fig9a_G',      [4.6305 0 179.4399 1.0735 0 0 0 190.5571],            230, true; ...
    'ddg',        'fig9b_DD',     [0.2941 0.4164 1.8954 19.3156 0.1019 0.2570 0.5435 24.4997], NaN, false; ...
    'g_pod_200',  'fig9c_G-POD',  [],                                                    230, true; ...
    'dd_pod_200', 'fig9d_DD-POD', [3.0429 1.7744 0.0077 0.0120 1.3805 0.1212 0.4098 7.2037],  NaN, false};
for q=1:4
  h=openfig([R P{q,1} '.fig'],'invisible'); h.Position=[1 1 1512 894]; ax=findobj(h,'Type','axes'); b=findobj(ax,'Type','bar');
  old=b.YData; if ~isempty(P{q,3}), b.YData=P{q,3}; end
  v=b.YData; 
  if isnan(P{q,4})
    raw=1.15*max(v); e=10^floor(log10(raw)); m=raw/e; nm=[1 1.2 1.5 2 2.5 3 4 5 6 8 10]; top=nm(find(nm>=m,1))*e; ax.YTickMode='auto';
  else, top=P{q,4}; end
  ax.YLim=[0 top];
  ax.Title.Position(2)=0.99*top;   % title is anchored in data units: keep it on top of the rescaled axes
  t=findobj(ax,'Type','text');
  for x=t'
    s=str2double(strjoin(string(x.String),'')); if isnan(s), continue; end
    xi=round(x.Position(1)); if ~isempty(P{q,3}), x.String=sprintf('%.2f',v(xi)); end
    x.Position(2)=v(xi)+0.012*top; x.VerticalAlignment='bottom';
  end
  fprintf('%-12s old=%s\n             new=%s  ylim=[0 %g]\n',P{q,2},mat2str(round(old,3)),mat2str(round(v,3)),top);
  h.PaperPositionMode='auto'; print(h,'-depsc','-vector','-loose',[OUT P{q,2} '.eps']); savefig(h,[OUT P{q,2} '.fig']); exportgraphics(h,[OUT P{q,2} '.png'],'Resolution',60); close(h);
end
