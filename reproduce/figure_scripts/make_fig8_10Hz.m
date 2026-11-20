% Fig. 8b (f = 10 Hz): three plotted values replaced by the logged values, legend label corrected.
%   G-POD total,   n_POD = 175: 14.3  -> 16.1070  (MT_G_ROM/outputs/svdtiming_nprocs175_..._TM_freq 10.0000.txt)
%   G-POD assembly,n_POD = 175: 12.29 -> 14.5306  (same file)
%   DD-POD assembly, n_POD = 200: 4.04 -> 3.7380  (DD_ROM_FREQ/outputs/svdtimingsnprocs200_..._TE_freq 10.0000.txt)
R='/Users/luis/Desktop/Code/Results/100x100x60(1airlayer)/'; OUT='/Users/luis/Desktop/2nd_revision/figures_rescaled/';
h=openfig([R 'computational_performance_10.fig'],'invisible'); h.Position=[1 1 1512 894]; ax=findobj(h,'Type','axes'); L=findobj(ax,'Type','line');
fix={'G-POD',175,16.1070; 'G - Assemble',175,14.5306; 'DD-POD - Assemble',200,3.7380};
for k=1:size(fix,1)
  for l=L'
    if strcmp(l.DisplayName,fix{k,1})
      i=find(l.XData==fix{k,2}); fprintf('%-18s n=%d: %.4f -> %.4f\n',fix{k,1},fix{k,2},l.YData(i),fix{k,3}); l.YData(i)=fix{k,3};
    end
  end
end
% the 'Iterative Scheme' arrow at n_POD = 200 ends just above the DD-POD assembly marker: move its lower end with the marker
dy = ax.Position(4)*(4.04-3.7380)/ax.YLim(2);
for a=findall(h,'Type','doubleendarrowshape')', if abs(a.X(1)-0.8671)<0.002, y=a.Y; y(2)=y(2)-dy; a.Y=y; fprintf('arrow lower end moved by %.4f (normalised units)\n',dy); end, end
for l=L', if strcmp(l.DisplayName,'G - Assemble'), l.DisplayName='G-POD - Assemble'; end, end
fprintf('YLim=%s\n',mat2str(ax.YLim));
h.PaperPositionMode='auto'; print(h,'-depsc','-vector','-loose',[OUT 'fig8b_10Hz.eps']); savefig(h,[OUT 'fig8b_10Hz.fig']); exportgraphics(h,[OUT 'fig8b_10Hz.png'],'Resolution',60); close(h);
