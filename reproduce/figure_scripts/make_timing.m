% Timing-bar figures (Figs. 6, 7, 9): right-hand panels get their own y-axis fitted to the data.
R='/Users/luis/Desktop/Code/Results/'; OUT='/Users/luis/Desktop/2nd_revision/figures_rescaled/';
P={'80x80x72 2/ddg_005','fig6b_DD',1; '80x80x72 2/ddg_pod_005','fig6d_DD-POD',1; '160x160x144 2/ddg_pod_005','fig7b_DD-POD',1; ...
   '100x100x60(1airlayer)/ddg','fig9b_DD',1; '100x100x60(1airlayer)/dd_pod_200','fig9d_DD-POD',1};
for q=1:size(P,1)
  h=openfig([R P{q,1} '.fig'],'invisible'); h.Position=[1 1 1512 894]; ax=findobj(h,'Type','axes'); b=findobj(ax,'Type','bar');
  mx=max(b.YData); raw=1.15*mx; e=10^floor(log10(raw)); m=raw/e; nm=[1 1.2 1.5 2 2.5 3 4 5 6 8 10]; top=nm(find(nm>=m,1))*e;
  % keep the original ylim for reference
  old=ax.YLim; ax.YLim=[0 top]; ax.YTickMode='auto';
  ax.Title.Position(2)=0.99*top;   % title is anchored in data units: keep it on top of the rescaled axes
  % value labels: put each just above its bar (original nudges were tuned for the shared axis)
  t=findobj(ax,'Type','text'); 
  for x=t'
    s=str2double(strjoin(string(x.String),'')); if isnan(s), continue; end
    xi=round(x.Position(1)); x.Position(2)=b.YData(xi)+0.012*top; x.VerticalAlignment='bottom';
  end
  fprintf('%-14s old YLim=[%g %g] -> [0 %g]  (max bar %.3g)  ticks=%s\n',P{q,2},old(1),old(2),top,mx,mat2str(ax.YTick));
  h.PaperPositionMode='auto'; print(h,'-depsc','-vector','-loose',[OUT P{q,2} '.eps']); savefig(h,[OUT P{q,2} '.fig']); exportgraphics(h,[OUT P{q,2} '.png'],'Resolution',60); close(h);
end
