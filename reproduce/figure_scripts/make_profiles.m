% Regenerates the 8 profile panels of Figs. 12 and 13 with data-fitted y-limits on the phase (right) panels.
D='/Users/luis/Desktop/Code/Results/100x100x60(1airlayer)/'; OUT='/Users/luis/Desktop/2nd_revision/figures_rescaled/';
offx=238667.43349/2000; offy=219608.976480/2000; coord=100;
pre.f001.G='rho_phase_1.00000E-03Hz_100x100x60_40.txt';
pre.f001.GPOD='rxy_y_fase_freq= 0.0010Hz_ngx100_ngy100_ngz60_np2_NS50_tol1.0e-02_relx0.000_penal1.000_maxit1_conv5.0e-04_svd.txt';
pre.f001.DD='rxy_y_fase_freq= 0.10000E-02Hz_ngx100_ngy100_ngz60_np40_no_svd.txt';
pre.f001.DDPOD='GUESSrxy_y_fase_freq= 0.10000E-02Hz_ngx100_ngy100_ngz60_np40_NS50_tol1.0e-02_relx0.050_penal0.050_maxit999_conv5.0e-04_svd.txt';
pre.f10.G='rho_phase_1.00000E+01Hz_100x100_60.txt';
pre.f10.GPOD='rxy_y_fase_freq=10.0000Hz_ngx100_ngy100_ngz60_np2_NS50_tol1.0e-02_relx0.000_penal1.000_maxit1_conv5.0e-04_svd.txt';
pre.f10.DD='rxy_y_fase_freq= 0.10000E+02Hz_ngx100_ngy100_ngz60_np40_no_svd.txt';
pre.f10.DDPOD='GUESSrxy_y_fase_freq= 0.10000E+02Hz_ngx100_ngy100_ngz60_np40_NS50_tol1.0e-02_relx0.050_penal0.050_maxit999_conv5.0e-04_svd.txt';
order={'DDPOD','DD','GPOD','G'};   % order of line objects in the .fig (top of the stack first)
for fq={'f001','f10'}, fq=fq{1};
  for pol=1:2
    for k=1:4
      f=pre.(fq).(order{k}); fn=f;
      if strcmp(order{k},'G'), raw=readmatrix([D fn]); if pol==1, r=raw(:,6);p=raw(:,7); else, r=raw(:,8);p=raw(:,9); end, X=raw;
      else, if pol==2, fn=strrep(f,'rxy','ryx'); end, X=readmatrix([D fn]); r=X(:,3); p=X(:,4); end
      x=flip(X(:,1))+offx; y=X(:,2)+offy; p=atan2d(sind(p),cosd(p));
      [~,i]=min(abs(y-coord)); m=abs(y-y(i))<1e-3; [xs,s]=sort(x(m)); rr=r(m); rr=rr(s); pp=p(m); pp=pp(s);
      C.(fq){pol,k}=struct('x',xs,'rho',rr,'phi',pp);
    end
  end
end
% panels: name, template fig, freq key, pol, quantity
P={'rho_xy','rho_xy','f001',1,'rho'; 'phi_xy','phi_xy','f001',1,'phi'; 'rho_yx','rho_yx','f001',2,'rho'; 'phi_yx','phi_yx','f001',2,'phi'; ...
   'rho_xy_10','rho_xy_10','f10',1,'rho'; 'phi_xy_10','phi_xy_10','f10',1,'phi'; 'rho_yx_10','rho_yx','f10',2,'rho'; 'phi_yx_10','phi_yx','f10',2,'phi'};
fprintf('panel | max|data-fig| per line (0 = fig data identical to files)\n');
for q=1:size(P,1)
  h=openfig([D P{q,2} '.fig'],'invisible'); h.Position=[1 1 1512 894]; ax=findobj(h,'Type','axes');
  L=flipud(findobj(ax,'Type','line')); 
  allv=[]; dif=zeros(1,4);
  for k=1:4
    s=C.(P{q,3}){P{q,4},k}; v=s.(P{q,5}); 
    if numel(L(k).XData)~=numel(s.x), error('size mismatch'); end
    dif(k)=max(abs(double(L(k).YData(:))-v)); L(k).YData=v'; L(k).XData=s.x'; allv=[allv;v]; %#ok
  end
  fprintf('%-10s | %s\n',P{q,1},sprintf('%.3g ',dif));
  if strcmp(P{q,5},'phi')
    lo=min(allv); hi=max(allv); st=30;
    a=[-180 -180 -30 0]; a=a(strcmp(P{q,1},{'phi_xy','phi_xy_10','phi_yx','phi_yx_10'}));   % 90-deg windows (true quadrants where the data allow it)
    b=a+90; assert(lo>=a && hi<=b); ax.YLim=[a b]; ax.YTick=a:st:b; ax.YTickLabelMode='auto';
    fprintf('           phase ylim=[%g %g] step=%g (data %.2f..%.2f)\n',a,b,st,lo,hi);
  end
  h.PaperPositionMode='auto'; print(h,'-depsc','-vector','-loose',[OUT P{q,1} '.eps']); savefig(h,[OUT P{q,1} '.fig']); exportgraphics(h,[OUT P{q,1} '.png'],'Resolution',60); close(h);
end
