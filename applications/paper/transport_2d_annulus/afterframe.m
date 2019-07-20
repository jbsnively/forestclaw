setviews;

% Geometry
[example,A,rinit,pstart, beta,theta,vcart, ...
        freq,cart_speed] = read_vars();
;
tfinal = 0.25;
vcart = cart_speed;
dlen = vcart*tfinal;

pend = pstart + [dlen, 0];

fprintf('%10s %24.16e\n','qmin',qmin);
fprintf('%10s %24.16e\n','qmax',qmax);


s = 1e-2;    
alim = [-1,1];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)


caxis([-1,1]);

set(gca,'fontsize',16);

axis([-0.707106781186547   0.707106781186548   0.282842712474619,1]);

%axis([-0.211745637207774, 0.211745637207774, 0.5, 0.85])


% showgridlines;

plot_path = true;
if (plot_path)    
    hold on;
    
    N = 128;
    th0 = linspace(0,2*pi,N+1);
    x0 = rinit*cos(th) + pstart(1);
    y0 = rinit*sin(th) + pstart(2);
        
    if (t > 0)
        Y0 = [x0(:); y0(:)];    
        [tout,yout] = ode45(@vel_ellipse,[0,t],Y0);
        xth = yout(end,1:N+1)';
        yth = yout(end,N+2:end)';
    else
        xth = x0;
        yth = y0;
    end
        
    plot(xth,yth,'k','linewidth',2);
    hold off
    
    if (example == 4)
        hold on;        
        t1 = theta(1);
        t2 = theta(2);
        t0 = atan2(y0,x0);
        m = t0 < 0;
        t0(m) = t0(m) + 2*pi;
        t0 = t0/(2*pi);
        xm = (t0-t1)/(t2-t1);
        rm = sqrt(x0.^2 + y0.^2);
        ym = (rm-beta)/(1-beta);
        % hq = ellipse_vel(4,1.5,xm,ym);
        % set(hq,'linewidth',2,'color','k');
        hold off
    end
end


%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
