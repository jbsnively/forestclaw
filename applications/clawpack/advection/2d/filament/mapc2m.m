function [xp,yp,zp] = mapc2m(xc,yc)

parms = read_vars();
if parms.example == 1
    map = 'cart';   % brick
elseif parms.example == 2
    map = 'fivepatch';
elseif parms.example == 3
    map = 'bilinear';
end


% This domain should be in [0,1],[0,1]

alpha = parms.alpha;

shift = [1,1,0];

switch map
    case 'cart'
        % (xc,yc) in [0,1]x[0,1]
        s = 0;
        [xc1,yc1,~] = mapc2m_brick(xc,yc,s);
        [xp,yp,~] = mapc2m_cart(xc1,yc1);
                
        xp = xp + shift(1);
        yp = yp + shift(2);
        
    case 'fivepatch'
        [xp,yp,~] = mapc2m_fivepatch(xc,yc,alpha);

        xp = xp + shift(1);
        yp = yp + shift(2);

    case 'bilinear'
        center = parms.center;
        [xp,yp,~] = mapc2m_bilinear(xc,yc,center);

        xp = xp + shift(1);
        yp = yp + shift(2);

end
zp = 0*xp;




end