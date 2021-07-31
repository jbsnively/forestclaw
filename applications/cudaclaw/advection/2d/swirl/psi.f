      double precision function psi(xp,yp)
      implicit none

      double precision xp,yp,pi
      common /compi/ pi

      psi = ((sin(pi*xp))**2 * (sin(pi*yp))**2) / pi

      return
      end
