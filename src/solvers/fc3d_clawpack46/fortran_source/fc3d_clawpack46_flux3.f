c
c
c     ==================================================================
      subroutine flux3(ixyz,maxm,meqn,maux,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 wave,s,amdq,apdq,cqxx,
     &                 bmamdq,bmapdq,bpamdq,bpapdq,
     &                 cmamdq,cmapdq,cpamdq,cpapdq,
     &                 cmamdq2,cmapdq2,cpamdq2,cpapdq2,
     &                 bmcqxxm,bmcqxxp,bpcqxxm,bpcqxxp,
     &                 cmcqxxm,cmcqxxp,cpcqxxm,cpcqxxp,
     &                 bmcmamdq,bmcmapdq,bpcmamdq,bpcmapdq,
     &                 bmcpamdq,bmcpapdq,bpcpamdq,bpcpapdq,
     &                 rpn3,rpt3, rptt3,mwaves,mcapa,method,mthlim,
     &                 use_fwaves)
c     ==================================================================
c
c     # clawpack routine ...  modified for AMRCLAW
c
c     # Compute the modification to fluxes f, g and h that are generated by
c     # all interfaces along a 1D slice of the 3D grid.
c     #    ixyz = 1  if it is a slice in x
c     #           2  if it is a slice in y
c     #           3  if it is a slice in z
c     # This value is passed into the Riemann solvers. The flux modifications
c     # go into the arrays fadd, gadd and hadd.  The notation is written
c     # assuming we are solving along a 1D slice in the x-direction.
c
c     # fadd(i,.) modifies F to the left of cell i
c     # gadd(i,.,1,slice) modifies G below cell i (in the z-direction)
c     # gadd(i,.,2,slice) modifies G above cell i
c     #                   The G flux in the surrounding slices may
c     #                   also be updated.
c     #                   slice  =  -1     The slice below in y-direction
c     #                   slice  =   0     The slice used in the 2D method
c     #                   slice  =   1     The slice above in y-direction
c     # hadd(i,.,1,slice) modifies H below cell i (in the y-direction)
c     # hadd(i,.,2,slice) modifies H above cell i
c     #                   The H flux in the surrounding slices may
c     #                   also be updated.
c     #                   slice  =  -1     The slice below in z-direction
c     #                   slice  =   0     The slice used in the 2D method
c     #                   slice  =   1     The slice above in z-direction
c     #
c     # The method used is specified by method(2) and method(3):
c
c        method(2) = 1 No correction waves
c                  = 2 if second order correction terms are to be added, with
c                      a flux limiter as specified by mthlim.  No transverse
c                      propagation of these waves.
c
c         method(3) specify how the transverse wave propagation
c         of the increment wave and the correction wave are performed.
c         Note that method(3) is given by a two digit number, in
c         contrast to what is the case for claw2. It is convenient
c         to define the scheme using the pair (method(2),method(3)).
c
c         method(3) <  0 Gives dimensional splitting using Godunov
c                        splitting, i.e. formally first order
c                        accurate.
c                      0 Gives the Donor cell method. No transverse
c                        propagation of neither the increment wave
c                        nor the correction wave.
c                   = 10 Transverse propagation of the increment wave
c                        as in 2D. Note that method (2,10) is
c                        unconditionally unstable.
c                   = 11 Corner transport upwind of the increment
c                        wave. Note that method (2,11) also is
c                        unconditionally unstable.
c                   = 20 Both the increment wave and the correction
c                        wave propagate as in the 2D case. Only to
c                        be used with method(2) = 2.
c                   = 21 Corner transport upwind of the increment wave,
c                        and the correction wave propagates as in 2D.
c                        Only to be used with method(2) = 2.
c                   = 22 3D propagation of both the increment wave and
c                        the correction wave. Only to be used with
c                        method(2) = 2.
c
c         Recommended settings:   First order schemes:
c                                       (1,10) Stable for CFL < 1/2
c                                       (1,11) Stable for CFL < 1
c                                 Second order schemes:
c                                        (2,20) Stable for CFL < 1/2
c                                        (2,22) Stable for CFL < 1
c
c         WARNING! The schemes (2,10), (2,11) are unconditionally
c                  unstable.
c
c                       ----------------------------------
c
c     Note that if method(6)=1 then the capa array comes into the second
c     order correction terms, and is already included in dtdx1d:
c     If ixyz = 1 then
c        dtdx1d(i) = dt/dx                      if method(6) = 0
c                  = dt/(dx*capa(i,jcom,kcom))  if method(6) = 1
c     If ixyz = 2 then
c        dtdx1d(j) = dt/dy                      if method(6) = 0
c                  = dt/(dy*capa(icom,j,kcom))  if method(6) = 1
c     If ixyz = 3 then
c        dtdx1d(k) = dt/dz                      if method(6) = 0
c                  = dt/(dz*capa(icom,jcom,k))  if method(6) = 1
c
c     Notation:
c        The jump in q (q1d(i,:)-q1d(i-1,:))  is split by rpn3 into
c            amdq =  the left-going flux difference  A^- Delta q
c            apdq = the right-going flux difference  A^+ Delta q
c        Each of these is split by rpt3 into
c            bmasdq = the down-going transverse flux difference B^- A^* Delta q
c            bpasdq =   the up-going transverse flux difference B^+ A^* Delta q
c        where A^* represents either A^- or A^+.
c
c        Finally, each bsasdq is split by rptt3 into :
c            cmbsasdq = C^- B^* A^* Dq
c            cpbsasdq = C^+ B^* A^* Dq
c
c
c      use amr_module
      implicit none

      external rpn3,rpt3, rptt3


      integer ixyz, maxm, meqn, maux, mbc, mx, mwaves, mcapa
      integer mthlim(mwaves), method(7)
      double precision cfl1d, dtdy, dtdz

      double precision     q1d(1-mbc:maxm+mbc, meqn)
      double precision    amdq(1-mbc:maxm+mbc, meqn)
      double precision    apdq(1-mbc:maxm+mbc, meqn)
      double precision  bmamdq(1-mbc:maxm+mbc, meqn)
      double precision  bmapdq(1-mbc:maxm+mbc, meqn)
      double precision  bpamdq(1-mbc:maxm+mbc, meqn)
      double precision  bpapdq(1-mbc:maxm+mbc, meqn)
      double precision    cqxx(1-mbc:maxm+mbc, meqn)
      double precision   faddm(1-mbc:maxm+mbc, meqn)
      double precision   faddp(1-mbc:maxm+mbc, meqn)
      double precision    gadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
      double precision    hadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
c
      double precision  cmamdq(1-mbc:maxm+mbc, meqn)
      double precision  cmapdq(1-mbc:maxm+mbc, meqn)
      double precision  cpamdq(1-mbc:maxm+mbc, meqn)
      double precision  cpapdq(1-mbc:maxm+mbc, meqn)
c
      double precision  cmamdq2(1-mbc:maxm+mbc, meqn)
      double precision  cmapdq2(1-mbc:maxm+mbc, meqn)
      double precision  cpamdq2(1-mbc:maxm+mbc, meqn)
      double precision  cpapdq2(1-mbc:maxm+mbc, meqn)
c
      double precision  bmcqxxm(1-mbc:maxm+mbc, meqn)
      double precision  bpcqxxm(1-mbc:maxm+mbc, meqn)
      double precision  bmcqxxp(1-mbc:maxm+mbc, meqn)
      double precision  bpcqxxp(1-mbc:maxm+mbc, meqn)
c
      double precision  cmcqxxm(1-mbc:maxm+mbc, meqn)
      double precision  cpcqxxm(1-mbc:maxm+mbc, meqn)
      double precision  cmcqxxp(1-mbc:maxm+mbc, meqn)
      double precision  cpcqxxp(1-mbc:maxm+mbc, meqn)

      double precision  bpcmamdq(1-mbc:maxm+mbc, meqn)
      double precision  bpcmapdq(1-mbc:maxm+mbc, meqn)
      double precision  bpcpamdq(1-mbc:maxm+mbc, meqn)
      double precision  bpcpapdq(1-mbc:maxm+mbc, meqn)
      double precision  bmcmamdq(1-mbc:maxm+mbc, meqn)
      double precision  bmcmapdq(1-mbc:maxm+mbc, meqn)
      double precision  bmcpamdq(1-mbc:maxm+mbc, meqn)
      double precision  bmcpapdq(1-mbc:maxm+mbc, meqn)
c
      double precision dtdx1d(1-mbc:maxm+mbc)
      double precision aux1(1-mbc:maxm+mbc, maux, 3)
      double precision aux2(1-mbc:maxm+mbc, maux, 3)
      double precision aux3(1-mbc:maxm+mbc, maux, 3)
c
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision  wave(1-mbc:maxm+mbc, meqn, mwaves)
c
      logical limit

      double precision dtcom, dxcom, dycom, dzcom, tcom
      integer icom,jcom,kcom
      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

      integer mw, m, i,m3, m4
      double precision dtdxave
      integer use_fwaves

      limit = .false.
      do mw = 1,mwaves
         if (mthlim(mw) .gt. 0) then 
            limit = .true.
         endif
      end do
c     
c     # initialize flux increments:
c     -----------------------------
c     
      do i = 1-mbc, mx+mbc
         do m = 1,meqn
            faddm(i,m) = 0.d0
            faddp(i,m) = 0.d0
            gadd(i,m,1,-1) = 0.d0
            gadd(i,m,1, 0) = 0.d0
            gadd(i,m,1, 1) = 0.d0
            hadd(i,m,1,-1) = 0.d0
            hadd(i,m,1, 0) = 0.d0
            hadd(i,m,1, 1) = 0.d0
            gadd(i,m,2,-1) = 0.d0
            gadd(i,m,2, 0) = 0.d0
            gadd(i,m,2, 1) = 0.d0
            hadd(i,m,2,-1) = 0.d0
            hadd(i,m,2, 0) = 0.d0
            hadd(i,m,2, 1) = 0.d0
         end do
      end do
     
c     # local method parameters
      if (method(3) .lt. 0) then
         m3 = -1
         m4 = 0
      else
         m3 = method(3)/10
         m4 = method(3) - 10*m3
      endif

c     -----------------------------------------------------------
c     # solve normal Riemann problem and compute Godunov updates
c     -----------------------------------------------------------
     
c     # aux2(1-mbc,1,2) is the start of a 1d array now used by rpn3
     
      if (maux > 0) then
         call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &           aux2(1-mbc,1,2),aux2(1-mbc,1,2),maux,
     &           wave,s,amdq,apdq)
      else
         call rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &           aux2,aux2,maux,
     &           wave,s,amdq,apdq)
      endif
     
c     # Set fadd for the donor-cell upwind method (Godunov)
      do i = 1,mx+1
         do m = 1,meqn
            faddp(i,m) = faddp(i,m) - apdq(i,m)
            faddm(i,m) = faddm(i,m) + amdq(i,m)
         end do
      end do
     
c     # compute maximum wave speed for checking Courant number:
      cfl1d = 0.d0
      do mw = 1,mwaves
         do i = 1,mx+1
c          !  cfl1d = dmax1(cfl1d,dtdx1d(i)*dabs(s(i,mw))) OLD WAY
c          # if s>0 use dtdx1d(i) to compute CFL,
c          # if s<0 use dtdx1d(i-1) to compute CFL:
           cfl1d = dmax1(cfl1d,  dtdx1d(i)*s(i,mw), 
     &                    -dtdx1d(i-1)*s(i,mw))
         end do
      end do

     
      if (method(2).eq.1) then
        go to 130
      endif
    
c         # modify F fluxes for second order q_{xx} correction terms:
c         -----------------------------------------------------------
     
c     # apply limiter to waves:
      if (limit) then 
         call clawpack46_inlinelimiter(maxm,meqn,mwaves,
     &                  mbc,mx,wave,s,mthlim)
      endif
     
      if (use_fwaves .eq. 0) then
          do i = 1, mx+1
     
              dtdxave = 0.5d0 * (dtdx1d(i-1) + dtdx1d(i))
     
              do m = 1,meqn
                  cqxx(i,m) = 0.d0
                  do mw = 1,mwaves
                      cqxx(i,m) = cqxx(i,m) + 0.5d0 * dabs(s(i,mw)) *
     &                     (1.d0 - dabs(s(i,mw))*dtdxave) * wave(i,m,mw)
                  end do
                  faddm(i,m) = faddm(i,m) + cqxx(i,m)
                  faddp(i,m) = faddp(i,m) + cqxx(i,m)
              enddo
          end do
      else
          do i = 1, mx+1

              dtdxave = 0.5d0 * (dtdx1d(i-1) + dtdx1d(i)) 
c
              do m = 1,meqn
                  cqxx(i,m) = 0.d0
                  do mw = 1,mwaves
                      cqxx(i,m) = cqxx(i,m) + 0.5d0 * 
     &                    dsign(1.d0, s(i,mw)) *
     &              (1.d0 - dabs(s(i,mw))*dtdxave) * wave(i,m,mw)
                  end do
                  faddm(i,m) = faddm(i,m) + cqxx(i,m)
                  faddp(i,m) = faddp(i,m) + cqxx(i,m)
              end do
          end do
      endif
c     
 130  continue
c     
      if ( m3 .le. 0) then
          return
      endif
c     
c     --------------------------------------------
c     # TRANSVERSE PROPAGATION
c     --------------------------------------------
c     
c     # split the left-going flux difference into down-going and up-going
c     # flux differences (in the y-direction).
c     
      call rpt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &          aux3,maux,1,amdq,bmamdq,bpamdq)
c
c     # split the right-going flux difference into down-going and up-going
c     # flux differences (in the y-direction).
c
      call rpt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &          aux3,maux,2,apdq,bmapdq,bpapdq)
c
c     # split the left-going flux difference into down-going and up-going
c     # flux differences (in the z-direction).
c
      call rpt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &          aux3,maux,1,amdq,cmamdq,cpamdq)
c
c     # split the right-going flux difference into down-going and up-going
c     # flux differences (in the y-direction).
c
      call rpt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &          aux3,maux,2,apdq,cmapdq,cpapdq)
c
c     
c     # Split the correction wave into transverse propagating waves
c     # in the y-direction and z-direction.
c     
      if (m3 .eq. 2) then
          if (maux .gt. 0) then
c            # The corrections cqxx affect both cell i-1 to left and cell i
c            # to right of interface.  Transverse splitting will affect
c            # fluxes on both sides.
c            # If there are aux arrays, then we must split cqxx twice in
c            # each transverse direction, once with imp=1 and once with imp=2:

c            # imp = 1 or 2 is used to indicate whether we are propagating
c            # amdq or apdq, i.e. cqxxm or cqxxp

c            # in the y-like direction with imp=1
              call rpt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &            aux1,aux2,aux3,maux,1,cqxx,bmcqxxm,bpcqxxm)

c            # in the y-like direction with imp=2
              call rpt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &            aux1,aux2,aux3,maux,2,cqxx,bmcqxxp,bpcqxxp)

c            # in the z-like direction with imp=1
              call rpt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &            aux1,aux2, aux3,maux,1,cqxx,cmcqxxm,cpcqxxm)

c            # in the z-like direction with imp=2
             call rpt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &            aux1,aux2,aux3,maux,2,cqxx,cmcqxxp,cpcqxxp)
           else
c            # aux arrays aren't being used, so we only need to split
c            # cqxx once in each transverse direction and the same result can
c            # presumably be used to left and right.  
c            # Set imp = 0 since this shouldn't be needed in rpt3 in this case.

c            # in the y-like direction 
              call rpt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &              aux1,aux2,aux3,maux,0,cqxx,bmcqxxm,bpcqxxm)

c            # in the z-like direction 
              call rpt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,
     &              aux1,aux2,aux3,maux,0,cqxx,cmcqxxm,cpcqxxm)

c             # use the same splitting to left and right:
              do i = 0,mx+2
                 do m = 1,meqn
                    bmcqxxp(i,m) = bmcqxxm(i,m)
                    bpcqxxp(i,m) = bpcqxxm(i,m)
                    cmcqxxp(i,m) = cmcqxxm(i,m)
                    cpcqxxp(i,m) = cpcqxxm(i,m)
                 enddo
              enddo
          endif
        endif
c
c     --------------------------------------------
c     # modify G fluxes in the y-like direction
c     --------------------------------------------
c     
c     # If the correction wave also propagates in a 3D sense, incorporate
c     # cpcqxx,... into cmamdq, cpamdq, ... so that it is split also.
c     
      if (m4 .eq. 1) then
          do m = 1, meqn
              do i = 0, mx+2
                  cpapdq2(i,m) = cpapdq(i,m)
                  cpamdq2(i,m) = cpamdq(i,m)
                  cmapdq2(i,m) = cmapdq(i,m)
                  cmamdq2(i,m) = cmamdq(i,m)
              end do
          end do
      else if (m4 .eq. 2) then
          do m = 1, meqn
              do i = 0, mx+2
                  cpapdq2(i,m) = cpapdq(i,m) - 3.d0*cpcqxxp(i,m)
                  cpamdq2(i,m) = cpamdq(i,m) + 3.d0*cpcqxxm(i,m)
                  cmapdq2(i,m) = cmapdq(i,m) - 3.d0*cmcqxxp(i,m)
                  cmamdq2(i,m) = cmamdq(i,m) + 3.d0*cmcqxxm(i,m)
              end do
          end do
      endif
c     
c     # The transverse flux differences in the z-direction are split
c     # into waves propagating in the y-direction. If m4 = 2,
c     # then the transverse propagating correction waves in the z-direction
c     # are also split. This yields terms of the form BCAu_{xzy} and
c     # BCAAu_{xxzy}.
c     
      if (m4 .gt. 0) then
          call rptt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,2,2,cpapdq2,bmcpapdq,bpcpapdq)
          call rptt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,1,2,cpamdq2,bmcpamdq,bpcpamdq)
          call rptt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,2,1,cmapdq2,bmcmapdq,bpcmapdq)
          call rptt3(ixyz,2,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,1,1,cmamdq2,bmcmamdq,bpcmamdq)
      endif
c     
c     -----------------------------
c     # The updates for G fluxes :
c     -----------------------------
c     
      meqng_loop : do m = 1,meqn
          ig_loop : do i = 1, mx+1
     
c             # Transverse propagation of the increment waves
c             # between cells sharing interfaces, i.e. the 2D approach.
c             # Yields BAu_{xy}.
c
            gadd(i-1,m,1,0) = gadd(i-1,m,1,0)
     &                      - 0.5d0*dtdx1d(i-1)*bmamdq(i,m)
            gadd(i-1,m,2,0) = gadd(i-1,m,2,0)
     &                      - 0.5d0*dtdx1d(i-1)*bpamdq(i,m)
            gadd(i,m,1,0)   = gadd(i,m,1,0)
     &                      - 0.5d0*dtdx1d(i-1)*bmapdq(i,m)
            gadd(i,m,2,0)   = gadd(i,m,2,0)
     &                      - 0.5d0*dtdx1d(i-1)*bpapdq(i,m)
c     
c             # Transverse propagation of the increment wave (and the
c             # correction wave if m4=2) between cells
c             # only having a corner or edge in common. Yields terms of the
c             # BCAu_{xzy} and BCAAu_{xxzy}.
c     
              if (m4 .gt. 0) then
c     
                gadd(i,m,2,0) = gadd(i,m,2,0)
     &                  + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                  * (bpcpapdq(i,m) - bpcmapdq(i,m))
                gadd(i,m,1,0) = gadd(i,m,1,0)
     &                  + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                  * (bmcpapdq(i,m) - bmcmapdq(i,m))


                gadd(i,m,2,1) = gadd(i,m,2,1)
     &                          - (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bpcpapdq(i,m)
                gadd(i,m,1,1) = gadd(i,m,1,1)
     &                          - (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bmcpapdq(i,m)
                gadd(i,m,2,-1) = gadd(i,m,2,-1)
     &                          + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bpcmapdq(i,m)
                gadd(i,m,1,-1) = gadd(i,m,1,-1)
     &                          + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bmcmapdq(i,m)

                gadd(i-1,m,2,0) = gadd(i-1,m,2,0)
     &                   + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                   * (bpcpamdq(i,m) - bpcmamdq(i,m))
                gadd(i-1,m,1,0) = gadd(i-1,m,1,0)
     &                   + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                   * (bmcpamdq(i,m) - bmcmamdq(i,m))


                gadd(i-1,m,2,1) = gadd(i-1,m,2,1)
     &                          - (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bpcpamdq(i,m)
                gadd(i-1,m,1,1) = gadd(i-1,m,1,1)
     &                          - (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bmcpamdq(i,m)
                gadd(i-1,m,2,-1) = gadd(i-1,m,2,-1)
     &                          + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bpcmamdq(i,m)
                gadd(i-1,m,1,-1) = gadd(i-1,m,1,-1)
     &                          + (1.d0/6.d0)*dtdx1d(i-1)*dtdz
     &                          * bmcmamdq(i,m)     
              endif
c     
c             # Transverse propagation of the correction wave between
c             # cells sharing faces. This gives BAAu_{xxy}.
c     
              if (m3 .ge. 2) then
                gadd(i,m,2,0)   = gadd(i,m,2,0)
     &                         + dtdx1d(i-1)*bpcqxxp(i,m)
                gadd(i,m,1,0)   = gadd(i,m,1,0)
     &                         + dtdx1d(i-1)*bmcqxxp(i,m)
                gadd(i-1,m,2,0) = gadd(i-1,m,2,0)
     &                         - dtdx1d(i-1)*bpcqxxm(i,m)
                gadd(i-1,m,1,0) = gadd(i-1,m,1,0)
     &                         - dtdx1d(i-1)*bmcqxxm(i,m)
              endif
          end do ig_loop
      end do meqng_loop

c     
c     
c     --------------------------------------------
c     # modify H fluxes in the z-like direction
c     --------------------------------------------
c     
c     # If the correction wave also propagates in a 3D sense, incorporate
c     # cqxx into bmamdq, bpamdq, ... so that is is split also.
c     
      if (m4 .eq. 2) then
          do m = 1, meqn
              do i = 0, mx+2
                  bpapdq(i,m) = bpapdq(i,m) - 3.d0*bpcqxxp(i,m)
                  bpamdq(i,m) = bpamdq(i,m) + 3.d0*bpcqxxm(i,m)
                  bmapdq(i,m) = bmapdq(i,m) - 3.d0*bmcqxxp(i,m)
                  bmamdq(i,m) = bmamdq(i,m) + 3.d0*bmcqxxm(i,m)
              end do
         end do
      endif
c     
c     # The transverse flux differences in the y-direction are split
c     # into waves propagating in the z-direction. If m4 = 2,
c     # then the transverse propagating correction waves in the y-direction
c     # are also split. This yields terms of the form BCAu_{xzy} and
c     # BCAAu_{xxzy}.
c     
c     # Note that the output to rptt3 below should logically be named
c     # cmbsasdq and cpbsasdq rather than bmcsasdq and bpcsasdq, but
c     # we are re-using the previous storage rather than requiring new arrays.
c     
      if (m4 .gt. 0) then
          call rptt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,2,2,bpapdq,bmcpapdq,bpcpapdq)
          call rptt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,1,2,bpamdq,bmcpamdq,bpcpamdq)
          call rptt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,2,1,bmapdq,bmcmapdq,bpcmapdq)
          call rptt3(ixyz,3,maxm,meqn,mwaves,mbc,mx,q1d,q1d,aux1,aux2,
     &              aux3,maux,1,1,bmamdq,bmcmamdq,bpcmamdq)
      endif
c     
c     -----------------------------
c     # The updates for H fluxes :
c     -----------------------------
c     
      meqnh_loop : do m = 1,meqn
          ih_loop : do i = 1, mx+1
c     
c             # Transverse propagation of the increment waves
c             # between cells sharing interfaces, i.e. the 2D approach.
c             # Yields CAu_{xy}.
c     
              hadd(i-1,m,1,0) = hadd(i-1,m,1,0)
     &                      - 0.5d0*dtdx1d(i-1)*cmamdq(i,m)
              hadd(i-1,m,2,0) = hadd(i-1,m,2,0)
     &                      - 0.5d0*dtdx1d(i-1)*cpamdq(i,m)
              hadd(i,m,1,0)   = hadd(i,m,1,0)
     &                      - 0.5d0*dtdx1d(i-1)*cmapdq(i,m)
              hadd(i,m,2,0)   = hadd(i,m,2,0)
     &                      - 0.5d0*dtdx1d(i-1)*cpapdq(i,m)
cc     
c             # Transverse propagation of the increment wave (and the
c             # correction wave if m4=2) between cells
c             # only having a corner or edge in common. Yields terms of the
c             # CBAu_{xzy} and CBAAu_{xxzy}.
c     
              if (m4 .gt. 0) then
c     
                hadd(i,m,2,0)  = hadd(i,m,2,0)
     &                  + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                  * (bpcpapdq(i,m) - bpcmapdq(i,m))
                hadd(i,m,1,0)  = hadd(i,m,1,0)
     &                  + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                  * (bmcpapdq(i,m) - bmcmapdq(i,m))


                hadd(i,m,2,1)  = hadd(i,m,2,1)
     &                         - (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                         * bpcpapdq(i,m)
                hadd(i,m,1,1)  = hadd(i,m,1,1)
     &                         - (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                         * bmcpapdq(i,m)
                hadd(i,m,2,-1) = hadd(i,m,2,-1)
     &                         + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                         * bpcmapdq(i,m)
                hadd(i,m,1,-1) = hadd(i,m,1,-1)
     &                         + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                         * bmcmapdq(i,m)
c
                hadd(i-1,m,2,0)  = hadd(i-1,m,2,0)
     &                   + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                   * (bpcpamdq(i,m) - bpcmamdq(i,m))
                hadd(i-1,m,1,0)  = hadd(i-1,m,1,0)
     &                   + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                   * (bmcpamdq(i,m) - bmcmamdq(i,m))


                hadd(i-1,m,2,1)  = hadd(i-1,m,2,1)
     &                           - (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                           * bpcpamdq(i,m)
                hadd(i-1,m,1,1)  = hadd(i-1,m,1,1)
     &                           - (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                           * bmcpamdq(i,m)
                hadd(i-1,m,2,-1) = hadd(i-1,m,2,-1)
     &                           + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                           * bpcmamdq(i,m)
                hadd(i-1,m,1,-1) = hadd(i-1,m,1,-1)
     &                           + (1.d0/6.d0)*dtdx1d(i-1)*dtdy
     &                           * bmcmamdq(i,m)
c
c     
              endif
c     
c             # Transverse propagation of the correction wave between
c             # cells sharing faces. This gives CAAu_{xxy}.
c     
              if (m3 .ge. 2) then
                hadd(i,m,2,0)   = hadd(i,m,2,0)
     &                         + dtdx1d(i-1)*cpcqxxp(i,m)
                hadd(i,m,1,0)   = hadd(i,m,1,0)
     &                         + dtdx1d(i-1)*cmcqxxp(i,m)
                hadd(i-1,m,2,0) = hadd(i-1,m,2,0)
     &                         - dtdx1d(i-1)*cpcqxxm(i,m)
                hadd(i-1,m,1,0) = hadd(i-1,m,1,0)
     &                         - dtdx1d(i-1)*cmcqxxm(i,m)
              endif
          end do ih_loop
      end do meqnh_loop  
c
      return
      end
