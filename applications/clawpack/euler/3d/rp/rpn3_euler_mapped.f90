subroutine clawpack46_rpn3_mapped(ixyz,maxm,meqn,mwaves,maux,mbc,mx,& 
    ql_cart,qr_cart, auxl,auxr,fwave_cart,s,amdq_cart,apdq_cart)
!!
!!     # Roe-solver for the Euler equations
!!      -----------------------------------------------------------
!!
!!     # solve Riemann problems along one slice of data.
!!     # This data is along a slice in the x-direction if ixyz=1
!!     #                               the y-direction if ixyz=2.
!!     #                               the z-direction if ixyz=3.
!!
!!     # On input, ql contains the state vector at the left edge of each cell
!!     #           qr contains the state vector at the right edge of each cell
!!
!!     # On output, wave contains the waves, s the speeds,
!!     # and amdq, apdq the left-going and right-going flux differences,
!!     # respectively.
!!
!!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!!     #                                    and right state ql(i,:)
!!     # From the basic clawpack routines, this routine is called with ql = qr
!!

    use setprob_mod, only : gamma, gamma1, mcapa
    implicit none

    double precision, external :: suppow
    double precision, external :: cflsupbee

    !! Riemann solvers all use (meqn,i,j,k) indexing (5.x indexing)
    integer ixyz, maxm, meqn, mwaves, mbc, mx, maux
    double precision fwave_cart(meqn,mwaves,1-mbc:maxm+mbc)
    double precision          s(mwaves,1-mbc:maxm+mbc)
    double precision    ql_cart(meqn,1-mbc:maxm+mbc)
    double precision    qr_cart(meqn,1-mbc:maxm+mbc)
    double precision  amdq_cart(meqn,1-mbc:maxm+mbc)
    double precision  apdq_cart(meqn,1-mbc:maxm+mbc)
    double precision       auxl(maux,1-mbc:maxm+mbc)
    double precision       auxr(maux,1-mbc:maxm+mbc)

    double precision dtcom, dxcom, dycom, dzcom, tcom
    integer icom, jcom, kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    double precision ql(5),qr(5), s_rot(3), fwave(5,3)
    double precision amdq(5), apdq(5), rot(9)

    !! Roe data
    double precision   u2v2w2, uvw2, u,v,w   ! uvw2 is square of flow, u2v2w2 is square of roe
    double precision    uvw(3,1-mbc:maxm+mbc)
    double precision   enth(  1-mbc:maxm+mbc)

    !! Beta and delta storage for limiting scheme
    double precision      b(5,1-mbc:maxm+mbc)
    double precision    bId(5)
    double precision  delta(5,1-mbc:maxm+mbc)

    !! Data for CFL-dependent limiting scheme
    double precision  dtdx1d( 1-mbc:maxm+mbc), dx
    integer              Id(3,1-mbc:maxm+mbc)
    double precision    CFL(3,1-mbc:maxm+mbc), CFLup
    double precision   wlim, wlimEu(5)

    integer i, j, ii, m, info, mws, locrot, locarea
    double precision a2, a, g1a2, euv
    double precision pr, pl, rhsq2, rhsqrtl, rhsqrtr, ur, ul
    double precision rho_im1, pim1, cim1, s0
    double precision rho1, rhou1, rhov1, rhow1, en1, p1, c1
    double precision s1, sfract, rhoi, pi, ci, s3
    double precision rho2, rhou2, rhov2, rhow2, en2, p2, c2
    double precision s2, df, area

    double precision theta, beta

    logical efix, nonlocal

    data efix /.true./    !# use entropy fix for transonic rarefactions

    data nonlocal /.false./ !# use nonlocal correction for CFL, maybe not helpful

    !! CFL Superbee beta-theta limiter
    !! (Theta = 1, Beta = 2/3 Recommended)
    theta = 0.95d0
    beta  = 0.666666666666666d0

    if (mwaves .ne. 3) then
        write(6,*) '*** Should have mwaves=3 for this Riemann solver'
        stop
    endif

    !! Locations in aux arrays
    call get_aux_locations_n(ixyz,mcapa,locrot,locarea)

    do i = 2-mbc, mx+mbc

        do m = 1,meqn
            ql(m) = ql_cart(m,i)
            qr(m) = qr_cart(m,i-1)
        enddo

        !!# Use Roe averaged values for normal solves
        if (qr(1) .lt. 0) then
            write(6,*) 'qr(1) < 0; ', qr(1)
            stop
        endif
         if (ql(1) .lt. 0) then
            write(6,*) 'ql(1) < 0; ', ql(1)
            stop
        endif

        !! We can use unrotated values here, since we only use 
        !! norm of velocity.
        rhsqrtl = sqrt(qr(1))
        rhsqrtr = sqrt(ql(1))
        rhsq2 = rhsqrtl + rhsqrtr

        uvw2 = qr(2)**2 + qr(3)**2 + qr(4)**2
        pl = gamma1*(qr(5) - 0.5d0*uvw2/qr(1))

        uvw2 = ql(2)**2 + ql(3)**2 + ql(4)**2
        pr = gamma1*(ql(5) - 0.5d0*uvw2/ql(1))

        enth(i) = (((qr(5)+pl)/rhsqrtl + (ql(5)+pr)/rhsqrtr)) / rhsq2


        !! --------------- Use rotated values in Riemann solve ------------
        do j = 1,9
            rot(j) = auxl(locrot+j-1,i)
        enddo

        do j = 1,3
            uvw(j,i) = (qr(j+1)/rhsqrtl + ql(j+1)/rhsqrtr) / rhsq2
        enddo

        call rotate3(rot,uvw(:,i))

        call rotate3(rot,ql(2))
        call rotate3(rot,qr(2))

        dx = dxcom*merge(1.d0,0.d0,(ixyz.EQ.1)) + &
             dycom*merge(1.d0,0.d0,(ixyz.EQ.2)) + &
             dzcom*merge(1.d0,0.d0,(ixyz.EQ.3))
        dtdx1d(i) = (dtcom/dx)/auxl(locarea,i)

        !! # Calculate flux differences (f-wave method)
        ur = ql(2)/ql(1)
        ul = qr(2)/qr(1)
        delta(1,i) = ql(2) - qr(2)
        delta(2,i) = (ql(2)*ur + pr) - (qr(2)*ul + pl)
        delta(3,i) = ur*ql(3) - ul*qr(3)
        delta(4,i) = ur*ql(4) - ul*qr(4)
        delta(5,i) = (ql(5)+pr)*ur - (qr(5)+pl)*ul   

        !! --------------- Use rotated values in Riemann solve ------------

        !! # Solve normal Riemann problem
        call solve_riemann(gamma1, uvw(:,i), enth(i), delta(:,i), fwave, s_rot, info)

        !! @ Output error diagnostics
        if (info > 0) then
            write(6,'(A)') 'In rpn3'
            write(6,1001) 'rhol = ',qr_cart(1,i-1)
            write(6,1001) 'rhor = ',ql_cart(1,i)
            write(6,1001) 'ul   = ',qr_cart(2,i-1)/qr_cart(1,i-1)
            write(6,1001) 'ur   = ',ql_cart(2,i)/ql_cart(1,i)
            write(6,1001) 'pl   = ',pl
            write(6,1001) 'pr   = ',pr
            do ii = 1,5
               write(6,1001) 'ql = ', ql(ii)
            enddo
            write(6,*) ' '
            do ii = 1,5
               write(6,1001) 'qr = ', qr(ii)
            enddo
            write(6,*) ' '
            do ii = 1,9
               write(6,1001) 'rot = ', rot(ii)
            enddo
            write(6,*) 'locarea = ', locarea
            write(6,*) 'locrot = ', locrot
            write(6,*) 'ixyz = ', ixyz
            write(6,'(A,I5)') 'i = ', i
            stop
        endif
1001    format(A,2E16.8)

        !! # Store betas from the waves
        b(1,i) = fwave(1,1)
        b(4,i) = fwave(1,2)
        b(2,i) = fwave(3,2)-b(4,i)*uvw(2,i)
        b(3,i) = fwave(4,2)-b(4,i)*uvw(3,i)
        b(5,i) = fwave(1,3)


        !! # Entropy Fix

        if (efix) then

            !! # check 1-wave:
            !! ---------------

            rho_im1 = qr(1)
            uvw2 = qr(2)**2 + qr(3)**2 + qr(4)**2
            pim1 = gamma1*(qr(5) - 0.5d0*uvw2/rho_im1)
            cim1 = sqrt(gamma*pim1/rho_im1)

            !! # u-c in left state (cell i-1)
            s0 = qr(2)/rho_im1 - cim1 


            !! # check for fully supersonic case:
            !! right_going = .false.
            if (s0 .ge. 0.d0 .and. s_rot(1) .gt. 0.d0) then
                !! # everything is right-going
                amdq = 0.d0
                !! right_going = .true.
            else
                rho1 = qr(1) + fwave(1,1)/s_rot(1)
                rhou1 = qr(2) + fwave(2,1)/s_rot(1)
                rhov1 = qr(3) + fwave(3,1)/s_rot(1)
                rhow1 = qr(4) + fwave(4,1)/s_rot(1)
                en1 = qr(5) + fwave(5,1)/s_rot(1)
                uvw2 = rhou1**2 + rhov1**2 + rhow1**2
                p1 = gamma1*(en1 - 0.5d0*uvw2/rho1)
                c1 = sqrt(gamma*p1/rho1)

                !! # u-c to right of 1-wave_rot
                s1 = rhou1/rho1 - c1     

                !!left_going = .false.
                if (s0 .lt. 0.d0 .and. s1 .gt. 0.d0) then
                    !! # transonic rarefaction in the 1-wave
                    sfract = s0 * (s1-s_rot(1)) / (s1-s0)
                else if (s_rot(1) .lt. 0.d0) then
                    !! # 1-wave is leftgoing
                    sfract = s_rot(1)
                else
                    !! # 1-wave is right-going
                    !! # this shouldn't happen since s0 < 0
                    sfract = 0.d0         
                endif
                do m = 1,meqn
                    amdq(m) = sfract*fwave(m,1)/s_rot(1)
                end do
            endif

            !! # check 2-wave:
            !! ---------------

            if (s_rot(2) .ge. 0.d0) then 
                !! # 2-,3- and 4- waves are rightgoing
                go to 200 
            endif

            do m=1,meqn
                amdq(m) = amdq(m) + fwave(m,2)
            enddo

           !! # check 3-wave:
           !! ---------------

            rhoi = ql(1)
            uvw2 = ql(2)**2 + ql(3)**2 + ql(4)**2
            pi = gamma1*(ql(5) - 0.5d0*(uvw2) / rhoi)
            ci = dsqrt(gamma*pi/rhoi)

            !! # u+c in right state  (cell i)
            s3 = ql(2)/rhoi + ci     

            rho2 = ql(1) - fwave(1,3)/s_rot(3)
            rhou2 = ql(2) - fwave(2,3)/s_rot(3)
            rhov2 = ql(3) - fwave(3,3)/s_rot(3)
            rhow2 = ql(4) - fwave(4,3)/s_rot(3)
            en2 = ql(5) - fwave(5,3)/s_rot(3)

            uvw2 = rhou2**2 + rhov2**2 + rhow2**2
            p2 = gamma1*(en2 - 0.5d0*(uvw2)/rho2)
            c2 = dsqrt(gamma*p2/rho2)
            s2 = rhou2/rho2 + c2        !# u+c to left of 3-wave
            if (s2 .lt. 0.d0 .and. s3 .gt. 0.d0 ) then
                !! # transonic rarefaction in the 3-wave
                sfract = s2 * (s3-s_rot(3)) / (s3-s2)
            else if (s_rot(3) .lt. 0.d0) then
                !! # 3-wave is leftgoing
                sfract = s_rot(3)
            else
                !! # 3-wave is rightgoing
                go to 200
            endif

            do m = 1,5
                amdq(m) = amdq(m) + sfract*fwave(m,3)/s_rot(3)
            enddo

  200       continue

           !! # compute the rightgoing flux differences:
           !! # df = SUM s*wave   is the total flux difference and apdq = df -
           !! # amdq

            do m = 1,meqn
               df = 0.d0
               do mws = 1,mwaves
                    df = df + fwave(m,mws)
               enddo
               apdq(m) = df - amdq(m)
            enddo

        else
            !! # No entropy fix
            do m=1,meqn
                amdq(m) = 0.d0
                apdq(m) = 0.d0
                do  mws = 1,mwaves
                    if (s_rot(mws).LT.0.d0) then
                        amdq(m) = amdq(m) + fwave(m,mws)
                    else     
                        apdq(m) = apdq(m) + fwave(m,mws)
                    endif
                enddo
            enddo
        endif  !! end of entropy fix

        !! # Store apdq and amdq

        area = auxl(locarea,i)
        call rotate3_tr(rot,apdq(2))
        call rotate3_tr(rot,amdq(2))
        do m = 1,meqn
            apdq_cart(m,i) = area*apdq(m)
            amdq_cart(m,i) = area*amdq(m)
        end do
        do mws = 1,mwaves
            s(mws,i) = area*s_rot(mws)
        enddo
    enddo  !! end of i loop over 1d sweep array

    !! # Limiting loops

    do i = 0, mx+1

        !! # Calculate upwind direction and local CFL
        do mws = 1, 3
            Id(mws,i)  = int(i-dsign(1.d0,s(mws,i)))
            CFL(mws,i) = dabs(s(mws,i)*dtdx1d(max(i,Id(mws,i))))
        enddo
    enddo !! end calculation

    wlimEu(:)=0.d0    
    do i = 1, mx+1 
        u = uvw(1,i)
        v = uvw(2,i)
        w = uvw(3,i)

        u2v2w2 = u**2 + v**2 + w**2
        a2 = gamma1*(enth(i) - 0.5d0*u2v2w2)
        a = sqrt(a2)
        g1a2 = gamma1 / a2
        euv = enth(i) - u2v2w2
        
     !! # Acoustic Wave 1
        if (b(1,i).NE.0.d0) then
            bId(4) = g1a2 * (euv*delta(1,Id(1,i))               &
                   + u*delta(2,Id(1,i)) + v*delta(3,Id(1,i))    &
                   + w*delta(4,Id(1,i)) - delta(5,Id(1,i)))
            bId(5) = (delta(2,Id(1,i)) + (a-u)*delta(1,Id(1,i)) &
                   - a*bId(4)) / (2.d0*a)
            bId(1) = delta(1,Id(1,i)) - bId(4) - bId(5)
            if (nonlocal) then
                CFLup=CFL(1,Id(1,i))
            else
                CFLup=CFL(1,i)
            endif
            wlimEu(1) = suppow(bId(1)/b(1,i),CFL(1,i),CFLup)
        endif

     !! # Shear Wave 2
        if (b(2,i).NE.0.d0) then
            bId(2) = delta(3,Id(2,i)) - v*delta(1,Id(2,i))
            if (nonlocal) then
                CFLup=CFL(2,Id(2,i))
            else
                CFLup=CFL(2,i)
            endif
            wlimEu(2) = cflsupbee(theta,beta,bId(2)/b(2,i),CFL(2,i),CFLup)
        endif

     !! # Shear Wave 3
        if (b(3,i).NE.0.d0) then
            bId(3) = delta(4,Id(2,i)) - w*delta(1,Id(2,i))
            if (nonlocal) then
                CFLup=CFL(2,Id(2,i))
            else
                CFLup=CFL(2,i)
            endif
            wlimEu(3) = cflsupbee(theta,beta,bId(3)/b(3,i),CFL(2,i),CFLup)
        endif

     !! # Entroppy Wave 4
        if (b(4,i).NE.0.d0) then
            bId(4) = g1a2 * (euv*delta(1,Id(2,i))             &
                + u*delta(2,Id(2,i)) + v*delta(3,Id(2,i))     &
                + w*delta(4,Id(2,i)) - delta(5,Id(2,i)))
            if (nonlocal) then
                CFLup=CFL(2,Id(2,i))
            else
                CFLup=CFL(2,i)
            endif
            wlimEu(4) = cflsupbee(theta,beta,bId(4)/b(4,i),CFL(2,i),CFLup)
        endif

     !! # Acoustic Wave 5
        if (b(5,i).NE.0.d0) then
            bId(4) = g1a2 * (euv*delta(1,Id(3,i))             &
                   + u*delta(2,Id(3,i)) + v*delta(3,Id(3,i))  &
                   + w*delta(4,Id(3,i)) - delta(5,Id(3,i)))
            bId(5) = (delta(2,Id(3,i)) + (a-u)*delta(1,Id(3,i))  &
                   - a*bId(4)) / (2.d0*a)
            if (nonlocal) then
                CFLup=CFL(3,Id(3,i))
            else
                CFLup=CFL(3,i)
            endif
            wlimEu(5) = suppow(bId(5)/b(5,i),CFL(3,i),CFLup)
        endif

        !! # Apply Limiter with rotations

        do j = 1,9
            rot(j) = auxl(locrot+j-1,i)
        enddo
        area = auxl(locarea,i)
  
        do m = 1, 5
            bId(m) = wlimEu(m)*b(m,i)
        enddo

        fwave(1,1) = bId(1)
        fwave(2,1) = bId(1)*(u-a)
        fwave(3,1) = bId(1)*v
        fwave(4,1) = bId(1)*w
        fwave(5,1) = bId(1)*(enth(i) - u*a)

        fwave(1,2) = bId(4)
        fwave(2,2) = bId(4)*u
        fwave(3,2) = bId(4)*v + bId(2)
        fwave(4,2) = bId(4)*w + bId(3)
        fwave(5,2) = bId(4)*0.5d0*u2v2w2  + bId(2)*v + bId(3)*w

        fwave(1,3) = bId(5)
        fwave(2,3) = bId(5)*(u+a)
        fwave(3,3) = bId(5)*v
        fwave(4,3) = bId(5)*w
        fwave(5,3) = bId(5)*(enth(i)+u*a)

        !! # Return f-waves
        do mws = 1, 3 
            call rotate3_tr(rot,fwave(2,mws))
            do m = 1, 5
                fwave_cart(m,mws,i) = area*fwave(m,mws)
            enddo 
        enddo


    enddo !! end of i loop over 1d sweep array for limiting


999 return
end subroutine clawpack46_rpn3_mapped


!! # Riemann Solver

subroutine solve_riemann(gamma1,uvw,enth,delta,wave,s,info)
    implicit none

    double precision enth, uvw(3), delta(5)
    double precision wave(5,3), s(3), gamma1
    double precision u2v2w2, a2, a, g1a2, euv
    double precision a1, a3, a4, a5,u,v,w
    integer info

    info = 0

    u = uvw(1)
    v = uvw(2)
    w = uvw(3)

    u2v2w2 = u**2 + v**2 + w**2
    a2 = gamma1*(enth - 0.5d0*u2v2w2)
    if (a2 .lt. 0.d0) then
        write(6,*) 'solve_riemann : '
        write(6,*) 'a2 .lt. 0; ', a2
        write(6,*) enth, u2v2w2
        info = 1
        return
    endif
    a = sqrt(a2)
    g1a2 = gamma1 / a2
    euv = enth - u2v2w2

    a4 = g1a2 * (euv*delta(1) + u*delta(2) + v*delta(3) + w*delta(4) &
                 - delta(5))
    a2 = delta(3) - v*delta(1)
    a3 = delta(4) - w*delta(1)
    a5 = (delta(2) + (a-u)*delta(1) - a*a4) / (2.d0*a)
    a1 = delta(1) - a4 - a5

    !! # Compute the waves.
    !! # Note that the 2-wave, 3-wave and 4-wave travel at the same speed
    !! # and are lumped together in wave(.,.,2).  The 5-wave is then stored
    !! # in wave(.,.,3).

    wave(1,1) = a1
    wave(2,1) = a1*(u-a)
    wave(3,1) = a1*v
    wave(4,1) = a1*w
    wave(5,1) = a1*(enth - u*a)
    s(1) = u - a

    wave(1,2) = a4
    wave(2,2) = a4*u
    wave(3,2) = a4*v + a2
    wave(4,2) = a4*w + a3
    wave(5,2) = a4*0.5d0*u2v2w2  + a2*v + a3*w
    s(2) = u

    wave(1,3) = a5
    wave(2,3) = a5*(u+a)
    wave(3,3) = a5*v
    wave(4,3) = a5*w
    wave(5,3) = a5*(enth+u*a)
    s(3) = u + a

end subroutine solve_riemann


!! ----------------------------------------------------
!!   Limiter functions for third-order MAGIC Schemes
!! ----------------------------------------------------

double precision function cflsupbee(theta,beta,r,CFL,CFLup)
      implicit double precision (a-h,o-z)

      !! # CFL-Superbee Limiter:
      CFLe = dmax1(0.001d0,dmin1(0.999d0,CFL))
      CFLupe = dmax1(0.001d0,dmin1(0.999d0,CFLup))
      weight = (1.d0 - CFLupe) / (1.d0 - CFLe)
      s1 = theta*weight*2.d0/CFLe
      s2 = (1.d0 + CFLe) / 3.d0
      phimax = theta * 2.d0 / (1.d0 - CFLe)
      ultra = dmax1(0.d0,dmin1(s1*r,phimax))
      c1 = 1.d0 + (s2 - beta/2.d0) * (r-1.d0)
      c2 = 1.d0 + (s2 + beta/2.d0) * (r-1.d0)
      cflsupbee = dmax1(0.d0, dmin1(ultra,dmax1(c1,c2)))

return
end


double precision function suppow(r,CFL,CFLup)
     implicit double precision (a-h,o-z)

      !! # Superpower Limiter:
      CFLe = dmax1(0.001d0,dmin1(0.999d0,CFL))
      CFLupe = dmax1(0.001d0,dmin1(0.999d0,CFLup))
      s2 = (1.d0 + CFLe) / 3.d0
      s3 = 1.d0 - s2
      weight = (1.d0 - CFLupe) / (1.d0 - CFLe)
      if (r.LE.1.d0) then
         pp = weight*(2.d0/CFLe)*2.d0*s3
      else
         pp = dabs(2.d0/(1.d0-CFLe))*2.d0*s2
      endif
      rabs = dabs(r)
      rfrac = dabs((1.d0-rabs)/(1.d0+rabs))
      signum = dmax1(0.d0,dsign(1.d0,r))
      suppow = signum * (s3+s2*r) * (1.d0-rfrac**pp)

return
end
