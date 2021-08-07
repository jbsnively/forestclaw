subroutine clawpack46_step3(maxm,meqn,maux,mbc,mx,my,mz, & 
    qold,aux,dx,dy,dz,dt,cflgrid, & 
    fm,fp,gm,gp,hm,hp, & 
    faddm,faddp,gadd,hadd, & 
    q1d,dtdx1d,dtdy1d,dtdz1d, & 
    aux1,aux2,aux3,work,mwork,rpn3,rpt3,rptt3, &
    mwaves,mcapa,method,mthlim,use_fwaves,ierror)

    !!  ==================================================================

    !!  # clawpack routine ...  modified for AMRCLAW

    !!
    !!  # Take one time step, updating q.
    !!  # On entry, qold gives
    !!  #    initial data for this step
    !!  #    and is unchanged in this version.
    !!
    !!  # fm, fp are fluxes to left and right of single cell edge
    !!
    !!  # See the flux3 documentation for more information.
    !!
    !!      
    implicit none
    external rpn3,rpt3,rptt3

    integer :: maxm, meqn, maux, mbc, use_fwaves
    integer :: mx,my,mz, mwaves, mcapa, ierror
    integer :: mthlim(mwaves), method(7)
    double precision :: dx,dy,dz,dt,cflgrid

    !! Claw 4 indexing
    double precision :: qold(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc, meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc, maux)

    !! Claw 5 indexing
    double precision :: q1d(1-mbc:maxm+mbc,meqn)
    double precision ::  fm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  fp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  gm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  gp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  hm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision ::  hp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    double precision :: faddm(1-mbc:maxm+mbc,meqn)
    double precision :: faddp(1-mbc:maxm+mbc,meqn)
    double precision ::  gadd(1-mbc:maxm+mbc,meqn, 2, -1:1)
    double precision ::  hadd(1-mbc:maxm+mbc,meqn, 2, -1:1)

    double precision :: aux1(1-mbc:maxm+mbc,maux, 3)
    double precision :: aux2(1-mbc:maxm+mbc,maux, 3)
    double precision :: aux3(1-mbc:maxm+mbc,maux, 3)
    double precision :: dtdx1d(1-mbc:maxm+mbc)
    double precision :: dtdy1d(1-mbc:maxm+mbc)
    double precision :: dtdz1d(1-mbc:maxm+mbc)
    double precision :: work(mwork)

    double precision :: dtcom, dxcom, dycom, dzcom, tcom
    integer :: icom,jcom,kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    integer :: i0wave, i0s, i0amdq, i0apdq, i0cqxx, i0bmamdq, i0bmapdq
    integer :: i0bpamdq, i0bpapdq, i0cmamdq, i0cmapdq, i0cpamdq
    integer :: i0cpapdq, i0cmamdq2, i0cmapdq2, i0cpamdq2, i0cpapdq2
    integer :: i0bmcqxxm, i0bmcqxxp, i0bpcqxxm, i0bpcqxxp, i0cmcqxxm
    integer :: i0cmcqxxp, i0cpcqxxm, i0cpcqxxp

    integer :: i0bmcmamdq, i0bmcmapdq, i0bpcmamdq, i0bpcmapdq
    integer :: i0bmcpamdq, i0bmcpapdq, i0bpcpamdq, i0bpcpapdq, iused
    integer :: mwork

    double precision :: dtdx, dtdy, dtdz, cfl1d
    integer :: i,j,k,m, ma


    !!  # store mesh parameters that may be needed in Riemann solver but not
    !!  # passed in...
    dxcom = dx
    dycom = dy
    dtcom = dt

    !!  # partition work array into pieces needed for local storage in
    !!  # flux3 routine.  Find starting index of each piece:

    i0wave     = 1
    i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
    i0amdq     = i0s        + (maxm+2*mbc)*mwaves
    i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
    i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
    i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
    i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
    i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
    i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
    i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
    i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
    i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
    i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
    i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
    i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
    i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
    i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn

    i0bmcqxxm   = i0cpapdq2  + (maxm+2*mbc)*meqn
    i0bmcqxxp   = i0bmcqxxm  + (maxm+2*mbc)*meqn
    i0bpcqxxm   = i0bmcqxxp  + (maxm+2*mbc)*meqn
    i0bpcqxxp   = i0bpcqxxm  + (maxm+2*mbc)*meqn
    i0cmcqxxm   = i0bpcqxxp  + (maxm+2*mbc)*meqn
    i0cmcqxxp   = i0cmcqxxm  + (maxm+2*mbc)*meqn
    i0cpcqxxm   = i0cmcqxxp  + (maxm+2*mbc)*meqn
    i0cpcqxxp   = i0cpcqxxm  + (maxm+2*mbc)*meqn

    work(i0bmcqxxm) = i0bmcqxxm
    work(i0bmcqxxp) = i0bmcqxxp
    work(i0bpcqxxm) = i0bpcqxxm
    work(i0bpcqxxp) = i0bpcqxxp
    work(i0cmcqxxm) = i0cmcqxxm
    work(i0cmcqxxp) = i0cmcqxxp
    work(i0cpcqxxm) = i0cpcqxxm
    work(i0cpcqxxp) = i0cpcqxxp

    i0bmcmamdq = i0cpcqxxp  + (maxm+2*mbc)*meqn
    i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
    i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
    i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
    i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
    i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
    i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
    i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
    iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1

    if (iused .gt. mwork) then
        !! # This shouldn't happen if mwork is set properly in stepgrid3
        write(6,*) '*** not enough work space in step3'
        write(6,*) '*** check parameter mwork in fc3d_clawpack46_step3'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif

    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    do m = 1,meqn
       do k = 1-mbc,mz+mbc
          do j = 1-mbc,my+mbc
             do i = 1-mbc,mx+mbc
                fm(i,j,k,m) = 0.d0
                fp(i,j,k,m) = 0.d0
                gm(i,j,k,m) = 0.d0
                gp(i,j,k,m) = 0.d0
                hm(i,j,k,m) = 0.d0
                hp(i,j,k,m) = 0.d0
             end do
          end do
       end do
    end do

    if (mcapa.eq.0) then
        !! # no capa array:
        do i = 1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
        end do
    end if

    !!  ==================
    !!  # perform x-sweeps
    !!  ==================


    kx_loop : do k = 0,mz+1
        jx_loop : do j = 0,my+1
            do i = 1-mbc, mx+mbc 
                do m = 1,meqn
                    !! # copy data along a slice into 1d array:
                    q1d(m,i) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do i = 1-mbc, mx+mbc
                    dtdx1d(i) = dtdx / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                 do i = 1-mbc, mx+mbc
                    do ma = 1,maux
                        aux1(ma,i,1) = aux(i,j-1,k-1,ma)
                        aux1(ma,i,2) = aux(i,j-1,k,ma)
                        aux1(ma,i,3) = aux(i,j-1,k+1,ma)
                        aux2(ma,i,1) = aux(i,j,k-1,ma)
                        aux2(ma,i,2) = aux(i,j,k,ma)
                        aux2(ma,i,3) = aux(i,j,k+1,ma)
                        aux3(ma,i,1) = aux(i,j+1,k-1,ma)
                        aux3(ma,i,2) = aux(i,j+1,k,ma)
                        aux3(ma,i,3) = aux(i,j+1,k+1,ma)
                    end do
                end do
            endif

            !! # Store the value of j and k along this slice in the common block
            !! # comxyt in case it is needed in the Riemann solver (for
            !!  # variable coefficient problems)

            jcom = j
            kcom = k

            !! # compute modifications fadd, gadd and hadd to fluxes along
            !! # this slice:

            call flux3(1,maxm,meqn,maux,mbc,mx,  & 
                    q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,  & 
                    faddm,faddp,gadd,hadd,cfl1d,  & 
                    work(i0wave),work(i0s),work(i0amdq),  & 
                    work(i0apdq),work(i0cqxx),  & 
                    work(i0bmamdq),work(i0bmapdq),   & 
                    work(i0bpamdq),work(i0bpapdq),   & 
                    work(i0cmamdq),work(i0cmapdq),   & 
                    work(i0cpamdq),work(i0cpapdq),   & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), & 

                    work(i0bmcqxxm),work(i0bmcqxxp), & 
                    work(i0bpcqxxm),work(i0bpcqxxp), &
                    work(i0cmcqxxm),work(i0cmcqxxp), &
                    work(i0cpcqxxm),work(i0cpcqxxp), &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim,use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)

            !! # update fluxes for use in AMR:

             do i = 1,mx+1
                do m = 1,meqn
                    fm(i,j,k,m) = fm(i,j,k,m) + faddm(m,i)
                    fp(i,j,k,m) = fp(i,j,k,m) + faddp(m,i)

                    gm(i,j  ,k-1,m) = gm(i,j  ,k-1,m) + gadd(i,m,1,-1)
                    gp(i,j  ,k-1,m) = gp(i,j  ,k-1,m) + gadd(i,m,1,-1)
                    gm(i,j  ,k  ,m) = gm(i,j  ,k  ,m) + gadd(i,m,1, 0)
                    gp(i,j  ,k  ,m) = gp(i,j  ,k  ,m) + gadd(i,m,1, 0)
                    gm(i,j  ,k+1,m) = gm(i,j  ,k+1,m) + gadd(i,m,1, 1)
                    gp(i,j  ,k+1,m) = gp(i,j  ,k+1,m) + gadd(i,m,1, 1)

                    gm(i,j+1,k-1,m) = gm(i,j+1,k-1,m) + gadd(i,m,2,-1)
                    gp(i,j+1,k-1,m) = gp(i,j+1,k-1,m) + gadd(i,m,2,-1)
                    gm(i,j+1,k  ,m) = gm(i,j+1,k  ,m) + gadd(i,m,2, 0)
                    gp(i,j+1,k  ,m) = gp(i,j+1,k  ,m) + gadd(i,m,2, 0)
                    gm(i,j+1,k+1,m) = gm(i,j+1,k+1,m) + gadd(i,m,2, 1)
                    gp(i,j+1,k+1,m) = gp(i,j+1,k+1,m) + gadd(i,m,2, 1)

                    hm(i,j-1,k,m) = hm(i,j-1,k,m) + hadd(i,m,1,-1)
                    hp(i,j-1,k,m) = hp(i,j-1,k,m) + hadd(i,m,1,-1)
                    hm(i,j  ,k,m) = hm(i,j,  k,m) + hadd(i,m,1, 0)
                    hp(i,j  ,k,m) = hp(i,j,  k,m) + hadd(i,m,1, 0)
                    hm(i,j+1,k,m) = hm(i,j+1,k,m) + hadd(i,m,1, 1)
                    hp(i,j+1,k,m) = hp(i,j+1,k,m) + hadd(i,m,1, 1)

                    hm(i,j-1,k+1,m) = hm(i,j-1,k+1,m) + hadd(i,m,2,-1)
                    hp(i,j-1,k+1,m) = hp(i,j-1,k+1,m) + hadd(i,m,2,-1)
                    hm(i,j  ,k+1,m) = hm(i,j,  k+1,m) + hadd(i,m,2, 0)
                    hp(i,j  ,k+1,m) = hp(i,j,  k+1,m) + hadd(i,m,2, 0)
                    hm(i,j+1,k+1,m) = hm(i,j+1,k+1,m) + hadd(i,m,2, 1)
                    hp(i,j+1,k+1,m) = hp(i,j+1,k+1,m) + hadd(i,m,2, 1)
                end do
            end do
        end do jx_loop
    end do kx_loop


    !!  ==================
    !!  # perform y sweeps
    !!  ==================

    ky_loop : do k = 0, mz+1
        iy_loop : do i = 0, mx+1

            do j = 1-mbc, my+mbc
                do m = 1,meqn
                    !! # copy data along a slice into 1d array:
                    q1d(m,j) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do j = 1-mbc, my+mbc
                    dtdy1d(j) = dtdy / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                do j = 1-mbc, my+mbc
                    do ma = 1,maux
                        aux1(ma,j,1) = aux(i-1,j,k-1,ma)
                        aux1(ma,j,2) = aux(i,  j,k-1,ma)
                        aux1(ma,j,3) = aux(i+1,j,k-1,ma)
                        aux2(ma,j,1) = aux(i-1,j,k,  ma)
                        aux2(ma,j,2) = aux(i,  j,k,  ma)
                        aux2(ma,j,3) = aux(i+1,j,k,  ma)
                        aux3(ma,j,1) = aux(i-1,j,k+1,ma)
                        aux3(ma,j,2) = aux(i,  j,k+1,ma)
                        aux3(ma,j,3) = aux(i+1,j,k+1,ma)
                    end do
                end do
            endif

            !! # Store the value of i and k along this slice in the common block
            !! # comxyzt in case it is needed in the Riemann solver (for
            !! # variable coefficient problems)

            icom = i
            kcom = k

            !! # compute modifications fadd, gadd and hadd to fluxes along this
            !! # slice:

            call flux3(2,maxm,meqn,maux,mbc,my, & 
                    q1d,dtdy1d,dtdz,dtdx,aux1,aux2,aux3, & 
                    faddm,faddp,gadd,hadd,cfl1d, & 
                    work(i0wave),work(i0s),work(i0amdq), & 
                    work(i0apdq),work(i0cqxx), & 
                    work(i0bmamdq),work(i0bmapdq), & 
                    work(i0bpamdq),work(i0bpapdq), & 
                    work(i0cmamdq),work(i0cmapdq), & 
                    work(i0cpamdq),work(i0cpapdq), & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), & 

                    work(i0bmcqxxm),work(i0bmcqxxp),  &
                    work(i0bpcqxxm),work(i0bpcqxxp),  &
                    work(i0cmcqxxm),work(i0cmcqxxp),  &
                    work(i0cpcqxxm),work(i0cpcqxxp),  &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim, & 
                    use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)

            !! # update fluxes for use in AMR:

            !! # Note that the roles of the flux updates are changed.
            !! # fadd - modifies the g-fluxes
            !! # gadd - modifies the h-fluxes
            !! # hadd - modifies the f-fluxes

            do  j = 1,my+1
                do  m = 1,meqn
                    gm(i,j,k,m) = gm(i,j,k,m) + faddm(m,j)
                    gp(i,j,k,m) = gp(i,j,k,m) + faddp(m,j)

                    hm(i-1,j,k,m) = hm(i-1,j,k,m) + gadd(j,m,1,-1)
                    hp(i-1,j,k,m) = hp(i-1,j,k,m) + gadd(j,m,1,-1)
                    hm(i  ,j,k,m) = hm(i  ,j,k,m) + gadd(j,m,1, 0)
                    hp(i  ,j,k,m) = hp(i  ,j,k,m) + gadd(j,m,1, 0)
                    hm(i+1,j,k,m) = hm(i+1,j,k,m) + gadd(j,m,1, 1)
                    hp(i+1,j,k,m) = hp(i+1,j,k,m) + gadd(j,m,1, 1)

                    hm(i-1,j,k+1,m) = hm(i-1,j,k+1,m) + gadd(j,m,2,-1)
                    hp(i-1,j,k+1,m) = hp(i-1,j,k+1,m) + gadd(j,m,2,-1)
                    hm(i  ,j,k+1,m) = hm(i  ,j,k+1,m) + gadd(j,m,2, 0)
                    hp(i  ,j,k+1,m) = hp(i  ,j,k+1,m) + gadd(j,m,2, 0)
                    hm(i+1,j,k+1,m) = hm(i+1,j,k+1,m) + gadd(j,m,2, 1)
                    hp(i+1,j,k+1,m) = hp(i+1,j,k+1,m) + gadd(j,m,2, 1)

                    fm(i  ,j,k-1,m) = fm(i  ,j,k-1,m) + hadd(j,m,1,-1)
                    fp(i  ,j,k-1,m) = fp(i  ,j,k-1,m) + hadd(j,m,1,-1)
                    fm(i  ,j,k  ,m) = fm(i  ,j,k  ,m) + hadd(j,m,1, 0)
                    fp(i  ,j,k  ,m) = fp(i  ,j,k  ,m) + hadd(j,m,1, 0)
                    fm(i  ,j,k+1,m) = fm(i  ,j,k+1,m) + hadd(j,m,1, 1)
                    fp(i  ,j,k+1,m) = fp(i  ,j,k+1,m) + hadd(j,m,1, 1)

                    fm(i+1,j,k-1,m) = fm(i+1,j,k-1,m) + hadd(j,m,2,-1)
                    fp(i+1,j,k-1,m) = fp(i+1,j,k-1,m) + hadd(j,m,2,-1)
                    fm(i+1,j,k  ,m) = fm(i+1,j,k  ,m) + hadd(j,m,2, 0)
                    fp(i+1,j,k  ,m) = fp(i+1,j,k  ,m) + hadd(j,m,2, 0)
                    fm(i+1,j,k+1,m) = fm(i+1,j,k+1,m) + hadd(j,m,2, 1)
                    fp(i+1,j,k+1,m) = fp(i+1,j,k+1,m) + hadd(j,m,2, 1)
                end do
            end do
        end do iy_loop
    end do ky_loop


    !! ==================
    !! # perform z sweeps
    !! ==================

    jz_loop : do j = 0, my+1
        iz_loop : do i = 0, mx+1

            do k = 1-mbc, mz+mbc
                do m = 1,meqn
                    !! # copy data along a slice into 1d array:
                    q1d(m,k) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do k = 1-mbc, mz+mbc
                    dtdz1d(k) = dtdz / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                do  k = 1-mbc, mz+mbc
                    do  ma = 1,maux
                        aux1(ma,k,1) = aux(i-1,j-1,k,ma)
                        aux1(ma,k,2) = aux(i-1,j,k,ma)
                        aux1(ma,k,3) = aux(i-1,j+1,k,ma)
                        aux2(ma,k,1) = aux(i,j-1,k,ma)
                        aux2(ma,k,2) = aux(i,j,k,ma)
                        aux2(ma,k,3) = aux(i,j+1,k,ma)
                        aux3(ma,k,1) = aux(i+1,j-1,k,ma)
                        aux3(ma,k,2) = aux(i+1,j,k,ma)
                        aux3(ma,k,3) = aux(i+1,j+1,k,ma)
                    end do
                end do
            endif

            !! # Store the value of i and j along this slice in the common block
            !! # comxyzt in case it is needed in the Riemann solver (for
            !! # variable coefficient problems)

            icom = i
            jcom = j

            !! # compute modifications fadd, gadd and hadd to fluxes along this
            !! # slice:

            call flux3(3,maxm,meqn,maux,mbc,mz,  & 
                    q1d,dtdz1d,dtdx,dtdy,aux1,aux2,aux3, & 
                    faddm,faddp,gadd,hadd,cfl1d, & 
                    work(i0wave),work(i0s),work(i0amdq), & 
                    work(i0apdq),work(i0cqxx), &
                    work(i0bmamdq),work(i0bmapdq), & 
                    work(i0bpamdq),work(i0bpapdq), &
                    work(i0cmamdq),work(i0cmapdq), & 
                    work(i0cpamdq),work(i0cpapdq), & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), &

                    work(i0bmcqxxm),work(i0bmcqxxp), &
                    work(i0bpcqxxm),work(i0bpcqxxp), & 
                    work(i0cmcqxxm),work(i0cmcqxxp), & 
                    work(i0cpcqxxm),work(i0cpcqxxp), &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim, & 
                    use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)


            !! # update fluxes for use in AMR:
            !! # Note that the roles of the flux updates are changed.
            !! # fadd - modifies the h-fluxes
            !! # gadd - modifies the f-fluxes
            !! # hadd - modifies the g-fluxes

             do k = 1,mz+1
                do m = 1,meqn
                    hm(i,j,k,m) = hm(m,i,j,k,m) + faddm(m,k)
                    hp(i,j,k,m) = hp(m,i,j,k,m) + faddp(m,k)

                    fm(i  ,j-1,k,m) = fm(i  ,j-1,k,m) + gadd(k,m,1,-1)
                    fp(i  ,j-1,k,m) = fp(i  ,j-1,k,m) + gadd(k,m,1,-1)
                    fm(i  ,j  ,k,m) = fm(i  ,j  ,k,m) + gadd(k,m,1, 0)
                    fp(i  ,j  ,k,m) = fp(i  ,j  ,k,m) + gadd(k,m,1, 0)
                    fm(i  ,j+1,k,m) = fm(i  ,j+1,k,m) + gadd(k,m,1, 1)
                    fp(i  ,j+1,k,m) = fp(i  ,j+1,k,m) + gadd(k,m,1, 1)

                    fm(i+1,j-1,k,m) = fm(i+1,j-1,k,m) + gadd(k,m,2,-1)
                    fp(i+1,j-1,k,m) = fp(i+1,j-1,k,m) + gadd(k,m,2,-1)
                    fm(i+1,j  ,k,m) = fm(i+1,j  ,k,m) + gadd(k,m,2, 0)
                    fp(i+1,j  ,k,m) = fp(i+1,j  ,k,m) + gadd(k,m,2, 0)
                    fm(i+1,j+1,k,m) = fm(i+1,j+1,k,m) + gadd(k,m,2, 1)
                    fp(i+1,j+1,k,m) = fp(i+1,j+1,k,m) + gadd(k,m,2, 1)

                    gm(i-1,j  ,k,m) = gm(i-1,j  ,k,m) + hadd(k,m,1,-1)
                    gp(i-1,j  ,k,m) = gp(i-1,j  ,k,m) + hadd(k,m,1,-1)
                    gm(i  ,j  ,k,m) = gm(i  ,j  ,k,m) + hadd(k,m,1, 0)
                    gp(i  ,j  ,k,m) = gp(i  ,j  ,k,m) + hadd(k,m,1, 0)
                    gm(i+1,j  ,k,m) = gm(i+1,j  ,k,m) + hadd(k,m,1, 1)
                    gp(i+1,j  ,k,m) = gp(i+1,j  ,k,m) + hadd(k,m,1, 1)

                    gm(i-1,j+1,k,m) = gm(i-1,j+1,k,m) + hadd(k,m,2,-1)
                    gp(i-1,j+1,k,m) = gp(i-1,j+1,k,m) + hadd(k,m,2,-1)
                    gm(i  ,j+1,k,m) = gm(i  ,j+1,k,m) + hadd(k,m,2, 0)
                    gp(i  ,j+1,k,m) = gp(i  ,j+1,k,m) + hadd(k,m,2, 0)
                    gm(i+1,j+1,k,m) = gm(i+1,j+1,k,m) + hadd(k,m,2, 1)
                    gp(i+1,j+1,k,m) = gp(i+1,j+1,k,m) + hadd(k,m,2, 1)                    
                end do 
            end do
        end do iz_loop
    end do jz_loop

    return
end subroutine clawpack46_step3
