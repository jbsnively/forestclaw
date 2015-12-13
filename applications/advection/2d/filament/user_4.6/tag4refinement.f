c     # Template function for setting refinement criteria.  The
c     # user can copy this file to their directory, and then set the
c     # vt.fort_tag4refinement = &tag4refinement.

      subroutine tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax
      double precision dq, dqi, dqj,xc,yc,rc

      tag_patch = 0

c     # Refine based only on first variable in system.
      qmin = q(1,1,1)
      qmax = q(1,1,1)
      dq = 0
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            if (init_flag .ne. 0) then
               xc = xlower + (i-0.5)*dx
               yc = ylower + (j-0.5)*dy
               rc = sqrt((xc-0.5)**2 + (yc-1.0)**2)
               if ((0.25-2*dx) < rc .and. rc < (0.25 + 2*dx)) then
                  tag_patch = 1
                  return
               endif
            else
               qmin = min(qmin,q(i,j,1))
               qmax = max(qmax,q(i,j,1))
               if (qmax-qmin .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            endif


c            do mq = 1,1
c                dqi = dabs(q(i+1,j,mq) - q(i-1,j,mq))
c                dqj = dabs(q(i,j+1,mq) - q(i,j-1,mq))
c                dq  = dmax1(dq, dqi, dqj)
c                if (dq .gt. tag_threshold) then
c                   tag_patch = 1
c                   return
c                endif
c            enddo
         enddo
      enddo


      end
