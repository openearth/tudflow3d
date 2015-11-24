
      subroutine advecu_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      implicit none
c
c********************************************************************
c
c     advecu calculates the advection of the u-velocity, which is
c     the velocity in the radial direction.
c
c     In formula:
c
c         1 d(ruu)     1 d(uv)     d(uw)     vv
c    - (  - ------  +  - -----  +  -----  -  --  )
c         r   dr       r  dphi      dz        r
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Utmp              : contains velocity at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      do k=kb,ke
      kp=k+1
      km=k-1
        do j=jb,je
        jp=j+1
        jm=j-1
          do  i=ib,ie
          ip=i+1
          im=i-1
          rhojp =0.25*(rho(i,j,k)+rho(i,jp,k)+rho(ip,j,k)+rho(ip,jp,k))
          rhojm =0.25*(rho(i,j,k)+rho(i,jm,k)+rho(ip,j,k)+rho(ip,jm,k))
          rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
          rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(ip,j,k)+rho(ip,j,km))
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 (Rp(ip)*(Uvel(i,j,k)+Uvel(ip,j,k))*(Uvel(i,j,k)+Uvel(ip,j,k))*rho(ip,j,k) -
     1  Rp(i )*(Uvel(i,j,k)+Uvel(im,j,k))*(Uvel(i,j,k)+Uvel(im,j,k))*rho(i ,j,k)  )
     1  / ( Ru(i) * ( Rp(ip)-Rp(i) ) )
     +                         + 
     2 (      (Vvel(i,j,k) +Vvel(ip,j,k) )*(Uvel(i,j,k)+Uvel(i,jp,k))*rhojp -
     2        (Vvel(i,jm,k)+Vvel(ip,jm,k))*(Uvel(i,j,k)+Uvel(i,jm,k))*rhojm  )
     2  / ( Ru(i) * dphi )
     +                         +
     3 (      (Wvel(i,j,k) +Wvel(ip,j,k) )*(Uvel(i,j,k)+Uvel(i,j,kp))*rhokp -
     3        (Wvel(i,j,km)+Wvel(ip,j,km))*(Uvel(i,j,k)+Uvel(i,j,km))*rhokm  )
     3  / ( dz )
     +                         -
     4 0.5*(rho(i,j,k)+rho(ip,j,k))*( Vvel(i,j,k) + Vvel(ip,j,k) + Vvel(i,jm,k) + Vvel(ip,jm,k) )**2
     4  / ( 4.0 * Ru(i) )
     +                         )
           enddo
        enddo
      enddo
      return
      end


      subroutine advecv_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      implicit none
c
c********************************************************************
c
c     advecv calculates the advection of the v-velocity, which is
c     the velocity in the tangential direction.
c
c     In formula:
c
c         1 d(ruv)     1 d(vv)     d(vw)     uv
c    - (  - ------  +  - -----  +  -----  +  --  )
c         r   dr       r  dphi      dz        r
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Vtmp              : contains velocity at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm,rhokp,rhokm
      do k=kb,ke
      kp=k+1
      km=k-1
        do j=jb,je
        jp=j+1
        jm=j-1
          do  i=ib,ie
          ip=i+1
          im=i-1
          rhoip =0.25*(rho(i,j,k)+rho(ip,j,k)+rho(i,jp,k)+rho(ip,jp,k))
          rhoim =0.25*(rho(i,j,k)+rho(im,j,k)+rho(i,jp,k)+rho(im,jp,k))
          rhojp =rho(i,jp,k)
          rhojm =rho(i,j,k )
          rhokp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
          rhokm =0.25*(rho(i,j,k)+rho(i,j,km)+rho(i,jp,k)+rho(i,jp,km))
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 ( Ru(i ) * Ru(i ) *rhoip*
     1     (Uvel(i,j,k) +Uvel(i,jp,k) )*(Vvel(i,j,k)+Vvel(ip,j,k)) -
     1   Ru(im) * Ru(im) *rhoim*
     1     (Uvel(im,j,k)+Uvel(im,jp,k))*(Vvel(i,j,k)+Vvel(im,j,k))  )
     1  / ( Rp(i) * Rp(i) * dr(i) )
     +                         + 
     2 (   (Vvel(i,j,k) +Vvel(i,jp,k) )*(Vvel(i,j,k)+Vvel(i,jp,k))*rhojp -
     2     (Vvel(i,jm,k)+Vvel(i,j,k)  )*(Vvel(i,j,k)+Vvel(i,jm,k))*rhojm  )
     2  / ( Rp(i) * dphi )
     +                         +
     3 (   (Wvel(i,j,k) +Wvel(i,jp,k) )*(Vvel(i,j,k)+Vvel(i,j,kp))*rhokp -
     3     (Wvel(i,j,km)+Wvel(i,jp,km))*(Vvel(i,j,k)+Vvel(i,j,km))*rhokm  )
     3  / ( dz )
     +                         )
           enddo
        enddo
      enddo
      return
      end

      subroutine advecw_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px)
      implicit none
c
c********************************************************************
c
c     advecw calculates the advection of the w-velocity, which is
c     the velocity in the axial direction.
c
c     In formula:
c
c         1 d(ruw)     1 d(wv)     d(ww)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Wtmp              : contains velocity at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke,rank,px
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm
      do k=kb,ke
      kp=k+1
      km=k-1
        do j=jb,je
        jp=j+1
        jm=j-1
          do i=ib,ie
          ip=i+1
          im=i-1
          rhoip =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
          rhoim =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
          rhojp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
          rhojm =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))
      putout(i,j,k) = 0.0
      putout(i,j,k) = - 0.25 * (
     1 (Ru(i )*(Uvel(i,j,k) +Uvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(ip,j,k))*rhoip -
     1  Ru(im)*(Uvel(im,j,k)+Uvel(im,j,kp))*(Wvel(i,j,k)+Wvel(im,j,k))*rhoim )
     1  / ( Rp(i) * dr(i) )
     +                         + 
     2 (   (Vvel(i,j,k) +Vvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,jp,k))*rhojp -
     2     (Vvel(i,jm,k)+Vvel(i,jm,kp))*(Wvel(i,j,k)+Wvel(i,jm,k))*rhojm  )
     2  / ( Rp(i) * dphi )
     +                         +
     3 (   (Wvel(i,j,k) +Wvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp) -
     3     (Wvel(i,j,k) +Wvel(i,j,km) )*(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k )  )
     3  / ( dz )
     +                         )
           enddo
         enddo
      enddo
      return
      end


      subroutine advecw_driftfluxCDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke)
      implicit none
c
c********************************************************************
c
c     advecw calculates the advection of the w-velocity, which is
c     the velocity in the axial direction.
c
c     In formula:
c
c         1 d(ruw)     1 d(wv)     d(ww)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          Wtmp              : contains velocity at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
      real rho(0:i1,0:j1,0:k1)
      real rhoip,rhoim,rhojp,rhojm,dzi

	dzi=1./dz
      do k=kb,ke
      kp=k+1
      km=k-1
        do j=jb,je
!        jp=j+1
!        jm=j-1
          do i=ib,ie
!          ip=i+1
!          im=i-1
!          rhoip =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(ip,j,k)+rho(ip,j,kp))
!          rhoim =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(im,j,k)+rho(im,j,kp))
!          rhojp =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jp,k)+rho(i,jp,kp))
!          rhojm =0.25*(rho(i,j,k)+rho(i,j,kp)+rho(i,jm,k)+rho(i,jm,kp))
      putout(i,j,k) = putout(i,j,k) - 0.25 * (
!     1 (Ru(i )*(Uvel(i,j,k) +Uvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(ip,j,k))*rhoip -
!     1  Ru(im)*(Uvel(im,j,k)+Uvel(im,j,kp))*(Wvel(i,j,k)+Wvel(im,j,k))*rhoim )
!     1  / ( Rp(i) * dr(i) )
!     +                         + 
!     2 (   (Vvel(i,j,k) +Vvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,jp,k))*rhojp -
!     2     (Vvel(i,jm,k)+Vvel(i,jm,kp))*(Wvel(i,j,k)+Wvel(i,jm,k))*rhojm  )
!     2  / ( Rp(i) * dphi )
!     +                         +
     3 (   (Wvel(i,j,k) +Wvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp) -
     3     (Wvel(i,j,k) +Wvel(i,j,km) )*(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k )  )
     3  *dzi
     +                         )
           enddo
         enddo
      enddo
      return
      end


      subroutine advecc(putout,putin,Uvel,Vvel,Wvel,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke)
      implicit none
c
c********************************************************************
c
c     advecc calculates the advection for a scalar variable, which is
c     situated in the center point of the grid cell.
c
c     In formula:
c
c         1 d(ruC)     1 d(vC)     d(wC)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1)
c
c     -------------------------------------------start k-loop
      do 100 k=kb,ke
c
      kp=k+1
      km=k-1
c
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1
c
c     -------------------------------------------start i-loop
          do 300 i=ib,ie
c
          ip=i+1
          im=i-1
c
c
      putout(i,j,k) = - 0.5 * (
     1 ( Ru(i)  * Uvel(i,j,k)  * ( putin(i,j,k)  + putin(ip,j,k)  ) -
     1   Ru(im) * Uvel(im,j,k) * ( putin(i,j,k)  + putin(im,j,k)  )   )
     1  / ( Rp(i) * dr(i) )
     +                        + 
     2 (      Vvel(i,j,k)  * ( putin(i,j,k) + putin(i,jp,k) ) -
     2        Vvel(i,jm,k) * ( putin(i,j,k) + putin(i,jm,k) )   )
     2  / ( Rp(i) * dphi )
     +                        +
     3 (      Wvel(i,j,k)  * ( putin(i,j,k)  + putin(i,j,kp)  ) -
     3        Wvel(i,j,km) * ( putin(i,j,k)  + putin(i,j,km)  )   )
     3  / ( dz )
     +                        )
c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop
c
      return
      end


	subroutine advecc_TVD(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
c
c********************************************************************
c
c     advecc calculates the advection for a scalar variable, which is
c     situated in the center point of the grid cell.
c
c     In formula:
c
c         1 d(ruC)     1 d(vC)     d(wC)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real 	dz_i,Rpdr_i,Rpdphi_i,var_grid_Rpos,var_grid_Rneg,var_grid_Lpos,var_grid_Lneg
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)

      dz_i=1./dz


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim

	putin2(0:i1,0:j1,0:k1)=putin
	Vvel2=Vvel  
	if (.false.) then !lateral inflow and outflow is specifically wanted sometimes, therefore lines below switched off... LdW 11-8-2015
	if ((periodicy.eq.0..or.periodicy.eq.2).and.rank.eq.0) then !! make vvel2 zero at lateral boundaries to keep sediment inside with waves
    	  j=0
	  do i=0,i1
	    do k=0,k1
	      Vvel2(i,j,k)=0.
	    enddo
	  enddo
	elseif ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
    	  j=je
	  do i=0,i1
	    do k=0,k1
	      Vvel2(i,j,k)=0.
	    enddo
	  enddo
	endif
	endif



	if (periodicx.eq.0) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
c get stuff from other CPU's
	  call shiftf2(putin,pbf)
	  call shiftb2(putin,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = putin(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =putin(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif


c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      var_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      var_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      var_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      var_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
      Rpdphi_i=1./(Rp2(i)*dphi)
    
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif

c
c     -------------------------------------------start k-loop
	  do 300 k=kb,ke
c
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
	IF (Uvel(i,j,k).ge.0.) THEN
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *var_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
          putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *var_grid_Rneg
	  cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
          putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
	ENDIF
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *var_grid_Lpos
	  cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im)))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *var_grid_Lneg
	  cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im)))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer
	  cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-dt*ABS(Vvel2(i,j ,k))*Rpdphi_i)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer 
	  cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-dt*ABS(Vvel2(i,j ,k))*Rpdphi_i)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRneg*Rpdphi_i 
	ENDIF		
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer
	  cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-dt*ABS(Vvel2(i,jm,k))*Rpdphi_i)
	  putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer 
	  cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-dt*ABS(Vvel2(i,jm,k))*Rpdphi_i)
	  putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel2(i,j ,k)+ABS(Vvel2(i,j ,k)))*cRpos + 0.5*(Vvel2(i,j ,k)-ABS(Vvel2(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel2(i,jm,k)+ABS(Vvel2(i,jm,k)))*cLpos + 0.5*(Vvel2(i,jm,k)-ABS(Vvel2(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRpos*dz_i 
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i 
	ENDIF
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &  ( 0.5*(Wvel(i,j,k )+ABS(Wvel(i,j,k )))*cRpos + 0.5*(Wvel(i,j,k )-ABS(Wvel(i,j,k )))*cRneg )  -
!     &  ( 0.5*(Wvel(i,j,km)+ABS(Wvel(i,j,km)))*cLpos + 0.5*(Wvel(i,j,km)-ABS(Wvel(i,j,km)))*cLneg ) )
!     &  *dz_i !/ ( dz ) 
c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop
c
      return
      end


	subroutine advecc_TVD_old(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
c
c********************************************************************
c
c     advecc calculates the advection for a scalar variable, which is
c     situated in the center point of the grid cell.
c
c     In formula:
c
c         1 d(ruC)     1 d(vC)     d(wC)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real 	dz_i,Rpdr_i,Rpdphi_i,var_grid_Rpos,var_grid_Rneg,var_grid_Lpos,var_grid_Lneg
	real  dt
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)

      dz_i=1./dz


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim

	putin2(0:i1,0:j1,0:k1)=putin

	if (periodicx.eq.0) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
c get stuff from other CPU's
	  call shiftf2(putin,pbf)
	  call shiftb2(putin,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = putin(i,0,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =putin(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-1,k) = pbf(i,k)
		   putin2(i,j1+1,k) =pbb(i,k)
		   enddo
		enddo
	endif


c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      var_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      var_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      var_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      var_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
      Rpdphi_i=1./(Rp2(i)*dphi)
    
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif

c
c     -------------------------------------------start k-loop
	  do 300 k=kb,ke
c
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
c
c	
	noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *var_grid_Rpos
  	rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *var_grid_Rneg
	noemer = putin2(i,j,k)-putin2(im,j,k)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *var_grid_Lneg
	rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *var_grid_Lpos

	cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
	cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im)))
	cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
	cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im)))

      putout(i,j,k) = - (
     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
		
	rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer
	rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer 

	noemer = putin2(i,j,k)-putin2(i,jm,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	
	rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer 
	rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer


	cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-dt*ABS(Vvel(i,j ,k))*Rpdphi_i)
	cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-dt*ABS(Vvel(i,jm,k))*Rpdphi_i)
	cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-dt*ABS(Vvel(i,j ,k))*Rpdphi_i)
	cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-dt*ABS(Vvel(i,jm,k))*Rpdphi_i)


	putout(i,j,k) = putout(i,j,k) - (
     &   (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
     &   (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 

	noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer 
	rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer

	cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)
	cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)

	putout(i,j,k) = putout(i,j,k) - (
     &  ( 0.5*(Wvel(i,j,k )+ABS(Wvel(i,j,k )))*cRpos + 0.5*(Wvel(i,j,k )-ABS(Wvel(i,j,k )))*cRneg )  -
     &  ( 0.5*(Wvel(i,j,km)+ABS(Wvel(i,j,km)))*cLpos + 0.5*(Wvel(i,j,km)-ABS(Wvel(i,j,km)))*cLneg ) )
     &  *dz_i !/ ( dz ) 

c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop
c
      return
      end

	subroutine adveccbot_TVD(putout,putin,Uvel,Vvel,Ru,Rp,dr,dphi,dz,i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
      implicit none
c
c********************************************************************
c
c     advecc calculates the advection for a scalar variable in i,j dir (in bot), which is
c     situated in the center point of the grid cell.
c
c     In formula:
c
c         1 d(ruC)     1 d(vC)     
c    - (  - ------  +  - -----  )
c         r   dr       r  dphi     
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,i1,j1,ib,ie,jb,je
      real     putout(0:i1,0:j1),putin(0:i1,0:j1),putin2(-1:i1+1,-1:j1+1),
     +         Uvel(0:j1),Vvel(0:j1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1),Rp2(-1:i1+1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m
      integer imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real 	dz_i,Rpdr_i,Rpdphi_i,var_grid_Rpos,var_grid_Rneg,var_grid_Lpos,var_grid_Lneg
	real  dt
      real pbf(0:i1)
      real pbb(0:i1)


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	putin2(0:i1,0:j1)=putin

	if (periodicx.eq.0) then
		putin2(-1,0:j1)=putin(0,0:j1)
		putin2(i1+1,0:j1)=putin(ie,0:j1)
	else 
		putin2(-1,0:j1)=putin(ie-1,0:j1)
		putin2(i1+1,0:j1)=putin(2,0:j1)
	endif
c get stuff from other CPU's
	  call shiftf_l2(putin,pbf)
	  call shiftb_l2(putin,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		   do i=1,ie
		   putin2(i,-1) = putin(i,0)
		   putin2(i,j1+1) =pbb(i)
		   enddo
	  elseif (rank.eq.px-1) then
		   do i=1,ie
		   putin2(i,-1) = pbf(i)
		   putin2(i,j1+1) =putin(i,j1)
		   enddo
	  else 
		   do i=1,ie
		   putin2(i,-1) = pbf(i)
		   putin2(i,j1+1) =pbb(i)
		   enddo
	  endif
	else 
		   do i=1,ie
		   putin2(i,-1) = pbf(i)
		   putin2(i,j1+1) =pbb(i)
		   enddo
	endif


      dz_i=1./dz

c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      var_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      var_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      var_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      var_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
      Rpdphi_i=1./(Rp2(i)*dphi)
    
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.0) then
	  if (j.eq.1) jmm=j
	endif
        if ((periodicy.eq.0.or.periodicy.eq.2).and.rank.eq.px-1) then
	  if (j.eq.je) jpp=j
	endif
	noemer = putin2(ip,j)-putin2(i,j)
        noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)


	rRpos = (putin2(i,j)-putin2(im,j))/noemer *var_grid_Rpos
  	rRneg = (putin2(ipp+1,j)-putin2(ip,j))/noemer *var_grid_Rneg
	noemer = putin2(i,j)-putin2(im,j)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(ip,j)-putin2(i,j))/noemer *var_grid_Lneg
	rLpos = (putin2(im,j)-putin2(imm-1,j))/noemer *var_grid_Lpos

	cRpos = putin2(i ,j) + 0.5*limiter(rRpos)*(putin2(ip,j) - putin2(i ,j))*(1.-dt*ABS(Uvel(j))/(Rp(ip)-Rp(i)))
	cLpos = putin2(im,j) + 0.5*limiter(rLpos)*(putin2( i,j) - putin2(im,j))*(1.-dt*ABS(Uvel(j))/(Rp(i)-Rp(im)))
	cRneg = putin2(ip,j) + 0.5*limiter(rRneg)*(putin2( i,j) - putin2(ip,j))*(1.-dt*ABS(Uvel(j))/(Rp(ip)-Rp(i)))
	cLneg = putin2(i ,j) + 0.5*limiter(rLneg)*(putin2(im,j) - putin2(i ,j))*(1.-dt*ABS(Uvel(j))/(Rp(i)-Rp(im)))

      putout(i,j) = - (
     &   Ru(i)  * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cRpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cRneg ) -
     &   Ru(im) * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cLpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cLneg ) )
     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	noemer = putin2(i,jp)-putin2(i,j)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
		
	rRpos = (putin2(i,j)-putin2(i,jm))/noemer
	rRneg = (putin2(i,jpp+1)-putin2(i,jp))/noemer 

	noemer = putin2(i,j)-putin2(i,jm)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	
	rLneg = (putin2(i,jp)-putin2(i,j))/noemer 
	rLpos = (putin2(i,jm)-putin2(i,jmm-1))/noemer


	cRpos = putin2(i ,j) + 0.5*limiter(rRpos)*(putin2(i,jp) - putin2(i ,j))*(1.-dt*ABS(Vvel(j ))*Rpdphi_i)
	cLpos = putin2(i,jm) + 0.5*limiter(rLpos)*(putin2( i,j) - putin2(i,jm))*(1.-dt*ABS(Vvel(jm))*Rpdphi_i)
	cRneg = putin2(i,jp) + 0.5*limiter(rRneg)*(putin2( i,j) - putin2(i,jp))*(1.-dt*ABS(Vvel(j ))*Rpdphi_i)
	cLneg = putin2(i,j ) + 0.5*limiter(rLneg)*(putin2(i,jm) - putin2(i,j ))*(1.-dt*ABS(Vvel(jm))*Rpdphi_i)


	putout(i,j) = putout(i,j) - (
     &   (  0.5*(Vvel(j )+ABS(Vvel(j )))*cRpos + 0.5*(Vvel(j )-ABS(Vvel(j )))*cRneg )  -
     &   (  0.5*(Vvel(jm)+ABS(Vvel(jm)))*cLpos + 0.5*(Vvel(jm)-ABS(Vvel(jm)))*cLneg ) )
     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end i-loop
c
      return
      end

	subroutine advecc_TVD_rho(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,dphi,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt)
      implicit none
c
c********************************************************************
c
c     advecc calculates the advection for a scalar variable, which is
c     situated in the center point of the grid cell.
c
c     In formula:
c
c         1 d(ruC)     1 d(vC)     d(wC)
c    - (  - ------  +  - -----  +  -----  )
c         r   dr       r  dphi      dz
c
c      on input :
c
c          putout            : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,dphi,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putout            : advection part
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),dphi,dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m
      integer kmm,kpp,imm,ipp,jmm,jpp
      real 	dz_i,Rpdr_i,Rpdphi_i,var_grid_Rpos,var_grid_Rneg,var_grid_Lpos,var_grid_Lneg
	real  dt

      dz_i=1./dz

c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (i.eq.1) imm=i
      if (i.eq.ie) ipp=i
      var_grid_Rpos=(Rp(ip)-Rp(i))/(Rp(i)-Rp(im))
      var_grid_Rneg=(Rp(ipp)-Rp(i))/(Rp(ipp+1)-Rp(ipp))
      var_grid_Lpos=(Rp(i)-Rp(imm))/(Rp(imm)-Rp(imm-1))
      var_grid_Lneg=(Rp(i)-Rp(im))/(Rp(ip)-Rp(i))
      Rpdr_i=1./(Rp(i)*dr(i))
      Rpdphi_i=1./(Rp(i)*dphi)
    
c     -------------------------------------------start j-loop
        do 200 j=jb,je
c
        jp=j+1
        jm=j-1
	jmm=jm
	jpp=jp
	if (j.eq.1) jmm=j
	if (j.eq.je) jpp=j

c
c     -------------------------------------------start k-loop
	  do 300 k=kb,ke
c
	  kp=k+1
	  km=k-1
	  kmm=km
	  kpp=kp
	  if (k.eq.1) kmm=k
	  if (k.eq.ke) kpp=k
c
c	
	noemer = putin(ip,j,k)-putin(i,j,k)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif
	


	rRpos = (putin(i,j,k)-putin(im,j,k))/noemer *var_grid_Rpos
!      &        (Rp(ip)-Rp(i))/(Rp(i)-Rp(im))

! 	if (i.ne.ie) then
! 	  rRneg = (putin(ip+1,j,k)-putin(ip,j,k))/noemer *
!      &        (Rp(ip)-Rp(i))/(Rp(ip+1)-Rp(ip))
! 	else
! 	  rRneg = 0.
! 	endif
	  rRneg = (putin(ipp+1,j,k)-putin(ip,j,k))/noemer *var_grid_Rneg
!      &        (Rp(ipp)-Rp(i))/(Rp(ipp+1)-Rp(ipp))

	noemer = putin(i,j,k)-putin(im,j,k)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif

	rLneg = (putin(ip,j,k)-putin(i,j,k))/noemer *var_grid_Lneg
!      &        (Rp(i)-Rp(im))/(Rp(ip)-Rp(i))

! 	if (im.ne.0) then
! 	  rLpos = (putin(im,j,k)-putin(im-1,j,k))/noemer *
!      &        (Rp(i)-Rp(im))/(Rp(im)-Rp(im-1))	
! 	else
! 	  rLpos = 0.
! 	endif
	  rLpos = (putin(im,j,k)-putin(imm-1,j,k))/noemer *var_grid_Lpos
!      &        (Rp(i)-Rp(imm))/(Rp(imm)-Rp(imm-1))

	cRpos = putin(i ,j,k) + 0.5*limiter(rRpos)*(putin(ip,j,k) - putin(i ,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp(ip)-Rp(i)))
	cLpos = putin(im,j,k) + 0.5*limiter(rLpos)*(putin( i,j,k) - putin(im,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp(i)-Rp(im)))
	cRneg = putin(ip,j,k) + 0.5*limiter(rRneg)*(putin( i,j,k) - putin(ip,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp(ip)-Rp(i)))
	cLneg = putin(i ,j,k) + 0.5*limiter(rLneg)*(putin(im,j,k) - putin(i ,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp(i)-Rp(im)))
	

	rho_p = 0.5*(rho(i,j,k)+rho(ip,j,k))
	rho_m = 0.5*(rho(i,j,k)+rho(im,j,k))
      putout(i,j,k) = - (
     &   rho_p*Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
     &   rho_m*Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	noemer = putin(i,jp,k)-putin(i,j,k)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif
		
	rRpos = (putin(i,j,k)-putin(i,jm,k))/noemer

! 	if (j.ne.je) then
! 	  rRneg = (putin(i,jp+1,k)-putin(i,jp,k))/noemer 
! 	else
! 	  rRneg = 0.
! 	endif
	rRneg = (putin(i,jpp+1,k)-putin(i,jp,k))/noemer 

	noemer = putin(i,j,k)-putin(i,jm,k)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif
	
	rLneg = (putin(i,jp,k)-putin(i,j,k))/noemer 

! 	if (jm.ne.0) then
! 	  rLpos = (putin(i,jm,k)-putin(i,jm-1,k))/noemer
! 	else
! 	  rLpos = 0.
! 	endif
	  rLpos = (putin(i,jm,k)-putin(i,jmm-1,k))/noemer


	cRpos = putin(i ,j,k) + 0.5*limiter(rRpos)*(putin(i,jp,k) - putin(i ,j,k))*(1.-dt*ABS(Vvel(i,j ,k))*Rpdphi_i)
	cLpos = putin(i,jm,k) + 0.5*limiter(rLpos)*(putin( i,j,k) - putin(i,jm,k))*(1.-dt*ABS(Vvel(i,jm,k))*Rpdphi_i)
	cRneg = putin(i,jp,k) + 0.5*limiter(rRneg)*(putin( i,j,k) - putin(i,jp,k))*(1.-dt*ABS(Vvel(i,j ,k))*Rpdphi_i)
	cLneg = putin(i,j ,k) + 0.5*limiter(rLneg)*(putin(i,jm,k) - putin(i,j ,k))*(1.-dt*ABS(Vvel(i,jm,k))*Rpdphi_i)	

	rho_p = 0.5*(rho(i,j,k)+rho(i,jp,k))
	rho_m = 0.5*(rho(i,j,k)+rho(i,jm,k))

	putout(i,j,k) = putout(i,j,k) - (
     &  rho_p * (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
     &  rho_m * (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	noemer = putin(i,j,kp)-putin(i,j,k)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif

	rRpos = (putin(i,j,k)-putin(i,j,km))/noemer

! 	if (k.ne.ke) then
! 	  rRneg = (putin(i,j,kp+1)-putin(i,j,kp))/noemer 
! 	else
! 	  rRneg = 0.
! 	endif
	rRneg = (putin(i,j,kpp+1)-putin(i,j,kp))/noemer 

	noemer = putin(i,j,k)-putin(i,j,km)
! 	if (abs(noemer).le.1.e-6) then
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
! 	endif

	rLneg = (putin(i,j,kp)-putin(i,j,k))/noemer 

! 	if (km.ne.0) then
! 	  rLpos = (putin(i,j,km)-putin(i,j,km-1))/noemer
! 	else
! 	  rLpos = 0.
! 	endif
	rLpos = (putin(i,j,km)-putin(i,j,kmm-1))/noemer



	cRpos = putin(i ,j,k) + 0.5*limiter(rRpos)*(putin(i,j,kp) - putin(i ,j,k))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	cLpos = putin(i,j,km) + 0.5*limiter(rLpos)*(putin(i,j, k) - putin(i,j,km))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)
	cRneg = putin(i,j,kp) + 0.5*limiter(rRneg)*(putin(i,j, k) - putin(i,j,kp))*(1.-dt*ABS(Wvel(i,j,k ))*dz_i)
	cLneg = putin(i,j,k ) + 0.5*limiter(rLneg)*(putin(i,j,km) - putin(i,j, k))*(1.-dt*ABS(Wvel(i,j,km))*dz_i)
	
	rho_p = 0.5*(rho(i,j,k)+rho(i,j,kp))
	rho_m = 0.5*(rho(i,j,k)+rho(i,j,km))

	putout(i,j,k) = putout(i,j,k) - (
     &  rho_p * ( 0.5*(Wvel(i,j,k )+ABS(Wvel(i,j,k )))*cRpos + 0.5*(Wvel(i,j,k )-ABS(Wvel(i,j,k )))*cRneg )  -
     &  rho_m * ( 0.5*(Wvel(i,j,km)+ABS(Wvel(i,j,km)))*cLpos + 0.5*(Wvel(i,j,km)-ABS(Wvel(i,j,km)))*cLneg ) )
     &  *dz_i !/ ( dz ) 

c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop
c
      return
      end



      real function limiter(r)
      real r

!	limiter=0.5
!	limiter=MAX(MIN(r,1.),0.) 		! MINMOD
!	limiter=1.				! CDS2
!	limiter=0.				!1st order upwind	
!	limiter=r 				!2nd order upwind
!	limiter=(3+r)/4. 			!Second order QUICK
	limiter=(r+ABS(r))/MAX(1.+r,1.) 	! Van Leer limiter
!	limiter=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	! Superbee limiter
!       limiter=MAX(0.,MIN(5.*r,1.),MIN(r,2.))  ! Superbee-4 limiter non-TVD, less dissipation than Superbee, in matlab test okay


      return
      end

