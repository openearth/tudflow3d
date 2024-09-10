
      subroutine advecu_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phiv,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
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
     2  / ( Ru(i) * (phiv(j)-phiv(jm)) )
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


      subroutine advecv_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phiv,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
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
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
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

      subroutine advecw_CDS2(putout,Uvel,Vvel,Wvel,RHO,Ru,Rp,dr,phiv,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
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
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
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


      subroutine advecw_driftfluxCDS2(putout,Wvel,RHO,Ru,Rp,dr,phiv,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
      real     putout(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
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
     3 (   (Wvel(i,j,k) +Wvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,j,kp)) -
     3     (Wvel(i,j,k) +Wvel(i,j,km) )*(Wvel(i,j,k)+Wvel(i,j,km))  )
     3  *dzi
     +                         )
!     3 (   (Wvel(i,j,k) +Wvel(i,j,kp) )*(Wvel(i,j,k)+Wvel(i,j,kp))*rho(i,j,kp) -
!     3     (Wvel(i,j,k) +Wvel(i,j,km) )*(Wvel(i,j,k)+Wvel(i,j,km))*rho(i,j,k )  )
!     3  *dzi
!     +                         )	 
           enddo
         enddo
      enddo
      return
      end


      subroutine advecc(putout,putin,Uvel,Vvel,Wvel,Ru,Rp,dr,phiv,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
     +         dr(0:i1),phiv(0:j1),dz,Ru(0:i1),Rp(0:i1)
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
     2  / ( Rp(i) * (phiv(j)-phiv(jm)) )
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


	subroutine advecc_TVDbackup(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! TVD scheme cell-wise, so each face is calculated twice
! this TVD routine uses limiter (instead of limiter2) with Van Leer
! can be called by nerd-option advec_conc='VLE'  
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)

      dz_i=1./dz


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

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



	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
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
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
          putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
          putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
	ENDIF
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k))*(1.-dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im)))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
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
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k))*
     & (1.-dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j))))
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k))*
     & (1.-dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j))))
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRneg*Rpdphi_i 
	ENDIF		
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k))
     & *(1.-dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm))))
	  putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k))
     & *(1.-dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm))))
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

	  	subroutine advecc_VLE(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
      implicit none
! TVD scheme called with nerd option advec_conc.eq.'VLE', efficiently face-based with limiter Van Leer incl (1-cfl) term
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
	  integer*8 kbed(0:i1,0:j1)
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	
	if (transporteq_fracs.eq.'massfrac') then
	  uvel2(0:ie,:,:)=uvel(0:ie,:,:)/(0.5*(rho(0:ie,:,:)+rho(1:i1,:,:)))
	  vvel2(:,0:je,:)=vvel(:,0:je,:)/(0.5*(rho(:,0:je,:)+rho(:,1:j1,:))) !vnew2 j1 not needed
	  wvel2(:,:,0:ke)=wvel(:,:,0:ke)/(0.5*(rho(:,:,0:ke)+rho(:,:,1:k1)))
	else
	  uvel2=uvel
	  vvel2=vvel
	  wvel2=wvel
	endif

	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
c
c     -------------------------------------------start k-loop
	  do 300 k=MAX(kb,MIN(kbed(i,j),kbed(i+1,j),kbed(i,j+1),kbed(i+1,j+1))),ke !kb,ke
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
	  cfl =dt*ABS(Uvel2(i ,j,k))/(Rp2(ip)-Rp2(i ))
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl =dt*ABS(Uvel2(i ,j,k))/(Rp2(ip)-Rp2(i ))
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cRneg = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel2(im,j,k))/(Rp2(i )-Rp2(im))
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cLpos = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel2(im,j,k))/(Rp2(i )-Rp2(im))
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLneg = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cRneg = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLpos = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLneg = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,k ))*dz_i
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,k ))*dz_i
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)
	  cRneg = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,km))*dz_i
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
	  cLpos = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,km))*dz_i
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)		  
	  cLneg = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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

      return
      end

	  	subroutine advecc_SBE(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy,transporteq_fracs,kbed)
      implicit none
! TVD scheme called with nerd option advec_conc.eq.'SBE', efficiently face-based with limiter Superbee incl (1-cfl) term
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
	  integer*8 kbed(0:i1,0:j1)
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	
	if (transporteq_fracs.eq.'massfrac') then
	  uvel2(0:ie,:,:)=uvel(0:ie,:,:)/(0.5*(rho(0:ie,:,:)+rho(1:i1,:,:)))
	  vvel2(:,0:je,:)=vvel(:,0:je,:)/(0.5*(rho(:,0:je,:)+rho(:,1:j1,:))) !vnew2 j1 not needed
	  wvel2(:,:,0:ke)=wvel(:,:,0:ke)/(0.5*(rho(:,:,0:ke)+rho(:,:,1:k1)))
	else
	  uvel2=uvel
	  vvel2=vvel
	  wvel2=wvel
	endif

	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
c
c     -------------------------------------------start k-loop
	  do 300 k=MAX(kb,MIN(kbed(i,j),kbed(i+1,j),kbed(i,j+1),kbed(i+1,j+1))),ke !k=kb,ke
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
	  cfl =dt*ABS(Uvel2(i ,j,k))/(Rp2(ip)-Rp2(i ))
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl =dt*ABS(Uvel2(i ,j,k))/(Rp2(ip)-Rp2(i ))
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cRneg = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel2(im,j,k))/(Rp2(i )-Rp2(im))
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cLpos = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel2(im,j,k))/(Rp2(i )-Rp2(im))
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLneg = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cRneg = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLpos = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLneg = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,k ))*dz_i
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,k ))*dz_i
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))
	  cRneg = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,km))*dz_i
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
	  cLpos = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel2(i,j,km))*dz_i
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))		  
	  cLneg = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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

      return
      end

	  

	  	subroutine advecc_VL2_nocfl(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! TVD scheme called with nerd option advec_conc.eq.'VL2', efficiently face-based with limiter Van Leer without (1-cfl) term --> mix upw1 and cds2
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim


	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cRneg = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cLpos = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLneg = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k))
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  cRneg = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k))
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLpos = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k))
	  putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  cLneg = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) 
	  putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)
	  cRneg = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  r=rLpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
	  cLpos = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  r=rLneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)		  
	  cLneg = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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

      return
      end


	  	subroutine advecc_SB2_nocfl(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! TVD scheme called with nerd option advec_conc.eq.'SB2', efficiently face-based with limiter Superbee without (1-cfl) term --> mix upw1 and cds2
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim


	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cRneg = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cLpos = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLneg = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k))
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k))
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  cRneg = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k))
	  putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLpos = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k))
	  putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  cLneg = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) 
	  putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel(i,j ,k)+ABS(Vvel(i,j ,k)))*cRpos + 0.5*(Vvel(i,j ,k)-ABS(Vvel(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel(i,jm,k)+ABS(Vvel(i,jm,k)))*cLpos + 0.5*(Vvel(i,jm,k)-ABS(Vvel(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
	  cRpos = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))
	  cRneg = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  r=rLpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
	  cLpos = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  r=rLneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))		  
	  cLneg = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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

      return
      end
	  
	  
	  	subroutine c_edges_VL2_nocfl(putoutU,putoutV,putoutW,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! determines c at cell edges with TVD van Leer scheme without (1-cfl) term
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
c          putoutU,putoutV,putoutW  : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putoutU            : variable interpolated TVD fashion at cell edge U-loc
c          putoutV            : variable interpolated TVD fashion at cell edge V-loc
c          putoutW            : variable interpolated TVD fashion at cell edge W-loc
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putoutU(0:i1,0:j1,0:k1),putoutV(0:i1,0:j1,0:k1),putoutW(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putoutU=0.
	  putoutV=0.
	  putoutW=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim


	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  putoutU(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))
          !putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		  !putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  putoutU(i,j,k) = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))
		  !putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		  !putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
	ENDIF
!	IF (i.eq.1) THEN  ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Uvel(im,j,k).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(im,j,k)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
!	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
!	  r=rLpos
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
!	  putoutU(im,j,k) = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k))
!          !putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(im,j,k)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
!	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
!	  r=rLneg
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
!	  putoutU(im,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k))
!          !putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!	ENDIF
!	ENDIF

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
	  putoutV(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k))
	  !putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 
	  putoutV(i,j,k) = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k))
	  !putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	ENDIF
!	IF (j.eq.1) THEN ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Vvel(i,jm,k).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(i,jm,k)
!  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
!	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
!	  r=rLpos
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
!	  putoutV(i,jm,k) = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k))
!	  !putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(i,jm,k)
!  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
!	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
!	  r=rLneg
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.) 	  
!	  putoutV(i,jm,k) = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) 
!	  !putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
!	ENDIF	
!	ENDIF
	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
	  putoutW(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))
	  !putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=(r+ABS(r))/MAX(1.+r,1.)
	  putoutW(i,j,k) = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))
	  !putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	ENDIF
!	IF (k.eq.1) THEN  ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Wvel(i,j,km).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(i,j,km)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
!	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
!	  r=rLpos
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.)	  
!	  putoutW(i,j,km) = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))
!	  !putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(i,j,km)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
!	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
!	  r=rLneg
!	  limiter99=(r+ABS(r))/MAX(1.+r,1.)		  
!	  putoutW(i,j,km) = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))
!	  !putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
!	ENDIF
!	ENDIF
c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop

      return
      end


	  	subroutine c_edges_SB2_nocfl(putoutU,putoutV,putoutW,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! determines c at cell edges with TVD Superbee scheme without (1-cfl) term
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
c          putoutU,putoutV,putoutW  : "empty" (initialised to zero)
c          putin             : variable for which the advection has
c                              to be calculated
c          Uvel,Vvel,Wvel    : contain velocities at former timestep
c          putinn            : contains subgrid energy at oldest timestep
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
c          i1,j1,k1          : parameters for array-dimensions
c          ib,ie,jb,je,kb,ke : range of gridpoints for which the
c                              advection has to be calculated
c          Ru,Rp             : radial positions of the U-velocity
c                              component and the pressure location
c                              respectively
c
c      on output :
c
c          putoutU            : variable interpolated TVD fashion at cell edge U-loc
c          putoutV            : variable interpolated TVD fashion at cell edge V-loc
c          putoutW            : variable interpolated TVD fashion at cell edge W-loc
c          other parameters  : all unchanged
c
c********************************************************************
      integer  i,j,k,im,ip,jm,jp,km,kp,i1,j1,k1,ib,ie,jb,je,kb,ke
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putoutU(0:i1,0:j1,0:k1),putoutV(0:i1,0:j1,0:k1),putoutW(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
     +         ,Uvel2(0:i1,0:j1,0:k1),Vvel2(0:i1,0:j1,0:k1),Wvel2(0:i1,0:j1,0:k1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt!,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip,limiter99,r
	character*8 transporteq_fracs

      dz_i=1./dz

	  putoutU=0.
	  putoutV=0.
	  putoutW=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim


	putin2(0:i1,0:j1,0:k1)=putin
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  putoutU(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(ip,j,k) - putin2(i ,j,k))
          !putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		  !putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  !cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  putoutU(i,j,k) = putin2(ip,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(ip,j,k))
		  !putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		  !putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
	ENDIF
!	IF (i.eq.1) THEN  ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Uvel(im,j,k).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(im,j,k)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
!	  !cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
!	  r=rLpos
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
!	  putoutU(im,j,k) = putin2(im,j,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(im,j,k))
!          !putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(im,j,k)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
!	  !cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
!	  r=rLneg
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
!	  putoutU(im,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(im,j,k) - putin2(i ,j,k))
!          !putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!	ENDIF
!	ENDIF

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
	  putoutV(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,jp,k) - putin2(i ,j,k))
	  !putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRpos*Rpdphi_i 
	  !putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRpos*Rpdphi_ip 
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  !cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 
	  putoutV(i,j,k) = putin2(i,jp,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jp,k))
	  !putout(i,j,k) = putout(i,j,k) - Vvel(i,j ,k)*cRneg*Rpdphi_i 
	  !putout(i,jp,k) = putout(i,jp,k) + Vvel(i,j ,k)*cRneg*Rpdphi_ip 
	ENDIF
!	IF (j.eq.1) THEN ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Vvel(i,jm,k).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(i,jm,k)
!  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
!	  !cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
!	  r=rLpos
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
!	  putoutV(i,jm,k) = putin2(i,jm,k) + 0.5*limiter99*(putin2( i,j,k) - putin2(i,jm,k))
!	  !putout(i,j,k)  = putout(i,j ,k) + Vvel(i,jm,k)*cLpos*Rpdphi_i 
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(i,jm,k)
!  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
!	  !cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
!	  r=rLneg
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) 	  
!	  putoutV(i,jm,k) = putin2(i,j ,k) + 0.5*limiter99*(putin2(i,jm,k) - putin2(i,j ,k)) 
!	  !putout(i,j,k) = putout(i,j,k) + Vvel(i,jm,k)*cLneg*Rpdphi_i 
!	ENDIF	
!	ENDIF
	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  !cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  r=rRpos
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
	  putoutW(i,j,k) = putin2(i ,j,k) + 0.5*limiter99*(putin2(i,j,kp) - putin2(i ,j,k))
	  !putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  !cRneg = putin2(i,j,kp) + 0.5*limiter(rRneg)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  r=rRneg
	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))
	  putoutW(i,j,k) = putin2(i,j,kp) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,kp))
	  !putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	ENDIF
!	IF (k.eq.1) THEN  ! skip all boundary rhoU,rhoV,rhoW determination --> simply use (rho_left+rho_right)/2 after this subroutine for all 0,i1,0,j1,0,k1 cells
!	IF (Wvel(i,j,km).ge.0.) THEN
!	  noemer = putin2(i,j,k)-putin2(i,j,km)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
!	  !cLpos = putin2(i,j,km) + 0.5*limiter(rLpos)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
!	  r=rLpos
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))	  
!	  putoutW(i,j,km) = putin2(i,j,km) + 0.5*limiter99*(putin2(i,j, k) - putin2(i,j,km))
!	  !putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
!	ELSE
!	  noemer = putin2(i,j,k)-putin2(i,j,km)
!	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
!	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
!	  !cLneg = putin2(i,j,k ) + 0.5*limiter(rLneg)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
!	  r=rLneg
!	  limiter99=MAX(0.,MIN(2.*r,1.),MIN(r,2.))		  
!	  putoutW(i,j,km) = putin2(i,j,k ) + 0.5*limiter99*(putin2(i,j,km) - putin2(i,j, k))
!	  !putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
!	ENDIF
!	ENDIF
c
300       continue
c     -------------------------------------------end i-loop
200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end k-loop

      return
      end
	  
	  
	  
	  	  
	subroutine advecc_ARO(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! TVD scheme called with nerd option advec_conc.eq.'ARO', efficiently face-based with limiter2 (3th order CFL dependent alpha limiter), LdW 25-11-15
! in test sim this scheme gave some values <0 and default advec_conc.eq.'VLE' did not... this shouldn't have happened but it did...
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter2,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

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



	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  cfl =dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i ))
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter2(rRpos,cfl)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl =dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i ))
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  cRneg = putin2(ip,j,k) + 0.5*limiter2(rRneg,cfl)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im))
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  cLpos = putin2(im,j,k) + 0.5*limiter2(rLpos,cfl)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im))
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  cLneg = putin2(i ,j,k) + 0.5*limiter2(rLneg,cfl)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u<0',i,Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!		endif		  
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter2(rRpos,cfl)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  cRneg = putin2(i,jp,k) + 0.5*limiter2(rRneg,cfl)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  cLpos = putin2(i,jm,k) + 0.5*limiter2(rLpos,cfl)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  putout(i,j,k)  = putout(i,j ,k) + Vvel2(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  cLneg = putin2(i,j ,k) + 0.5*limiter2(rLneg,cfl)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel2(i,j ,k)+ABS(Vvel2(i,j ,k)))*cRpos + 0.5*(Vvel2(i,j ,k)-ABS(Vvel2(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel2(i,jm,k)+ABS(Vvel2(i,jm,k)))*cLpos + 0.5*(Vvel2(i,jm,k)-ABS(Vvel2(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,k ))*dz_i
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  cRpos = putin2(i ,j,k) + 0.5*limiter2(rRpos,cfl)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,k ))*dz_i
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  cRneg = putin2(i,j,kp) + 0.5*limiter2(rRneg,cfl)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,km))*dz_i
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  cLpos = putin2(i,j,km) + 0.5*limiter2(rLpos,cfl)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,km))*dz_i
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  cLneg = putin2(i,j,k ) + 0.5*limiter2(rLneg,cfl)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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

	  
	subroutine advecc_VOF(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
! TVD scheme called with frac(.)%type=-1 for VOF fraction, efficiently face-based with limiter4 (variant of hyperBEE) 
! keeps steep gradients steep, some over/undershoot is fixed by simple clipping to force values between 0-1 after time update 
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter4,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real cfl,flux,Rpdr_ip,Rpdphi_ip

      dz_i=1./dz

	  putout(ib:ie,jb:je,kb:ke)=0.
	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

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



	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  cfl =dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i ))
	  rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter4(rRpos,cfl)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cfl)
	    !  if (i.eq.ib) then
		!  putout(i,j,k)  = - flux
		!  else
          putout(i,j,k)  = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_i
		!  endif
		!  if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRpos*Rpdr_ip
		!  endif
		!if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
		!write(*,*),'i,fluxR u>0',i,flux
		!endif		
	ELSE
	  noemer = putin2(ip,j,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl =dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i ))
  	  rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	  cRneg = putin2(ip,j,k) + 0.5*limiter4(rRneg,cfl)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cfl)
	     ! if (i.eq.ib) then
         ! putout(i,j,k)  = - flux
		 ! else
		  putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_i
		 ! endif
		 ! if (ip.le.ie) then
		  putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cRneg*Rpdr_ip
		 ! endif
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-1-12) then
!		write(*,*),'i,fluxR u<0',i,flux
!		endif		 
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im))
	  rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	  cLpos = putin2(im,j,k) + 0.5*limiter4(rLpos,cfl)*(putin2( i,j,k) - putin2(im,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u>0',i,Ru(im)*Uvel(im,j,k)*cLpos*Rpdr_i
!		endif		  
	ELSE
	  noemer = putin2(i,j,k)-putin2(im,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Uvel(im,j,k))/(Rp2(i )-Rp2(im))
	  rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	  cLneg = putin2(i ,j,k) + 0.5*limiter4(rLneg,cfl)*(putin2(im,j,k) - putin2(i ,j,k)) *(1.-cfl)
          putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!		if (rank.eq.4.and.j.eq.1.and.k.eq.26.and.i.eq.ie-12) then
!		write(*,*),'i,fluxL u<0',i,Ru(im)*Uvel(im,j,k)*cLneg*Rpdr_i
!		endif		  
	ENDIF
	ENDIF
!      putout(i,j,k) = - (
!     &   Ru(i)  * ( 0.5*(Uvel(i ,j,k)+ABS(Uvel(i ,j,k)))*cRpos + 0.5*(Uvel(i ,j,k)-ABS(Uvel(i ,j,k)))*cRneg ) -
!     &   Ru(im) * ( 0.5*(Uvel(im,j,k)+ABS(Uvel(im,j,k)))*cLpos + 0.5*(Uvel(im,j,k)-ABS(Uvel(im,j,k)))*cLneg ) )
!     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	IF (Vvel(i,j ,k).ge.0.) THEN
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	  cRpos = putin2(i ,j,k) + 0.5*limiter4(rRpos,cfl)*(putin2(i,jp,k) - putin2(i ,j,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRpos*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cRpos*Rpdphi_ip 
	  !endif
	ELSE
	  noemer = putin2(i,jp,k)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	  rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 
	  cRneg = putin2(i,jp,k) + 0.5*limiter4(rRneg,cfl)*(putin2( i,j,k) - putin2(i,jp,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cRneg*Rpdphi_i 
	  !if (jp.le.je) then
	  putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cRneg*Rpdphi_ip 
	  !endif
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	  cLpos = putin2(i,jm,k) + 0.5*limiter4(rLpos,cfl)*(putin2( i,j,k) - putin2(i,jm,k)) *(1.-cfl)
	  putout(i,j,k)  = putout(i,j ,k) + Vvel2(i,jm,k)*cLpos*Rpdphi_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,jm,k)
  	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	  rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	  cLneg = putin2(i,j ,k) + 0.5*limiter4(rLneg,cfl)*(putin2(i,jm,k) - putin2(i,j ,k)) *(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cLneg*Rpdphi_i 
	ENDIF	
	ENDIF
!	putout(i,j,k) = putout(i,j,k) - (
!     &   (  0.5*(Vvel2(i,j ,k)+ABS(Vvel2(i,j ,k)))*cRpos + 0.5*(Vvel2(i,j ,k)-ABS(Vvel2(i,j ,k)))*cRneg )  -
!     &   (  0.5*(Vvel2(i,jm,k)+ABS(Vvel2(i,jm,k)))*cLpos + 0.5*(Vvel2(i,jm,k)-ABS(Vvel2(i,jm,k)))*cLneg ) )
!     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 

	IF (Wvel(i,j,k ).ge.0.) THEN
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,k ))*dz_i
	  rRpos = (putin2(i,j,k    )-putin2(i,j,km))/noemer
	  cRpos = putin2(i ,j,k) + 0.5*limiter4(rRpos,cfl)*(putin2(i,j,kp) - putin2(i ,j,k))*(1.-cfl)
	  putout(i,j,k)  = putout(i,j,k ) - Wvel(i,j,k )*cRpos*dz_i 
	  !if (kp.le.ke) then
	  putout(i,j,kp) = putout(i,j,kp) + Wvel(i,j,k )*cRpos*dz_i 
	  !endif
	ELSE
	  noemer = putin2(i,j,kp)-putin2(i,j,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,k ))*dz_i
	  rRneg = (putin2(i,j,kpp+1)-putin2(i,j,kp))/noemer 
 	  cRneg = putin2(i,j,kp) + 0.5*limiter4(rRneg,cfl)*(putin2(i,j, k) - putin2(i,j,kp))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cRneg*dz_i
	  !if (kp.le.ke) then
	  putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cRneg*dz_i
	  !endif	  
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,km))*dz_i
	  rLpos = (putin2(i,j,km)-putin2(i,j,kmm-1))/noemer
	  cLpos = putin2(i,j,km) + 0.5*limiter4(rLpos,cfl)*(putin2(i,j, k) - putin2(i,j,km))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLpos*dz_i 
	ELSE
	  noemer = putin2(i,j,k)-putin2(i,j,km)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	  cfl=dt*ABS(Wvel(i,j,km))*dz_i
	  rLneg = (putin2(i,j,kp)-putin2(i,j,k))/noemer
	  cLneg = putin2(i,j,k ) + 0.5*limiter4(rLneg,cfl)*(putin2(i,j,km) - putin2(i,j, k))*(1.-cfl)
	  putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km)*cLneg*dz_i 
	ENDIF
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
	  

	  

	subroutine advecc_NVD(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
!	Face-wise advection scheme according to NVD is similar to TVD, but manner of expression is different
!	See for example:
!	# Hirsch Chap 8.4.
!	# Darwish 1994, Normalized Variable and Space Formulation Methodology for High-Resolution Schemes
!	# Leonard 1991, The ULTIMATE conservative difference scheme applied to unsteady one-dimensional advection
!	NVD combined with UPW5 high order scheme to minimise numerical diffusion and avoid under/overshoot -> that works, however in diagonal test waviness of the boundary occurs
!   nerd-option advec_conc='NVD' calls NVD Fromm which gave best results in diagonal edge and diagonal square tests (<diff than vLeer and Arora and same small waviness edge as vLeer)
!   this part is kept for future testing perhaps 
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-2:i1+2,-2:j1+2,-2:k1+2)
	real xD,xC,xU,xF,cD,cC,cU,cF,phic,phif,eps,xxf,xxc,cfl,alpha,del,Rpdr_ip,Rpdphi_ip,eps2,numdif
	integer kkp,kkm,kkpp,kkppp,kkmm,kkmmm,jjp,jjm,jjpp,jjppp,jjmm,jjmmm,iip,iim,iipp,iippp,iimm,iimmm
	real Ax,Bx,Cx,DDx
	real Ay,By,Cy,DDy
	real Az,Bz,Cz,DDz	

      dz_i=1./dz
	  eps=1.e-6
	  eps2=1.e-2
	  
	  putout(ib:ie,jb:je,kb:ke)=0.

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

	numdif=0. !1./60. ! upwind5
		  Az=37.
		  Bz=-8.
		  Cz=1.
		  DDz=1./60.
		  Ay=37.
		  By=-8.
		  Cy=1.
		  DDy=1./60.
		  Ax=37.
		  Bx=-8.
		  Cx=1.
		  DDx=1./60.	
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


!c get stuff from other CPU's
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
	  call shiftf3(putin,pbf)
	  call shiftb3(putin,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = putin(i,0,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =putin(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	endif
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-2,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
		putin2(i1+2,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-2,0:j1,0:k1)=putin(ie-2,0:j1,0:k1)
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
		putin2(i1+2,0:j1,0:k1)=putin(3,0:j1,0:k1)
	endif

	putin2(0:i1,0:j1,-2)=putin(0:i1,0:j1,0)
	putin2(0:i1,0:j1,-1)=putin(0:i1,0:j1,0)
	putin2(0:i1,0:j1,k1+1)=putin(0:i1,0:j1,ke)
	putin2(0:i1,0:j1,k1+2)=putin(0:i1,0:j1,ke)



c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))
      iip=i+1
      iim=i-1
	  iipp=i+2
	  iippp=i+3
	  iimm=i-2
	  iimmm=i-3	  
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
		Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
		Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
		jjp=j+1
        jjm=j-1
		jjpp=j+2
		jjppp=j+3
		jjmm=j-2
	    jjmmm=j-3
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
      kkp=k+1
      kkm=k-1
	  kkpp=k+2
	  kkppp=k+3
	  kkmm=k-2
	  kkmmm=k-3	  
	IF (Uvel(i,j,k).ge.0.) THEN
	  cD=putin2(ip,j,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(im,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i !UPW 
		putout(ip,j,k)= putout(ip,j,k)+ Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip !UPW 
	  ELSE
		cfl = (dt*(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
		cF = 
     1 DDx*(Ax*(putin2(i,j,k)+putin2(iip,j,k))+Bx*(putin2(iim,j,k)+putin2(iipp,j,k))+Cx*(putin2(iimm,j,k)+putin2(iippp,j,k))) -		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(iip,j,k))-5.*(-putin2(iim,j,k)+putin2(iipp,j,k))+(-putin2(iimm,j,k)+putin2(iippp,j,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
        putout(i,j,k)  = putout(i,j,k)  - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF
	ELSE
	  cD=putin2(i    ,j,k)
	  cC=putin2(ip   ,j,k)
	  cU=putin2(ipp+1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i  !UPW 
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip	 !UPW 	
	  ELSE
		cfl = (dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
		cF = 
     1 DDx*(Ax*(putin2(i,j,k)+putin2(iip,j,k))+Bx*(putin2(iim,j,k)+putin2(iipp,j,k))+Cx*(putin2(iimm,j,k)+putin2(iippp,j,k))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(iip,j,k))-5.*(-putin2(iim,j,k)+putin2(iipp,j,k))+(-putin2(iimm,j,k)+putin2(iippp,j,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF	
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  cD=putin2(i,j,k)
	  cC=putin2(im ,j,k)
	  cU=putin2(imm-1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW
	  ELSE
		cfl = (dt*(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
		cF = 
     1 DDx*(Ax*(putin2(iim,j,k)+putin2(i,j,k))+Bx*(putin2(iimm,j,k)+putin2(iip,j,k))+Cx*(putin2(iimmm,j,k)+putin2(iipp,j,k))) -		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(iim,j,k)+putin2(i,j,k))-5.*(-putin2(iimm,j,k)+putin2(iip,j,k))+(-putin2(iimmm,j,k)+putin2(iipp,j,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF
	ELSE
	  cD=putin2(im,j,k)
	  cC=putin2(i   ,j,k)
	  cU=putin2(ip,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW 
	  ELSE
		cfl = (dt*ABS(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
		cF = 
     1 DDx*(Ax*(putin2(iim,j,k)+putin2(i,j,k))+Bx*(putin2(iimm,j,k)+putin2(iip,j,k))+Cx*(putin2(iimmm,j,k)+putin2(iipp,j,k))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(iim,j,k)+putin2(i,j,k))-5.*(-putin2(iimm,j,k)+putin2(iip,j,k))+(-putin2(iimmm,j,k)+putin2(iipp,j,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF		
	ENDIF
	ENDIF

	IF (Vvel(i,j ,k).ge.0.) THEN
	  cD=putin2(i,jp,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(i,jm,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
		cfl = dt*(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
		cF = 
     1 DDy*(Ay*(putin2(i,j,k)+putin2(i,jjp,k))+By*(putin2(i,jjm,k)+putin2(i,jjpp,k))+Cy*(putin2(i,jjmm,k)+putin2(i,jjppp,k))) -		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,jjp,k))-5.*(-putin2(i,jjm,k)+putin2(i,jjpp,k))+(-putin2(i,jjmm,k)+putin2(i,jjppp,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cF*Rpdphi_ip 
	  ENDIF
	ELSE
	  cD=putin2(i, j,k)
	  cC=putin2(i ,jp,k)
	  cU=putin2(i,jpp+1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
		cfl = dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
		cF = 
     1 DDy*(Ay*(putin2(i,j,k)+putin2(i,jjp,k))+By*(putin2(i,jjm,k)+putin2(i,jjpp,k))+Cy*(putin2(i,jjmm,k)+putin2(i,jjppp,k))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,jjp,k))-5.*(-putin2(i,jjm,k)+putin2(i,jjpp,k))+(-putin2(i,jjmm,k)+putin2(i,jjppp,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic
		
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
	    putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cF*Rpdphi_ip 	  
	  ENDIF
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  cD=putin2(i ,j,k)
	  cC=putin2(i ,jm,k)
	  cU=putin2(i ,jmm-1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
		cfl = dt*(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
		cF = 
     1 DDy*(Ay*(putin2(i,jjm,k)+putin2(i,j,k))+By*(putin2(i,jjmm,k)+putin2(i,jjp,k))+Cy*(putin2(i,jjmmm,k)+putin2(i,jjpp,k))) -		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(i,jjm,k)+putin2(i,j,k))-5.*(-putin2(i,jjmm,k)+putin2(i,jjp,k))+(-putin2(i,jjmmm,k)+putin2(i,jjpp,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ELSE
	  cD=putin2(i, jm,k)
	  cC=putin2(i ,j ,k)
	  cU=putin2(i ,jp,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
		cfl = dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
		cF = 
     1 DDy*(Ay*(putin2(i,jjm,k)+putin2(i,j,k))+By*(putin2(i,jjmm,k)+putin2(i,jjp,k))+Cy*(putin2(i,jjmmm,k)+putin2(i,jjpp,k))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(i,jjm,k)+putin2(i,j,k))-5.*(-putin2(i,jjmm,k)+putin2(i,jjp,k))+(-putin2(i,jjmmm,k)+putin2(i,jjpp,k)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ENDIF	
	ENDIF

	IF (Wvel(i,j,k ).ge.0.) THEN
	  cD=putin2(i,j ,kp)
	  cC=putin2(i ,j,k )
	  cU=putin2(i,j ,km)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,k ))*dz_i
		cF = 
     1 DDz*(Az*(putin2(i,j,k)+putin2(i,j,kkp))+Bz*(putin2(i,j,kkm)+putin2(i,j,kkpp))+Cz*(putin2(i,j,kkmm)+putin2(i,j,kkppp)))-		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,j,kkp))-5.*(-putin2(i,j,kkm)+putin2(i,j,kkpp))+(-putin2(i,j,kkmm)+putin2(i,j,kkppp)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i 
	  ENDIF
	  
	ELSE
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,kp)
	  cU=putin2(i,j ,kpp+1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,k ))*dz_i
		cF = 
     1 DDz*(Az*(putin2(i,j,k)+putin2(i,j,kkp))+Bz*(putin2(i,j,kkm)+putin2(i,j,kkpp))+Cz*(putin2(i,j,kkmm)+putin2(i,j,kkppp))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,j,kkp))-5.*(-putin2(i,j,kkm)+putin2(i,j,kkpp))+(-putin2(i,j,kkmm)+putin2(i,j,kkppp)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i
	  ENDIF
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,km)
	  cU=putin2(i,j ,kmm-1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,km ))*dz_i
		cF = 
     1 DDz*(Az*(putin2(i,j,kkm)+putin2(i,j,k))+Bz*(putin2(i,j,kkmm)+putin2(i,j,kkp))+Cz*(putin2(i,j,kkmmm)+putin2(i,j,kkpp)))-		! pos velocity minus numdif
     1 numdif*(10.*(-putin2(i,j,kkm)+putin2(i,j,k))-5.*(-putin2(i,j,kkmm)+putin2(i,j,kkp))+(-putin2(i,j,kkmmm)+putin2(i,j,kkpp)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ELSE
	  cD=putin2(i,j ,km )
	  cC=putin2(i ,j,k)
	  cU=putin2(i,j ,kp)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,km ))*dz_i
		cF = 
     1 DDz*(Az*(putin2(i,j,kkm)+putin2(i,j,k))+Bz*(putin2(i,j,kkmm)+putin2(i,j,kkp))+Cz*(putin2(i,j,kkmmm)+putin2(i,j,kkpp))) +		! neg velocity plus numdif
     1 numdif*(10.*(-putin2(i,j,kkm)+putin2(i,j,k))-5.*(-putin2(i,j,kkmm)+putin2(i,j,kkp))+(-putin2(i,j,kkmmm)+putin2(i,j,kkpp)))
		phif = (cF-cU)/del
		phif = (1.-cfl)*phif+cfl*phic

		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ENDIF
	ENDIF

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

	subroutine advecc_NVD_VOF(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
!	Face-wise advection scheme according to NVD is similar to TVD, but manner of expression is different
!	See for example:
!	# Hirsch Chap 8.4.
!	# Darwish 1994, Normalized Variable and Space Formulation Methodology for High-Resolution Schemes
!	# Leonard 1991, The ULTIMATE conservative difference scheme applied to unsteady one-dimensional advection
!	NVD-HRIC scheme; for VOF HRIC NVD scheme is used often (for example Fluent and StarCCM+)
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1),putin2(-2:i1+2,-2:j1+2,-2:k1+2)
	real xD,xC,xU,xF,cD,cC,cU,cF,phic,phif,eps,xxf,xxc,cfl,alpha,del,Rpdr_ip,Rpdphi_ip,eps2,numdif
	integer kkp,kkm,kkpp,kkppp,kkmm,kkmmm,jjp,jjm,jjpp,jjppp,jjmm,jjmmm,iip,iim,iipp,iippp,iimm,iimmm
	real Ax,Bx,Cx,DDx
	real Ay,By,Cy,DDy
	real Az,Bz,Cz,DDz	

      dz_i=1./dz
	  eps=1.e-6
	  eps2=1.e-2
	  
	  putout(ib:ie,jb:je,kb:ke)=0.

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

	numdif=0. !1./60. ! upwind5
		  Az=37.
		  Bz=-8.
		  Cz=1.
		  DDz=1./60.
		  Ay=37.
		  By=-8.
		  Cy=1.
		  DDy=1./60.
		  Ax=37.
		  Bx=-8.
		  Cx=1.
		  DDx=1./60.	
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


!c get stuff from other CPU's
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
	  call shiftf3(putin,pbf)
	  call shiftb3(putin,pbb) 

	if (periodicy.eq.0.or.periodicy.eq.2) then
	  if (rank.eq.0) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = putin(i,0,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	  elseif (rank.eq.px-1) then
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =putin(i,j1,k)
		   enddo
		enddo
	  else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	  endif
	else 
		do k=1,ke
		   do i=1,ie
		   putin2(i,-2,k) = pbf(i,k)
		   putin2(i,j1+2,k) =pbb(i,k)
		   enddo
		enddo
	endif
	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-2,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
		putin2(i1+2,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-2,0:j1,0:k1)=putin(ie-2,0:j1,0:k1)
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
		putin2(i1+2,0:j1,0:k1)=putin(3,0:j1,0:k1)
	endif

	putin2(0:i1,0:j1,-2)=putin(0:i1,0:j1,0)
	putin2(0:i1,0:j1,-1)=putin(0:i1,0:j1,0)
	putin2(0:i1,0:j1,k1+1)=putin(0:i1,0:j1,ke)
	putin2(0:i1,0:j1,k1+2)=putin(0:i1,0:j1,ke)



c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))
      iip=i+1
      iim=i-1
	  iipp=i+2
	  iippp=i+3
	  iimm=i-2
	  iimmm=i-3	  
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
		Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
		Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
		jjp=j+1
        jjm=j-1
		jjpp=j+2
		jjppp=j+3
		jjmm=j-2
	    jjmmm=j-3
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
      kkp=k+1
      kkm=k-1
	  kkpp=k+2
	  kkppp=k+3
	  kkmm=k-2
	  kkmmm=k-3	  
	IF (Uvel(i,j,k).ge.0.) THEN
	  cD=putin2(ip,j,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(im,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i !UPW 
		putout(ip,j,k)= putout(ip,j,k)+ Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip !UPW 
	  ELSE
		cfl = (dt*(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
!		cF = 
!     1 DDx*(Ax*(putin2(i,j,k)+putin2(iip,j,k))+Bx*(putin2(iim,j,k)+putin2(iipp,j,k))+Cx*(putin2(iimm,j,k)+putin2(iippp,j,k))) -		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(iip,j,k))-5.*(-putin2(iim,j,k)+putin2(iipp,j,k))+(-putin2(iimm,j,k)+putin2(iippp,j,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
        putout(i,j,k)  = putout(i,j,k)  - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF
	ELSE
	  cD=putin2(i    ,j,k)
	  cC=putin2(ip   ,j,k)
	  cU=putin2(ipp+1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i  !UPW 
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip	 !UPW 	
	  ELSE
		cfl = (dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
!		cF = 
!     1 DDx*(Ax*(putin2(i,j,k)+putin2(iip,j,k))+Bx*(putin2(iim,j,k)+putin2(iipp,j,k))+Cx*(putin2(iimm,j,k)+putin2(iippp,j,k))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(iip,j,k))-5.*(-putin2(iim,j,k)+putin2(iipp,j,k))+(-putin2(iimm,j,k)+putin2(iippp,j,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF	
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  cD=putin2(i,j,k)
	  cC=putin2(im ,j,k)
	  cU=putin2(imm-1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW
	  ELSE
		cfl = (dt*(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
!		cF = 
!     1 DDx*(Ax*(putin2(iim,j,k)+putin2(i,j,k))+Bx*(putin2(iimm,j,k)+putin2(iip,j,k))+Cx*(putin2(iimmm,j,k)+putin2(iipp,j,k))) -		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(iim,j,k)+putin2(i,j,k))-5.*(-putin2(iimm,j,k)+putin2(iip,j,k))+(-putin2(iimmm,j,k)+putin2(iipp,j,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF
	ELSE
	  cD=putin2(im,j,k)
	  cC=putin2(i   ,j,k)
	  cU=putin2(ip,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW 
	  ELSE
		cfl = (dt*ABS(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
!		cF = 
!     1 DDx*(Ax*(putin2(iim,j,k)+putin2(i,j,k))+Bx*(putin2(iimm,j,k)+putin2(iip,j,k))+Cx*(putin2(iimmm,j,k)+putin2(iipp,j,k))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(iim,j,k)+putin2(i,j,k))-5.*(-putin2(iimm,j,k)+putin2(iip,j,k))+(-putin2(iimmm,j,k)+putin2(iipp,j,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF		
	ENDIF
	ENDIF

	IF (Vvel(i,j ,k).ge.0.) THEN
	  cD=putin2(i,jp,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(i,jm,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
		cfl = dt*(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
!		cF = 
!     1 DDy*(Ay*(putin2(i,j,k)+putin2(i,jjp,k))+By*(putin2(i,jjm,k)+putin2(i,jjpp,k))+Cy*(putin2(i,jjmm,k)+putin2(i,jjppp,k))) -		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,jjp,k))-5.*(-putin2(i,jjm,k)+putin2(i,jjpp,k))+(-putin2(i,jjmm,k)+putin2(i,jjppp,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cF*Rpdphi_ip 
	  ENDIF
	ELSE
	  cD=putin2(i, j,k)
	  cC=putin2(i ,jp,k)
	  cU=putin2(i,jpp+1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
		cfl = dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
!		cF = 
!     1 DDy*(Ay*(putin2(i,j,k)+putin2(i,jjp,k))+By*(putin2(i,jjm,k)+putin2(i,jjpp,k))+Cy*(putin2(i,jjmm,k)+putin2(i,jjppp,k))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,jjp,k))-5.*(-putin2(i,jjm,k)+putin2(i,jjpp,k))+(-putin2(i,jjmm,k)+putin2(i,jjppp,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!		
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
	    putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cF*Rpdphi_ip 	  
	  ENDIF
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  cD=putin2(i ,j,k)
	  cC=putin2(i ,jm,k)
	  cU=putin2(i ,jmm-1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
		cfl = dt*(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
!		cF = 
!     1 DDy*(Ay*(putin2(i,jjm,k)+putin2(i,j,k))+By*(putin2(i,jjmm,k)+putin2(i,jjp,k))+Cy*(putin2(i,jjmmm,k)+putin2(i,jjpp,k))) -		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(i,jjm,k)+putin2(i,j,k))-5.*(-putin2(i,jjmm,k)+putin2(i,jjp,k))+(-putin2(i,jjmmm,k)+putin2(i,jjpp,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ELSE
	  cD=putin2(i, jm,k)
	  cC=putin2(i ,j ,k)
	  cU=putin2(i ,jp,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
		cfl = dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
!		cF = 
!     1 DDy*(Ay*(putin2(i,jjm,k)+putin2(i,j,k))+By*(putin2(i,jjmm,k)+putin2(i,jjp,k))+Cy*(putin2(i,jjmmm,k)+putin2(i,jjpp,k))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(i,jjm,k)+putin2(i,j,k))-5.*(-putin2(i,jjmm,k)+putin2(i,jjp,k))+(-putin2(i,jjmmm,k)+putin2(i,jjpp,k)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ENDIF	
	ENDIF

	IF (Wvel(i,j,k ).ge.0.) THEN
	  cD=putin2(i,j ,kp)
	  cC=putin2(i ,j,k )
	  cU=putin2(i,j ,km)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,k ))*dz_i
!		cF = 
!     1 DDz*(Az*(putin2(i,j,k)+putin2(i,j,kkp))+Bz*(putin2(i,j,kkm)+putin2(i,j,kkpp))+Cz*(putin2(i,j,kkmm)+putin2(i,j,kkppp)))-		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,j,kkp))-5.*(-putin2(i,j,kkm)+putin2(i,j,kkpp))+(-putin2(i,j,kkmm)+putin2(i,j,kkppp)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i 
	  ENDIF
	  
	ELSE
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,kp)
	  cU=putin2(i,j ,kpp+1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,k ))*dz_i
!		cF = 
!     1 DDz*(Az*(putin2(i,j,k)+putin2(i,j,kkp))+Bz*(putin2(i,j,kkm)+putin2(i,j,kkpp))+Cz*(putin2(i,j,kkmm)+putin2(i,j,kkppp))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(i,j,k)+putin2(i,j,kkp))-5.*(-putin2(i,j,kkm)+putin2(i,j,kkpp))+(-putin2(i,j,kkmm)+putin2(i,j,kkppp)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i
	  ENDIF
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,km)
	  cU=putin2(i,j ,kmm-1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,km ))*dz_i
!		cF = 
!     1 DDz*(Az*(putin2(i,j,kkm)+putin2(i,j,k))+Bz*(putin2(i,j,kkmm)+putin2(i,j,kkp))+Cz*(putin2(i,j,kkmmm)+putin2(i,j,kkpp)))-		! pos velocity minus numdif
!     1 numdif*(10.*(-putin2(i,j,kkm)+putin2(i,j,k))-5.*(-putin2(i,j,kkmm)+putin2(i,j,kkp))+(-putin2(i,j,kkmmm)+putin2(i,j,kkpp)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ELSE
	  cD=putin2(i,j ,km )
	  cC=putin2(i ,j,k)
	  cU=putin2(i,j ,kp)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,km ))*dz_i
!		cF = 
!     1 DDz*(Az*(putin2(i,j,kkm)+putin2(i,j,k))+Bz*(putin2(i,j,kkmm)+putin2(i,j,kkp))+Cz*(putin2(i,j,kkmmm)+putin2(i,j,kkpp))) +		! neg velocity plus numdif
!     1 numdif*(10.*(-putin2(i,j,kkm)+putin2(i,j,k))-5.*(-putin2(i,j,kkmm)+putin2(i,j,kkp))+(-putin2(i,j,kkmmm)+putin2(i,j,kkpp)))
!		phif = (cF-cU)/del
!		phif = (1.-cfl)*phif+cfl*phic
!
!		phif = max(phif,phic)
!		phif = min(phif,1.)
!	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
!		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
!		phif = min((1.-cfl)*(1-phic)+phic,phif) ! this limiter corresponds to Sweby TVD 2 (for TVD needed in practice)
		phif = min(1.,2.*phic)
		phif=(1-cfl)*phif+cfl*phic
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ENDIF
	ENDIF

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
	  
	subroutine advecc_NVD_fromm(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
!	Face-wise advection scheme according to NVD is similar to TVD, but manner of expression is different
!	See for example:
!	# Hirsch Chap 8.4.
!	# Darwish 1994, Normalized Variable and Space Formulation Methodology for High-Resolution Schemes
!	# Leonard 1991, The ULTIMATE conservative difference scheme applied to unsteady one-dimensional advection
!	The only real difference I could find is the cfl dependent limiter phif = min(phif,phic/(cfl+eps2)) 
!	which I cannot find a TVD equivalent for and which is crucial to the NVD succes to give no under/overshoot
!   it turned out to be necessary to include an extra limiter phif = min(phif,2.*phic) which corresponds to Sweby TVD 2r 
!	in practice TUDflow3d needs 2r limiter in TVD (and not 2r/cfl as theory claims) and also NVD needs this extra limiter min(phif,2.*phic) in practice!
!	20-11-2017 in matlab NVD FROMM or QUICK scheme with proper NVD limiters give stable results with EE1 time integration 
!   with less diffusion then TVD, and no under or overshoot. FROMM shape is slightly better than QUICK, and NVD Fromm gave best results in diagonal edge and diagonal square tests TUDflow3d
!	(<diff than Superbee, vLeer, Arora, NVD-upw5 and smallest waviness edge of all just like vLeer) therefore: 22-11-2017 NVD scheme is switched to NVD-FROMM
! can be called by nerd-option advec_conc='NVD'
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real xD,xC,xU,xF,cD,cC,cU,cF,phic,phif,eps,xxf,xxc,cfl,alpha,del,Rpdr_ip,Rpdphi_ip,eps2

      dz_i=1./dz
	  eps=1.e-6
	  eps2=1.e-2
	  
	  putout(ib:ie,jb:je,kb:ke)=0.

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

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



	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      Rpdr_i=1./(Rp2(i)*dr(i))
	  Rpdr_ip=1./(Rp2(ip)*dr(ip))

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
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  Rpdphi_ip=1./(Rp2(i)*(phiv(jp)-phiv(j)))
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
	  cD=putin2(ip,j,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(im,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i !UPW 
		putout(ip,j,k)= putout(ip,j,k)+ Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip !UPW 
	  ELSE
	    xD=Rp2(ip)
	    xC=Rp2(i )
	    xU=Rp2(im)
		xF=Ru (i)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
        putout(i,j,k)  = putout(i,j,k)  - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF
	ELSE
	  cD=putin2(i    ,j,k)
	  cC=putin2(ip   ,j,k)
	  cU=putin2(ipp+1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i  !UPW 
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cC*Rpdr_ip	 !UPW 	
	  ELSE
	    xD=Rp2(i )
	    xC=Rp2(ip)
	    xU=Rp2(ipp+1)
		xF=Ru (i)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
	    putout(i,j,k)  = putout(i ,j,k) - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	    putout(ip,j,k) = putout(ip,j,k) + Ru(i)*Uvel(i ,j,k)*cF*Rpdr_ip		
	  ENDIF	
	ENDIF
	IF (i.eq.1) THEN
	IF (Uvel(im,j,k).ge.0.) THEN
	  cD=putin2(i,j,k)
	  cC=putin2(im ,j,k)
	  cU=putin2(imm-1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW
	  ELSE
	    xD=Rp2(i)
	    xC=Rp2(im )
	    xU=Rp2(imm-1)
		xF=Ru (im)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF
	ELSE
	  cD=putin2(im,j,k)
	  cC=putin2(i   ,j,k)
	  cU=putin2(ip,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW 
	  ELSE
	    xD=Rp2(im)
	    xC=Rp2(i )
	    xU=Rp2(ip)
		xF=Ru (im)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*ABS(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF		
	ENDIF
	ENDIF

	IF (Vvel(i,j ,k).ge.0.) THEN
	  cD=putin2(i,jp,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(i,jm,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+jp)
	    xC=Rp2(i)*phipt2(rank*je+j )
	    xU=Rp2(i)*phipt2(rank*je+jm)
		xF=Rp2(i)*phiv(j )
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cF*Rpdphi_ip 
	  ENDIF
	ELSE
	  cD=putin2(i, j,k)
	  cC=putin2(i ,jp,k)
	  cU=putin2(i,jpp+1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
		putout(i,jp,k)= putout(i,jp,k)+ Vvel2(i,j ,k)*cC*Rpdphi_ip !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+j)
	    xC=Rp2(i)*phipt2(rank*je+jp )
	    xU=Rp2(i)*phipt2(rank*je+jpp+1)
		xF=Rp2(i)*phiv(j )
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
		putout(i,jp,k) = putout(i,jp,k) + Vvel2(i,j ,k)*cF*Rpdphi_ip 	  
	  ENDIF
	ENDIF
	IF (j.eq.1) THEN
	IF (Vvel(i,jm,k).ge.0.) THEN
	  cD=putin2(i ,j,k)
	  cC=putin2(i ,jm,k)
	  cU=putin2(i ,jmm-1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+j)
	    xC=Rp2(i)*phipt2(rank*je+jm )
	    xU=Rp2(i)*phipt2(rank*je+jmm-1)
		xF=Rp2(i)*phiv(jm)
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ELSE
	  cD=putin2(i, jm,k)
	  cC=putin2(i ,j ,k)
	  cU=putin2(i ,jp,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+jm)
	    xC=Rp2(i)*phipt2(rank*je+j )
	    xU=Rp2(i)*phipt2(rank*je+jp)
		xF=Rp2(i)*phiv(jm)
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ENDIF	
	ENDIF

	IF (Wvel(i,j,k ).ge.0.) THEN
	  cD=putin2(i,j ,kp)
	  cC=putin2(i ,j,k )
	  cU=putin2(i,j ,km)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,k ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = 0.375+0.75*phic	!QUICK on equidistant grid --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i 
	  ENDIF
	  
	ELSE
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,kp)
	  cU=putin2(i,j ,kpp+1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,k ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = 0.375+0.75*phic	!QUICK on equidistant grid --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
		putout(i,j,kp)= putout(i,j,kp)+ Wvel(i,j,k )*cF*dz_i
	  ENDIF
	ENDIF
	IF (k.eq.1) THEN
	IF (Wvel(i,j,km).ge.0.) THEN
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,km)
	  cU=putin2(i,j ,kmm-1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,km ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = 0.375+0.75*phic	!QUICK on equidistant grid --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ELSE
	  cD=putin2(i,j ,km )
	  cC=putin2(i ,j,k)
	  cU=putin2(i,j ,kp)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,km ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) for time dependent problems
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = 0.375+0.75*phic	!QUICK on equidistant grid --> need to be combined with phif=(1-cfl)*phif+cfl*phic
		!phif = (1.-cfl)*phif+cfl*phic
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  ! this limiter corresponds to Sweby TVD 2r/cfl (for TVD scheme also not stable)
		phif = min(phif,2.*phic)  ! this limiter corresponds to Sweby TVD 2r (for NVD and TVD needed in practice to avoid <0)
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ENDIF
	ENDIF

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
	  
	  	subroutine advecc_NVDold2(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
     +                  i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,periodicx,periodicy)
      implicit none
!	Cell-wise advection scheme according to NVD is similar to TVD, but manner of expression is different
!	See for example:
!	# Hirsch Chap 8.4.
!	# Darwish 1994, Normalized Variable and Space Formulation Methodology for High-Resolution Schemes
!	# Leonard 1991, The ULTIMATE conservative difference scheme applied to unsteady one-dimensional advection
!	The only real difference I could find is the cfl dependent
!	limiter phif = min(phif,phic/(cfl+eps2)) which I cannot find a TVD equivalent for and which is crucial to the NVD succes to give no under/overshoot!!
!	test results for NVD alpha scheme are same for TVD alpha scheme, but TVD has advantage to
!	deal with non-regular grids in proper manner; therefore I stay on the TVD track, LdW 25-11-15
!	20-11-2017 in matlab NVD FROMM or QUICK scheme with proper NVD limiters give stable results with EE1 time integration without need for (1-cfl) LW term
!   hence much less diffusion then TVD, no under or overshoot in matlab, but artificial steepening of sine wave is visible
!   20-11-2017 NVD scheme is switched to NVD-QUICK
! can be called by nerd-option advec_conc='NVD'
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,Vvel2(0:i1,0:j1,0:k1)
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)
	real xD,xC,xU,xF,cD,cC,cU,cF,phic,phif,eps,xxf,xxc,cfl,alpha,del,eps2

      dz_i=1./dz
	  eps=1.e-6
	  eps2=1.e-4

	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	

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



	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      Rpdr_i=1./(Rp2(i)*dr(i))

    
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
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))		
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
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
	  cD=putin2(ip,j,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(im,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i !UPW 
	  ELSE
	    xD=Rp2(ip)
	    xC=Rp2(i )
	    xU=Rp2(im)
		xF=Ru (i)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  
		cF = phif*del+cU
		putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	  ENDIF
	ELSE
	  cD=putin2(i    ,j,k)
	  cC=putin2(ip   ,j,k)
	  cU=putin2(ipp+1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cC*Rpdr_i !UPW 
	  ELSE
	    xD=Rp2(i )
	    xC=Rp2(ip)
	    xU=Rp2(ipp+1)
		xF=Ru (i)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*ABS(Uvel(i ,j,k))/(Rp2(ip)-Rp2(i )))
		phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = - Ru(i)*Uvel(i ,j,k)*cF*Rpdr_i
	  ENDIF		  
	ENDIF
	IF (Uvel(im,j,k).ge.0.) THEN
	  cD=putin2(i,j,k)
	  cC=putin2(im ,j,k)
	  cU=putin2(imm-1,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW
	  ELSE
	    xD=Rp2(i)
	    xC=Rp2(im )
	    xU=Rp2(imm-1)
		xF=Ru (im)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
		phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF		  
		  
	ELSE
	  cD=putin2(im,j,k)
	  cC=putin2(i   ,j,k)
	  cU=putin2(ip,j,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cC*Rpdr_i !UPW 
	  ELSE
	    xD=Rp2(im)
	    xC=Rp2(i )
	    xU=Rp2(ip)
		xF=Ru (im)		
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = (dt*ABS(Uvel(im,j,k))/(Rp2(i)-Rp2(im)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Ru(im)*Uvel(im,j,k)*cF*Rpdr_i
	  ENDIF		  
	  
	ENDIF

	IF (Vvel2(i,j ,k).ge.0.) THEN
	  cD=putin2(i,jp,k)
	  cC=putin2(i ,j,k)
	  cU=putin2(i,jm,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+jp)
	    xC=Rp2(i)*phipt2(rank*je+j )
	    xU=Rp2(i)*phipt2(rank*je+jm)
		xF=Rp2(i)*phiv(j )
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
	  ENDIF
	  
	ELSE
	  cD=putin2(i, j,k)
	  cC=putin2(i ,jp,k)
	  cU=putin2(i,jpp+1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+j)
	    xC=Rp2(i)*phipt2(rank*je+jp )
	    xU=Rp2(i)*phipt2(rank*je+jpp+1)
		xF=Rp2(i)*phiv(j )
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*ABS(Vvel2(i,j ,k))/(Rp(i)*(phipt2(rank*je+jp)-phipt2(rank*je+j)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU		
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2))  
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Vvel2(i,j ,k)*cF*Rpdphi_i 
	  ENDIF

	  
	ENDIF		
	IF (Vvel2(i,jm,k).ge.0.) THEN
	  cD=putin2(i ,j,k)
	  cC=putin2(i ,jm,k)
	  cU=putin2(i ,jmm-1,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+j)
	    xC=Rp2(i)*phipt2(rank*je+jm )
	    xU=Rp2(i)*phipt2(rank*je+jmm-1)
		xF=Rp2(i)*phiv(jm)
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	ELSE
	  cD=putin2(i, jm,k)
	  cC=putin2(i ,j ,k)
	  cU=putin2(i ,jp,k)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cC*Rpdphi_i !UPW
	  ELSE
	    xD=Rp2(i)*phipt2(rank*je+jm)
	    xC=Rp2(i)*phipt2(rank*je+j )
	    xU=Rp2(i)*phipt2(rank*je+jp)
		xF=Rp2(i)*phiv(jm)
	    xxf  = (xF-xU)/(xD-xU) 
		xxc  = (xC-xU)/(xD-xU)
		cfl = dt*ABS(Vvel2(i,jm,k))/(Rp(i)*(phipt2(rank*je+j)-phipt2(rank*je+jm)))
	    phif = phic + (1.-cfl)*(xxf-xxc) !FROMM on non-equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + xxf-xxc !FROMM on non-equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Vvel2(i,jm,k)*cF*Rpdphi_i 
	  ENDIF
	  
	ENDIF	

	IF (Wvel(i,j,k ).ge.0.) THEN
	  cD=putin2(i,j ,kp)
	  cC=putin2(i ,j,k )
	  cU=putin2(i,j ,km)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,k ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!phif = 0.375+0.75*phic	! QUICK on equidistant grid without CFL correction
		!phif = max(phif,phic)
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
	  ENDIF
	  
	ELSE
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,kp)
	  cU=putin2(i,j ,kpp+1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,k ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!phif = 0.375+0.75*phic	! QUICK on equidistant grid without CFL correction
		!phif = max(phif,phic)
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) - Wvel(i,j,k )*cF*dz_i 
	  ENDIF
	  
	ENDIF
	IF (Wvel(i,j,km).ge.0.) THEN
	  cD=putin2(i,j ,k )
	  cC=putin2(i ,j,km)
	  cU=putin2(i,j ,kmm-1)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*(Wvel(i,j,km ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!phif = 0.375+0.75*phic	! QUICK on equidistant grid without CFL correction
		!phif = max(phif,phic)
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ELSE
	  cD=putin2(i,j ,km )
	  cC=putin2(i ,j,k)
	  cU=putin2(i,j ,kp)
	  del = (cD-cU)
	  phic=(cC-cU)/del
	  IF (ABS(del)<eps.or.phic.gt.1.or.phic.lt.0.) THEN
	    putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cC*dz_i 
	  ELSE
		cfl = dt*ABS(Wvel(i,j,km ))*dz_i
	    phif = phic + (1.-cfl)*0.25 !FROMM on equidistant grid (Darwish 1994) with CFL correction --> this gave best results in matlab tests and less CPU
		!phif = phic + 0.25 !FROMM on equidistant grid (Darwish 1994) and without CFL correction --> this gave best results in matlab tests and less CPU
		!phif = xxf + (xxf*(xxf-1.))/(xxc*(xxc-1.))*(phic-xxc) !QUICK on non-equidistant grid (Darwish 1994) and without CFL correction
		!phif = 0.375+0.75*phic	! QUICK on equidistant grid without CFL correction
		!phif = max(phif,phic)
		!alpha=0.6667-0.3333*cfl
		!phif=max(0.,min(min((1.5-alpha)*phic+0.5*alpha,2.*phic),1.)) !Third order alpha scheme with CFL dependence; no correction for variable grid
		phif = max(phif,phic)
		phif = min(phif,1.)
	    phif = min(phif,phic/(cfl+eps2)) 
		cF = phif*del+cU
		putout(i,j,k) = putout(i,j,k) + Wvel(i,j,km )*cF*dz_i 
	  ENDIF	  
	ENDIF
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

	  

	subroutine advecc_TVD_old2(putout,putin,Uvel,Vvel,Wvel,rho,Ru,Rp,dr,phiv,phipt,dz,
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer kmm,kpp,imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),
     +         Uvel(0:i1,0:j1,0:k1),rho(0:i1,0:j1,0:k1),
     +         Vvel(0:i1,0:j1,0:k1),Wvel(0:i1,0:j1,0:k1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m,putin2(-1:i1+1,-1:j1+1,0:k1),Rp2(-1:i1+1)
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt
	real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)

      dz_i=1./dz


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim

	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	
	
	putin2(0:i1,0:j1,0:k1)=putin

	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
    
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

      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
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

	rRpos = (putin2(i    ,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
  	rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	noemer = putin2(i,j,k)-putin2(im,j,k)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(ip,j,k)-putin2(i    ,j,k))/noemer *varx_grid_Lneg
	rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos

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
		
	rRpos = (putin2(i,j    ,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 

	noemer = putin2(i,j,k)-putin2(i,jm,k)
	  noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	
	rLneg = (putin2(i,jp,k)-putin2(i,j    ,k))/noemer*vary_grid_Lneg 
	rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos


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

	subroutine adveccbot_TVD(putout,putin,Uvel,Vvel,Ru,Rp,dr,phiv,phipt,dz,i1,j1,ib,ie,jb,je,dt,rank,px,periodicx,periodicy)
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
	  integer imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1),putin(0:i1,0:j1),putin2(-1:i1+1,-1:j1+1),
     +         Uvel(0:j1),Vvel(0:j1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Rp2(-1:i1+1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,cflL,cflR
      real pbf(0:i1)
      real pbb(0:i1)


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim
	
	putin2(0:i1,0:j1)=putin

	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1)=putin(0,0:j1)
		putin2(i1+1,0:j1)=putin(ie,0:j1)
	else 
		putin2(-1,0:j1)=putin(ie-1,0:j1)
		putin2(i1+1,0:j1)=putin(2,0:j1)
	endif
!c get stuff from other CPU's
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
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
    
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
	
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  
	noemer = putin2(ip,j)-putin2(i,j)
        noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)


	rRpos = (putin2(i,j)-putin2(im,j))/noemer *varx_grid_Rpos
  	rRneg = (putin2(ipp+1,j)-putin2(ip,j))/noemer *varx_grid_Rneg
	noemer = putin2(i,j)-putin2(im,j)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(ip,j)-putin2(i,j))/noemer *varx_grid_Lneg
	rLpos = (putin2(im,j)-putin2(imm-1,j))/noemer *varx_grid_Lpos
	cflR = dt*ABS(Uvel(j))/(Rp(ip)-Rp(i))
	cflL = dt*ABS(Uvel(j))/(Rp(i)-Rp(im))
!	cRpos = putin2(i ,j) + 0.5*limiter2(rRpos,cflR)*(putin2(ip,j) - putin2(i ,j))*(1.-cflR)
!	cLpos = putin2(im,j) + 0.5*limiter2(rLpos,cflL)*(putin2( i,j) - putin2(im,j))*(1.-cflL)
!	cRneg = putin2(ip,j) + 0.5*limiter2(rRneg,cflR)*(putin2( i,j) - putin2(ip,j))*(1.-cflR)
!	cLneg = putin2(i ,j) + 0.5*limiter2(rLneg,cflL)*(putin2(im,j) - putin2(i ,j))*(1.-cflL)
	cRpos = putin2(i ,j) + 0.5*limiter(rRpos)*(putin2(ip,j) - putin2(i ,j))*(1.-cflR)
	cLpos = putin2(im,j) + 0.5*limiter(rLpos)*(putin2( i,j) - putin2(im,j))*(1.-cflL)
	cRneg = putin2(ip,j) + 0.5*limiter(rRneg)*(putin2( i,j) - putin2(ip,j))*(1.-cflR)
	cLneg = putin2(i ,j) + 0.5*limiter(rLneg)*(putin2(im,j) - putin2(i ,j))*(1.-cflL)
      putout(i,j) = - (
     &   Ru(i)  * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cRpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cRneg ) -
     &   Ru(im) * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cLpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cLneg ) )
     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	noemer = putin2(i,jp)-putin2(i,j)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
		
	rRpos = (putin2(i,j)-putin2(i,jm))/noemer*vary_grid_Rpos
	rRneg = (putin2(i,jpp+1)-putin2(i,jp))/noemer*vary_grid_Rneg 

	noemer = putin2(i,j)-putin2(i,jm)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	
	rLneg = (putin2(i,jp)-putin2(i,j))/noemer *vary_grid_Lneg
	rLpos = (putin2(i,jm)-putin2(i,jmm-1))/noemer*vary_grid_Lpos
	cflL = dt*ABS(Vvel(jm))*Rpdphi_i
	cflR = dt*ABS(Vvel(j ))*Rpdphi_i

!	cRpos = putin2(i ,j) + 0.5*limiter2(rRpos,cflR)*(putin2(i,jp) - putin2(i ,j))*(1.-cflR)
!	cLpos = putin2(i,jm) + 0.5*limiter2(rLpos,cflL)*(putin2( i,j) - putin2(i,jm))*(1.-cflL)
!	cRneg = putin2(i,jp) + 0.5*limiter2(rRneg,cflR)*(putin2( i,j) - putin2(i,jp))*(1.-cflR)
!	cLneg = putin2(i,j ) + 0.5*limiter2(rLneg,cflL)*(putin2(i,jm) - putin2(i,j ))*(1.-cflL)
	cRpos = putin2(i ,j) + 0.5*limiter(rRpos)*(putin2(i,jp) - putin2(i ,j))*(1.-cflR)
	cLpos = putin2(i,jm) + 0.5*limiter(rLpos)*(putin2( i,j) - putin2(i,jm))*(1.-cflL)
	cRneg = putin2(i,jp) + 0.5*limiter(rRneg)*(putin2( i,j) - putin2(i,jp))*(1.-cflR)
	cLneg = putin2(i,j ) + 0.5*limiter(rLneg)*(putin2(i,jm) - putin2(i,j ))*(1.-cflL)

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

	subroutine adveccbot3d_TVD(putout,putin,Uvel,Vvel,Ru,Rp,dr,phiv,phipt,dz,i1,j1,k1,ib,ie,jb,je,kb,ke,dt,rank,px,
     & periodicx,periodicy)
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
c          dr,phiv,dz        : grid spacing in r, phi and z-direction
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
      integer  i,j,k,im,ip,jm,jp,km,kp,kb,ke,k1,i1,j1,ib,ie,jb,je
	  integer imm,ipp,jmm,jpp,rank,px,periodicx,periodicy
      real     putout(0:i1,0:j1,0:k1),putin(0:i1,0:j1,0:k1),putin2(-1:i1+1,-1:j1+1,0:k1),
     +         Uvel(0:j1),Vvel(0:j1),
     +         dr(0:i1),phiv(0:j1),phipt(0:je*px+1),dz,Ru(0:i1),Rp(0:i1),Rp2(-1:i1+1),
     +         noemer,rRpos,rRneg,rLpos,rLneg,cRpos,cRneg,cLpos,cLneg,limiter,
     +         rho_p,rho_m
      
      real 	dz_i,Rpdr_i,Rpdphi_i,varx_grid_Rpos,varx_grid_Rneg,varx_grid_Lpos,varx_grid_Lneg
	  real  vary_grid_Rpos,vary_grid_Rneg,vary_grid_Lpos,vary_grid_Lneg,phipt2(-1:je*px+2)
	real  dt,cflL,cflR
	  real*8 pbb(0:i1,0:k1),pbf(0:i1,0:k1)


	Rp2(0:i1)=Rp
	Rp2(-1)=Rp(0)-(Rp(2)-Rp(1)) !needed for periodicx sim
 	Rp2(i1+1)=Rp(i1)+(Rp(i1)-Rp(i1-1)) !needed for periodicx sim
	
	phipt2(0:je*px+1)=phipt
	phipt2(-1)=phipt(0)-(phipt(1)-phipt(0)) !needed for periodicy sim
 	phipt2(je*px+1+1)=phipt(je*px+1)+(phipt(je*px+1)-phipt(je*px+1-1)) !needed for periodicy sim


	putin2(0:i1,0:j1,0:k1)=putin

	if (periodicx.eq.0.or.periodicx.eq.2) then
		putin2(-1,0:j1,0:k1)=putin(0,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(ie,0:j1,0:k1)
	else 
		putin2(-1,0:j1,0:k1)=putin(ie-1,0:j1,0:k1)
		putin2(i1+1,0:j1,0:k1)=putin(2,0:j1,0:k1)
	endif
!c get stuff from other CPU's
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



      dz_i=1./dz

c
c     -------------------------------------------start i-loop
      do 100 i=ib,ie
c
      ip=i+1
      im=i-1
      imm=im
      ipp=ip
      if (periodicx.eq.0.or.periodicx.eq.2) then
        if (i.eq.1) imm=i
        if (i.eq.ie) ipp=i
      endif
      varx_grid_Rpos=(Rp2(ip)-Rp2(i))/(Rp2(i)-Rp2(im))
      varx_grid_Rneg=(Rp2(ipp)-Rp2(i))/(Rp2(ipp+1)-Rp2(ipp))
      varx_grid_Lpos=(Rp2(i)-Rp2(imm))/(Rp2(imm)-Rp2(imm-1))
      varx_grid_Lneg=(Rp2(i)-Rp2(im))/(Rp2(ip)-Rp2(i))
      Rpdr_i=1./(Rp2(i)*dr(i))
    
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
	
      vary_grid_Rpos=(phipt2(rank*je+jp) -phipt2(rank*je+j))  /(phipt2(rank*je+j)    -phipt2(rank*je+jm))
      vary_grid_Rneg=(phipt2(rank*je+jpp)-phipt2(rank*je+j))  /(phipt2(rank*je+jpp+1)-phipt2(rank*je+jpp))
      vary_grid_Lpos=(phipt2(rank*je+j)  -phipt2(rank*je+jmm))/(phipt2(rank*je+jmm)  -phipt2(rank*je+jmm-1))
      vary_grid_Lneg=(phipt2(rank*je+j)  -phipt2(rank*je+jm)) /(phipt2(rank*je+jp)   -phipt2(rank*je+j))	
      Rpdphi_i=1./(Rp2(i)*(phiv(j)-phiv(jm)))
	  
	  do k=kb,ke
	noemer = putin2(ip,j,k)-putin2(i,j,k)
        noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)


	rRpos = (putin2(i,j,k)-putin2(im,j,k))/noemer *varx_grid_Rpos
  	rRneg = (putin2(ipp+1,j,k)-putin2(ip,j,k))/noemer *varx_grid_Rneg
	noemer = putin2(i,j,k)-putin2(im,j,k)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)

	rLneg = (putin2(ip,j,k)-putin2(i,j,k))/noemer *varx_grid_Lneg
	rLpos = (putin2(im,j,k)-putin2(imm-1,j,k))/noemer *varx_grid_Lpos
	cflR = dt*ABS(Uvel(j))/(Rp(ip)-Rp(i))
	cflL = dt*ABS(Uvel(j))/(Rp(i)-Rp(im))
!	cRpos = putin2(i ,j) + 0.5*limiter2(rRpos,cflR)*(putin2(ip,j) - putin2(i ,j))*(1.-cflR)
!	cLpos = putin2(im,j) + 0.5*limiter2(rLpos,cflL)*(putin2( i,j) - putin2(im,j))*(1.-cflL)
!	cRneg = putin2(ip,j) + 0.5*limiter2(rRneg,cflR)*(putin2( i,j) - putin2(ip,j))*(1.-cflR)
!	cLneg = putin2(i ,j) + 0.5*limiter2(rLneg,cflL)*(putin2(im,j) - putin2(i ,j))*(1.-cflL)
	cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(ip,j,k) - putin2(i ,j,k))*(1.-cflR)
	cLpos = putin2(im,j,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(im,j,k))*(1.-cflL)
	cRneg = putin2(ip,j,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(ip,j,k))*(1.-cflR)
	cLneg = putin2(i ,j,k) + 0.5*limiter(rLneg)*(putin2(im,j,k) - putin2(i ,j,k))*(1.-cflL)
      putout(i,j,k) = - (
     &   Ru(i)  * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cRpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cRneg ) -
     &   Ru(im) * ( 0.5*(Uvel(j)+ABS(Uvel(j)))*cLpos + 0.5*(Uvel(j)-ABS(Uvel(j)))*cLneg ) )
     &  *Rpdr_i !/ ( Rp(i) * dr(i) ) 

	noemer = putin2(i,jp,k)-putin2(i,j,k)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
		
	rRpos = (putin2(i,j,k)-putin2(i,jm,k))/noemer*vary_grid_Rpos
	rRneg = (putin2(i,jpp+1,k)-putin2(i,jp,k))/noemer*vary_grid_Rneg 

	noemer = putin2(i,j,k)-putin2(i,jm,k)
	noemer=MAX(ABS(noemer),1.e-6)*sign(1.,noemer)
	
	rLneg = (putin2(i,jp,k)-putin2(i,j,k))/noemer *vary_grid_Lneg
	rLpos = (putin2(i,jm,k)-putin2(i,jmm-1,k))/noemer*vary_grid_Lpos
	cflL = dt*ABS(Vvel(jm))*Rpdphi_i
	cflR = dt*ABS(Vvel(j ))*Rpdphi_i

!	cRpos = putin2(i ,j) + 0.5*limiter2(rRpos,cflR)*(putin2(i,jp) - putin2(i ,j))*(1.-cflR)
!	cLpos = putin2(i,jm) + 0.5*limiter2(rLpos,cflL)*(putin2( i,j) - putin2(i,jm))*(1.-cflL)
!	cRneg = putin2(i,jp) + 0.5*limiter2(rRneg,cflR)*(putin2( i,j) - putin2(i,jp))*(1.-cflR)
!	cLneg = putin2(i,j ) + 0.5*limiter2(rLneg,cflL)*(putin2(i,jm) - putin2(i,j ))*(1.-cflL)
	cRpos = putin2(i ,j,k) + 0.5*limiter(rRpos)*(putin2(i,jp,k) - putin2(i ,j,k))*(1.-cflR)
	cLpos = putin2(i,jm,k) + 0.5*limiter(rLpos)*(putin2( i,j,k) - putin2(i,jm,k))*(1.-cflL)
	cRneg = putin2(i,jp,k) + 0.5*limiter(rRneg)*(putin2( i,j,k) - putin2(i,jp,k))*(1.-cflR)
	cLneg = putin2(i,j ,k) + 0.5*limiter(rLneg)*(putin2(i,jm,k) - putin2(i,j ,k))*(1.-cflL)

	putout(i,j,k) = putout(i,j,k) - (
     &   (  0.5*(Vvel(j )+ABS(Vvel(j )))*cRpos + 0.5*(Vvel(j )-ABS(Vvel(j )))*cRneg )  -
     &   (  0.5*(Vvel(jm)+ABS(Vvel(jm)))*cLpos + 0.5*(Vvel(jm)-ABS(Vvel(jm)))*cLneg ) )
     &  *Rpdphi_i !/ ( Rp(i) * dphi ) 
	 
		enddo

200     continue
c     -------------------------------------------end j-loop
100   continue
c     -------------------------------------------end i-loop
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

      real function limiter2(r,cfl)
      real r,cfl,alpha
! limiter = third order accurate (with cfl=0 equal to Koren limiter) CFL dependent limiter from Arora and Roe 1997
	alpha=(1.+cfl)/3.
	limiter2=MAX(0.,MIN(MIN(2./MAX(cfl,1e-12)*r,1.+alpha*(r-1.)),1.99/(1.-cfl))) !looser TVD contstraints gives in TUDflow3d e-7 undershoot, therefore not advised...

	!limiter2=MAX(0.,MIN(2.*r,1.),MIN(r,1.99/(1.-cfl)))	! Superbee limiter loosened CFL-dependent restrictions --> almost equal to HyperBEE (difference 2*r vs 2/cfl*r) in matlab good; not stable here...
      return
      end
	  
      real function limiter4(r,cfl)
      real r,cfl,alpha

	 limiter4=MAX(0.,MIN(2.*r,1.),MIN(r,2.)) ! Superbee
	 !limiter4=MAX(0.,MIN(2.*r,1.),MIN(r,1.99/(1.-cfl)))	! Superbee limiter loosened CFL-dependent restrictions --> almost equal to HyperBEE (difference 2*r vs 2/cfl*r) because less strange sharp edges
	 !limiter4=MAX(0.,MIN(2.*r/MAX(cfl,1e-12),1.),MIN(r,1.99/(1.-cfl)))	! Hyperbee limiter loosened CFL-dependent restrictions --> very aggresively sharpening edges 
	 !limiter4=(r+ABS(r))/MAX(1.+r,1.)  !Van Leer limiter
	 !limiter4=0. !upw1 for testing purpose
      return
      end	  
	  
