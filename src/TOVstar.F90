!-*-f90-*-
module tov_mod
  real*8 :: tov_R, tov_A, tov_tkick, tov_Akick, tov_alpBH
  integer :: tov_ptype = 0
  integer :: tov_ndiag
  logical :: kicked = .false.
  character(len=256) :: tov_file
end module tov_mod

subroutine TOVstar
  use GR1D_module
  use ideal_eos_module
  use tov_mod
  use atmos
  implicit none
  integer np,i,j
  real*8, allocatable :: pr(:),prho(:),peps(:)
  real*8 xmax, wt, rr, rhoatm, ptmp, d1, d2, d3, d4
  open(77,file='tov_params.in',status='old')
  read(77,'(a)') tov_file
  read(77,*) tov_A
  read(77,*) tov_tkick
  read(77,*) tov_Akick
  read(77,*) tov_ndiag
  read(77,*) tov_alpBH
  read(77,*) xmax
  read(77,*,end=91) tov_ptype
91 continue
  close(77)
  idealgamma = 2.0d0
  inquire(file='pp_eos.in',exist=pp_on)
  if (pp_on) then
     open(76,file='pp_eos.in',status='old')
     read(76,*) npp, pp_gth
     do i=1,npp
        read(76,*) pp_rlo(i), pp_K(i), pp_G(i), pp_a(i)
     enddo
     close(76)
     write(*,*) 'Using piecewise-polytrope EoS with',npp,' pieces, Gamma_th=',pp_gth
  endif
  inquire(file='conv.in',exist=conv_on)
  if (conv_on .and. pp_on) then
     open(75,file='conv.in',status='old')
     read(75,*) conv_mode, conv_par, conv_w
     close(75)
     open(76,file='pp_eos_had.in',status='old')
     read(76,*) npq
     do i=1,npq
        read(76,*) pq_rlo(i), pq_K(i), pq_G(i), pq_a(i)
     enddo
     close(76)
     write(*,*) 'Finite-speed conversion: mode',conv_mode,' par',conv_par,' width',conv_w,' hadronic pieces',npq
  else
     conv_on = .false.
  endif
  dynamic_output_control = .false.
  ntinfo = 100000000
  ntout = -1
  ntout_scalar = -1
  eoskey = 4
  open(78,file=trim(tov_file),status='old')
  read(78,*) np, tov_R
  allocate(pr(np),prho(np),peps(np))
  do i=1,np
     read(78,*) pr(i),prho(i),peps(i)
  enddo
  close(78)
  dt = 1.0d-3
  call grid(0.0d0,xmax,0.0d0)
  rhoatm = atmo_rho_abs_min*rho_gf
  call atmosphere_init
  do i=1,n1
     rr = x1(i)
     if (rr.lt.pr(np)) then
        j = 1
        do while (pr(j+1).lt.rr)
           j = j+1
        enddo
        wt = (rr-pr(j))/(pr(j+1)-pr(j))
        rho(i) = (1-wt)*prho(j)+wt*prho(j+1)
        eps(i) = (1-wt)*peps(j)+wt*peps(j+1)
     else
        rho(i) = 0.0d0; eps(i)=0.0d0
     endif
     if (rho(i).lt.rhoatm) then
        call atmos_eos(i,rho(i),press(i),eps(i),cs2(i))
     else
        if (pp_on) then
           ptmp = 0.0d0
           call ideal_eos(rho(i)/rho_gf,ptmp,d1,d2,d3,d4,1)
           eps(i) = ptmp*eps_gf
           press(i) = d1*press_gf
        else
           call ideal_eos(rho(i)/rho_gf,eps(i)/eps_gf,ptmp,d1,d2,d3,0)
           press(i) = ptmp*press_gf
        endif
     endif
     if (rr.lt.tov_R) then
        v(i) = tov_A*sin(pi*rr/tov_R)
     else
        v(i) = 0.0d0
     endif
     v1(i) = v(i)
     W(i) = 1.0d0/sqrt(1.0d0-v(i)**2)
     ye(i) = 0.5d0
  enddo
  call mass_interior
  call boundaries(0,0)
  mgrav(:) = mass(:)
  call GR_terms_initial
  call GR_boundaries
  call mass_interior
  if (tov_ptype.eq.1 .or. tov_ptype.eq.2) then
     ! uniform coordinate velocity U = dr/dt = V0 inside the star: v = a U / alpha, metric re-solved iteratively
     do j=1,6
        do i=1,n1
           if (x1(i).lt.tov_R) then
              if (tov_ptype.eq.1) v(i) = X(i)*tov_A/alp(i)
              if (tov_ptype.eq.2) v(i) = X(i)*tov_A*(x1(i)/tov_R)/alp(i)
           else
              v(i) = 0.0d0
           endif
           v1(i) = v(i); W(i) = 1.0d0/sqrt(1.0d0-v(i)**2)
        enddo
        call prim2con
        call boundaries(0,0)
        call con2GR
        call GR_alp
        call GR_boundaries
     enddo
     call prim2con
     call mass_interior
     write(*,*) "coordinate-velocity perturbation type",tov_ptype," V0=",tov_A,"  v_c=",v(ghosts1+1)," v_max=",minval(v)
  endif
  write(*,*) "TOV: R=",tov_R," A=",tov_A," Mgrav=",totalmass," Mb=",mass(n1-ghosts1+1)
end subroutine TOVstar

subroutine tov_step_hook
  use GR1D_module
  use tov_mod
  use ideal_eos_module
  implicit none
  integer i, ic, k
  real*8 munb, munbB, h, maxrho, emax, rq_c, rt_c, mw, mr, mq, mm
  ic = ghosts1+1
  ! mid-run kick
  if (tov_tkick.gt.0.0d0 .and. .not.kicked .and. time*1.0d3.ge.tov_tkick) then
     kicked = .true.
     do i=1,n1
        if (x1(i).lt.tov_R) then
           v(i) = v(i) + tov_Akick*sin(pi*x1(i)/tov_R)
           v1(i) = v(i)
           W(i) = 1.0d0/sqrt(1.0d0-v(i)**2)
        endif
     enddo
     do k=1,3
        call prim2con
        call boundaries(0,0)
        call con2GR
        call GR_alp
        call GR_boundaries
     enddo
     call prim2con
     call mass_interior
     write(*,*) "KICK applied at t(ms)=",time*1.0d3
  endif
  if (mod(nt,tov_ndiag).eq.0) then
     munb=0; munbB=0
     do i=ghosts1+1,n1-ghosts1
        h = 1.0d0+eps(i)+press(i)/rho(i)
        if (v(i).gt.0.0d0 .and. alp(i)*W(i).gt.1.0d0) munb = munb + q(i,1)*volume(i)
        if (v(i).gt.0.0d0 .and. h*alp(i)*W(i).gt.1.0d0) munbB = munbB + q(i,1)*volume(i)
     enddo
     ! smoothed diagnostics: mass-averaged rho for r<0.5 km, baryon mass in quark / mixed phase
     rq_c = 1.0d99; rt_c = 1.0d99
     if (pp_on) then
        do k=1,npp-1
           if (pp_rlo(k).gt.1.0d14 .and. pp_G(k).lt.1.0d0) then
              rt_c = pp_rlo(k)*rho_gf; rq_c = pp_rlo(k+1)*rho_gf
           endif
        enddo
     endif
     mw=0; mr=0; mq=0; mm=0
     do i=ghosts1+1,n1-ghosts1
        if (x1(i).lt.0.5d0/1.476625d0) then
           mw = mw + q(i,1)*volume(i); mr = mr + rho(i)*q(i,1)*volume(i)
        endif
        if (rho(i).ge.rq_c) mq = mq + q(i,1)*volume(i)
        if (rho(i).ge.rt_c .and. rho(i).lt.rq_c) mm = mm + q(i,1)*volume(i)
     enddo
     open(80,file=trim(adjustl(outdir))//'/tovdiag2.dat',position='append')
     write(80,'(7es16.8)') time*1.0d3, rho(ic), mr/mw, mq, mm, minval(v(ghosts1+1:n1-ghosts1)), maxval(v(ghosts1+1:n1-ghosts1))
     close(80)
     open(79,file=trim(adjustl(outdir))//'/tovdiag.dat',position='append')
     write(79,'(8es16.8)') time*1.0d3, rho(ic), alp(ic), minval(alp(ghosts1+1:n1-ghosts1)), &
          maxval(rho(ghosts1+1:n1-ghosts1)), munb, munbB, mass(n1-ghosts1+1)
     close(79)
  endif
  if (alp(ic).lt.tov_alpBH) then
     open(79,file=trim(adjustl(outdir))//'/BH',status='unknown')
     write(79,*) time*1.0d3
     close(79)
     write(*,*) "BH: central lapse below threshold at t(ms)=",time*1.0d3
     stop
  endif
end subroutine tov_step_hook
