program pisano_udsm_material_point_driver
  use, intrinsic :: ieee_arithmetic
  implicit none

  double precision, parameter :: initial_pressure_kpa = 294d0
  double precision, parameter :: initial_void_ratio = 0.808d0
  double precision, parameter :: cyclic_axial_amplitude = 0.0025d0
  integer, parameter :: steps_per_cycle = 400
  integer, parameter :: cycles = 8
  integer, parameter :: literature_cycles = 4
  double precision, parameter :: literature_q_amplitude_kpa = 114.2d0
  integer, parameter :: monotonic_steps = 1000
  double precision, parameter :: monotonic_final_axial_strain = -0.02d0

  write(*,'(A)') 'implementation,path_id,tolerance,step,cycle_time,axial_strain,' // &
       'p_effective_pa,q_signed_pa,ru,void_ratio,plastic_strain,alpha_norm,' // &
       'alpha_in_norm,alpha_memory_norm,memory_radius,fabric_norm,' // &
       'controller_target_q_pa,controller_residual_pa,status'
  call run_monotonic(1d-5)
  call run_cyclic(1d-5)
  call run_literature_q_controlled(1d-5)
  call run_monotonic(1d-7)
  call run_cyclic(1d-7)
  call run_literature_q_controlled(1d-7)

contains

  subroutine initialise_properties(props, tolerance)
    double precision, intent(out) :: props(25)
    double precision, intent(in) :: tolerance
    props = 0d0
    ! Liu et al. (2019), Table A1; cyclic memory parameters from Fig. A2b.
    props(1:16) = (/125d0, 0.05d0, 1.25d0, 0.712d0, 0.019d0, &
         0.934d0, 0.7d0, 0.01d0, 7.05d0, 0.968d0, 1.1d0, &
         0.704d0, 3.5d0, 45d0, 1d-5, 16.5d0/)
    ! emax/emin only bound void-ratio evolution; constant-volume paths retain
    ! e=0.808.  A high bounding-surface cap prevents an artificial cap here.
    props(17:25) = (/tolerance, tolerance, initial_void_ratio, 0.977d0, &
         0.597d0, 4d0, 2d0, 1d0, 0d0/)
  end subroutine initialise_properties

  subroutine initialise_point(tolerance, props, sig0, st0, swp0)
    double precision, intent(in) :: tolerance
    double precision, intent(out) :: props(25), sig0(6), st0(31), swp0
    double precision :: sig(6), st(31), deps(12), dmat(6,6), bulkw, swp
    double precision :: x, y, z, time0, dtime
    integer :: idtask, imod, isundr, istep, iter, iel, intpt
    integer :: ipl, nstat, nonsym, istrsdep, itimedep, itang, iabort

    call initialise_properties(props, tolerance)
    sig0 = 0d0
    sig0(1:3) = -initial_pressure_kpa
    st0 = 0d0
    deps = 0d0
    dmat = 0d0
    swp0 = 0d0
    x=0d0; y=0d0; z=0d0; time0=0d0; dtime=1d0
    imod=1; isundr=0; iter=1; iel=1; intpt=1
    idtask=1; istep=0
    call sanisandms(idtask,imod,isundr,istep,iter,iel,intpt,x,y,z, &
         time0,dtime,props,sig0,swp0,st0,deps,dmat,bulkw,sig,swp,st, &
         ipl,nstat,nonsym,istrsdep,itimedep,itang,iabort)
    if (iabort /= 0) error stop 'Pisano UDSM initialization aborted'
    st0 = st
  end subroutine initialise_point

  subroutine update_point(props, sig0, st0, swp0, total_strain, &
                          axial_increment, step, sig, st)
    double precision, intent(in) :: props(25), axial_increment
    double precision, intent(inout) :: sig0(6), st0(31), swp0
    double precision, intent(inout) :: total_strain(6)
    integer, intent(in) :: step
    double precision, intent(out) :: sig(6), st(31)
    double precision :: deps(12), dmat(6,6), bulkw, swp
    double precision :: x, y, z, time0, dtime
    integer :: idtask, imod, isundr, iter, iel, intpt
    integer :: ipl, nstat, nonsym, istrsdep, itimedep, itang, iabort

    deps = 0d0
    deps(1) = -0.5d0 * axial_increment
    deps(2) = -0.5d0 * axial_increment
    deps(3) = axial_increment
    total_strain = total_strain + deps(1:6)
    deps(7:12) = total_strain
    dmat = 0d0
    x=0d0; y=0d0; z=0d0; time0=0d0; dtime=1d0
    idtask=2; imod=1; isundr=0; iter=1; iel=1; intpt=1
    call sanisandms(idtask,imod,isundr,step,iter,iel,intpt,x,y,z,time0, &
         dtime,props,sig0,swp0,st0,deps,dmat,bulkw,sig,swp,st,ipl,nstat, &
         nonsym,istrsdep,itimedep,itang,iabort)
    if (iabort /= 0) error stop 'Pisano UDSM stress integration aborted'
    sig0 = sig
    st0 = st
    swp0 = swp
  end subroutine update_point

  subroutine emit(path_id, tolerance, step, cycle_time, axial_strain, sig, st, &
                  controller_target, controller_residual, status)
    character(len=*), intent(in) :: path_id
    double precision, intent(in) :: tolerance, cycle_time, axial_strain
    double precision, intent(in) :: sig(6), st(31)
    integer, intent(in) :: step
    double precision, intent(in) :: controller_target, controller_residual
    character(len=*), intent(in) :: status
    double precision :: p, q, ru, alpha_norm, alpha_in_norm
    double precision :: alpha_memory_norm, nan_value

    p = -(sig(1)+sig(2)+sig(3))/3d0
    q = sig(1)-sig(3)
    if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(q) .or. p <= 0d0) &
         error stop 'Pisano UDSM produced an invalid material-point state'
    ru = 1d0-p/initial_pressure_kpa
    alpha_norm = voigt_norm(st(1:6))
    alpha_memory_norm = voigt_norm(st(7:12))
    alpha_in_norm = voigt_norm(st(15:20))
    nan_value = ieee_value(0d0, ieee_quiet_nan)
    write(*,'(A,",",A,",",ES24.16,",",I0,14(",",ES24.16),",",A)') &
         'pisano_sanisand_ms_udsm', trim(path_id), tolerance, step, &
         cycle_time, axial_strain, p*1000d0, q*1000d0, ru, st(21), &
         st(23), alpha_norm, alpha_in_norm, alpha_memory_norm, st(31), &
         nan_value, controller_target, controller_residual, trim(status)
  end subroutine emit

  subroutine run_monotonic(tolerance)
    double precision, intent(in) :: tolerance
    double precision :: props(25), sig0(6), sig(6), st0(31), st(31), swp0
    double precision :: total_strain(6), axial_strain, increment
    integer :: step
    call initialise_point(tolerance, props, sig0, st0, swp0)
    total_strain = 0d0
    sig = sig0
    st = st0
    axial_strain = 0d0
    call emit('toyoura_monotonic_constant_volume', tolerance, 0, 0d0, &
         axial_strain, sig, st, ieee_value(0d0,ieee_quiet_nan), &
         ieee_value(0d0,ieee_quiet_nan),'ok')
    increment = monotonic_final_axial_strain / dble(monotonic_steps)
    do step=1,monotonic_steps
      call update_point(props,sig0,st0,swp0,total_strain,increment,step,sig,st)
      axial_strain = axial_strain + increment
      if (mod(step,5)==0) call emit('toyoura_monotonic_constant_volume', &
           tolerance,step,0d0,axial_strain,sig,st, &
           ieee_value(0d0,ieee_quiet_nan), &
           ieee_value(0d0,ieee_quiet_nan),'ok')
    enddo
  end subroutine run_monotonic

  subroutine run_cyclic(tolerance)
    double precision, intent(in) :: tolerance
    double precision :: props(25), sig0(6), sig(6), st0(31), st(31), swp0
    double precision :: total_strain(6), previous_target, target, increment
    double precision :: cycle_time
    integer :: step
    call initialise_point(tolerance, props, sig0, st0, swp0)
    total_strain = 0d0
    sig = sig0
    st = st0
    previous_target = 0d0
    call emit('toyoura_cyclic_constant_volume', tolerance,0,0d0,0d0,sig,st, &
         ieee_value(0d0,ieee_quiet_nan), &
         ieee_value(0d0,ieee_quiet_nan),'ok')
    do step=1,cycles*steps_per_cycle
      cycle_time = dble(step)/dble(steps_per_cycle)
      target = -cyclic_axial_amplitude*triangle(cycle_time)
      increment = target-previous_target
      previous_target = target
      call update_point(props,sig0,st0,swp0,total_strain,increment,step,sig,st)
      if (mod(step,5)==0) call emit('toyoura_cyclic_constant_volume', &
           tolerance,step,cycle_time,target,sig,st, &
           ieee_value(0d0,ieee_quiet_nan), &
           ieee_value(0d0,ieee_quiet_nan),'ok')
    enddo
  end subroutine run_cyclic

  subroutine evaluate_increment(props, sig_source, st_source, swp_source, &
       total_source, axial_increment, step, sig_trial, st_trial, swp_trial, &
       total_trial, iabort)
    double precision, intent(in) :: props(25), sig_source(6), st_source(31)
    double precision, intent(in) :: swp_source, total_source(6), axial_increment
    integer, intent(in) :: step
    double precision, intent(out) :: sig_trial(6), st_trial(31), swp_trial
    double precision, intent(out) :: total_trial(6)
    integer, intent(out) :: iabort
    double precision :: sig0(6), st0(31), deps(12), dmat(6,6), bulkw
    double precision :: x, y, z, time0, dtime
    integer :: idtask, imod, isundr, iter, iel, intpt
    integer :: ipl, nstat, nonsym, istrsdep, itimedep, itang

    sig0 = sig_source
    st0 = st_source
    swp_trial = swp_source
    deps = 0d0
    deps(1) = -0.5d0*axial_increment
    deps(2) = -0.5d0*axial_increment
    deps(3) = axial_increment
    total_trial = total_source+deps(1:6)
    deps(7:12) = total_trial
    dmat = 0d0
    x=0d0; y=0d0; z=0d0; time0=0d0; dtime=1d0
    idtask=2; imod=1; isundr=0; iter=1; iel=1; intpt=1
    call sanisandms(idtask,imod,isundr,step,iter,iel,intpt,x,y,z,time0, &
         dtime,props,sig0,swp_source,st0,deps,dmat,bulkw,sig_trial,swp_trial, &
         st_trial,ipl,nstat,nonsym,istrsdep,itimedep,itang,iabort)
  end subroutine evaluate_increment

  subroutine solve_q_target(props, sig_source, st_source, swp_source, &
       total_source, target_q, step, sig_best, st_best, swp_best, total_best, &
       axial_increment, residual, controller_ok)
    double precision, intent(in) :: props(25), sig_source(6), st_source(31)
    double precision, intent(in) :: swp_source, total_source(6), target_q
    integer, intent(in) :: step
    double precision, intent(out) :: sig_best(6), st_best(31), swp_best
    double precision, intent(out) :: total_best(6), axial_increment, residual
    logical, intent(out) :: controller_ok
    double precision :: sig_lower(6), sig_upper(6), sig_mid(6)
    double precision :: st_lower(31), st_upper(31), st_mid(31)
    double precision :: total_lower(6), total_upper(6), total_mid(6)
    double precision :: swp_lower, swp_upper, swp_mid
    double precision :: span, lower_inc, upper_inc, mid_inc
    double precision :: lower_res, upper_res, mid_res, tolerance
    integer :: expansion, iteration, abort_lower, abort_upper, abort_mid
    logical :: bracketed

    tolerance = 0.005d0*literature_q_amplitude_kpa
    span = 2d-5
    bracketed = .false.
    do expansion=0,11
      lower_inc = -span
      upper_inc = span
      call evaluate_increment(props,sig_source,st_source,swp_source, &
           total_source,lower_inc,step,sig_lower,st_lower,swp_lower, &
           total_lower,abort_lower)
      call evaluate_increment(props,sig_source,st_source,swp_source, &
           total_source,upper_inc,step,sig_upper,st_upper,swp_upper, &
           total_upper,abort_upper)
      if (abort_lower==0 .and. abort_upper==0) then
        lower_res = sig_lower(1)-sig_lower(3)-target_q
        upper_res = sig_upper(1)-sig_upper(3)-target_q
        if (lower_res==0d0 .or. upper_res==0d0 .or. &
             sign(1d0,lower_res)/=sign(1d0,upper_res)) then
          bracketed=.true.
          exit
        endif
      endif
      span=2d0*span
    enddo
    if (.not.bracketed) then
      sig_best=sig_source
      st_best=st_source
      swp_best=swp_source
      total_best=total_source
      axial_increment=0d0
      residual=sig_source(1)-sig_source(3)-target_q
      controller_ok=.false.
      return
    endif

    do iteration=1,60
      mid_inc=0.5d0*(lower_inc+upper_inc)
      call evaluate_increment(props,sig_source,st_source,swp_source, &
           total_source,mid_inc,step,sig_mid,st_mid,swp_mid,total_mid,abort_mid)
      if (abort_mid/=0) error stop 'Pisano UDSM q-controller trial aborted'
      mid_res=sig_mid(1)-sig_mid(3)-target_q
      if (abs(mid_res)<=tolerance) exit
      if (sign(1d0,mid_res)==sign(1d0,lower_res)) then
        lower_inc=mid_inc
        lower_res=mid_res
        sig_lower=sig_mid; st_lower=st_mid; swp_lower=swp_mid
        total_lower=total_mid
      else
        upper_inc=mid_inc
        upper_res=mid_res
        sig_upper=sig_mid; st_upper=st_mid; swp_upper=swp_mid
        total_upper=total_mid
      endif
    enddo
    sig_best=sig_mid; st_best=st_mid; swp_best=swp_mid; total_best=total_mid
    axial_increment=mid_inc
    residual=mid_res
    controller_ok=abs(residual)<=1.002d0*tolerance
  end subroutine solve_q_target

  subroutine run_literature_q_controlled(tolerance)
    double precision, intent(in) :: tolerance
    double precision :: props(25), sig0(6), sig(6), st0(31), st(31), swp0
    double precision :: total_strain(6), total_next(6), target_q, increment
    double precision :: residual, cycle_time, axial_strain, swp_next
    integer :: step
    logical :: controller_ok
    call initialise_point(tolerance,props,sig0,st0,swp0)
    total_strain=0d0
    sig=sig0; st=st0; axial_strain=0d0
    call emit('toyoura_literature_q_controlled',tolerance,0,0d0,0d0,sig,st, &
         0d0,0d0,'ok')
    do step=1,literature_cycles*steps_per_cycle
      cycle_time=dble(step)/dble(steps_per_cycle)
      target_q=literature_q_amplitude_kpa*triangle(cycle_time)
      call solve_q_target(props,sig0,st0,swp0,total_strain,target_q,step, &
           sig,st,swp_next,total_next,increment,residual,controller_ok)
      axial_strain=axial_strain+increment
      sig0=sig; st0=st; swp0=swp_next; total_strain=total_next
      if (controller_ok) then
        call emit('toyoura_literature_q_controlled',tolerance,step,cycle_time, &
             axial_strain,sig,st,target_q*1000d0,residual*1000d0,'ok')
      else
        call emit('toyoura_literature_q_controlled',tolerance,step,cycle_time, &
             axial_strain,sig,st,target_q*1000d0,residual*1000d0, &
             'controller_limit')
        exit
      endif
    enddo
  end subroutine run_literature_q_controlled

  double precision function triangle(cycle_time)
    double precision, intent(in) :: cycle_time
    double precision :: fraction
    fraction = modulo(cycle_time,1d0)
    if (fraction < 0.25d0) then
      triangle = 4d0*fraction
    else if (fraction < 0.75d0) then
      triangle = 2d0-4d0*fraction
    else
      triangle = -4d0+4d0*fraction
    endif
  end function triangle

  double precision function voigt_norm(tensor)
    double precision, intent(in) :: tensor(6)
    voigt_norm = sqrt(tensor(1)**2+tensor(2)**2+tensor(3)**2 + &
         2d0*(tensor(4)**2+tensor(5)**2+tensor(6)**2))
  end function voigt_norm

end program pisano_udsm_material_point_driver
