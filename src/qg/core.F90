module core

  use qg, only: qg_type
  use vorticity, only: vort_step
  use pressure, only: solve_pressure
  use vorsubs, only: ocqbdy

  implicit none

  private

  public step_qg

contains

  subroutine step_qg(tdt, ent, wekp, txis, txin, qg)

    type(qg_type), intent(inout) :: qg
    double precision, intent(in) :: tdt
    double precision, intent(in) :: ent(qg%b%nxp,qg%b%nyp,qg%b%nl)
    double precision, intent(in) :: wekp(qg%b%nxp,qg%b%nyp)
    double precision, intent(in) :: txis, txin

    ! Step qg model
    call vort_step(qg, tdt, ent, wekp, txis, txin)

    ! Invert pv into pressures
    call solve_pressure(qg, qg%b, qg%topo, tdt, ent)

    ! Compute pv on boundaries (mixed condition)
    call ocqbdy (qg%p, qg%mod%amat, qg%b, qg%bcco, qg%topo, qg%q)

  end subroutine step_qg

end module core
