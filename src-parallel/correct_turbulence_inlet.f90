subroutine correct_turbulence_inlet()
  use parameters, only: TurbModel
  use k_epsilon_std
  use k_epsilon_rng
  use k_omega_sst
  use k_eqn_eddy
  use spalart_allmaras
  implicit none

  select case (TurbModel)
    case (1)
      call correct_turbulence_inlet_k_epsilon_std()
    case (2)
      call correct_turbulence_inlet_k_epsilon_rng()
    case (3)
      call correct_turbulence_inlet_k_omega_sst()
    case (4)
      call correct_turbulence_inlet_k_omega_sst()
    case (5)
      call correct_turbulence_inlet_spalart_allmaras()
    case (6)
      call correct_turbulence_inlet_k_eqn_eddy()
    case default
  end select

end subroutine