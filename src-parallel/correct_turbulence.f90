subroutine correct_turbulence()
  use parameters, only: TurbModel
  use k_epsilon_std
  use k_epsilon_rng
  use k_omega_sst
  use k_eqn_eddy
  use spalart_allmaras
  implicit none

  ! Necessary for almost any turb model are fresh velocity gradients,
  ! and strain magnitude calculated here:
  call find_strain_rate

  select case (TurbModel)

    case (1)
      call correct_turbulence_k_epsilon_std()

    case (2)
      call correct_turbulence_k_epsilon_rng()

    case (3)
      LowRe = .false.
      call correct_turbulence_k_omega_sst()

    case (4)
      LowRe = .true.
      call correct_turbulence_k_omega_sst()

    case (5)
      call correct_turbulence_spalart_allmaras()

    case (6)
      call correct_turbulence_k_eqn_eddy()

    case default
      
  end select

end subroutine

