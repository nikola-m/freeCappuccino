subroutine correct_turbulence_inlet()
  use parameters, only: TurbModel
  use k_epsilon_std
  implicit none

  select case (TurbModel)
    case (1)
      call correct_turbulence_inlet_k_epsilon_std()
    case default
      call correct_turbulence_inlet_k_epsilon_std()
  end select

end subroutine