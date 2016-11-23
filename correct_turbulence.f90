subroutine correct_turbulence()
  use parameters, only: TurbModel
  use k_epsilon_std
  implicit none

  ! Necessary for almost any turb model are fresh velocity gradients,
  ! and strain magnitude calculated here:
  call find_strain_rate

  select case (TurbModel)
    case (1)
      call correct_turbulence_k_epsilon_std()
    case default
      call correct_turbulence_k_epsilon_std()
  end select

end subroutine