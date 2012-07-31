program unit_tests

  use test_linalg_mod, only: test_linalg
  use test_inhomog_mod, only: test_inhomog

  implicit none

  call main()

contains

  subroutine main()

    call test_inhomog()
    call test_linalg()

  end subroutine main
  
end program unit_tests
