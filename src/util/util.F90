module util

  implicit none

  private

  public streq
  public strlen
  public s2s

contains

  integer pure function strlen(s)
    character, allocatable, intent(in) :: s(:)
    integer :: i
    strlen = -1
    do i=1,size(s)
       if ((ichar(s(i)) == 0 .or. ichar(s(i)) == 32) .and. strlen == -1) then
          strlen = i-1
       endif
    enddo

  end function strlen

  logical pure function streq(s1, s2)
    character, allocatable, intent(in) :: s1(:)
    character (len=*), intent(in) :: s2
    
    streq = s2s(s1) == s2

  end function streq

  character (len=strlen(s)) pure function s2s(s)

    character, allocatable, intent(in) :: s(:)

    s2s = transfer(s(:), s2s)

  end function s2s

end module util
