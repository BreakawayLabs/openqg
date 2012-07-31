module testlib

  implicit none

  private

  public start_suite
  public end_suite
  public start_test
  public end_test
  public check_threshold

contains

  subroutine start_suite(suite_name)
    character (len=*), intent(in) :: suite_name
    print *, "##teamcity[testSuiteStarted name='", suite_name, "']"
  end subroutine start_suite

  subroutine end_suite(suite_name)
    character (len=*), intent(in) :: suite_name
    print *, "##teamcity[testSuiteFinished name='", suite_name, "']"
  end subroutine end_suite

  subroutine start_test(test_name)
    character (len=*), intent(in) :: test_name
    print *, "##teamcity[testStarted name='", test_name, "' captureStandardOutput='true']"
  end subroutine start_test

  subroutine end_test(test_name)
    character (len=*), intent(in) :: test_name
    print *, "##teamcity[testFinished name='", test_name, "']"
  end subroutine end_test

  subroutine check_threshold(result, threshold, test_name)
    double precision, intent(in) :: result
    double precision, intent(in) :: threshold
    character (len=*), intent(in) :: test_name

    if (result > threshold) then
       print *, "##teamcity[testFailed type='comparisonFailure' name='", test_name, &
            "' message='result > threshold' expected='", threshold , &
            "' actual='", result , "']]"
    end if
  end subroutine check_threshold


end module testlib
