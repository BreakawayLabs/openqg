module units

  implicit none

  private

  public m3_to_Sv
  public m_to_km

  interface m_to_km
     module procedure m_to_km
     module procedure m_to_km_1d
  end interface m_to_km

contains

  double precision function m_to_km(m)
    double precision :: m

    m_to_km = 1.0d-3*m

  end function m_to_km

  function m_to_km_1d(m)
    double precision :: m(:)
    double precision :: m_to_km_1d(size(m))

    m_to_km_1d = 1.0d-3*m

  end function m_to_km_1d
    
  double precision function m3_to_Sv(m3)
    double precision :: m3

    m3_to_Sv = 1.0d-6*m3

  end function m3_to_Sv

end module units
