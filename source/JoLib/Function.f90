
  function fgauss(x,sigw)
    implicit none
    real fgauss,x,sigw
    fgauss=exp(-x**2/2/sigw**2)
  end function fgauss

