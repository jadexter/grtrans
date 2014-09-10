  program test_interpemis

    use emissivity, only: interpemis

    implicit none

    double precision, dimension(:), allocatable :: nu, freqarr, jnu, anu
    double precision, dimension(:,:), allocatable :: K
    integer :: nr, nf, i

    nr=5; nf=2

    allocate(nu(nr)); allocate(freqarr(nf))
    allocate(jnu(nr*nf)); allocate(anu(nr*nf))
    allocate(K(nr,11))

    nu=(/(1.,i=1,nr)/)
    freqarr=(/.8,1.2/)
    jnu=(/0.2, 0.4, 0.8, 0.9, 1.0, 0.05, 0.08, 0.12, 0.15, 0.17/)
    anu=0.

    call interpemis(nu,freqarr,jnu,anu,K)

    write(6,*) 'K1: ',K(:,1)
    write(6,*) 'K2: ',K(:,5)

  end program
