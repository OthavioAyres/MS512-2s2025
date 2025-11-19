      program poisson
      implicit none
      integer, parameter :: m = 21 ! numero de nos
      integer :: i,im,ip,j,jm,jp,k
      integer :: n,nitmax,iconv
      real(8) :: x,y,f,h,erro,eps
      real(8), dimension(m,m) :: u,uk,b,r,p,g
!
      nitmax=1000000
      eps=1.d-8
!
      h=1.d0/(m-1)
!
      u=0.d0
      b=0.d0
      r=0.d0
!
!     termo de contorno
!
      call gbc(g,m,h)
!
!     u inicial
!
      u=g
!
      do k=1,nitmax
!
!     Jacobi
!
      call jacobi(b,uk,u,m,h)
!
!     Checa a convergencia
!
      call convcheck(b,u,uk,r,m,h,eps,erro,iconv)
!
!
      if(iconv.eq.1) then
      write(*,"('Convergiu com ',i8,' iteracoes, Erro=',e15.10)")  k,erro
      exit
      end if
!
      end do ! k
! !
! !     imprime a solucao, descomente abaixo
! !
!       do j=1,m
!       y=h*(j-1)
!       do i=1,m
!       x=h*(i-1)
!       write(100,"(2(f13.6,2x),e15.8)") x,y,u(i,j)
!       end do ! i
!       end do ! j
! !
! !     para visualizar: gnuplot 'g-sol.plt'
! !
      end program
!
!-------------------------------------------------------
!
      subroutine gbc(g,m,h)
      implicit none
!
      integer :: m
      real(8) :: h
      real(8), dimension(m,m) :: g
!
      integer :: i,im,ip,j,jm,jp
      real(8) :: x,y,f
!     
      g=0.d0
!
      do j=1,m
      y=h*(j-1)
      do i=1,m
      x=h*(i-1)
!
      if(i.eq.1) g(i,j)=0.d0
      if(i.eq.m) g(i,j)=y
      if(j.eq.1) g(i,j)=(x-1.d0)*dsin(x)
      if(j.eq.m) g(i,j)=x*(2-x)
!
      end do ! i
      end do ! j
!
!      g=0.d0
!
      end subroutine
!
!-------------------------------------------------------
!
      subroutine jacobi(b,uk,u,m,h)
      implicit none
!
      integer :: m
      real(8) :: h
      real(8), dimension(m,m) :: b,u,uk
!
      integer :: i,im,ip,j,jm,jp
      real(8) :: x,y,f
!
      uk=u
!
      do j=2,m-1
      jm=j-1
      jp=j+1
      y=h*(j-1)
      do i=2,m-1
      im=i-1
      ip=i+1
      x=h*(i-1)
!
!     jacobi
!
      u(i,j)=(uk(im,j)+uk(ip,j)+uk(i,jm)+uk(i,jp)+h**2.d0*b(i,j))/4.d0
!
      end do ! i
      end do ! j
!
      end subroutine
!
!-------------------------------------------------------
!
      subroutine convcheck(b,u,uk,r,m,h,eps,erro,iconv)
      implicit none
!
      integer :: m,iconv
      real(8) :: h,eps
      real(8), dimension(m,m) :: b,u,uk,r
!
      real(8) :: erro,eu,ed
      integer :: i,im,ip,j,jm,jp
      real(8) :: x,y,f
!
      iconv=0
!
!     checa convergencia  
!
      erro=0.d0
      eu=0.d0
      ed=0.d0
!
      do j=1,m
      do i=1,m
!
      eu=eu+(u(i,j))**2.d0
      ed=ed+(u(i,j)-uk(i,j))**2.d0
!
      end do ! i
      end do ! j
!
!     usando o erro relativo
!
      erro=dsqrt(ed/eu)
!
      if(erro.lt.eps) iconv=1
!
      end subroutine
!