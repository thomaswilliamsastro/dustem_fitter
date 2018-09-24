subroutine covariance_matrix(a,b,c,chisq,m,n)

   implicit none

   integer, intent(in) :: m
   integer, intent(in) :: n

   real(8), dimension(n,m), intent(in) :: a
   real(8), dimension(m,m), intent(in) :: b
   real(8), dimension(m,n), intent(in) :: c
   real(8), dimension(m,m) :: b_inv

   real(8), dimension(m) :: work
   real(8), dimension(m) ::ipiv

   integer :: info

   real(8), dimension(n,m) :: term_1
   real(8), dimension(n,n), intent(out) :: chisq

   external DGETRF
   external DGETRI

   b_inv = b

   call DGETRF(m,m,b_inv,m,ipiv,info)
   call DGETRI(m,b_inv,m,ipiv,work,m,info)

   term_1 = matmul(a,b_inv)
   chisq = matmul(term_1,c)

end subroutine covariance_matrix

subroutine trapz(y,x,integrand,m)

   implicit none

   integer, intent(in) :: m

   real(8), dimension(m), intent(in) :: y
   real(8), dimension(m), intent(in) :: x
   real(8), intent(out) :: integrand

   associate(n => size(x))
      integrand = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
   end associate

end subroutine trapz
