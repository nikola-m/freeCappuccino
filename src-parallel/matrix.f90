   module matrix_module
!
!  Functions for basic dense matrix calculations
!
   use types
 
   public

   contains

!=======================================================================
   function eye(n) 
!
!  Make identity matrix
!

!  Result
   real(dp), dimension(n,n) :: eye

!  Input
   integer, intent(in) :: n

   do i=1,n
      do j=1,n
         eye(i,j) = 0.0_dp
      enddo
      eye(i,i) = 1.0_dp
   enddo
   end function eye

!=======================================================================
   function matmult (a,b,n) result(c)
!
!  Compute matrix product c=ab
!
   implicit none

!  Result
   real(dp),dimension(n,n) :: c(n,n)

!  Input
   integer :: n   
   real(dp),dimension(n,n), intent(in) :: a,b

!  Locals
   integer :: i,j,k
   real(dp) :: x

      do i=1,n
         do j=1,n
            x=0.0_dp
            do k=1,n
             x=x+a(i,k)*b(k,j) 
            enddo
            c(i,j)=x
          enddo
      enddo
   end function matmult

!=======================================================================
   function rank_one_update(A,m,n,y,x,alpha) result(c)
!
!  rank-1 update A = A + alpha*y*xT
!
   implicit none
!  Result
   real(dp), dimension(m,n) :: c
!  Input
   real(dp), dimension(m,n), intent(in) :: A
   real(dp), dimension(m), intent(in) :: y
   real(dp), dimension(n), intent(in) :: x
   real(dp), intent(in) :: alpha
   integer, intent(in) :: m, n

!  Locals
   integer :: i,j

   do j = 1,n
     do i =1,m
       c(i,j) = a(i,j)+alpha*y(i)*x(j)
     enddo
   enddo

   end function rank_one_update

!=======================================================================
   function matrix_update(c,m,n,a,b,alpha) result(d)
!
!  Compute C  = C + alpha*A*B
!
   use types
   implicit none

!  Result
   real(dp),dimension(m,m) :: d
   

!  Input
   integer, intent(in) :: m,n   
   real(dp),dimension(m,m), intent(in) :: c
   real(dp),dimension(m,n), intent(in) :: a(m,n)
   real(dp),dimension(n,m), intent(in) :: b(n,m)
   real(dp), intent(in) :: alpha

!  Locals
   integer :: i,j,k
   real(dp) :: x

      do i=1,m
         do j=1,n
            x=0.0_dp
            do k=1,n
             x=x+a(i,k)*b(k,j) 
            enddo
            d(i,j) = c(i,j)+ alpha*x
          enddo
      enddo
   end function matrix_update
!=======================================================================
   function det(a,n) 
!  Result
   real(dp) :: det
!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a

   det = 0.0_dp
   if (n==3) then  
      det =  a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
            -a(1,3)*a(2,2)*a(3,1)-a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2) 
   else if (n==2) then
      det = a(1,1)*a(2,2)-a(1,2)*a(2,1)
   endif

   end function
!=======================================================================
   function inv(a)
!
!  An inverse of a 3x3 matrix.
!
!  Result
   real(dp), dimension(3,3) :: inv
!  Input
   real(dp), dimension(3,3), intent(in) :: a

     inv(1,1)= (a(2,2)*a(3,3) - a(2,3)*a(3,2)) /(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(1,2)= -(a(1,2)*a(3,3) - a(1,3)*a(3,2))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))  

     inv(1,3)= (a(1,2)*a(2,3) - a(1,3)*a(2,2))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2)  &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(2,1)= -(a(2,1)*a(3,3) - a(2,3)*a(3,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
    - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(2,2)=  (a(1,1)*a(3,3) - a(1,3)*a(3,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
    - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(2,3)= -(a(1,1)*a(2,3) - a(1,3)*a(2,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
    - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(3,1)=  (a(2,1)*a(3,2) - a(2,2)*a(3,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))

     inv(3,2)= -(a(1,1)*a(3,2) - a(1,2)*a(3,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))
     
     inv(3,3)=  (a(1,1)*a(2,2) - a(1,2)*a(2,1))/(a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) &
     - a(1,2)*a(2,1)*a(3,3) + a(1,2)*a(2,3)*a(3,1) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1))
 
   end function

!=======================================================================
   function solve(A,b,n) result(x)
!
!  Solve system with 3x3 system matrix
!

!  Result
   real(dp), dimension(n) :: x

!  Input
   real(dp), dimension(n,n), intent(in) :: A
   real(dp), dimension(n), intent(in) :: b
   integer, intent(in) :: n

!  Locals
   real(dp), dimension(n,n) :: invA

   invA = inv(A)
   x = matmul(invA,b)
   end function

!=======================================================================
   function trace(a)
!
!  Trace of a 3x3 matrix.
!

!  Result
   real(dp), dimension(3) :: trace
!  Input
   real(dp), dimension(3,3), intent(in) :: a

   trace = a(1,1)+a(2,2)+a(3,3)

   end function

!=======================================================================
   function solve_leastsq(A,b,m,n) result(x)
!
!  Solve system with m x n system matrix in least square sense (minimizing Euclidean norm).
!  System is overdetermined so we solve A'A * x = A'b, where A' is transpose of A.
!  The solution is given by x = (A'A)^(-1)A'*b
!

!  Result
   real(dp), dimension(n) :: x

!  Input
   integer, intent(in) :: n,m
   real(dp), dimension(m,n), intent(in) :: A
   real(dp), dimension(m), intent(in) :: b

!  Locals
   real(dp), dimension(n,m) ::  At
   real(dp), dimension(n,n) ::  AtA
   real(dp), dimension(n) ::  Atb

   At = transpose(A)      ! transpose of A : A'
   AtA = matmul(At,A)     ! left-hand side multiplication of A by A': A'A
   Atb = matmul(At,b)     ! A'b
   x = matmul(inv(AtA),Atb) ! x = (A'A)^(-1) * A'b

   end function

!=======================================================================
   function solve_leastsq_cholesky(A,b,m,n) result(x)
!
!  Solve system with m x n system matrix in least square sense (minimizing Euclidean norm).
!  System is overdetermined so we solve A'A * x = A'b, where A' is transpose of A.
!

!  Result
   real(dp), dimension(n) :: x

!  Input
   real(dp), dimension(m,n), intent(in) :: A
   real(dp), dimension(m), intent(in) :: b

!  Locals
   real(dp), dimension(n,m) ::  At
   real(dp), dimension(n,n) ::  AtA
   real(dp), dimension(n) ::  Atb

   At = transpose(A)    ! transpose of A : A'
   AtA = matmul(At,A)   ! left-hand side multiplication of A by A': A'A
   Atb = matmul(At,b)   ! A'b
   call solve_cholesky(AtA, n, Atb, x)

   end function

!=======================================================================
   function solve_leastsq_mgs(A,b,m,n) result(x)
!
!  Solve overdetermined linear system, using Modified Gram-Schmidt for QR decomposition
!

!  Result
   real(dp), dimension(n) :: x

!  Input
   real(dp), dimension(m,n), intent(in) :: A
   real(dp), dimension(m), intent(in) :: b

!  Locals
   real(dp), dimension(m,n) ::  Q
   real(dp), dimension(n,n) ::  R
   real(dp), dimension(n,m) ::  Qt
   real(dp), dimension(n) ::  Qtb

!  1. Decompose A=QR using Gram-Schmidt
   call mgs_qr(A, m, n, Q, R)

!  2. Multiply b with QT to get RHS of Rx=QTb   
   Qt = transpose(Q)
   Qtb = matmul(Qt,b)

!  3. Backsolve to get x=R^(-1)QTb
   call rsolve(R,n,Qtb,x) 

   end function

!=======================================================================
   function solve_lstsq_householder(A,b,m,n) result(x)
!
!  Solve overdetermined linear system, using Householder reflections for QR decomposition
!

!  Result
   real(dp), dimension(n) :: x

!  Input
   real(dp), dimension(m,n), intent(in) :: A
   real(dp), dimension(m), intent(in) :: b

!  Locals
   real(dp), dimension(m,m) ::  Q
   real(dp), dimension(m,n) ::  R
   real(dp), dimension(m,m) ::  Qt
   real(dp), dimension(m) ::  Qtb

!  1. Decompose A=QR using Householder
   call householder_qr(A, m, n, Q, R)

!  2. Multiply b with QT to get RHS of Rx=QTb   
   Qt = transpose(Q)
   Qtb = matmul(Qt,b)

!  3. Backsolve to get x=R**-1*QTb
   call rsolve(R,n,Qtb,x) 

   end function

!=======================================================================
   SUBROUTINE householder_qr(a,m,n,q,r)
!
!  QR decomposition via Householder reflections for a generic m x n matrix.
!
   use types
   implicit none
   real(dp), dimension(m,n), intent(in) :: a
   real(dp), dimension(m,n), intent(inout) :: R
   real(dp), dimension(m,m), intent(inout) :: Q
   integer :: m,n
!  Locals
   real(dp), dimension(m) :: w
   real(dp), dimension(n) :: u
   real(dp), dimension(m,m) :: WWT
   real(dp) :: s,g
   integer :: i,k
   
   R = a
   Q = eye(m)
   DO k = 1,N
!  Uzmi sve elemente u koloni koji se nalaze ispod glavne dijagonale:
   w = 0.0_dp
   DO i=k,m
     w(i) = R(i,k) 
   END DO
!  Nadji 2 normu tog vektora:
   s = sqrt(sum(w*w))!(dot_product(w,w))  
!  Dodaj tu vrednost prvom elementu u w-vektoru:
   w(k) = w(k) - s  
   g = sqrt(dot_product(w,w))
   DO i=k,m
   w(i) = w(i)/g
   END DO
   u = 2.0_dp*matmul(w,R)
   R = rank_one_update(R,m,n,W,u,-1.0_dp)    ! R = R – W*uT 
!  Sad pravimo Q:
   WWT = 0.0_dp                         ! Initialize WWT
   WWT = rank_one_update(WWT,m,m,W,W,1.0_dp) ! WWT = WWT + W*WT
   Q = matrix_update(Q,m,n,Q,WWT,-2.0_dp)  ! Q = Q – 2*Q*WWT)
   END DO
   RETURN
   END SUBROUTINE householder_qr

!=======================================================================
   subroutine mgs_qr(a, m, n, q, r)
!
!  QR decomposition via modified Gram-Schmidt procedure for a generic m x n matrix.
!
   use types
   implicit none

!  Result
   real(dp), dimension(m,n), intent(inout) :: q  ! Orthogonal matrix
   real(dp), dimension(n,n), intent(inout) :: r  ! Upper triangular matrix

!  Input
   integer, intent(in) :: m, n               ! matrix type m x n
   real(dp), dimension(m,n), intent(in) :: a ! matrix

!  Locals
   integer i,j,k
   real(dp) :: z


!  Initialize Q and R
   q = 0.0_dp
   r = 0.0_dp

!  Copy A to Q
   q = a

   do j=1,n ! For all columns of A

     z = 0.0_dp
     do i=1,m
     z = z + q(i,j)**2
     enddo
     r(j,j) = sqrt(z) 
     
     do i=1,m
     q(i,j) = q(i,j)/r(j,j)
     enddo

     do k = j+1,n !........................
       z = 0.0_dp                         !
       do i=1,m                           !
        z = z + q(i,j)*q(i,k)             !
       enddo                              ! 
       r(j,k)  = z                        !
       do i=1,m                           !
       q(i,k) = q(i,k) - r(j,k)*q(i,j)    !
       enddo                              !
     enddo !..............................!

   enddo

   return
   end subroutine mgs_qr

!=======================================================================
   function pnorm(w,p)
!
!  Calculates p-norm of a vector (2norm moze kao: sqrt(dot_product(w,w))
!

!  Result
   real(dp) :: pnorm

!  Input
   integer :: p 
   real(dp),dimension(:), intent(in) :: w

!  Locals
   integer :: i, n
   real(dp) :: pr

   n = size(w)
   pr = 1/p

   pnorm = 0.0_dp
   do i=1,n
     pnorm = pnorm + w(i)**p
   enddo
   pnorm = pnorm**(pr)

   end function

!=======================================================================
   subroutine rsolve(a, n, b, x)
!
!  Solve linear system R x = b, where R is upper triangular matrix.
!
   use types
   implicit none

!  Input
   integer :: n
   real(dp), dimension(n,n), intent(in) :: a ! system matrix
   real(dp), dimension(n), intent(in) :: b   ! rhs vector 
   real(dp), dimension(n), intent(out) :: x  ! result

!  Local
   integer :: i, j
   real(dp) :: suma

   do i= n,1,-1
     suma = b(i) 
     do j=i+1,n
       suma = suma - a(i,j)*x(j)
     enddo
     x(i) = suma/a(i,i)
   end do
 
   return
   end subroutine rsolve

!=======================================================================
   subroutine cholesky(a, n, l)
!
!  Cholesky decomposition of a symmetric positive definite matrix
!
   use types
   implicit none

!  Result
   real(dp), dimension(n,n), intent(out) :: l

!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a

!  Locals
   integer :: i, k, m
   real(dp) :: sum1, sum2

   do k=1,n
     sum1 = a(k,k)
     do m=1,k-1
        sum1 = sum1 - l(k,m)**2
     enddo
     l(k,k) = sqrt(sum1)
     do i=k+1,n
        sum2 = a(i,k)
        do m=1,k-1
           sum2 = sum2 - l(i,m)*l(k,m)
        enddo
        l(i,k) = sum2/l(k,k)
        l(k,i) = l(i,k)  ! making LT
     enddo
   enddo
   return
   end subroutine cholesky

!=======================================================================
   function solve_spd_cholesky(a, n, b) result(x)
!
!  Solves Ax=b, where A is symmetric pos-def, by using Cholesky decomposition. 
!  
   use types
   implicit none

!  Result
   real(dp), dimension(n) :: x
   
!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a
   real(dp), dimension(n), intent(in) :: b

!  Locals
   integer :: i, k
   real(dp) :: suma
   real(dp), dimension(n,n) :: l

   call cholesky(a, n, l)

!  Solve L*y=b storing y in x
   do i=1,n
     suma = b(i)
     do k=i-1,1,-1
       suma = suma-l(i,k)*x(k)
     enddo
     x(i) = suma/l(i,i)
   enddo

!  Solve LT*x = y
   do i=n,1,-1
     suma = x(i)
     do k=i+1,n
       suma = suma-l(k,i)*x(k)
     enddo
     x(i) = suma/l(i,i)
   enddo
 
   return
   end function solve_spd_cholesky

!=======================================================================
   function solve_already_choleskyed(l, n, b) result(x)
!
!  Solves L*L^Tx=b, where L is lower triangular with sq-rooted main diagonal,
!  and is found by using Cholesky decomposition. 
!  NOTE: "allready_cholesked" is a tongue-in-cheek way to say that Cholesky 
!  decomposition has already been perfomed on that system matrix.
!  
   use types
   implicit none

!  Result
   real(dp), dimension(n) :: x
   
!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: l
   real(dp), dimension(n), intent(in) :: b

!  Locals
   integer :: i, k
   real(dp) :: suma

!  Solve L*y=b storing y in x
   do i=1,n
     suma = b(i)
     do k=i-1,1,-1
       suma = suma-l(i,k)*x(k)
     enddo
     x(i) = suma/l(i,i) !sqrt(l(i,i))
   enddo

!  Solve LT*x = y
   do i=n,1,-1
     suma = x(i)
     do k=i+1,n
       suma = suma-l(k,i)*x(k)
     enddo
     x(i) = suma/l(i,i) !sqrt(l(i,i))
   enddo
 
   return
   end function solve_already_choleskyed

!=======================================================================
   subroutine solve_cholesky(a, n, b, x)
!
!  Solves Ax=b, where A is pos-def, by using Cholesky decomposition. 
!  
   use types
   implicit none

!  Result
   real(dp), dimension(n), intent(out) :: x
   
!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a
   real(dp), dimension(n), intent(in) :: b

!  Locals
   integer :: i, k
   real(dp) :: suma
   real(dp), dimension(n,n) :: l

   call cholesky(a, n, l)

!  Solve L*y=b storing y in x
   do i=1,n
     suma = b(i)
     do k=i-1,1,-1
       suma = suma-l(i,k)*x(k)
     enddo
     x(i) = suma/l(i,i)
   enddo

!  Solve LT*x = y
   do i=n,1,-1
     suma = x(i)
     do k=i+1,n
       suma = suma-l(k,i)*x(k)
     enddo
     x(i) = suma/l(i,i)
   enddo
 
   return
   end subroutine solve_cholesky

!=======================================================================
   function frobenius_product(a,b,n) result(c)
!
!  Returns matrix C - a result of Frobenius inner product 
!  C= A:B or in component form Cij = AijBij
!
   
!  Result
   real(dp), dimension(n,n) :: c

!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a,b

!  Locals
   integer :: i,j

   do i=1,n
     do j=1,n
       c(i,j) = a(i,j)*b(i,j)
     enddo
   enddo

   end function frobenius_product

!=======================================================================
   subroutine eig(a,n,lam,v) !result(lam,v)
!
!  Finds eigenvalues of symmetric matrix using QR decomposition iteratively 
!  The eignevalues are stored in array 'lam'
!
   use types
   implicit none
!  Result
   real(dp), dimension(n,n) :: lam
   real(dp), dimension(n,n) :: v

!  Input
   integer, intent(in) :: n
   real(dp), dimension(n,n), intent(in) :: a

!  Locals
   integer :: i
   real(dp), dimension(n,n) :: q,r

   lam=a
   call householder_qr(lam,n,n,q,r)
   v=q
   lam = matmul(r,q)
   do i=1,11
     call householder_qr(lam,n,n,q,r)
     v = matmul(v,q)
     lam = matmul(r,q)
!     if (maxval(lam)<1e-1) exit
   end do
!   print*,'eig converged in', i ,'iterations.'
   end subroutine eig

!=======================================================================
   function norm2_offdiag(a,n) result(norm)
!
!  Calculates p-norm of a vector (2norm moze kao: sqrt(dot_product(w,w))
!

!  Result
   real(dp) :: norm

!  Input
   integer, intent(in) :: n
   real(dp),dimension(n,n), intent(in) :: a

!  Locals
   integer :: i, j

   norm = 0.0_dp
   do i=1,n
     do j=1,n
     if (j==i) cycle
     norm = norm + a(i,j)*a(i,j)
     end do
   enddo
   norm = sqrt(norm)

   end function

   end module matrix_module


