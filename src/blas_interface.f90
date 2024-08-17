module blas_types
  use iso_fortran_env
  implicit none
  integer, parameter :: sp = real32        ! Single precision
  integer, parameter :: dp = real64        ! Double precision
  integer, parameter :: cp = real32        ! Complex single precision
  integer, parameter :: zp = real64        ! Complex double precision
  integer, parameter :: blasint = int32    ! BLAS integer type
end module blas_types

module blas_level1
  use blas_types
  implicit none
  private

  ! Public interface for BLAS Level-1 functions
  public :: srotg, srotmg, srot, srotm, sswap, sscal, scopy, saxpy, sdot, snrm2, sasum, isamax
  public :: drotg, drotmg, drot, drotm, dswap, dscal, dcopy, daxpy, ddot, dnrm2, dasum, idamax
  public :: dsdot, sdsdot
  public :: cswap, cscal, csscal, ccopy, caxpy, cdotu, cdotc, cnrm2, casum, icamax
  public :: zswap, zscal, zdscal, zcopy, zaxpy, zdotu, zdotc, znrm2, zasum, izamax

  interface
      ! Givens rotation generation
      subroutine srotg(sa, sb, c, s)
          import :: sp
          real(sp), intent(inout) :: sa, sb
          real(sp), intent(out) :: c, s
      end subroutine srotg

      subroutine drotg(da, db, c, s)
          import :: dp
          real(dp), intent(inout) :: da, db
          real(dp), intent(out) :: c, s
      end subroutine drotg

      ! Modified Givens rotation generation
      subroutine srotmg(d1, d2, x1, y1, param)
          import :: sp, blasint
          real(sp), intent(inout) :: d1, d2, x1
          real(sp), intent(in) :: y1
          real(sp), intent(out) :: param(5)
      end subroutine srotmg

      subroutine drotmg(d1, d2, x1, y1, param)
          import :: dp, blasint
          real(dp), intent(inout) :: d1, d2, x1
          real(dp), intent(in) :: y1
          real(dp), intent(out) :: param(5)
      end subroutine drotmg

      ! Givens rotation
      subroutine srot(n, dx, incx, dy, incy, c, s)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: c, s
          real(sp), intent(inout) :: dx(*), dy(*)
      end subroutine srot

      subroutine drot(n, dx, incx, dy, incy, c, s)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: c, s
          real(dp), intent(inout) :: dx(*), dy(*)
      end subroutine drot

      ! Modified Givens rotation
      subroutine srotm(n, dx, incx, dy, incy, param)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: param(5)
          real(sp), intent(inout) :: dx(*), dy(*)
      end subroutine srotm

      subroutine drotm(n, dx, incx, dy, incy, param)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: param(5)
          real(dp), intent(inout) :: dx(*), dy(*)
      end subroutine drotm

      ! swap (x <-> y)
      subroutine sswap(n, dx, incx, dy, incy)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(inout) :: dx(*), dy(*)
      end subroutine sswap

      subroutine dswap(n, dx, incx, dy, incy)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(inout) :: dx(*), dy(*)
      end subroutine dswap

      subroutine cswap(n, cx, incx, cy, incy)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(inout) :: cx(*), cy(*)
      end subroutine cswap

      subroutine zswap(n, cx, incx, cy, incy)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(inout) :: cx(*), cy(*)
      end subroutine zswap

      ! scale (x = alpha*x)
      subroutine sscal(n, sa, dx, incx)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: sa
          real(sp), intent(inout) :: dx(*)
      end subroutine sscal

      subroutine dscal(n, da, dx, incx)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: da
          real(dp), intent(inout) :: dx(*)
      end subroutine dscal

      subroutine cscal(n, ca, cx, incx)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: ca
          complex(cp), intent(inout) :: cx(*)
      end subroutine cscal

      subroutine zscal(n, ca, cx, incx)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: ca
          complex(zp), intent(inout) :: cx(*)
      end subroutine zscal

      subroutine csscal(n, sa, cx, incx)
          import :: sp, cp, blasint
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: sa
          complex(cp), intent(inout) :: cx(*)
      end subroutine csscal

      subroutine zdscal(n, da, cx, incx)
          import :: dp, zp, blasint
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: da
          complex(zp), intent(inout) :: cx(*)
      end subroutine zdscal

      ! copy (y = x)
      subroutine scopy(n, dx, incx, dy, incy)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: dx(*)
          real(sp), intent(out) :: dy(*)
      end subroutine scopy

      subroutine dcopy(n, dx, incx, dy, incy)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: dx(*)
          real(dp), intent(out) :: dy(*)
      end subroutine dcopy

      subroutine ccopy(n, cx, incx, cy, incy)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: cx(*)
          complex(cp), intent(out) :: cy(*)
      end subroutine ccopy

      subroutine zcopy(n, cx, incx, cy, incy)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: cx(*)
          complex(zp), intent(out) :: cy(*)
      end subroutine zcopy

      ! axpy (y = alpha*x + y)
      subroutine saxpy(n, sa, dx, incx, dy, incy)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: sa, dx(*)
          real(sp), intent(inout) :: dy(*)
      end subroutine saxpy

      subroutine daxpy(n, da, dx, incx, dy, incy)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: da, dx(*)
          real(dp), intent(inout) :: dy(*)
      end subroutine daxpy

      subroutine caxpy(n, ca, cx, incx, cy, incy)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: ca, cx(*)
          complex(cp), intent(inout) :: cy(*)
      end subroutine caxpy

      subroutine zaxpy(n, ca, cx, incx, cy, incy)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: ca, cx(*)
          complex(zp), intent(inout) :: cy(*)
      end subroutine zaxpy

      ! Dot product
      function sdot(n, dx, incx, dy, incy) result(dot)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: dx(*), dy(*)
          real(sp) :: dot
      end function sdot

      function ddot(n, dx, incx, dy, incy) result(dot)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: dx(*), dy(*)
          real(dp) :: dot
      end function ddot

      function dsdot(n, sx, incx, sy, incy) result(dot)
          import :: sp, dp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: sx(*), sy(*)
          real(dp) :: dot
      end function dsdot

      function sdsdot(n, sb, sx, incx, sy, incy) result(dot)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: sb, sx(*), sy(*)
          real(sp) :: dot
      end function sdsdot

      function cdotu(n, cx, incx, cy, incy) result(dot)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: cx(*), cy(*)
          complex(cp) :: dot
      end function cdotu

      function zdotu(n, cx, incx, cy, incy) result(dot)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: cx(*), cy(*)
          complex(zp) :: dot
      end function zdotu

      function cdotc(n, cx, incx, cy, incy) result(dot)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: cx(*), cy(*)
          complex(cp) :: dot
      end function cdotc

      function zdotc(n, cx, incx, cy, incy) result(dot)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: cx(*), cy(*)
          complex(zp) :: dot
      end function zdotc

      ! Euclidean norm
      function snrm2(n, x, incx) result(norm)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: x(*)
          real(sp) :: norm
      end function snrm2

      function dnrm2(n, x, incx) result(norm)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: x(*)
          real(dp) :: norm
      end function dnrm2

      function cnrm2(n, x, incx) result(norm)
          import :: sp, cp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: x(*)
          real(sp) :: norm
      end function cnrm2

      function znrm2(n, x, incx) result(norm)
          import :: dp, zp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: x(*)
          real(dp) :: norm
      end function znrm2

      ! Sum of absolute values
      function sasum(n, x, incx) result(sum)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: x(*)
          real(sp) :: sum
      end function sasum

      function dasum(n, x, incx) result(sum)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: x(*)
          real(dp) :: sum
      end function dasum

      function casum(n, x, incx) result(sum)
          import :: sp, cp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: x(*)
          real(sp) :: sum
      end function casum

      function zasum(n, x, incx) result(sum)
          import :: dp, zp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: x(*)
          real(dp) :: sum
      end function zasum

      ! Index of max absolute value
      function isamax(n, dx, incx) result(index)
          import :: sp, blasint
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: dx(*)
          integer(blasint) :: index
      end function isamax

      function idamax(n, dx, incx) result(index)
          import :: dp, blasint
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: dx(*)
          integer(blasint) :: index
      end function idamax

      function icamax(n, cx, incx) result(index)
          import :: cp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: cx(*)
          integer(blasint) :: index
      end function icamax

      function izamax(n, cx, incx) result(index)
          import :: zp, blasint
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: cx(*)
          integer(blasint) :: index
      end function izamax
  end interface
end module blas_level1

module blas_level2
  use blas_types
  implicit none
  private

  ! Public interface for BLAS Level-2 functions
  public :: sgemv, sgbmv, ssymv, ssbmv, sspmv, strmv, stbmv, stpmv, strsv, stbsv, stpsv, sger, ssyr, sspr, ssyr2, sspr2
  public :: dgemv, dgbmv, dsymv, dsbmv, dspmv, dtrmv, dtbmv, dtpmv, dtrsv, dtbsv, dtpsv, dger, dsyr, dspr, dsyr2, dspr2
  public :: cgemv, cgbmv, chemv, chbmv, chpmv, ctrmv, ctbmv, ctpmv, ctrsv, ctbsv, ctpsv, cgeru, cgerc, cher, chpr, cher2, chpr2
  public :: zgemv, zgbmv, zhemv, zhbmv, zhpmv, ztrmv, ztbmv, ztpmv, ztrsv, ztbsv, ztpsv, zgeru, zgerc, zher, zhpr, zher2, zhpr2

  interface
      ! General matrix-vector multiplication
      subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: sp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, lda, incx, incy
          real(sp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(sp), intent(inout) :: y(*)
      end subroutine sgemv

      subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: dp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, lda, incx, incy
          real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(dp), intent(inout) :: y(*)
      end subroutine dgemv

      subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: cp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, lda, incx, incy
          complex(cp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(cp), intent(inout) :: y(*)
      end subroutine cgemv

      subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: zp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, lda, incx, incy
          complex(zp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(zp), intent(inout) :: y(*)
      end subroutine zgemv

      ! General band matrix-vector multiplication
      subroutine sgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
          import :: sp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, kl, ku, lda, incx, incy
          real(sp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(sp), intent(inout) :: y(*)
      end subroutine sgbmv

      subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
          import :: dp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, kl, ku, lda, incx, incy
          real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(dp), intent(inout) :: y(*)
      end subroutine dgbmv

      subroutine cgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
          import :: cp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, kl, ku, lda, incx, incy
          complex(cp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(cp), intent(inout) :: y(*)
      end subroutine cgbmv

      subroutine zgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
          import :: zp, blasint
          character, intent(in) :: trans
          integer(blasint), intent(in) :: m, n, kl, ku, lda, incx, incy
          complex(zp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(zp), intent(inout) :: y(*)
      end subroutine zgbmv

      ! Symmetric/Hermitian matrix-vector multiplication
      subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, lda, incx, incy
          real(sp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(sp), intent(inout) :: y(*)
      end subroutine ssymv

      subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, lda, incx, incy
          real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(dp), intent(inout) :: y(*)
      end subroutine dsymv

      subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, lda, incx, incy
          complex(cp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(cp), intent(inout) :: y(*)
      end subroutine chemv

      subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
          import :: zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, lda, incx, incy
          complex(zp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(zp), intent(inout) :: y(*)
      end subroutine zhemv

      ! Symmetric/Hermitian band matrix-vector multiplication
      subroutine ssbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, k, lda, incx, incy
          real(sp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(sp), intent(inout) :: y(*)
      end subroutine ssbmv

      subroutine dsbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, k, lda, incx, incy
          real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
          real(dp), intent(inout) :: y(*)
      end subroutine dsbmv

      subroutine chbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
          import :: cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, k, lda, incx, incy
          complex(cp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(cp), intent(inout) :: y(*)
      end subroutine chbmv

      subroutine zhbmv(uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)
          import :: zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, k, lda, incx, incy
          complex(zp), intent(in) :: alpha, beta, a(lda, *), x(*)
          complex(zp), intent(inout) :: y(*)
      end subroutine zhbmv

      ! Symmetric/Hermitian packed matrix-vector multiplication
      subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: alpha, beta, ap(*), x(*)
          real(sp), intent(inout) :: y(*)
      end subroutine sspmv

      subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: alpha, beta, ap(*), x(*)
          real(dp), intent(inout) :: y(*)
      end subroutine dspmv

      subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
          import :: cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: alpha, beta, ap(*), x(*)
          complex(cp), intent(inout) :: y(*)
      end subroutine chpmv

      subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
          import :: zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: alpha, beta, ap(*), x(*)
          complex(zp), intent(inout) :: y(*)
      end subroutine zhpmv

      ! Triangular matrix-vector multiplication
      subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          real(sp), intent(in) :: a(lda, *)
          real(sp), intent(inout) :: x(*)
      end subroutine strmv

      subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          real(dp), intent(in) :: a(lda, *)
          real(dp), intent(inout) :: x(*)
      end subroutine dtrmv

      subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          complex(cp), intent(in) :: a(lda, *)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctrmv

      subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          complex(zp), intent(in) :: a(lda, *)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztrmv

      ! Triangular band matrix-vector multiplication
      subroutine stbmv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          real(sp), intent(in) :: a(lda, *)
          real(sp), intent(inout) :: x(*)
      end subroutine stbmv

      subroutine dtbmv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          real(dp), intent(in) :: a(lda, *)
          real(dp), intent(inout) :: x(*)
      end subroutine dtbmv

      subroutine ctbmv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          complex(cp), intent(in) :: a(lda, *)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctbmv

      subroutine ztbmv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          complex(zp), intent(in) :: a(lda, *)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztbmv

      ! Triangular packed matrix-vector multiplication
      subroutine stpmv(uplo, trans, diag, n, ap, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: ap(*)
          real(sp), intent(inout) :: x(*)
      end subroutine stpmv

      subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: ap(*)
          real(dp), intent(inout) :: x(*)
      end subroutine dtpmv

      subroutine ctpmv(uplo, trans, diag, n, ap, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: ap(*)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctpmv

      subroutine ztpmv(uplo, trans, diag, n, ap, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: ap(*)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztpmv

      ! Triangular solve
      subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          real(sp), intent(in) :: a(lda, *)
          real(sp), intent(inout) :: x(*)
      end subroutine strsv

      subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          real(dp), intent(in) :: a(lda, *)
          real(dp), intent(inout) :: x(*)
      end subroutine dtrsv

      subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          complex(cp), intent(in) :: a(lda, *)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctrsv

      subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, lda, incx
          complex(zp), intent(in) :: a(lda, *)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztrsv

      ! Triangular band solve
      subroutine stbsv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          real(sp), intent(in) :: a(lda, *)
          real(sp), intent(inout) :: x(*)
      end subroutine stbsv

      subroutine dtbsv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          real(dp), intent(in) :: a(lda, *)
          real(dp), intent(inout) :: x(*)
      end subroutine dtbsv

      subroutine ctbsv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          complex(cp), intent(in) :: a(lda, *)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctbsv

      subroutine ztbsv(uplo, trans, diag, n, k, a, lda, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, k, lda, incx
          complex(zp), intent(in) :: a(lda, *)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztbsv

      ! Triangular packed solve
      subroutine stpsv(uplo, trans, diag, n, ap, x, incx)
          import :: sp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: ap(*)
          real(sp), intent(inout) :: x(*)
      end subroutine stpsv

      subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
          import :: dp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: ap(*)
          real(dp), intent(inout) :: x(*)
      end subroutine dtpsv

      subroutine ctpsv(uplo, trans, diag, n, ap, x, incx)
          import :: cp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          complex(cp), intent(in) :: ap(*)
          complex(cp), intent(inout) :: x(*)
      end subroutine ctpsv

      subroutine ztpsv(uplo, trans, diag, n, ap, x, incx)
          import :: zp, blasint
          character, intent(in) :: uplo, trans, diag
          integer(blasint), intent(in) :: n, incx
          complex(zp), intent(in) :: ap(*)
          complex(zp), intent(inout) :: x(*)
      end subroutine ztpsv

      ! Rank-1 update
      subroutine sger(m, n, alpha, x, incx, y, incy, a, lda)
          import :: sp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          real(sp), intent(in) :: alpha, x(*), y(*)
          real(sp), intent(inout) :: a(lda, *)
      end subroutine sger

      subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
          import :: dp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          real(dp), intent(in) :: alpha, x(*), y(*)
          real(dp), intent(inout) :: a(lda, *)
      end subroutine dger

      subroutine cgeru(m, n, alpha, x, incx, y, incy, a, lda)
          import :: cp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          complex(cp), intent(in) :: alpha, x(*), y(*)
          complex(cp), intent(inout) :: a(lda, *)
      end subroutine cgeru

      subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
          import :: zp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          complex(zp), intent(in) :: alpha, x(*), y(*)
          complex(zp), intent(inout) :: a(lda, *)
      end subroutine zgeru

      subroutine cgerc(m, n, alpha, x, incx, y, incy, a, lda)
          import :: cp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          complex(cp), intent(in) :: alpha, x(*), y(*)
          complex(cp), intent(inout) :: a(lda, *)
      end subroutine cgerc

      subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
          import :: zp, blasint
          integer(blasint), intent(in) :: m, n, incx, incy, lda
          complex(zp), intent(in) :: alpha, x(*), y(*)
          complex(zp), intent(inout) :: a(lda, *)
      end subroutine zgerc

      ! Symmetric/Hermitian rank-1 update
      subroutine ssyr(uplo, n, alpha, x, incx, a, lda)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, lda
          real(sp), intent(in) :: alpha, x(*)
          real(sp), intent(inout) :: a(lda, *)
      end subroutine ssyr

      subroutine dsyr(uplo, n, alpha, x, incx, a, lda)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, lda
          real(dp), intent(in) :: alpha, x(*)
          real(dp), intent(inout) :: a(lda, *)
      end subroutine dsyr

      subroutine cher(uplo, n, alpha, x, incx, a, lda)
          import :: sp, cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, lda
          real(sp), intent(in) :: alpha
          complex(cp), intent(in) :: x(*)
          complex(cp), intent(inout) :: a(lda, *)
      end subroutine cher

      subroutine zher(uplo, n, alpha, x, incx, a, lda)
          import :: dp, zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, lda
          real(dp), intent(in) :: alpha
          complex(zp), intent(in) :: x(*)
          complex(zp), intent(inout) :: a(lda, *)
      end subroutine zher

      ! Symmetric/Hermitian packed rank-1 update
      subroutine sspr(uplo, n, alpha, x, incx, ap)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: alpha, x(*)
          real(sp), intent(inout) :: ap(*)
      end subroutine sspr

      subroutine dspr(uplo, n, alpha, x, incx, ap)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: alpha, x(*)
          real(dp), intent(inout) :: ap(*)
      end subroutine dspr

      subroutine chpr(uplo, n, alpha, x, incx, ap)
          import :: sp, cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx
          real(sp), intent(in) :: alpha
          complex(cp), intent(in) :: x(*)
          complex(cp), intent(inout) :: ap(*)
      end subroutine chpr

      subroutine zhpr(uplo, n, alpha, x, incx, ap)
          import :: dp, zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx
          real(dp), intent(in) :: alpha
          complex(zp), intent(in) :: x(*)
          complex(zp), intent(inout) :: ap(*)
      end subroutine zhpr

      ! Symmetric/Hermitian rank-2 update
      subroutine ssyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy, lda
          real(sp), intent(in) :: alpha, x(*), y(*)
          real(sp), intent(inout) :: a(lda, *)
      end subroutine ssyr2

      subroutine dsyr2(uplo, n, alpha, x, incx, y, incy, a, lda)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy, lda
          real(dp), intent(in) :: alpha, x(*), y(*)
          real(dp), intent(inout) :: a(lda, *)
      end subroutine dsyr2

      subroutine cher2(uplo, n, alpha, x, incx, y, incy, a, lda)
          import :: cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy, lda
          complex(cp), intent(in) :: alpha, x(*), y(*)
          complex(cp), intent(inout) :: a(lda, *)
      end subroutine cher2

      subroutine zher2(uplo, n, alpha, x, incx, y, incy, a, lda)
          import :: zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy, lda
          complex(zp), intent(in) :: alpha, x(*), y(*)
          complex(zp), intent(inout) :: a(lda, *)
      end subroutine zher2

      ! Symmetric/Hermitian packed rank-2 update
      subroutine sspr2(uplo, n, alpha, x, incx, y, incy, ap)
          import :: sp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          real(sp), intent(in) :: alpha, x(*), y(*)
          real(sp), intent(inout) :: ap(*)
      end subroutine sspr2

      subroutine dspr2(uplo, n, alpha, x, incx, y, incy, ap)
          import :: dp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          real(dp), intent(in) :: alpha, x(*), y(*)
          real(dp), intent(inout) :: ap(*)
      end subroutine dspr2

      subroutine chpr2(uplo, n, alpha, x, incx, y, incy, ap)
          import :: cp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          complex(cp), intent(in) :: alpha, x(*), y(*)
          complex(cp), intent(inout) :: ap(*)
      end subroutine chpr2

      subroutine zhpr2(uplo, n, alpha, x, incx, y, incy, ap)
          import :: zp, blasint
          character, intent(in) :: uplo
          integer(blasint), intent(in) :: n, incx, incy
          complex(zp), intent(in) :: alpha, x(*), y(*)
          complex(zp), intent(inout) :: ap(*)
      end subroutine zhpr2
  end interface
end module blas_level2

module blas_level3
  use blas_types
  implicit none
  private

  ! Public interface for BLAS Level-3 functions
  public :: sgemm, ssymm, ssyrk, ssyr2k, strmm, strsm
  public :: dgemm, dsymm, dsyrk, dsyr2k, dtrmm, dtrsm
  public :: cgemm, csymm, chemm, csyrk, cherk, csyr2k, cher2k, ctrmm, ctrsm
  public :: zgemm, zsymm, zhemm, zsyrk, zherk, zsyr2k, zher2k, ztrmm, ztrsm

  interface
      ! General matrix-matrix multiplication
      subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: sp, blasint
          character, intent(in) :: transa, transb
          integer(blasint), intent(in) :: m, n, k, lda, ldb, ldc
          real(sp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(sp), intent(inout) :: c(ldc, *)
      end subroutine sgemm

      subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: dp, blasint
          character, intent(in) :: transa, transb
          integer(blasint), intent(in) :: m, n, k, lda, ldb, ldc
          real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(dp), intent(inout) :: c(ldc, *)
      end subroutine dgemm

      subroutine cgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: cp, blasint
          character, intent(in) :: transa, transb
          integer(blasint), intent(in) :: m, n, k, lda, ldb, ldc
          complex(cp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine cgemm

      subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: zp, blasint
          character, intent(in) :: transa, transb
          integer(blasint), intent(in) :: m, n, k, lda, ldb, ldc
          complex(zp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zgemm

      ! Symmetric matrix-matrix multiplication
      subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: sp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          real(sp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssymm

      subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: dp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsymm

      subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: cp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          complex(cp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine csymm

      subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: zp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          complex(zp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zsymm

      subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: cp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          complex(cp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine chemm

      subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: zp, blasint
          character, intent(in) :: side, uplo
          integer(blasint), intent(in) :: m, n, lda, ldb, ldc
          complex(zp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zhemm

      ! Symmetric rank-k update
      subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: sp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          real(sp), intent(in) :: alpha, beta, a(lda, *)
          real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssyrk

      subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: dp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          real(dp), intent(in) :: alpha, beta, a(lda, *)
          real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsyrk

      subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: cp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          complex(cp), intent(in) :: alpha, beta, a(lda, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine csyrk

      subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: zp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          complex(zp), intent(in) :: alpha, beta, a(lda, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zsyrk

      ! Hermitian rank-k update
      subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: sp, cp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          real(sp), intent(in) :: alpha, beta
          complex(cp), intent(in) :: a(lda, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine cherk

      subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
          import :: dp, zp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldc
          real(dp), intent(in) :: alpha, beta
          complex(zp), intent(in) :: a(lda, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zherk

      ! Symmetric rank-2k update
      subroutine ssyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: sp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          real(sp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(sp), intent(inout) :: c(ldc, *)
      end subroutine ssyr2k

      subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: dp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          real(dp), intent(inout) :: c(ldc, *)
      end subroutine dsyr2k

      subroutine csyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: cp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          complex(cp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine csyr2k

      subroutine zsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: zp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          complex(zp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zsyr2k

      ! Hermitian rank-2k update
      subroutine cher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: sp, cp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          complex(cp), intent(in) :: alpha, a(lda, *), b(ldb, *)
          real(sp), intent(in) :: beta
          complex(cp), intent(inout) :: c(ldc, *)
      end subroutine cher2k

      subroutine zher2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
          import :: dp, zp, blasint
          character, intent(in) :: uplo, trans
          integer(blasint), intent(in) :: n, k, lda, ldb, ldc
          complex(zp), intent(in) :: alpha, a(lda, *), b(ldb, *)
          real(dp), intent(in) :: beta
          complex(zp), intent(inout) :: c(ldc, *)
      end subroutine zher2k

      ! Triangular matrix-matrix multiplication
      subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: sp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          real(sp), intent(in) :: alpha, a(lda, *)
          real(sp), intent(inout) :: b(ldb, *)
      end subroutine strmm

      subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: dp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          real(dp), intent(in) :: alpha, a(lda, *)
          real(dp), intent(inout) :: b(ldb, *)
      end subroutine dtrmm

      subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: cp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          complex(cp), intent(in) :: alpha, a(lda, *)
          complex(cp), intent(inout) :: b(ldb, *)
      end subroutine ctrmm

      subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: zp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          complex(zp), intent(in) :: alpha, a(lda, *)
          complex(zp), intent(inout) :: b(ldb, *)
      end subroutine ztrmm

      ! Triangular solve
      subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: sp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          real(sp), intent(in) :: alpha, a(lda, *)
          real(sp), intent(inout) :: b(ldb, *)
      end subroutine strsm

      subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: dp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          real(dp), intent(in) :: alpha, a(lda, *)
          real(dp), intent(inout) :: b(ldb, *)
      end subroutine dtrsm

      subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: cp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          complex(cp), intent(in) :: alpha, a(lda, *)
          complex(cp), intent(inout) :: b(ldb, *)
      end subroutine ctrsm

      subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
          import :: zp, blasint
          character, intent(in) :: side, uplo, transa, diag
          integer(blasint), intent(in) :: m, n, lda, ldb
          complex(zp), intent(in) :: alpha, a(lda, *)
          complex(zp), intent(inout) :: b(ldb, *)
      end subroutine ztrsm
  end interface
end module blas_level3

module blas_interface
  use blas_level1
  use blas_level2
  use blas_level3
  implicit none
  private

  ! Public generic interfaces
  public :: blas_rotg, blas_rotmg, blas_rot, blas_rotm
  public ::  blas_swap, blas_scal, blas_copy, blas_axpy, blas_dot, blas_dotu, blas_dotc, blas_nrm2, blas_asum, blas_iamax
  public :: blas_gemv, blas_gbmv, blas_symv, blas_sbmv, blas_spmv
  public :: blas_trmv, blas_tbmv, blas_tpmv, blas_trsv, blas_tbsv, blas_tpsv
  public :: blas_ger, blas_geru, blas_gerc, blas_syr, blas_spr, blas_syr2, blas_spr2
  public :: blas_hemv, blas_hbmv, blas_hpmv, blas_her, blas_hpr, blas_her2, blas_hpr2
  public :: blas_gemm, blas_symm, blas_hemm, blas_syrk, blas_herk, blas_syr2k, blas_her2k, blas_trmm, blas_trsm

  ! Generic interfaces for BLAS Level-1 functions
  interface blas_rotg
      procedure srotg, drotg
  end interface blas_rotg

  interface blas_rotmg
      procedure srotmg, drotmg
  end interface blas_rotmg

  interface blas_rot
      procedure srot, drot
  end interface blas_rot

  interface blas_rotm
      procedure srotm, drotm
  end interface blas_rotm

  interface blas_swap
      procedure sswap, dswap, cswap, zswap
  end interface blas_swap

  interface blas_scal
      procedure sscal, dscal, cscal, zscal
  end interface blas_scal

  interface blas_copy
      procedure scopy, dcopy, ccopy, zcopy
  end interface blas_copy

  interface blas_axpy
      procedure saxpy, daxpy, caxpy, zaxpy
  end interface blas_axpy

  interface blas_dot
      procedure sdot, ddot, cdotc, zdotc
  end interface blas_dot

  interface blas_dotu
      procedure sdot, ddot, cdotu, zdotu
  end interface blas_dotu

  interface blas_dotc
      procedure sdot, ddot, cdotc, zdotc
  end interface blas_dotc

  interface blas_nrm2
      procedure snrm2, dnrm2, cnrm2, znrm2
  end interface blas_nrm2

  interface blas_asum
      procedure sasum, dasum, casum, zasum
  end interface blas_asum

  interface blas_iamax
      procedure isamax, idamax, icamax, izamax
  end interface blas_iamax

  ! Generic interfaces for BLAS Level-2 functions
  interface blas_gemv
      procedure sgemv, dgemv, cgemv, zgemv
  end interface blas_gemv

  interface blas_gbmv
      procedure sgbmv, dgbmv, cgbmv, zgbmv
  end interface blas_gbmv

  interface blas_symv
      procedure ssymv, dsymv
  end interface blas_symv

  interface blas_sbmv
      procedure ssbmv, dsbmv
  end interface blas_sbmv

  interface blas_spmv
      procedure sspmv, dspmv
  end interface blas_spmv

  interface blas_trmv
      procedure strmv, dtrmv, ctrmv, ztrmv
  end interface blas_trmv

  interface blas_tbmv
      procedure stbmv, dtbmv, ctbmv, ztbmv
  end interface blas_tbmv

  interface blas_tpmv
      procedure stpmv, dtpmv, ctpmv, ztpmv
  end interface blas_tpmv

  interface blas_trsv
      procedure strsv, dtrsv, ctrsv, ztrsv
  end interface blas_trsv

  interface blas_tbsv
      procedure stbsv, dtbsv, ctbsv, ztbsv
  end interface blas_tbsv

  interface blas_tpsv
      procedure stpsv, dtpsv, ctpsv, ztpsv
  end interface blas_tpsv

  interface blas_ger
      procedure sger, dger, cgerc, zgerc
  end interface blas_ger

  interface blas_geru
      procedure sger, dger, cgeru, zgeru
  end interface blas_geru

  interface blas_gerc
      procedure sger, dger, cgerc, zgerc
  end interface blas_gerc

  interface blas_syr
      procedure ssyr, dsyr
  end interface blas_syr

  interface blas_spr
      procedure sspr, dspr
  end interface blas_spr

  interface blas_syr2
      procedure ssyr2, dsyr2
  end interface blas_syr2

  interface blas_spr2
      procedure sspr2, dspr2
  end interface blas_spr2

  interface blas_hemv
      procedure ssymv, dsymv, chemv, zhemv
  end interface blas_hemv

  interface blas_hbmv
      procedure ssbmv, dsbmv, chbmv, zhbmv
  end interface blas_hbmv

  interface blas_hpmv
      procedure sspmv, dspmv, chpmv, zhpmv
  end interface blas_hpmv

  interface blas_her
      procedure ssyr, dsyr, cher, zher
  end interface blas_her

  interface blas_hpr
      procedure sspr, dspr, chpr, zhpr
  end interface blas_hpr

  interface blas_her2
      procedure ssyr2, dsyr2, cher2, zher2
  end interface blas_her2

  interface blas_hpr2
      procedure sspr2, dspr2, chpr2, zhpr2
  end interface blas_hpr2

  ! Generic interfaces for BLAS Level-3 functions
  interface blas_gemm
      procedure sgemm, dgemm, cgemm, zgemm
  end interface blas_gemm

  interface blas_symm
      procedure ssymm, dsymm, csymm, zsymm
  end interface blas_symm

  interface blas_hemm
      procedure ssymm, dsymm, chemm, zhemm
  end interface blas_hemm

  interface blas_syrk
      procedure ssyrk, dsyrk, csyrk, zsyrk
  end interface blas_syrk

  interface blas_herk
      procedure ssyrk, dsyrk, cherk, zherk
  end interface blas_herk

  interface blas_syr2k
      procedure ssyr2k, dsyr2k, csyr2k, zsyr2k
  end interface blas_syr2k

  interface blas_her2k
      procedure ssyr2k, dsyr2k, cher2k, zher2k
  end interface blas_her2k

  interface blas_trmm
      procedure strmm, dtrmm, ctrmm, ztrmm
  end interface blas_trmm

  interface blas_trsm
      procedure strsm, dtrsm, ctrsm, ztrsm
  end interface blas_trsm
end module blas_interface
