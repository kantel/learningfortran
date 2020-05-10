!-------------------------------------------------------------------------------
!
! Fortran 2008 example that solves Lotka-Volterra ODE using RK4.
!
! $ gfortran -O2 -std=f2008 -o lotka lotka.f90
! $ ./lotka > output.txt
! $ gnuplot -c plot.plt -p
!
!-------------------------------------------------------------------------------
!
! Functions:
!
!     u(t) >= 0 -   Population size of prey at time t.
!     v(t) >= 0 -   Population size of predator at time t.
!
! Equations:
!
!     d/dt u = u * (alpha - beta * v)
!     d/dt v = -v * (gamma - delta * u)
!
! Model Parameters:
!
!     u         -   Prey population.
!     v         -   Predator population.
!     alpha     -   Reproduction rate of prey.
!     beta      -   Death rate of prey by predator.
!     gamma     -   Death rate of predator.
!     delta     -   Reproduction rate of predator by prey.
!
!-------------------------------------------------------------------------------
program main
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    real(kind=wp), parameter :: t_max = 30_wp
    real(kind=wp), parameter :: h     = 0.001_wp
    integer,       parameter :: n     = t_max / h
    integer                  :: i
    real(kind=wp)            :: t(n)  = [ (h * i, i = 1,  n) ]
    real(kind=wp)            :: r(2)
    real(kind=wp)            :: x(n), y(n)

    r = [ 20.0_wp, & ! Initial prey population.
           5.0_wp ]  ! Initial predator population.

    do i = 1, n
        x(i) = r(1)
        y(i) = r(2)

        r = r + rk4(r, t(i), h)

        print '(f15.8, 2(" ", f15.8))', t(i), x(i), y(i)
    end do
contains
    function rk4(r, t, h)
        !! Runge-Kutta 4th order solver.
        real(kind=wp), intent(in) :: r(2)   ! Initial values.
        real(kind=wp), intent(in) :: t      ! Step.
        real(kind=wp), intent(in) :: h      ! Step size.
        real(kind=wp)             :: rk4(2)
        real(kind=wp)             :: k1(2), k2(2), k3(2), k4(2)

        k1 = h * f(r,            t)
        k2 = h * f(r + 0.5 * k1, t + 0.5 * h)
        k3 = h * f(r + 0.5 * k2, t + 0.5 * h)
        k4 = h * f(r + k3,       t + h)

        rk4 = (k1 + (2 * k2) + (2 * k3) + k4) / 6
    end function rk4

    function f(r, t)
        !! Lotka-Volterra ODE.
        !!
        !! Point estimates for model parameters taken from:
        !!     https://www.math.tamu.edu/~phoward/m442/modbasics.pdf
        real(kind=wp), parameter  :: alpha = 0.47_wp
        real(kind=wp), parameter  :: beta  = 0.024_wp
        real(kind=wp), parameter  :: delta = 0.023_wp
        real(kind=wp), parameter  :: gamma = 0.76_wp
        real(kind=wp), intent(in) :: r(2)
        real(kind=wp), intent(in) :: t
        real(kind=wp)             :: f(2)
        real(kind=wp)             :: u, v

        u = r(1)
        v = r(2)

        f(1) =  u * (alpha - beta * v)
        f(2) = -v * (gamma - delta * u)
    end function f
end program main
