program generate_lambert
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    implicit none

    integer, parameter :: dp = real64

    ! Constants and parameters
    integer(kind=int64), parameter :: Q      = 10000000_int64
    real(dp),            parameter :: mu     = 1.0_dp
    real(dp),            parameter :: a_min  = 1.0_dp
    real(dp),            parameter :: a_max  = 10.0_dp
    real(dp),            parameter :: e_max  = 0.9_dp
    real(dp),            parameter :: minSep = 1.0e-6_dp
    real(dp),            parameter :: PI     = 3.141592653589793238462643383279502884_dp
    real(dp),            parameter :: TWO_PI = 2.0_dp * PI

    integer(kind=int64), parameter :: block_size = 1000000_int64

    real(dp), allocatable :: r1(:,:), v1(:,:), r2(:,:), v2(:,:)
    real(dp), allocatable :: tof(:), direction(:)

    integer(kind=int64) :: i0, i1, n, i
    integer :: file_unit, io_stat
    integer(kind=int64) :: rng_s0, rng_s1

    real(dp) :: a, e, inc, RAAN, AOP, nu1, nu2, dnu
    real(dp) :: p, sqrtMP, rm1, rm2, rp1x, rp1y, rp2x, rp2y, vp1x, vp1y, vp2x, vp2y
    real(dp) :: cO, sO, ci, si, cw, sw
    real(dp) :: R11, R12, R21, R22, R31, R32
    real(dp) :: E1, E2, M1, M2, dM, nMean, dt, d

    print '(A, I0, A)', 'Setting up ', Q, ' Lambert problems...'

    allocate(r1(3, block_size), v1(3, block_size))
    allocate(r2(3, block_size), v2(3, block_size))
    allocate(tof(block_size), direction(block_size))

    open(newunit=file_unit, file='datalambert_fortran.bin', &
         access='stream', form='unformatted', status='replace', action='write', iostat=io_stat)
    if (io_stat /= 0) error stop 'Error opening output file.'

    ! File layout expected by benchmark_lambert.f90 (little-endian stream):
    !   int64 Q | direction(Q) | r1(3,Q) | r2(3,Q) | tof(Q) | v1(3,Q) | v2(3,Q)
    ! Each block is written at its offset inside every array section.
    write(file_unit, pos=1) Q

    i0 = 1_int64
    do while (i0 <= Q)
        i1 = min(i0 + block_size - 1_int64, Q)
        n  = i1 - i0 + 1_int64

        !$omp parallel do default(shared) private(i, rng_s0, rng_s1, a, e, inc, RAAN, AOP, &
        !$omp nu1, nu2, dnu, p, sqrtMP, rm1, rm2, rp1x, rp1y, rp2x, rp2y, vp1x, vp1y, &
        !$omp vp2x, vp2y, cO, sO, ci, si, cw, sw, R11, R12, R21, R22, R31, R32, &
        !$omp E1, E2, M1, M2, dM, nMean, dt, d) schedule(static)
        do i = 1, n
            ! Thread-safe local PRNG state initialized uniquely per index
            rng_s0 = (i0 + i - 1_int64) * 6364136223846793005_int64 + 1_int64
            rng_s1 = rng_s0 * 1442695040888963407_int64 + 1_int64

            ! 1. Classical orbital elements
            a    = a_min + (a_max - a_min) * get_uniform_rand(rng_s0, rng_s1)
            e    = e_max * get_uniform_rand(rng_s0, rng_s1)
            inc  = PI * get_uniform_rand(rng_s0, rng_s1)
            RAAN = TWO_PI * get_uniform_rand(rng_s0, rng_s1)
            AOP  = TWO_PI * get_uniform_rand(rng_s0, rng_s1)

            nu1  = TWO_PI * get_uniform_rand(rng_s0, rng_s1)
            nu2  = TWO_PI * get_uniform_rand(rng_s0, rng_s1)

            ! Collinear check
            dnu = modulo(nu2 - nu1, TWO_PI)
            do while (dnu < minSep .or. dnu > (TWO_PI - minSep))
                nu2 = TWO_PI * get_uniform_rand(rng_s0, rng_s1)
                dnu = modulo(nu2 - nu1, TWO_PI)
            end do

            ! 2. Perifocal states
            p      = a * (1.0_dp - e**2)
            sqrtMP = sqrt(mu / p)

            rm1 = p / (1.0_dp + e * cos(nu1))
            rm2 = p / (1.0_dp + e * cos(nu2))

            rp1x = rm1 * cos(nu1)
            rp1y = rm1 * sin(nu1)
            rp2x = rm2 * cos(nu2)
            rp2y = rm2 * sin(nu2)

            vp1x = -sqrtMP * sin(nu1)
            vp1y =  sqrtMP * (e + cos(nu1))
            vp2x = -sqrtMP * sin(nu2)
            vp2y =  sqrtMP * (e + cos(nu2))

            ! 3. Rotation to inertial
            cO = cos(RAAN);  sO = sin(RAAN)
            ci = cos(inc);   si = sin(inc)
            cw = cos(AOP);   sw = sin(AOP)

            R11 =  cO * cw - sO * sw * ci
            R12 = -cO * sw - sO * cw * ci
            R21 =  sO * cw + cO * sw * ci
            R22 = -sO * sw + cO * cw * ci
            R31 =  sw * si
            R32 =  cw * si

            r1(1, i) = R11 * rp1x + R12 * rp1y
            r1(2, i) = R21 * rp1x + R22 * rp1y
            r1(3, i) = R31 * rp1x + R32 * rp1y

            v1(1, i) = R11 * vp1x + R12 * vp1y
            v1(2, i) = R21 * vp1x + R22 * vp1y
            v1(3, i) = R31 * vp1x + R32 * vp1y

            r2(1, i) = R11 * rp2x + R12 * rp2y
            r2(2, i) = R21 * rp2x + R22 * rp2y
            r2(3, i) = R31 * rp2x + R32 * rp2y

            v2(1, i) = R11 * vp2x + R12 * vp2y
            v2(2, i) = R21 * vp2x + R22 * vp2y
            v2(3, i) = R31 * vp2x + R32 * vp2y

            ! 4. Anomalies
            E1 = 2.0_dp * atan2(sqrt(1.0_dp - e) * sin(nu1 * 0.5_dp), sqrt(1.0_dp + e) * cos(nu1 * 0.5_dp))
            E2 = 2.0_dp * atan2(sqrt(1.0_dp - e) * sin(nu2 * 0.5_dp), sqrt(1.0_dp + e) * cos(nu2 * 0.5_dp))

            M1 = modulo(E1 - e * sin(E1), TWO_PI)
            M2 = modulo(E2 - e * sin(E2), TWO_PI)

            nMean = sqrt(mu / (a**3))

            ! 5. Time of flight
            dM = modulo(M2 - M1, TWO_PI)
            dt = dM / nMean
            tof(i) = dt

            ! 6. Transfer direction
            d = 1.0_dp
            if (dnu > PI) d = -1.0_dp
            direction(i) = d
        end do
        !$omp end parallel do

        ! Stream block into each array section (8-byte header, 8-byte reals)
        write(file_unit, pos=section_pos(0_int64,  1_int64)) direction(1:n)
        write(file_unit, pos=section_pos(1_int64,  3_int64)) r1(:, 1:n)
        write(file_unit, pos=section_pos(4_int64,  3_int64)) r2(:, 1:n)
        write(file_unit, pos=section_pos(7_int64,  1_int64)) tof(1:n)
        write(file_unit, pos=section_pos(8_int64,  3_int64)) v1(:, 1:n)
        write(file_unit, pos=section_pos(11_int64, 3_int64)) v2(:, 1:n)

        print '(A, I0, A, I0, A)', 'Processed block: ', i0, ' to ', i1, '...'
        i0 = i0 + block_size
    end do

    close(file_unit)
    print *, "Data saved to datalambert_fortran.bin"

    deallocate(r1, v1, r2, v2, tof, direction)

contains

    ! Byte position of case i0 within a section that starts after `offset`
    ! values per case and stores `width` values per case.
    integer(kind=int64) function section_pos(offset, width)
        integer(kind=int64), intent(in) :: offset, width
        section_pos = 9_int64 + 8_int64 * (offset * Q + width * (i0 - 1_int64))
    end function section_pos

    ! Note: Removed 'pure' keyword so dummy arguments s0 and s1 can be mutated
    function get_uniform_rand(s0, s1) result(u)
        integer(kind=int64), intent(inout) :: s0, s1
        real(dp) :: u
        integer(kind=int64) :: x, y
        real(dp), parameter :: norm = 1.0_dp / 9007199254740992.0_dp  ! 2^-53

        x  = s0
        y  = s1
        s0 = y
        x  = ieor(x, ishft(x, 23))
        s1 = ieor(ieor(x, y), ieor(ishft(x, -17), ishft(y, -26)))

        u = real(ishft(s0 + s1, -11), kind=dp) * norm
    end function get_uniform_rand

end program generate_lambert