program benchmark_lambert
    use iso_fortran_env, only: int64, real64, output_unit
    use ieee_arithmetic, only: ieee_is_nan
    use omp_lib, only: omp_set_dynamic, omp_set_num_threads, omp_get_max_threads
    use exact_mod_fast, only: dp, exactsolution
    use ivLamIOmod, only: iu, ru
    implicit none

    integer, parameter :: data_count = 10000000
    integer, parameter :: default_tests = 5000000
    integer, parameter :: histogram_edge_count = 60
    real(dp), parameter :: mu = 1.0_dp

    integer :: test_count, selected_count, i, j, unit, ios, valid_count
    integer :: warm_count, arg_status
    integer(int64) :: file_count, start_tick, stop_tick, clock_rate
    integer(kind=iu) :: init_info, unload_info
    integer(kind=iu), allocatable :: direction(:), info_status(:), half_rev_status(:)
    integer, allocatable :: selected_index(:)
    real(dp), allocatable :: direction_all(:), r1_all(:,:), r2_all(:,:), tof_all(:)
    real(dp), allocatable :: truth_all(:,:)
    real(dp), allocatable :: r1(:,:), r2(:,:), tof(:), v1_true(:,:), v2_true(:,:)
    real(dp), allocatable :: v1_exact(:,:), v2_exact(:,:), v1_ivlam(:,:), v2_ivlam(:,:)
    real(dp), allocatable :: warm_v1(:,:), warm_v2(:,:)
    real(dp), allocatable :: e_a1(:), e_i1(:), e_a2(:), e_i2(:)
    real(dp), allocatable :: e_a_axis(:), e_i_axis(:), e_a_h(:), e_i_h(:)
    real(dp) :: h_vec(3), h_true(3), h_exact(3), h_ivlam(3)
    real(dp) :: r_mag, inv_a_true, inv_a_exact, inv_a_ivlam
    real(dp) :: t_exact, t_ivlam, median_a1, median_i1
    real(dp) :: h(3), zero_rev_tof(1), zero_rev_r1(3,1), zero_rev_r2(3,1)
    real(dp) :: zero_rev_v1(3,1), zero_rev_v2(3,1)
    integer(kind=iu) :: zero_rev_direction(1), zero_rev_info(1), zero_rev_half(1)
    character(len=512) :: data_file, coefficient_file, arg
    logical :: valid

    test_count = default_tests
    call get_command_argument(1, arg, status=arg_status)
    if (arg_status == 0 .and. len_trim(arg) > 0) then
        read(arg, *, iostat=ios) test_count
        if (ios /= 0 .or. test_count < 1 .or. test_count > data_count) then
            error stop 'First argument must be a test count from 1 to 10000000.'
        end if
    end if

    data_file = 'datalambert_fortran.bin'
    coefficient_file = 'ivLamTree_20210202_160219_i2d8.bin'
    call get_command_argument(2, arg, status=arg_status)
    if (arg_status == 0 .and. len_trim(arg) > 0) coefficient_file = trim(arg)

    open(newunit=unit, file=trim(data_file), form='unformatted', access='stream', &
         status='old', action='read', convert='little_endian', iostat=ios)
    if (ios /= 0) error stop 'Could not open datalambert_fortran.bin.'
    read(unit, iostat=ios) file_count
    if (ios /= 0 .or. file_count /= int(data_count, int64)) then
        error stop 'Dataset header does not contain the expected 10000000 cases.'
    end if

    allocate(direction_all(data_count), r1_all(3,data_count), r2_all(3,data_count))
    read(unit, iostat=ios) direction_all
    if (ios /= 0) error stop 'Could not read dataset direction values.'
    read(unit, iostat=ios) r1_all
    if (ios /= 0) error stop 'Could not read dataset r1 values.'
    read(unit, iostat=ios) r2_all
    if (ios /= 0) error stop 'Could not read dataset r2 values.'

    allocate(direction(test_count), selected_index(test_count))
    allocate(r1(3,test_count), r2(3,test_count), tof(test_count))
    allocate(v1_true(3,test_count), v2_true(3,test_count))
    selected_count = 0
    do i = 1, data_count
        h_vec(3) = r1_all(1,i)*r2_all(2,i) - r1_all(2,i)*r2_all(1,i)
        if (h_vec(3) == 0.0_dp) h_vec(3) = 1.0_dp
        if (sign(1.0_dp, h_vec(3))*direction_all(i) > 0.0_dp) then
            selected_count = selected_count + 1
            if (selected_count > test_count) exit
            selected_index(selected_count) = i
            direction(selected_count) = nint(direction_all(i), kind=iu)
            r1(:,selected_count) = r1_all(:,i)
            r2(:,selected_count) = r2_all(:,i)
        end if
    end do
    if (selected_count < test_count) error stop 'Not enough cases satisfy the requested geometry filter.'
    deallocate(direction_all, r1_all, r2_all)

    allocate(tof_all(data_count))
    read(unit, iostat=ios) tof_all
    if (ios /= 0) error stop 'Could not read dataset tof values.'
    do i = 1, test_count
        tof(i) = tof_all(selected_index(i))
    end do
    deallocate(tof_all)

    allocate(truth_all(3,data_count))
    read(unit, iostat=ios) truth_all
    if (ios /= 0) error stop 'Could not read dataset v1 values.'
    do i = 1, test_count
        v1_true(:,i) = truth_all(:,selected_index(i))
    end do
    deallocate(truth_all)

    allocate(truth_all(3,data_count))
    read(unit, iostat=ios) truth_all
    if (ios /= 0) error stop 'Could not read dataset v2 values.'
    do i = 1, test_count
        v2_true(:,i) = truth_all(:,selected_index(i))
    end do
    deallocate(truth_all, selected_index)
    close(unit)

    call omp_set_dynamic(.false.)
    call ivLam_initialize(-1_iu, trim(coefficient_file), init_info)
    if (init_info /= 0_iu) error stop 'ivLam initialization failed.'

    ! -------------------------------------------------------------------------
    ! WARM-UP SECTION (100 Cases)
    ! -------------------------------------------------------------------------
    h = [0.0_dp, 0.0_dp, 1.0_dp]
    warm_count = min(test_count, 100)
    allocate(warm_v1(3,warm_count), warm_v2(3,warm_count))
    call exactsolution(r1(:,1:warm_count), r2(:,1:warm_count), &
                          tof(1:warm_count), mu, h, warm_v1, warm_v2, warm_count)
    deallocate(warm_v1, warm_v2)

    zero_rev_r1(:,1) = r1(:,1)
    zero_rev_r2(:,1) = r2(:,1)
    zero_rev_tof(1) = tof(1)
    zero_rev_direction(1) = direction(1)
    call omp_set_num_threads(1)
    call ivLam_zeroRev_multipleInput(1_iu, zero_rev_r1, zero_rev_r2, zero_rev_tof, &
                                     zero_rev_direction, zero_rev_v1, zero_rev_v2, &
                                     zero_rev_info, zero_rev_half)

    write(output_unit,'(A, I0)') '>>> Active OpenMP thread limit: ', omp_get_max_threads()

    ! -------------------------------------------------------------------------
    ! FULL BENCHMARK SECTION (Single Allocate, Pre-Allocated Timing)
    ! -------------------------------------------------------------------------
    allocate(v1_exact(3,test_count), v2_exact(3,test_count))
    allocate(v1_ivlam(3,test_count), v2_ivlam(3,test_count))
    allocate(info_status(test_count), half_rev_status(test_count))

    ! 1. Run Exact Solution
    call omp_set_num_threads(1)
    call system_clock(start_tick, clock_rate)
    call exactsolution(r1, r2, tof, mu, h, v1_exact, v2_exact, test_count)
    call system_clock(stop_tick)
    t_exact = real(stop_tick - start_tick, dp)/real(clock_rate, dp)

    ! 2. Run ivLam
    call omp_set_num_threads(1)
    call system_clock(start_tick, clock_rate)
    call ivLam_zeroRev_multipleInput(int(test_count,kind=iu), r1, r2, tof, direction, &
                                     v1_ivlam, v2_ivlam, info_status, half_rev_status)
    call system_clock(stop_tick)
    t_ivlam = real(stop_tick - start_tick, dp)/real(clock_rate, dp)

    call finish_benchmark()

contains

    subroutine finish_benchmark()
        integer :: k, histogram_unit
        real(dp) :: err1, err2
        real(dp) :: axis_true, axis_a, axis_i

        valid_count = 0
        do k = 1, test_count
            valid = .not. any(ieee_is_nan(v1_exact(:,k))) .and. &
                    .not. any(ieee_is_nan(v2_exact(:,k))) .and. info_status(k) == 0_iu
            if (valid) valid_count = valid_count + 1
        end do
        if (valid_count == 0) error stop 'No valid solutions to summarize.'

        allocate(e_a1(valid_count), e_i1(valid_count), e_a2(valid_count), e_i2(valid_count))
        allocate(e_a_axis(valid_count), e_i_axis(valid_count), e_a_h(valid_count), e_i_h(valid_count))
        j = 0
        do k = 1, test_count
            valid = .not. any(ieee_is_nan(v1_exact(:,k))) .and. &
                    .not. any(ieee_is_nan(v2_exact(:,k))) .and. info_status(k) == 0_iu
            if (.not. valid) cycle
            j = j + 1
            err1 = sqrt(sum((v1_exact(:,k) - v1_true(:,k))**2))
            err2 = sqrt(sum((v1_ivlam(:,k) - v1_true(:,k))**2))
            e_a1(j) = err1
            e_i1(j) = err2
            e_a2(j) = sqrt(sum((v2_exact(:,k) - v2_true(:,k))**2))
            e_i2(j) = sqrt(sum((v2_ivlam(:,k) - v2_true(:,k))**2))

            r_mag = sqrt(sum(r1(:,k)**2))
            axis_true = 2.0_dp/r_mag - sum(v1_true(:,k)**2)/mu
            axis_a = 2.0_dp/r_mag - sum(v1_exact(:,k)**2)/mu
            axis_i = 2.0_dp/r_mag - sum(v1_ivlam(:,k)**2)/mu
            e_a_axis(j) = abs(axis_a-axis_true)/abs(axis_true)
            e_i_axis(j) = abs(axis_i-axis_true)/abs(axis_true)

            call cross_product(r1(:,k), v1_true(:,k), h_true)
            call cross_product(r1(:,k), v1_exact(:,k), h_exact)
            call cross_product(r1(:,k), v1_ivlam(:,k), h_ivlam)
            e_a_h(j) = sqrt(sum((h_exact-h_true)**2))/sqrt(sum(h_true**2))
            e_i_h(j) = sqrt(sum((h_ivlam-h_true)**2))/sqrt(sum(h_true**2))
        end do

        call sort_ascending(e_a1)
        call sort_ascending(e_i1)
        call sort_ascending(e_a2)
        call sort_ascending(e_i2)
        call sort_ascending(e_a_axis)
        call sort_ascending(e_i_axis)
        call sort_ascending(e_a_h)
        call sort_ascending(e_i_h)
        median_a1 = percentile_sorted(e_a1, 50.0_dp)
        median_i1 = percentile_sorted(e_i1, 50.0_dp)

        write(output_unit,'(/,a)') '========================================================================'
        write(output_unit,'(a,i0,a)') '                           BENCHMARK RESULTS (', test_count, ' Orbits)'
        write(output_unit,'(a)') '========================================================================'
        write(output_unit,'(a)') '  Solver   | Total Time (s) | Per Solution (ns) | Median Error (|dV1|)'
        write(output_unit,'(a)') '-----------+----------------+-------------------+-----------------------'
        write(output_unit,'(a,f8.4,a,f6.1,a,a)') &
            '  Exact Solution    |   ', t_exact, ' s   |     ', t_exact/real(test_count,dp)*1.0e9_dp, &
            ' ns      |      ', sci(median_a1, 4)
        write(output_unit,'(a,f8.4,a,f6.1,a,a)') &
            '  ivLam    |   ', t_ivlam, ' s   |     ', t_ivlam/real(test_count,dp)*1.0e9_dp, &
            ' ns      |      ', sci(median_i1, 4)
        write(output_unit,'(a)') '========================================================================'
        if (t_exact < t_ivlam) then
            write(output_unit,'(a,f0.2,a,a,a)') ' Exact Solution  is ', t_ivlam/t_exact, 'x FASTER '
        else
            write(output_unit,'(a,f0.2,a)') ' ivLam is ', t_exact/t_ivlam, 'x faster.'
        end if
        write(output_unit,'(a,/)') '========================================================================'
        write(output_unit,'(a,i0,a,i0)') 'Valid solutions: ', valid_count, &
            ' / ', test_count

        open(newunit=histogram_unit, file='benchmark_histograms.csv', status='replace', action='write')
        write(histogram_unit,'(a)') 'figure,metric,solver,bin_left,bin_right,count'
        call write_histogram_pair(histogram_unit, 'velocity', 'dV1', e_a1, e_i1)
        call write_histogram_pair(histogram_unit, 'velocity', 'dV2', e_a2, e_i2)
        call write_histogram_pair(histogram_unit, 'orbit', 'relative_inverse_a', e_a_axis, e_i_axis)
        call write_histogram_pair(histogram_unit, 'orbit', 'relative_angular_momentum', e_a_h, e_i_h)
        close(histogram_unit)
        write(output_unit,'(a)') 'Histogram counts written to benchmark_histograms.csv.'

        call ivLam_unloadData(unload_info, .true.)
        if (unload_info /= 0_iu) write(output_unit,'(a,i0)') 'ivLam unload status: ', unload_info
    end subroutine finish_benchmark

    subroutine cross_product(a, b, c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp), intent(out) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end subroutine cross_product

    subroutine sort_ascending(values)
        real(dp), intent(inout) :: values(:)
        integer :: n, start, last, root, child
        real(dp) :: temp
        n = size(values)
        do start = n/2, 1, -1
            root = start
            do while (2*root <= n)
                child = 2*root
                if (child < n) then
                    if (values(child) < values(child+1)) child = child + 1
                end if
                if (values(root) >= values(child)) exit
                temp = values(root); values(root) = values(child); values(child) = temp
                root = child
            end do
        end do
        do last = n, 2, -1
            temp = values(1); values(1) = values(last); values(last) = temp
            root = 1
            do while (2*root < last)
                child = 2*root
                if (child+1 < last) then
                    if (values(child) < values(child+1)) child = child + 1
                end if
                if (values(root) >= values(child)) exit
                temp = values(root); values(root) = values(child); values(child) = temp
                root = child
            end do
        end do
    end subroutine sort_ascending

    function sci(x, digits) result(text)
        real(dp), intent(in) :: x
        integer, intent(in) :: digits
        character(len=:), allocatable :: text
        character(len=32) :: buffer, fmt
        integer :: k
        write(fmt,'(a,i0,a,i0,a)') '(es', digits+8, '.', digits, 'e2)'
        write(buffer, fmt) x
        k = index(buffer, 'E')
        if (k > 0) buffer(k:k) = 'e'
        text = trim(adjustl(buffer))
    end function sci

    real(dp) function percentile_sorted(values, percentile)
        real(dp), intent(in) :: values(:), percentile
        real(dp) :: position, fraction
        integer :: lower, upper
        position = 1.0_dp + real(size(values)-1,dp)*percentile/100.0_dp
        lower = int(position)
        upper = min(lower+1, size(values))
        fraction = position-real(lower,dp)
        percentile_sorted = values(lower) + fraction*(values(upper)-values(lower))
    end function percentile_sorted

    subroutine write_histogram_pair(file_unit, figure_name, metric_name, values_a, values_i)
        integer, intent(in) :: file_unit
        character(len=*), intent(in) :: figure_name, metric_name
        real(dp), intent(in) :: values_a(:), values_i(:)
        real(dp) :: cutoff_a, cutoff_i, min_log, max_log
        real(dp) :: edges(histogram_edge_count)
        integer :: counts_a(histogram_edge_count-1), counts_i(histogram_edge_count-1)
        integer :: k

        cutoff_a = percentile_sorted(values_a, 99.9_dp)
        cutoff_i = percentile_sorted(values_i, 99.9_dp)
        min_log = huge(1.0_dp)
        max_log = -huge(1.0_dp)
        call find_log_bounds(values_a, cutoff_a, min_log, max_log)
        call find_log_bounds(values_i, cutoff_i, min_log, max_log)
        if (min_log > max_log) return
        if (min_log == max_log) max_log = min_log + 1.0e-12_dp
        do k = 1, histogram_edge_count
            edges(k) = min_log + real(k-1,dp)*(max_log-min_log)/ &
                       real(histogram_edge_count-1,dp)
        end do
        counts_a = 0
        counts_i = 0
        call count_log_bins(values_a, cutoff_a, min_log, max_log, counts_a)
        call count_log_bins(values_i, cutoff_i, min_log, max_log, counts_i)
        do k = 1, histogram_edge_count-1
            write(file_unit,'(a,",",a,",",a,",",es16.8,",",es16.8,",",i0)') &
                trim(figure_name), trim(metric_name), 'Exact Solution', edges(k), edges(k+1), counts_a(k)
            write(file_unit,'(a,",",a,",",a,",",es16.8,",",es16.8,",",i0)') &
                trim(figure_name), trim(metric_name), 'ivLam', edges(k), edges(k+1), counts_i(k)
        end do
    end subroutine write_histogram_pair

    subroutine find_log_bounds(values, cutoff, lower, upper)
        real(dp), intent(in) :: values(:), cutoff
        real(dp), intent(inout) :: lower, upper
        integer :: k
        do k = 1, size(values)
            if (values(k) <= 0.0_dp .or. values(k) >= cutoff) cycle
            lower = min(lower, log10(values(k)))
            upper = max(upper, log10(values(k)))
        end do
    end subroutine find_log_bounds

    subroutine count_log_bins(values, cutoff, lower, upper, counts)
        real(dp), intent(in) :: values(:), cutoff, lower, upper
        integer, intent(inout) :: counts(:)
        integer :: k, bin
        real(dp) :: log_value
        do k = 1, size(values)
            if (values(k) <= 0.0_dp .or. values(k) >= cutoff) cycle
            log_value = log10(values(k))
            bin = int((log_value-lower)/(upper-lower)*real(size(counts),dp)) + 1
            bin = max(1, min(size(counts), bin))
            counts(bin) = counts(bin) + 1
        end do
    end subroutine count_log_bins

end program benchmark_lambert