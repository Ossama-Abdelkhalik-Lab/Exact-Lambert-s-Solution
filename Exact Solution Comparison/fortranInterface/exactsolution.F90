#ifndef EXACT_STANDALONE
#include "fintrf.h"
#endif
module exact_mod_fast
    use omp_lib
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: sqrt2 = 1.4142135623730950488016887_dp
    real(dp), parameter :: c_tp = 0.4714045207910316829338962_dp  ! sqrt(2)/3
    logical :: initialized = .false.
    
    integer, parameter :: N1 = 16, N2 = 8
    complex(dp) :: cp_pts1(0:N1), cp_z_pts1(0:N1)
    complex(dp) :: K1(0:N1), C3_Sz1(0:N1)
    real(dp)    :: u2_re(0:N2), u2_im(0:N2)
    real(dp)    :: du2_re(0:N2), du2_im(0:N2)
    real(dp)    :: nan_val
contains

    subroutine init_cache()
        integer :: j
        real(dp) :: th1, th2, R1_circ, z0, W, zero
        complex(dp) :: z_pts1, sqrt_z1, Sz1_tmp, C_inv1_tmp
        
        z0 = 2.0_dp * pi**2
        R1_circ = 19.5_dp
        
        do j = 0, N1
            th1 = real(j, dp) * pi / real(N1, dp)
            z_pts1  = z0 + R1_circ * exp(cmplx(0.0_dp, th1, dp))
            cp_pts1(j) = cmplx(0.0_dp, R1_circ, dp) * exp(cmplx(0.0_dp, th1, dp))
            
            sqrt_z1 = sqrt(z_pts1)
            Sz1_tmp = (sqrt_z1 - sin(sqrt_z1)) / (z_pts1 * sqrt_z1)
            C_inv1_tmp = 1.0_dp / sqrt((1.0_dp - cos(sqrt_z1)) / z_pts1)
            K1(j) = (z_pts1 * Sz1_tmp - 1.0_dp) * C_inv1_tmp
            C3_Sz1(j) = (C_inv1_tmp * C_inv1_tmp * C_inv1_tmp) * Sz1_tmp
            
            W = 2.0_dp
            if (j == 0 .or. j == N1) W = 1.0_dp
            
            cp_pts1(j) = cp_pts1(j) * W
            cp_z_pts1(j) = cp_pts1(j) * z_pts1
        end do
        
        do j = 0, N2
            th2 = real(j, dp) * pi / real(N2, dp)
            u2_re(j) = cos(th2)
            u2_im(j) = sin(th2)
            
            W = 2.0_dp
            if (j == 0 .or. j == N2) W = 1.0_dp
            du2_re(j) = -sin(th2) * W
            du2_im(j) =  cos(th2) * W
        end do
        
        zero = 0.0_dp
        nan_val = zero / zero
        initialized = .true.
    end subroutine init_cache

    pure subroutine fast_cpx_sqrt(xr, xi, ur, ui)
        real(dp), intent(in)  :: xr, xi
        real(dp), intent(out) :: ur, ui
        real(dp) :: mag
        mag = sqrt(xr*xr + xi*xi)
        ur = sqrt(0.5_dp * (mag + xr))
        if (abs(ur) > 1.0e-18_dp) then
            ui = 0.5_dp * xi / ur
        else
            ui = 0.0_dp
        end if
    end subroutine fast_cpx_sqrt

    pure subroutine fast_cpx_trig(u, v, sin_w, cos_w)
        real(dp), intent(in)  :: u, v
        complex(dp), intent(out) :: sin_w, cos_w
        real(dp) :: sh, ch, sn, cs
        sn = sin(u); cs = cos(u)
        sh = sinh(v); ch = cosh(v)
        sin_w = cmplx(sn * ch, cs * sh, dp)
        cos_w = cmplx(cs * ch, -sn * sh, dp)
    end subroutine fast_cpx_trig

    subroutine exactsolution(R1, R2, dt, mu, H, V1, V2, Q)
        integer, intent(in) :: Q
        real(dp), intent(in) :: R1(3,Q), R2(3,Q), dt(Q), mu, H(3)
        real(dp), intent(out) :: V1(3,Q), V2(3,Q)
        
        integer :: k, j, num_procs
        real(dp) :: r1_mag, r2_mag, inv_r1, inv_r2, r12, inv_sqrt_mu, sqrt_mu, mu_t
        real(dp) :: c1, c2, c3, h_dot, sign_H, dot12, A
        real(dp) :: ry_z0, tp, z_coarse, z_fine, R2_rad, z
        real(dp) :: f, g, dg, inv_g, y_real, ry_real
        real(dp) :: t_eval, dt_dz, h_cpx, sz, Sz_real, rCz
        real(dp) :: den1, den2, re_d, im_d
        real(dp) :: zr, zi, sqr_r, sqr_i
        
        complex(dp) :: i1_1, i2_1, i1_2, i2_2, ryz1, Fz1, denom_cpx
        complex(dp) :: z_pts2, cp_pts2, sqrt_z2, Sz2_val, rCz2, ryz2, Fz2
        complex(dp) :: sin_w, cos_w, z_c, sqrt_zc, S_c, C_c, y_c, t_c
        complex(dp) :: ry_rC, y_over_C
        
        if (.not. initialized) call init_cache()
        
        sqrt_mu = sqrt(mu)
        inv_sqrt_mu = 1.0_dp / sqrt_mu
        h_cpx = 1.0e-8_dp
        
        !$omp parallel do default(shared) schedule(static, 4096) &
        !$omp private(k, j, r1_mag, r2_mag, inv_r1, inv_r2, r12, mu_t, c1, c2, c3, h_dot, sign_H, dot12, A, &
        !$omp         ry_z0, tp, i1_1, i2_1, ryz1, Fz1, z_coarse, R2_rad, i1_2, i2_2, &
        !$omp         zr, zi, sqr_r, sqr_i, z_pts2, cp_pts2, sqrt_z2, Sz2_val, rCz2, ryz2, Fz2, z_fine, &
        !$omp         sin_w, cos_w, z_c, sqrt_zc, S_c, C_c, y_c, t_c, ry_rC, y_over_C, &
        !$omp         z, t_eval, dt_dz, sz, Sz_real, rCz, y_real, ry_real, &
        !$omp         f, g, dg, inv_g, denom_cpx, re_d, im_d, den1, den2)
        do k = 1, Q
            r1_mag = sqrt(R1(1,k)*R1(1,k) + R1(2,k)*R1(2,k) + R1(3,k)*R1(3,k))
            r2_mag = sqrt(R2(1,k)*R2(1,k) + R2(2,k)*R2(2,k) + R2(3,k)*R2(3,k))
            inv_r1 = 1.0_dp / r1_mag
            inv_r2 = 1.0_dp / r2_mag
            r12 = r1_mag + r2_mag
            mu_t = dt(k) * sqrt_mu
            
            c1 = R1(2,k)*R2(3,k) - R1(3,k)*R2(2,k)
            c2 = R1(3,k)*R2(1,k) - R1(1,k)*R2(3,k)
            c3 = R1(1,k)*R2(2,k) - R1(2,k)*R2(1,k)
            h_dot = c1*H(1) + c2*H(2) + c3*H(3)
            
            sign_H = 1.0_dp
            if (h_dot < 0.0_dp) sign_H = -1.0_dp
            
            dot12 = R1(1,k)*R2(1,k) + R1(2,k)*R2(2,k) + R1(3,k)*R2(3,k)
            A = sign_H * sqrt(max(0.0_dp, r1_mag*r2_mag + dot12))
            
            ry_z0 = sqrt(max(0.0_dp, r12 - A * sqrt2))
            tp = ((ry_z0 * ry_z0 * ry_z0) * c_tp + A * ry_z0) * inv_sqrt_mu
            if (dt(k) < tp .or. abs(A) < 1.0e-12_dp) then
                V1(1:3, k) = nan_val
                V2(1:3, k) = nan_val
                cycle
            end if
            
            ! Contour 1: N1 = 16
            i1_1 = (0.0_dp, 0.0_dp)
            i2_1 = (0.0_dp, 0.0_dp)
            do j = 0, N1
                ryz1 = sqrt(r12 + K1(j)*A)
                denom_cpx = (ryz1 * ryz1 * ryz1) * C3_Sz1(j) + A * ryz1 - mu_t
                re_d = real(denom_cpx, dp)
                im_d = aimag(denom_cpx)
                den1 = 1.0_dp / (re_d*re_d + im_d*im_d + 1.0e-30_dp)
                Fz1 = cmplx(re_d * den1, -im_d * den1, dp)
                i1_1 = i1_1 + cp_z_pts1(j) * Fz1
                i2_1 = i2_1 + cp_pts1(j) * Fz1
            end do
            
            z_coarse = aimag(i1_1) / aimag(i2_1)
            z_coarse = max(0.1_dp, min(39.3784_dp, z_coarse))
            
            ! Contour 2: Fine Cauchy Circle (N2 = 8)
            R2_rad = max(1.0e-3_dp, min(2.0_dp, min(z_coarse - 1.0e-3_dp, 39.4774_dp - z_coarse)))
            i1_2 = (0.0_dp, 0.0_dp)
            i2_2 = (0.0_dp, 0.0_dp)
            do j = 0, N2
                zr = z_coarse + R2_rad * u2_re(j)
                zi = R2_rad * u2_im(j)
                z_pts2 = cmplx(zr, zi, dp)
                cp_pts2 = cmplx(R2_rad * du2_re(j), R2_rad * du2_im(j), dp)
                
                call fast_cpx_sqrt(zr, zi, sqr_r, sqr_i)
                sqrt_z2 = cmplx(sqr_r, sqr_i, dp)
                
                call fast_cpx_trig(sqr_r, sqr_i, sin_w, cos_w)
                Sz2_val = (sqrt_z2 - sin_w) / (z_pts2 * sqrt_z2)
                rCz2 = sqrt((1.0_dp - cos_w) / z_pts2)
                ryz2 = sqrt(r12 + A * ((z_pts2 * Sz2_val - 1.0_dp) / rCz2))
                
                ry_rC = ryz2 / rCz2
                denom_cpx = (ry_rC * ry_rC * ry_rC) * Sz2_val + A * ryz2 - mu_t
                re_d = real(denom_cpx, dp)
                im_d = aimag(denom_cpx)
                den2 = 1.0_dp / (re_d*re_d + im_d*im_d + 1.0e-30_dp)
                Fz2 = cmplx(re_d * den2, -im_d * den2, dp)
                
                i1_2 = i1_2 + cp_pts2 * z_pts2 * Fz2
                i2_2 = i2_2 + cp_pts2 * Fz2
            end do
            
            if (abs(aimag(i2_2)) < 1.0e-12_dp) then
                z_fine = z_coarse
            else
                z_fine = aimag(i1_2) / aimag(i2_2)
            end if
            if (z_fine /= z_fine) z_fine = z_coarse
            
            ! Complex-Step Polish Stage
            z_c = cmplx(z_fine, h_cpx, dp)
            call fast_cpx_sqrt(z_fine, h_cpx, sqr_r, sqr_i)
            sqrt_zc = cmplx(sqr_r, sqr_i, dp)
            
            call fast_cpx_trig(sqr_r, sqr_i, sin_w, cos_w)
            S_c = (sqrt_zc - sin_w) / (z_c * sqrt_zc)
            C_c = (1.0_dp - cos_w) / z_c
            y_c = r12 + A * (z_c * S_c - 1.0_dp) / sqrt(C_c)
            
            y_over_C = y_c / C_c
            t_c = (y_over_C * sqrt(y_over_C) * S_c + A * sqrt(y_c)) * inv_sqrt_mu
            
            t_eval = real(t_c, dp)
            dt_dz = aimag(t_c) / h_cpx
            if (abs(dt_dz) > 1.0e-14_dp) then
                z = z_fine - (t_eval - dt(k)) / dt_dz
            else
                z = z_fine
            end if
            z = max(1.0e-8_dp, z)
            
            ! Velocity Extraction
            sz = sqrt(z)
            Sz_real = (sz - sin(sz)) / (z * sz)
            rCz = sqrt(z / max(1.0e-14_dp, 1.0_dp - cos(sz)))
            y_real = max(0.0_dp, r12 + A * (z * Sz_real - 1.0_dp) * rCz)
            ry_real = sqrt(y_real)
            
            f = 1.0_dp - y_real * inv_r1
            g = A * ry_real * inv_sqrt_mu
            dg = 1.0_dp - y_real * inv_r2
            inv_g = 1.0_dp / g
            
            V1(1,k) = (R2(1,k) - f*R1(1,k)) * inv_g
            V1(2,k) = (R2(2,k) - f*R1(2,k)) * inv_g
            V1(3,k) = (R2(3,k) - f*R1(3,k)) * inv_g
            V2(1,k) = (dg*R2(1,k) - R1(1,k)) * inv_g
            V2(2,k) = (dg*R2(2,k) - R1(2,k)) * inv_g
            V2(3,k) = (dg*R2(3,k) - R1(3,k)) * inv_g
        end do
    end subroutine exactsolution
end module exact_mod_fast

#ifndef EXACT_STANDALONE
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use exact_mod_fast
    use iso_c_binding
    implicit none
    
    integer, intent(in) :: nlhs, nrhs
    mwPointer, intent(in) :: prhs(*)
    mwPointer, intent(out) :: plhs(*)
    mwPointer :: mxGetPr, mxCreateDoubleMatrix
    mwSize :: mxGetN, Q
    real(dp), pointer :: R1(:,:), R2(:,:), dt(:), mu_array(:), H(:), V1(:,:), V2(:,:)
    
    if (nrhs /= 5) call mexErrMsgTxt("5 inputs required")
    Q = mxGetN(prhs(1))
    
    plhs(1) = mxCreateDoubleMatrix(3, Q, 0)
    plhs(2) = mxCreateDoubleMatrix(3, Q, 0)
    
    call c_f_pointer(transfer(mxGetPr(prhs(1)), c_null_ptr), R1, [3, Q])
    call c_f_pointer(transfer(mxGetPr(prhs(2)), c_null_ptr), R2, [3, Q])
    call c_f_pointer(transfer(mxGetPr(prhs(3)), c_null_ptr), dt, [Q])
    call c_f_pointer(transfer(mxGetPr(prhs(4)), c_null_ptr), mu_array, [1])
    call c_f_pointer(transfer(mxGetPr(prhs(5)), c_null_ptr), H, [3])
    call c_f_pointer(transfer(mxGetPr(plhs(1)), c_null_ptr), V1, [3, Q])
    call c_f_pointer(transfer(mxGetPr(plhs(2)), c_null_ptr), V2, [3, Q])
    
    call exactsolution(R1, R2, dt, mu_array(1), H, V1, V2, Q)
end subroutine mexFunction
#endif