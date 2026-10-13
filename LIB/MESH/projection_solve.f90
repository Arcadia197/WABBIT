!> \brief Poisson solve of a projection, exact on adaptive grids
!> \details A projection solves L phi = b with b = D u and then sets u' = u - G phi. On adaptive grids the ghost-node
!! syncs do not commute with the stencils at coarse-fine interfaces, so L differs from the composition A = D G in a strip
!! of a few points along the interfaces and the divergence D u' = (b - L phi) + (L - A) phi stalls at (L - A) phi, no
!! matter how well the V-cycles converge. This routine solves A phi = b instead, so that D u' = b - A phi is exactly
!! the residual of the solver.
!!
!! The unchanged multigrid V-cycle on L stays the work horse:
!!   - poisson_exact_projection=.false., or an equidistant grid with the caller using the composite Poisson pair
!!     (then A = L exactly): plain multigrid_solve, the current solver
!!   - otherwise: classic V-cycles on L as long as each reduces max|b - A phi| by at least the factor
!!     poisson_exact_switch (they stall once phi reaches the incompatibility), then FGMRES(m) on A with one V-cycle
!!     per iteration as preconditioner, each followed by poisson_exact_sweeps Jacobi sweeps on A (FGMRES(1) is the
!!     minimal residual step). The V-cycle returns corrections of the right shape but A sees only part of them near
!!     block corners; FGMRES supplies the step length, the sweeps correct the interface rows locally.
!!
!! No memory is allocated: all vectors are slots of hvy_v (hvy_tmp in the callers), which is also the array the
!! V-cycle moves when it balances the load. Slots: phi, b (overwritten by the residual), work (V-cycle work and the
!! gradient scratch of A) and the free slots i_free. FGMRES(m) needs 2m free slots (v_2..v_{m+1}, z_1..z_m), so
!! m = min(poisson_exact_m, n_free/2). With less than two free slots, the current solver is used (with a warning).
!!
!! Divergence and gradient of the caller: D = D_L S_v (sync_ghosts_tree, or sync_ghosts_RHS_tree with g_RHS layers if
!! sync_rhs) and G = D_R S_phi (sync with poisson_order_predictor), with the FD1 stencil pair fd_order. On return, phi
!! is synced with S_phi on the leaf grid, as after multigrid_solve, and hvy_v(:,:,:,i_b,:) holds the divergence that
!! u - G phi will have. Details and test cases: report "Exact Helmholtz projection on adaptive grids" (10/2026).
!!
!! Variants for analysis (wabbit-post --proto-projection-krylov): poisson_exact_switch < 0 skips the classic V-cycles
!! (pure FGMRES), poisson_exact_m = 0 replaces FGMRES by defect correction (step length 1). With the optional hook
!! hist, the solver calls hist after every V-cycle with the current iterate and its residual b0 - A x (b0 = the
!! original b, kept by the caller in slot i_b0), and at the end of each FGMRES cycle with its Hessenberg matrix. The
!! iterate inside an FGMRES cycle is formed in slot i_hist, the residual in slot i_hres. timing returns the time
!! spent in the V-cycles and in A on this rank.
subroutine projection_solve(params, hvy_v, i_phi, i_b, i_work, i_free, tree_ID, fd_order, sync_rhs, init_0, time, &
    hist, i_b0, i_hist, i_hres, timing)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    !> physically allocated array holding all vectors of the solve (one slot each), also used for load balancing
    real(kind=rk), intent(inout)       :: hvy_v(:, :, :, :, :)
    !> slots: solution phi (initial guess unless init_0), right-hand side b = D u (destroyed), work
    integer(kind=ik), intent(in)       :: i_phi, i_b, i_work
    !> slots that are free during the solve, used for the Krylov vectors
    integer(kind=ik), intent(in)       :: i_free(:)
    integer(kind=ik), intent(in)       :: tree_ID
    !> FD1 stencil pair of the caller's divergence D and gradient G
    character(len=*), intent(in)       :: fd_order
    !> the caller syncs the vector field with sync_ghosts_RHS_tree (RK stage) instead of sync_ghosts_tree
    logical, intent(in)                :: sync_rhs
    logical, intent(in)                :: init_0     !< start from phi = 0
    real(kind=rk), intent(in)          :: time       !< only for logging
    !> analysis only: called after every V-cycle with the slots of the iterate x and of b0 - A x, and with Hbar at the
    !> end of each FGMRES cycle (then x is the iterate after its last step)
    interface
        subroutine hist(params, hvy_v, i_x, i_r, n_vc, Hbar)
            import :: type_params, rk, ik
            type (type_params), intent(inout)     :: params
            real(kind=rk), intent(inout)          :: hvy_v(:, :, :, :, :)
            integer(kind=ik), intent(in)          :: i_x, i_r, n_vc
            real(kind=rk), intent(in), optional   :: Hbar(:, :)
        end subroutine hist
    end interface
    optional                           :: hist
    !> analysis only, required with hist: slots of the original b (kept), the iterate and its residual (scratch)
    integer(kind=ik), intent(in), optional :: i_b0, i_hist, i_hres
    !> analysis only: time spent in the V-cycles, in A and in the history calls (incl. hist) on this rank
    real(kind=rk), intent(out), optional   :: timing(1:3)

    real(kind=rk), allocatable         :: FD1_l(:), FD1_r(:)
    integer(kind=ik)                   :: FD1_ls, FD1_le, FD1_rs, FD1_re
    integer(kind=ik)                   :: m, n_vc, n_vc_max, i, k, hvy_id, lgt_id, is(1:3), ie(1:3)
    integer(kind=ik), allocatable      :: i_v(:), i_z(:)
    integer(kind=ik), allocatable, save :: ref_backup(:)
    real(kind=rk)                      :: c0, b_norm(1:1), r_Linf, r_L2, r_L2_old, d_prev, mean_b, x0(1:3), dx(1:3), mean_phi
    real(kind=rk)                      :: resid(1:3), t_start, t_block, t_vc, t_A, t_hist
    logical                            :: is_equidistant, done, use_tolerance, dc
    logical, save                      :: scheme_printed = .false.
    character(len=cshort)              :: fname

    !---------------------------------------------------------------------------------------------------------------
    ! current solver where it is exact: A = L on equidistant grids if the caller uses the Poisson stencil pair and L
    ! is its composite stencil (the *_comp_* orders; for central orders L is the narrow FD2 stencil, not D_L * D_R).
    ! With less than two free slots there is no room for FGMRES, the current solver is used as well
    if (present(timing)) timing = 0.0_rk
    is_equidistant = (minActiveLevel_tree(tree_ID) == maxActiveLevel_tree(tree_ID))
    if (params%poisson_exact_projection .and. size(i_free) < 2 .and. params%rank == 0 .and. .not. scheme_printed) then
        write(*, '(A)') "WARNING projection_solve: less than 2 free slots in hvy_tmp, using the current solver (not exact)"
        scheme_printed = .true.
    endif
    if (.not. params%poisson_exact_projection .or. size(i_free) < 2 .or. (is_equidistant .and. .not. present(hist) .and. &
        trim(fd_order) == trim(params%poisson_order) .and. index(params%poisson_order, "_comp_") > 0)) then
        call multigrid_solve(params, hvy_v(:,:,:,i_phi:i_phi,:), hvy_v(:,:,:,i_b:i_b,:), hvy_v(:,:,:,i_work:i_work,:), &
            tree_ID, init_0=init_0, verbose=.false., hvy_full=hvy_v, time=time)
        return
    endif

    if (params%poisson_order == "FD_6th_mehrstellen") then
        call abort(261009, "projection_solve: Mehrstellen has no D/G stencil pair, use poisson_exact_projection=0")
    endif
    if (.not. all(params%periodic_BC(1:params%dim))) then
        ! the gradient components of A are synced one at a time, i.e. symmetry conditions would treat all of them as
        ! the first component of a vector
        call abort(261009, "projection_solve: the exact projection is implemented for periodic boundaries only")
    endif

    t_start = MPI_Wtime()
    t_vc = 0.0_rk
    t_A = 0.0_rk
    t_hist = 0.0_rk

    ! Krylov vectors: v_1 is the residual in slot i_b, v_2..v_{m+1} and z_1..z_m come from the free slots.
    ! Defect correction (m=0) uses z_1 and v_2 as sweep scratch, the slots of FGMRES(1)
    dc = (params%poisson_exact_m == 0)
    m = max(1, min(params%poisson_exact_m, size(i_free) / 2))
    allocate(i_v(1:m+1), i_z(1:m))
    i_v(1) = i_b
    do i = 1, m
        i_v(i+1) = i_free(i)
        i_z(i) = i_free(m+i)
    enddo

    if (params%rank == 0 .and. .not. scheme_printed) then
        if (dc) then
            write(fname, '(A)') "defect correction"
        else
            write(fname, '(A, i0, A)') "FGMRES(", m, ")"
        endif
        write(*, '(A, f5.2, A, A, A, i0, A, i0, A)') "projection_solve: classic V-cycles while they reduce the divergence by ", &
            params%poisson_exact_switch, ", then ", trim(fname), " + ", params%poisson_exact_sweeps, " Jacobi sweeps on A (", &
            size(i_free), " free slots)"
        scheme_printed = .true.
    endif

    ! stencils of D and G, and the diagonal of the composite stencil D_L * D_R per direction (for the Jacobi sweeps)
    call setup_FD1_left_stencil(fd_order, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(fd_order, FD1_r, FD1_rs, FD1_re)
    c0 = 0.0_rk
    do k = FD1_ls, FD1_le
        if (-k >= FD1_rs .and. -k <= FD1_re) c0 = c0 + FD1_l(k) * FD1_r(-k)
    enddo

    ! interior points of a block
    is = 1
    ie = 1
    is(1:params%dim) = params%g + 1
    ie(1:params%dim) = params%g + params%Bs(1:params%dim)

    use_tolerance = (params%poisson_cycle_end_criteria == "tolerance")
    if (params%poisson_cycle_end_criteria == "fixed_iterations") then
        n_vc_max = params%poisson_cycle_it
    elseif (use_tolerance) then
        n_vc_max = params%poisson_cycle_max_it
    else
        call abort(250903, "Don't know how to stop! Please choose between 'fixed_iterations' and 'tolerance', thank you :)")
    endif

    ! norm of b for relative tolerances, as in multigrid_solve
    call componentWiseNorm_tree(params, hvy_v(:,:,:,i_b:i_b,:), tree_ID, "Linfty", b_norm, threshold_state_vector=.false.)

    !---------------------------------------------------------------------------------------------------------------
    ! full tree for the V-cycles, as in multigrid_solve. The V-cycle uses the refinement status for its own
    ! purposes, the one of the caller is restored at the end
    if (allocated(ref_backup)) then
        if (size(ref_backup) /= size(lgt_block,1)) deallocate(ref_backup)
    endif
    if (.not. allocated(ref_backup)) allocate(ref_backup(1:size(lgt_block,1)))
    ref_backup(:) = 0
    do k = 1, lgt_n(tree_ID)
        ref_backup(lgt_active(k, tree_ID)) = lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS)
    enddo

    call init_full_tree(params, tree_ID, Jmin_set=params%poisson_Jmin, set_ref=0)
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        if (.not. block_is_leaf(params, hvy_id) .or. init_0) hvy_v(:,:,:,i_phi,hvy_id) = 0.0_rk
    enddo

    !---------------------------------------------------------------------------------------------------------------
    ! classic V-cycles on L phi = b, as long as they pay: they reduce the divergence b - A phi quickly until the
    ! iterate reaches the incompatibility (L - A) phi, then not at all. After each cycle, d = b - A phi goes into the
    ! first free slot, and the iteration switches to A once a cycle reduced max|d| by less than poisson_exact_switch.
    ! The V-cycle removes the mean of its right-hand side (harmless for L). On adaptive grids D is not conservative
    ! and a constant shift changes the problem for A, so the shift is undone. poisson_exact_switch < 0: no classic cycles
    call leaf_copy(i_b, i_free(1))
    call apply_A(i_phi, i_free(1), -1.0_rk, accumulate=.true.)
    call leaf_norms(i_free(1), r_Linf, r_L2)
    n_vc = 0
    done = converged(r_Linf)
    do while (.not. done .and. n_vc < n_vc_max .and. params%poisson_exact_switch >= 0.0_rk)
        mean_b = leaf_mean(i_b)
        t_block = MPI_Wtime()
        call multigrid_vcycle(params, hvy_v(:,:,:,i_phi:i_phi,:), hvy_v(:,:,:,i_b:i_b,:), hvy_v(:,:,:,i_work:i_work,:), &
            tree_ID, residual_out=resid, verbose=.false., hvy_full=hvy_v)
        t_vc = t_vc + MPI_Wtime() - t_block
        call leaf_add_const(i_b, mean_b - leaf_mean(i_b))
        n_vc = n_vc + 1
        if (present(hist)) call call_hist(i_phi)

        d_prev = r_Linf
        call leaf_copy(i_b, i_free(1))
        call apply_A(i_phi, i_free(1), -1.0_rk, accumulate=.true.)
        call leaf_norms(i_free(1), r_Linf, r_L2)
        done = converged(r_Linf)
        if (r_Linf > params%poisson_exact_switch * d_prev) exit
    enddo

    ! from here on, slot i_b holds the residual r = b - A phi, i.e. the divergence u - G phi would have
    call leaf_copy(i_free(1), i_b)
    done = done .or. (n_vc >= n_vc_max)

    !---------------------------------------------------------------------------------------------------------------
    ! outer iteration on A phi = b, one preconditioning V-cycle per step
    do while (.not. done)
        r_L2_old = r_L2
        if (dc) then
            call dc_step()
        else
            call fgmres_cycle()
        endif
        call leaf_norms(i_b, r_Linf, r_L2)

        if (.not. (r_Linf <= LIM_DIVERGED)) then   ! also catches NaN
            write(fname, '(A, i0, A, es10.3)') "Poisson solver diverged, aborting at it ", n_vc, ", Max Res Linfty: ", r_Linf
            open (77, file='NSPP_diverged', status='replace')
            close(77)
            call abort(260121, fname)
        endif

        ! FGMRES minimizes the L2 norm, so it cannot increase: once it stagnates, the least-squares floor (the part of
        ! b that no gradient can remove, as D is not conservative at interfaces) is reached
        done = (n_vc >= n_vc_max) .or. converged(r_Linf) .or. (use_tolerance .and. r_L2 > 0.9_rk * r_L2_old)
    enddo

    ! constants are in the kernel of G, the mean of phi is arbitrary. Never remove the mean of the residual
    mean_phi = leaf_mean(i_phi)
    call leaf_add_const(i_phi, -mean_phi)

    if (params%rank == 0) then
        write(*, '(A, i0, A, es10.3, A, es10.3)') "   Final Residual after ", n_vc, " it (exact projection), Linfty: ", &
            r_Linf, " , rel to RHS: ", r_Linf / max(b_norm(1), tiny(1.0_rk))
        call append_t_file('poisson_solver.t', (/ time, dble(n_vc), r_Linf, r_Linf / max(b_norm(1), tiny(1.0_rk)) /))
    endif

    call prune_fulltree2leafs(params, tree_ID)
    do k = 1, lgt_n(tree_ID)
        lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) = ref_backup(lgt_active(k, tree_ID))
    enddo

    ! the caller computes G phi with the ghost values of phi
    call sync_poisson(i_phi)

    call toc("projection_solve: apply A", 10060, t_A)
    call toc("projection_solve", 10061, MPI_Wtime() - t_start)
    if (present(timing)) timing = (/ t_vc, t_A, t_hist /)

contains

    !===============================================================================================================
    !> stop criterion of the parameters, applied to max|b - A phi|
    logical function converged(r_max)
        real(kind=rk), intent(in) :: r_max
        converged = .false.
        if (use_tolerance) then
            converged = (r_max < params%poisson_cycle_tol_abs) .or. &
                (r_max / b_norm(1) < params%poisson_cycle_tol_rel .and. b_norm(1) > 1.0e-15_rk)
        endif
    end function converged

    !> one restart cycle of FGMRES(m), the residual is in v_1 = slot i_b on entry and on exit.
    !> A Z_j = V_{j+1} Hbar_j, y minimizes |beta e_1 - Hbar_j y|, phi += Z_j y. The new residual is formed from the
    !> basis, r = V_{j+1} (beta e_1 - Hbar_j y), so b is never needed again
    subroutine fgmres_cycle()
        real(kind=rk) :: beta, H(1:m+1, 1:m), Hbar(1:m+1, 1:m), cs(1:m), sn(1:m), gv(1:m+1), y(1:m), h1(1:m+1), h2(1:m+1)
        real(kind=rk) :: tmp, denom, c(1:m+1)
        integer(kind=ik) :: j, jend, l
        logical :: breakdown, last

        beta = r_L2
        call leaf_scale(i_v(1), 1.0_rk / beta)
        gv = 0.0_rk
        gv(1) = beta
        H = 0.0_rk
        Hbar = 0.0_rk
        jend = 0
        do j = 1, m
            ! z_j = V(v_j) + sweeps, v_{j+1} is free until A z_j is written into it and serves as scratch
            call precondition(i_z(j), i_v(j), i_v(j+1))
            n_vc = n_vc + 1
            call apply_A(i_z(j), i_v(j+1), 1.0_rk, accumulate=.false.)

            ! classical Gram-Schmidt with one re-orthogonalization (CGS2)
            call leaf_dots(i_v(j+1), i_v(1:j), h1(1:j))
            do l = 1, j
                call leaf_axpy(-h1(l), i_v(l), i_v(j+1))
            enddo
            call leaf_dots(i_v(j+1), i_v(1:j), h2(1:j))
            do l = 1, j
                call leaf_axpy(-h2(l), i_v(l), i_v(j+1))
            enddo
            h1(1:j) = h1(1:j) + h2(1:j)
            call componentWiseNorm_tree(params, hvy_v(:,:,:,i_v(j+1):i_v(j+1),:), tree_ID, "L2", h1(j+1:j+1), &
                threshold_state_vector=.false.)
            jend = j
            Hbar(1:j+1, j) = h1(1:j+1)
            H(1:j+1, j) = h1(1:j+1)
            ! breakdown (h1(j+1) ~ 0): the residual is in the Krylov space, v_{j+1} is not needed
            breakdown = (h1(j+1) <= 1.0e-14_rk * beta)
            if (.not. breakdown) call leaf_scale(i_v(j+1), 1.0_rk / h1(j+1))

            ! Givens rotations: H becomes upper triangular, |gv(j+1)| is the residual norm after this step
            do l = 1, j-1
                tmp      =  cs(l)*H(l,j) + sn(l)*H(l+1,j)
                H(l+1,j) = -sn(l)*H(l,j) + cs(l)*H(l+1,j)
                H(l,j)   = tmp
            enddo
            denom = max(sqrt(H(j,j)**2 + H(j+1,j)**2), tiny(1.0_rk))
            cs(j) = H(j,j) / denom
            sn(j) = H(j+1,j) / denom
            H(j,j) = denom
            H(j+1,j) = 0.0_rk
            gv(j+1) = -sn(j)*gv(j)
            gv(j)   =  cs(j)*gv(j)

            last = (j == m .or. breakdown .or. n_vc >= n_vc_max)
            if (present(hist)) then
                ! iterate after this step, x = phi + Z_j y_j
                call back_substitute(H, gv, j, y)
                call leaf_copy(i_phi, i_hist)
                do l = 1, j
                    call leaf_axpy(y(l), i_z(l), i_hist)
                enddo
                if (last) then
                    call call_hist(i_hist, Hbar(1:j+1, 1:j))
                else
                    call call_hist(i_hist)
                endif
            endif
            if (last) exit
        enddo

        ! back substitution H y = gv and update phi += Z y
        call back_substitute(H, gv, jend, y)
        do l = 1, jend
            call leaf_axpy(y(l), i_z(l), i_phi)
        enddo

        ! r = V_{jend+1} c with c = beta e_1 - Hbar y, built in place in v_1
        c = 0.0_rk
        c(1) = beta
        do l = 1, jend+1
            c(l) = c(l) - sum(Hbar(l, 1:jend) * y(1:jend))
        enddo
        if (breakdown) c(jend+1) = 0.0_rk   ! v_{jend+1} was not normalized
        call leaf_scale(i_v(1), c(1))
        do l = 2, jend+1
            call leaf_axpy(c(l), i_v(l), i_v(1))
        enddo
    end subroutine fgmres_cycle

    !> back substitution R y = g for the leading j x j block of the rotated Hessenberg matrix
    subroutine back_substitute(R, g, j, y)
        real(kind=rk), intent(in)    :: R(:,:), g(:)
        integer(kind=ik), intent(in) :: j
        real(kind=rk), intent(inout) :: y(:)
        integer(kind=ik) :: l
        do l = j, 1, -1
            y(l) = (g(l) - sum(R(l, l+1:j) * y(l+1:j))) / R(l, l)
        enddo
    end subroutine back_substitute

    !> one step of defect correction (analysis only, poisson_exact_m=0): phi += z, r -= A z with z = V(r) + sweeps
    subroutine dc_step()
        call precondition(i_z(1), i_b, i_v(2))
        n_vc = n_vc + 1
        call leaf_axpy(1.0_rk, i_z(1), i_phi)
        call apply_A(i_z(1), i_b, -1.0_rk, accumulate=.true.)
        if (present(hist)) call call_hist(i_phi)
    end subroutine dc_step

    !> analysis only: residual b0 - A x into slot i_hres, then the caller's hook. Its time is not counted in t_A
    subroutine call_hist(x, Hbar)
        integer(kind=ik), intent(in)        :: x
        real(kind=rk), intent(in), optional :: Hbar(:,:)
        real(kind=rk) :: t0, t_A_backup
        t0 = MPI_Wtime()
        t_A_backup = t_A
        call leaf_copy(i_b0, i_hres)
        call apply_A(x, i_hres, -1.0_rk, accumulate=.true.)
        t_A = t_A_backup
        call hist(params, hvy_v, x, i_hres, n_vc, Hbar)
        t_hist = t_hist + MPI_Wtime() - t0
    end subroutine call_hist

    !> z = V(v): one unchanged V-cycle on L with zero initial guess and v as right-hand side. The V-cycle gives v back
    !> on the leaf blocks up to a constant (its mean removal), which is undone. Then Jacobi sweeps on A with the
    !> diagonal of L, z += (v - A z) / diag(L), with s as scratch: the V-cycle inverted L, the sweeps correct the
    !> rows next to interfaces where A /= L (for smooth residuals they change almost nothing, so no mask is needed)
    subroutine precondition(z, v, s)
        integer(kind=ik), intent(in) :: z, v, s
        integer(kind=ik) :: kk, ks
        real(kind=rk) :: mean_v, dinv

        mean_v = leaf_mean(v)
        do kk = 1, hvy_n(tree_ID)
            hvy_v(:,:,:,z,hvy_active(kk, tree_ID)) = 0.0_rk
        enddo
        t_block = MPI_Wtime()
        call multigrid_vcycle(params, hvy_v(:,:,:,z:z,:), hvy_v(:,:,:,v:v,:), hvy_v(:,:,:,i_work:i_work,:), &
            tree_ID, residual_out=resid, verbose=.false., hvy_full=hvy_v)
        t_vc = t_vc + MPI_Wtime() - t_block
        call leaf_add_const(v, mean_v - leaf_mean(v))

        do ks = 1, params%poisson_exact_sweeps
            call leaf_copy(v, s)
            call apply_A(z, s, -1.0_rk, accumulate=.true.)
            ! z = z + s / diag(L) on the leaf blocks, diag(L) = c0 * sum_i 1/dx_i^2
            do kk = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(kk, tree_ID)
                if (.not. block_is_leaf(params, hvy_id)) cycle
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)
                dinv = 1.0_rk / (c0 * sum(1.0_rk / dx(1:params%dim)**2))
                hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),z,hvy_id) = hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),z,hvy_id) &
                    + dinv * hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),s,hvy_id)
            enddo
        enddo
    end subroutine precondition

    !===============================================================================================================
    !> out = factor * A x, or out = out + factor * A x (accumulate), with A = sum_i D_L,i S_v D_R,i S_phi.
    !> One direction at a time: D_R,i x into the work slot, sync it as the caller syncs the velocity, add D_L,i of it
    subroutine apply_A(x, out, factor, accumulate)
        integer(kind=ik), intent(in) :: x, out
        real(kind=rk), intent(in)    :: factor
        logical, intent(in)          :: accumulate
        integer(kind=ik) :: idir, kk
        real(kind=rk) :: t0

        t0 = MPI_Wtime()
        ! on the full tree only leaf blocks may send and receive ghost values: leaves get refinement status 0,
        ! mother blocks REF_TMP_EMPTY
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            if (block_is_leaf(params, hvy_id)) then
                lgt_block(lgt_id, IDX_REFINE_STS) = 0
            else
                lgt_block(lgt_id, IDX_REFINE_STS) = REF_TMP_EMPTY
            endif
        enddo
        call synchronize_lgt_data(params, refinement_status_only=.true.)
        call sync_poisson(x)

        if (.not. accumulate) then
            do kk = 1, hvy_n(tree_ID)
                hvy_v(:,:,:,out,hvy_active(kk, tree_ID)) = 0.0_rk
            enddo
        endif

        do idir = 1, params%dim
            do kk = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(kk, tree_ID)
                if (.not. block_is_leaf(params, hvy_id)) cycle
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)
                hvy_v(:,:,:,i_work,hvy_id) = 0.0_rk
                call add_fd1(hvy_v(:,:,:,x,hvy_id), hvy_v(:,:,:,i_work,hvy_id), FD1_r, FD1_rs, FD1_re, idir, 1.0_rk/dx(idir))
            enddo

            if (sync_rhs) then
                call sync_ghosts_RHS_tree(params, hvy_v(:,:,:,i_work:i_work,:), tree_ID, g_minus=params%g_RHS, g_plus=params%g_RHS)
            else
                call sync_ghosts_tree(params, hvy_v(:,:,:,i_work:i_work,:), tree_ID)
            endif

            do kk = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(kk, tree_ID)
                if (.not. block_is_leaf(params, hvy_id)) cycle
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)
                call add_fd1(hvy_v(:,:,:,i_work,hvy_id), hvy_v(:,:,:,out,hvy_id), FD1_l, FD1_ls, FD1_le, idir, factor/dx(idir))
            enddo
        enddo
        t_A = t_A + MPI_Wtime() - t0
    end subroutine apply_A

    !> out(interior) = out(interior) + scale * sum_k c(k) f(. + k e_idir)
    subroutine add_fd1(f, out, c, cs_, ce_, idir, scale)
        real(kind=rk), intent(in)    :: f(:,:,:)
        real(kind=rk), intent(inout) :: out(:,:,:)
        integer(kind=ik), intent(in) :: cs_, ce_, idir
        real(kind=rk), intent(in)    :: c(cs_:)
        real(kind=rk), intent(in)    :: scale
        integer(kind=ik) :: kk, sh(1:3)
        do kk = cs_, ce_
            sh = 0
            sh(idir) = kk
            out(is(1):ie(1),is(2):ie(2),is(3):ie(3)) = out(is(1):ie(1),is(2):ie(2),is(3):ie(3)) &
                + scale * c(kk) * f(is(1)+sh(1):ie(1)+sh(1), is(2)+sh(2):ie(2)+sh(2), is(3)+sh(3):ie(3)+sh(3))
        enddo
    end subroutine add_fd1

    !> sync phi with the Poisson predictor (S_phi), as the multigrid solver does
    subroutine sync_poisson(x)
        integer(kind=ik), intent(in) :: x
        character(len=cshort) :: order_backup
        order_backup = params%order_predictor
        params%order_predictor = params%poisson_order_predictor
        call sync_ghosts_tree(params, hvy_v(:,:,:,x:x,:), tree_ID)
        params%order_predictor = order_backup
    end subroutine sync_poisson

    !===============================================================================================================
    ! vector operations on the interior points of the leaf blocks, inner products weighted with the cell volume

    !> h(l) = <w, v_l> for all slots in v(:), one allreduce
    subroutine leaf_dots(w, v, h)
        integer(kind=ik), intent(in) :: w, v(:)
        real(kind=rk), intent(out)   :: h(:)
        integer(kind=ik) :: kk, l, mpierr
        h = 0.0_rk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (.not. block_is_leaf(params, hvy_id)) cycle
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            do l = 1, size(v)
                h(l) = h(l) + product(dx(1:params%dim)) * sum(hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),w,hvy_id) &
                    * hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),v(l),hvy_id))
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, h, size(v), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
    end subroutine leaf_dots

    !> maximum norm and (volume weighted) L2 norm on the leaf blocks
    subroutine leaf_norms(a, n_Linf, n_L2)
        integer(kind=ik), intent(in) :: a
        real(kind=rk), intent(out)   :: n_Linf, n_L2
        real(kind=rk) :: nrm(1:1)
        call componentWiseNorm_tree(params, hvy_v(:,:,:,a:a,:), tree_ID, "Linfty", nrm, threshold_state_vector=.false.)
        n_Linf = nrm(1)
        call componentWiseNorm_tree(params, hvy_v(:,:,:,a:a,:), tree_ID, "L2", nrm, threshold_state_vector=.false.)
        n_L2 = nrm(1)
    end subroutine leaf_norms

    !> volume weighted mean on the leaf blocks
    function leaf_mean(a) result(mean)
        integer(kind=ik), intent(in) :: a
        real(kind=rk) :: mean, nrm(1:1)
        call componentWiseNorm_tree(params, hvy_v(:,:,:,a:a,:), tree_ID, "Mean", nrm, threshold_state_vector=.false.)
        mean = nrm(1)
    end function leaf_mean

    !> y = y + alpha x
    subroutine leaf_axpy(alpha, x, y)
        real(kind=rk), intent(in)    :: alpha
        integer(kind=ik), intent(in) :: x, y
        integer(kind=ik) :: kk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (block_is_leaf(params, hvy_id)) hvy_v(:,:,:,y,hvy_id) = hvy_v(:,:,:,y,hvy_id) + alpha * hvy_v(:,:,:,x,hvy_id)
        enddo
    end subroutine leaf_axpy

    !> x = alpha x
    subroutine leaf_scale(x, alpha)
        integer(kind=ik), intent(in) :: x
        real(kind=rk), intent(in)    :: alpha
        integer(kind=ik) :: kk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (block_is_leaf(params, hvy_id)) hvy_v(:,:,:,x,hvy_id) = alpha * hvy_v(:,:,:,x,hvy_id)
        enddo
    end subroutine leaf_scale

    !> y = x
    subroutine leaf_copy(x, y)
        integer(kind=ik), intent(in) :: x, y
        integer(kind=ik) :: kk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (block_is_leaf(params, hvy_id)) hvy_v(:,:,:,y,hvy_id) = hvy_v(:,:,:,x,hvy_id)
        enddo
    end subroutine leaf_copy

    !> x = x + const
    subroutine leaf_add_const(x, const)
        integer(kind=ik), intent(in) :: x
        real(kind=rk), intent(in)    :: const
        integer(kind=ik) :: kk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (block_is_leaf(params, hvy_id)) hvy_v(:,:,:,x,hvy_id) = hvy_v(:,:,:,x,hvy_id) + const
        enddo
    end subroutine leaf_add_const

end subroutine projection_solve
