!> \brief Prototype driver: convergence history of the exact Helmholtz projection on adaptive grids
!> \details
!! On adaptive grids the Laplacian L used by the multigrid solver differs from the composition A = D*S*G
!! (divergence of the ghost-synced gradient) in a few points next to coarse-fine interfaces. Solving
!! L phi = D u therefore leaves a divergence D(u - G phi) = (L - A) phi behind. projection_solve solves A phi = D u;
!! this driver runs it with the settings of each method and records the divergence left after every V-cycle:
!!
!!   L           - classic V-cycles on L phi = b only (poisson_exact_switch = 1e30)
!!   richardson  - defect correction on A phi = b, one V-cycle (on L) as preconditioner (poisson_exact_m = 0)
!!   fgmres      - restarted FGMRES(m) on A phi = b, one V-cycle per iteration (poisson_exact_switch < 0)
!!   hybrid      - classic V-cycles while they pay, then FGMRES(m) (the production scheme): hybrid_rule=stall
!!                 switches once a V-cycle reduces max|b - A phi| by less than hybrid_stall, hybrid_rule=kappa
!!                 (hybrid_kappa=4) after the first V-cycle
!!
!! ismooth=k: k Jacobi sweeps on A after every preconditioning V-cycle (poisson_exact_sweeps).
!! The Hessenberg matrices of FGMRES are written to projection_hessenberg.csv (Ritz values of A*V).
!!
!! Slots are scalar leaf fields (one hvy slot each) needed in production, the V-cycle work slot included.
!! Diagnostic slots (the original b, the iterate inside an FGMRES cycle and its residual) are not counted: they exist
!! only to report the true divergence after each V-cycle. Since D(u - G phi) = b - A phi exactly, the A residual is
!! the divergence left after the projection. After each method, the projection is repeated with the standard
!! routines on the pruned leaf grid (subtract_scalar_gradient, sync, compute_divergence) to verify that identity
!! independently.
!-----------------------------------------------------------------------------------------------------
module module_proto_projection
    use mpi
    use module_globals
    use module_params
    use module_mesh
    use module_poisson
    use module_forestMetaData

    implicit none

    ! slot layout in hvy_tmp
    integer(kind=ik), parameter :: S_B = 1   !< b = D u (diagnostic, kept)
    integer(kind=ik), parameter :: S_XD = 2  !< diagnostic iterate (FGMRES) / verification divergence
    integer(kind=ik), parameter :: S_RD = 3  !< diagnostic residual
    integer(kind=ik), parameter :: S_W = 4   !< V-cycle work slot and one-component gradient scratch of A
    integer(kind=ik), parameter :: S_X = 5   !< solution phi
    integer(kind=ik), parameter :: S_K = 6   !< b, overwritten by the residual; then the free slots of projection_solve

    ! state of the current method, shared with the history hook
    integer(kind=ik)            :: tree_ID = 1, h_m, h_slots, h_restart, h_n_vc
    character(len=cshort)       :: h_tag
    real(kind=rk)               :: b_L2, t_start, t_hook, h_div(1:3)
    real(kind=rk), allocatable  :: blk(:, :, :)

contains

    !> history hook of projection_solve: x is the iterate after V-cycle n_vc, r = b - A x (the divergence u - G x
    !> would have). Writes the true divergence and the L-residual b - L x to projection_convergence.csv, and Hbar (end of
    !> an FGMRES cycle) to projection_hessenberg.csv. x is synced with the Poisson predictor (by A in projection_solve)
    subroutine proto_projection_hist(params, hvy_v, i_x, i_r, n_vc, Hbar)
        type (type_params), intent(inout)     :: params
        real(kind=rk), intent(inout)          :: hvy_v(:, :, :, :, :)
        integer(kind=ik), intent(in)          :: i_x, i_r, n_vc
        real(kind=rk), intent(in), optional   :: Hbar(:, :)

        real(kind=rk) :: dLinf, dL2, dmean, dLinf0, dL20, lLinf, lL2, nrm(1:1), x0(1:3), dx(1:3), tt, tsolve
        integer(kind=ik) :: kk, hvy_id, lgt_id, i, j, is(1:3), ie(1:3)

        tt = MPI_Wtime()
        tsolve = tt - t_start - t_hook
        is = 1
        ie = 1
        is(1:params%dim) = params%g + 1
        ie(1:params%dim) = params%g + params%Bs(1:params%dim)

        dLinf = slot_norm(i_r, "Linfty")
        dL2 = slot_norm(i_r, "L2")
        dmean = slot_norm(i_r, "Mean")
        ! without its mean, the constant part is not in the range of A on adaptive grids
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (block_is_leaf(params, hvy_id)) hvy_v(:,:,:,i_r,hvy_id) = hvy_v(:,:,:,i_r,hvy_id) - dmean
        enddo
        dLinf0 = slot_norm(i_r, "Linfty")
        dL20 = slot_norm(i_r, "L2")

        ! b - L x with the composite stencil, as the multigrid solver sees it
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (.not. block_is_leaf(params, hvy_id)) cycle
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            call GS_compute_Ax(params, hvy_v(:,:,:,i_x,hvy_id), blk, dx, apply_B_RHS=.false.)
            hvy_v(:,:,:,i_r,hvy_id) = hvy_v(:,:,:,S_B,hvy_id)
            hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),i_r,hvy_id) = hvy_v(is(1):ie(1),is(2):ie(2),is(3):ie(3),i_r,hvy_id) &
                - blk(is(1):ie(1),is(2):ie(2),is(3):ie(3))
        enddo
        lLinf = slot_norm(i_r, "Linfty")
        lL2 = slot_norm(i_r, "L2")

        h_n_vc = n_vc
        h_div = (/ dLinf, dL2, dmean /)
        if (params%rank == 0) then
            write(*, '(A, A, A, i0, A, i3, 4(A, es10.3))') "  ", trim(h_tag), " m=", h_m, " vc=", n_vc, &
                "  div Linf=", dLinf, "  div L2 rel=", dL2/b_L2, "  div mean=", dmean, "  div-mean Linf=", dLinf0
            open(unit=91, file="projection_convergence.csv", status="old", position="append")
            write(91, '(A, ",", i0, ",", i0, ",", i0, 10(",", es14.6))') trim(h_tag), h_m, h_slots, n_vc, dLinf, dL2, &
                dL2/b_L2, dmean, dLinf0, dL20/b_L2, lLinf, lL2/b_L2, -1.0_rk, tsolve
            close(91)

            ! Hessenberg matrix of this cycle, A Z_j = V_{j+1} Hbar: its eigenvalues are Ritz values of A*V
            if (present(Hbar)) then
                h_restart = h_restart + 1
                open(unit=93, file="projection_hessenberg.csv", status="old", position="append")
                do j = 1, size(Hbar, 2)
                    do i = 1, size(Hbar, 1)
                        write(93, '(A, 4(",", i0), ",", es22.14)') trim(h_tag), h_m, h_restart, i, j, Hbar(i, j)
                    enddo
                enddo
                close(93)
            endif
        endif
        t_hook = t_hook + MPI_Wtime() - tt

    contains

        real(kind=rk) function slot_norm(a, which)
            integer(kind=ik), intent(in) :: a
            character(len=*), intent(in) :: which
            call componentWiseNorm_tree(params, hvy_v(:,:,:,a:a,:), tree_ID, which, nrm, threshold_state_vector=.false.)
            slot_norm = nrm(1)
        end function slot_norm

    end subroutine proto_projection_hist

end module module_proto_projection


subroutine proto_projection_krylov(params)
    use mpi
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_fft
    use module_poisson
    use module_helpers, only: split_string
    use module_proto_projection

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: FD1_l(:), FD1_r(:)
    integer(kind=ik)                   :: FD1_ls, FD1_le, FD1_rs, FD1_re

    character(len=chuge)               :: arg, key, val, fname
    character(len=chuge)               :: files(1:6), methods_str, m_str, sol_str
    character(len=cshort), allocatable :: methods(:), m_items(:), sol_files(:)
    character(len=cshort)              :: wavelet_in, order_in, hybrid_rule
    integer(kind=ik)                   :: n_files, i_arg, i_meth, i_m, ipos, mmax

    integer(kind=ik)                   :: k, hvy_id, lgt_id, it_max, it, tc_length, Bs(1:3), nslot
    integer(kind=ik)                   :: is(1:3), ie(1:3), save_fields, save_u, ismooth, ismooth_width, mpierr
    integer(kind=ik), allocatable      :: m_list(:)
    real(kind=rk)                      :: time, domain(1:3), x0(1:3), dx(1:3), omega, ismooth_omega
    real(kind=rk)                      :: hybrid_kappa, hybrid_stall
    real(kind=rk)                      :: div0_Linf, div0_L2, b_mean, slot_GB, err0_Linf, err0_L2
    logical                            :: has_sol

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    ! command line: files (*.h5) and key=value pairs
    call get_command_argument(2, arg)
    if (arg=='--help' .or. arg=='--h' .or. arg=='-h' .or. len_trim(arg) == 0) then
        if (params%rank==0) then
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Wabbit prototype: exact Helmholtz projection (projection_solve), history per V-cycle"
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') " Reads a velocity field, computes b = div(u) and solves for phi with the chosen methods,"
            write(*, '(A)') " reporting the divergence left after u - grad(phi) after every V-cycle."
            write(*, '(A)') "-----------------------------------------------------------"
            write(*, '(A)') "./wabbit-post --proto-projection-krylov ux.h5 uy.h5 [uz.h5] [key=value ...]"
            write(*, '(A)') "    wavelet=CDF42            wavelet used for ghost node syncing and multigrid"
            write(*, '(A)') "    order=FD_4th_comp_1_3    D/G stencil pair (poisson_order), Laplacian is its composite"
            write(*, '(A)') "    methods=L,richardson,fgmres   methods to run (comma separated), also: hybrid"
            write(*, '(A)') "    m=1,2,3                  memory depths for fgmres and hybrid (comma separated), FGMRES(1) = minimal residual"
            write(*, '(A)') "    hybrid_rule=kappa        hybrid switch: kappa (only hybrid_kappa=4: switch after the first V-cycle), or stall"
            write(*, '(A)') "    hybrid_stall=0.3         stall: switch once a V-cycle reduced max|b - DG phi| by less than this factor"
            write(*, '(A)') "    ismooth=0                Jacobi sweeps on A after every preconditioning V-cycle (ismooth_width=-1,"
            write(*, '(A)') "                             ismooth_omega=1.0 are accepted for compatibility, other values are not)"
            write(*, '(A)') "    it=20                    number of V-cycles per method"
            write(*, '(A)') "    GS_it=8 sync_it=2 Jmin=0 multigrid settings (poisson_GS_it, poisson_Sync_it, poisson_Jmin)"
            write(*, '(A)') "    fft=FD                   accuracy of the coarsest FFT solve (FD or spectral)"
            write(*, '(A)') "    save=0                   save phi and div(u - grad phi) for every method"
            write(*, '(A)') "    save_u=0                 save the projected velocity u - grad phi of the last method as ux/uy/uz_000000000000.h5"
            write(*, '(A)') "    sol=usx.h5,usy.h5        exact divergence-free field on the same grid: error of u - grad(phi) to it"
            write(*, '(A)') " Output: projection_convergence.csv (per V-cycle), projection_summary.csv (per method),"
            write(*, '(A)') "         projection_hessenberg.csv (FGMRES cycles)"
            write(*, '(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    ! defaults
    wavelet_in = "CDF42"
    order_in = "FD_4th_comp_1_3"
    methods_str = "L,richardson,fgmres"
    m_str = "1,2,3"
    it_max = 20
    omega = 1.0_rk
    ismooth = 0
    ismooth_width = -1
    ismooth_omega = 1.0_rk
    hybrid_rule = "kappa"
    hybrid_kappa = 4.0_rk
    hybrid_stall = 0.3_rk
    params%poisson_GS_it = 8
    params%poisson_Sync_it = 2
    params%poisson_Jmin = 0
    params%FFT_accuracy = "FD"
    save_fields = 0
    save_u = 0
    n_files = 0
    sol_str = ""

    do i_arg = 2, command_argument_count()
        call get_command_argument(i_arg, arg)
        ipos = index(arg, "=")
        if (ipos == 0) then
            n_files = n_files + 1
            if (n_files > 3) call abort(260929, "proto_projection_krylov: at most 3 velocity files")
            files(n_files) = arg
            call check_file_exists(trim(arg))
            cycle
        endif
        key = arg(1:ipos-1)
        val = arg(ipos+1:)
        select case (trim(key))
        case ("wavelet")
            wavelet_in = trim(val)
        case ("order")
            order_in = trim(val)
        case ("methods")
            methods_str = trim(val)
        case ("m")
            m_str = trim(val)
        case ("it")
            read(val, *) it_max
        case ("tol")
            ! the history always runs it V-cycles (fixed_iterations); scripts pass tol=1e-13, i.e. never stop
            if (params%rank == 0) write(*, '(A)') " note: tol= is ignored, every method runs it= V-cycles"
        case ("omega")
            read(val, *) omega
        case ("ismooth")
            read(val, *) ismooth
        case ("ismooth_width")
            read(val, *) ismooth_width
        case ("ismooth_omega")
            read(val, *) ismooth_omega
        case ("hybrid_kappa")
            read(val, *) hybrid_kappa
        case ("hybrid_rule")
            hybrid_rule = trim(val)
        case ("hybrid_stall")
            read(val, *) hybrid_stall
        case ("GS_it")
            read(val, *) params%poisson_GS_it
        case ("sync_it")
            read(val, *) params%poisson_Sync_it
        case ("Jmin")
            read(val, *) params%poisson_Jmin
        case ("fft")
            params%FFT_accuracy = trim(val)
        case ("save")
            read(val, *) save_fields
        case ("save_u")
            read(val, *) save_u
        case ("sol")
            sol_str = trim(val)
        case default
            ! removed 2026-10-13 (development only): chebyshev, cheb_a, cheb_b, hybrid_rule=predict, hybrid_rho0,
            ! masked interface sweeps, ritz_shifts, diag
            call abort(260929, "proto_projection_krylov: unknown or removed option "//trim(key))
        end select
    enddo

    ! the variants of the development phase are gone, only what projection_solve implements is accepted
    if (omega /= 1.0_rk) call abort(261013, "proto_projection_krylov: omega/=1 removed (defect correction has step length 1)")
    if (ismooth_omega /= 1.0_rk) call abort(261013, "proto_projection_krylov: ismooth_omega/=1 removed")
    if (ismooth_width /= -1) call abort(261013, "proto_projection_krylov: masked sweeps removed, use ismooth_width=-1")
    if (trim(hybrid_rule) == "kappa" .and. hybrid_kappa /= 4.0_rk) &
        call abort(261013, "proto_projection_krylov: hybrid_rule=kappa only with hybrid_kappa=4 (switch after one V-cycle)")
    if (trim(hybrid_rule) /= "kappa" .and. trim(hybrid_rule) /= "stall") &
        call abort(261013, "proto_projection_krylov: hybrid_rule must be stall or kappa")

    call split_string(methods_str, methods, ",")
    call split_string(m_str, m_items, ",")
    allocate(m_list(1:size(m_items)))
    do i_m = 1, size(m_items)
        read(m_items(i_m), *) m_list(i_m)
    enddo
    has_sol = len_trim(sol_str) > 0
    if (has_sol) then
        call split_string(sol_str, sol_files, ",")
        do i_arg = 1, size(sol_files)
            call check_file_exists(trim(sol_files(i_arg)))
            files(n_files + i_arg) = sol_files(i_arg)
        enddo
    endif
    mmax = 1
    do i_meth = 1, size(methods)
        if (trim(methods(i_meth)) == "fgmres" .or. trim(methods(i_meth)) == "hybrid") mmax = max(mmax, maxval(m_list))
    enddo

    !-----------------------------------------------------------------------------------------------------
    ! setup parameters, in the style of proto_GS_multigrid
    call read_attributes(files(1), lgt_n(tree_ID), time, it, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
    if (n_files /= params%dim) call abort(260929, "proto_projection_krylov: need one velocity file per dimension")
    if (has_sol) then
        if (size(sol_files) /= params%dim) call abort(260929, "proto_projection_krylov: need one solution file per dimension")
    endif

    params%domainSizeCropped = domain
    params%Jmax = tc_length
    params%Jmin = 0
    params%domain_size(1:3) = domain(1:3)
    params%Bs = Bs
    params%n_eqn = params%dim
    if (has_sol) params%n_eqn = 2*params%dim
    params%wavelet = wavelet_in
    params%poisson_order = order_in
    params%poisson_coarsest = "FFT"
    params%poisson_balanceLoad = .false.
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component(1) = "x"
    params%symmetry_vector_component(2) = "y"
    if (params%dim == 3) params%symmetry_vector_component(3) = "z"
    if (has_sol) params%symmetry_vector_component(params%dim+1:2*params%dim) = params%symmetry_vector_component(1:params%dim)
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1
    select case(params%poisson_order)
    case ("FD_4th_comp_2_2")
        params%order_discretization = "FD_4th_central"
    case ("FD_6th_comp_3_3")
        params%order_discretization = "FD_6th_central"
    case default
        params%order_discretization = params%poisson_order
    end select
    if (params%poisson_order == "FD_6th_mehrstellen") then
        call abort(260929, "proto_projection_krylov: Mehrstellen has no D/G pair, the exact projection does not apply")
    endif

    call setup_wavelet(params, params%g)
    params%useCoarseExtension = params%isLiftedWavelet
    ! also initializes the FFT for the coarsest level
    call setup_laplacian_stencils(params, params%g)

    read(params%poisson_order_predictor(17:17), *) k
    if (k == 1) read(params%poisson_order_predictor(17:18), *) k
    params%g = max(params%g, k/2)

    call setup_FD1_left_stencil(params%poisson_order, FD1_l, FD1_ls, FD1_le)
    call setup_FD1_right_stencil(params%poisson_order, FD1_r, FD1_rs, FD1_re)
    if (params%g < max(abs(FD1_ls), FD1_le, abs(FD1_rs), FD1_re)) call abort(260929, "proto_projection_krylov: g too small for D/G stencils")

    ! projection_solve: exact projection, every method runs it_max V-cycles
    params%poisson_exact_projection = .true.
    params%poisson_cycle_end_criteria = "fixed_iterations"
    params%poisson_cycle_it = it_max

    ! we need to decompose values, so we need atleast 2**dim/2**(dim-1) blocks
    params%number_blocks = ceiling( real(lgt_n(tree_ID))/real(params%number_procs) * 2**params%dim/(2**params%dim-1) + 8 )

    nslot = max(S_K + 2*mmax, S_K - 1 + params%dim)
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, nrhs_slots1=1, neqn_hvy_tmp=nslot)
    allocate(blk(size(hvy_tmp,1), size(hvy_tmp,2), size(hvy_tmp,3)))
    slot_GB = real(size(hvy_tmp,1),rk)*real(size(hvy_tmp,2),rk)*real(size(hvy_tmp,3),rk)*real(size(hvy_tmp,5),rk)*8.0e-9_rk

    ! interior index range
    is = 1
    ie = 1
    is(1:params%dim) = params%g + 1
    ie(1:params%dim) = params%g + Bs(1:params%dim)

    call readHDF5vct_tree( files(1:params%n_eqn), params, hvy_block, tree_ID)
    call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID )
    call balanceLoad_tree( params, hvy_block, tree_ID )

    !-----------------------------------------------------------------------------------------------------
    ! b = D u on the leaf grid, exactly as in helmholtz_projection
    call sync_ghosts_tree( params, hvy_block(:,:,:,1:params%dim,:), tree_ID )
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        call compute_divergence(hvy_block(:,:,:,1:params%dim,hvy_id), dx, Bs, params%g, params%poisson_order, hvy_tmp(:,:,:,S_B,hvy_id))
    enddo
    call leaf_norms(S_B, div0_Linf, div0_L2)
    b_L2 = div0_L2
    ! on adaptive grids D is not conservative, so the mean of D u is not zero. A constant is not in the range
    ! of DG (its left null vector is close to constant), this part of the divergence cannot be removed by a gradient
    b_mean = leaf_mean(S_B)

    ! error of u to the exact divergence-free field before the projection
    err0_Linf = -1.0_rk
    err0_L2 = -1.0_rk
    if (has_sol) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            hvy_tmp(:,:,:,S_K:S_K+params%dim-1,hvy_id) = hvy_block(:,:,:,1:params%dim,hvy_id)
        enddo
        call vector_error(S_K, err0_Linf, err0_L2)
    endif

    if (params%rank == 0) then
        write(*, '(A)') "-----------------------------------------------------------"
        write(*, '(A, A, A, A, A, i0)') " wavelet=", trim(params%wavelet), "  order=", trim(params%poisson_order), "  g=", params%g
        write(*, '(A, A, A, A)') " predictor=", trim(params%order_predictor), "  poisson predictor=", trim(params%poisson_order_predictor)
        write(*, '(A, i0, A, i0, A, i0, A, i0)') " leaf blocks=", lgt_n(tree_ID), "  levels ", minActiveLevel_tree(tree_ID), "-", maxActiveLevel_tree(tree_ID), "  GS_it=", params%poisson_GS_it
        write(*, '(A, es12.4, A, es12.4, A, es12.4)') " before projection: max|div u|=", div0_Linf, "  L2=", div0_L2, "  mean=", b_mean
        write(*, '(A, es12.4, A)') " one slot = ", slot_GB, " GB per rank (allocated for number_blocks)"
        if (has_sol) write(*, '(A, es12.4, A, es12.4)') " error to divergence-free solution before projection: Linf=", err0_Linf, "  L2=", err0_L2
        write(*, '(A, i0, A, A, A, es12.4)') " ismooth=", ismooth, "  hybrid_rule=", trim(hybrid_rule), "  hybrid_stall=", hybrid_stall
        write(*, '(A)') "-----------------------------------------------------------"

        open(unit=91, file="projection_convergence.csv", status="replace")
        write(91, '(A)') "method,m,slots,vcycle,div_Linf,div_L2,div_L2_rel,div_mean,div_Linf_nomean,div_L2_rel_nomean,resL_Linf,resL_L2_rel,est_L2_rel,t_solver"
        close(91)
        open(unit=92, file="projection_summary.csv", status="replace")
        write(92, '(A)') "method,m,slots,vcycles,div0_Linf,div0_L2,div0_mean,div_Linf,div_L2,div_L2_rel,div_mean,verif_div_Linf,verif_div_L2,err0_Linf,err0_L2,err_Linf,err_L2,t_total,t_vcycle,t_op,t_krylov,t_diag"
        close(92)
        open(unit=93, file="projection_hessenberg.csv", status="replace")
        write(93, '(A)') "method,m,restart,i,j,H"
        close(93)
    endif

    !-----------------------------------------------------------------------------------------------------
    ! warm-up: the first V-cycles pay for first-touch allocation of the full-tree blocks and of the saved work
    ! arrays. Without this, the method that runs first looks slower, so these V-cycles are not timed
    params%poisson_cycle_it = 2
    call set_scheme(-1.0_rk, 1, 0)
    call solve(.false.)
    params%poisson_cycle_it = it_max

    !-----------------------------------------------------------------------------------------------------
    ! run all methods
    do i_meth = 1, size(methods)
        select case (trim(methods(i_meth)))
        case ("L")
            call run_method("L", 0)
        case ("richardson")
            call run_method("richardson", 0)
        case ("fgmres", "hybrid")
            do i_m = 1, size(m_list)
                call run_method(trim(methods(i_meth)), m_list(i_m))
            enddo
        case default
            call abort(260929, "proto_projection_krylov: unknown or removed method "//trim(methods(i_meth)))
        end select
    enddo

    call fft_destroy(params)
    call summarize_profiling( WABBIT_COMM )

contains

    !=====================================================================================================
    !> run one method: solve with projection_solve (history through the hook), verify, report
    subroutine run_method(method, m)
        character(len=*), intent(in) :: method
        integer(kind=ik), intent(in) :: m

        real(kind=rk) :: v_Linf, v_L2, e_Linf, e_L2, t_sum(1:5)

        ! projection_solve settings of the method; slots needed in production as in the report
        select case (method)
        case ("L")
            call set_scheme(1.0e30_rk, 1, 0)
            h_slots = 3
        case ("richardson")
            call set_scheme(-1.0_rk, 0, ismooth)
            h_slots = 4
            if (ismooth > 0) h_slots = 5
        case ("fgmres")
            call set_scheme(-1.0_rk, m, ismooth)
            h_slots = 2*m + 3
        case ("hybrid")
            if (trim(hybrid_rule) == "kappa") then
                call set_scheme(0.0_rk, m, ismooth)
            else
                call set_scheme(hybrid_stall, m, ismooth)
            endif
            h_slots = 2*m + 3
        end select
        h_tag = method
        h_m = m
        h_restart = 0
        h_n_vc = 0
        h_div = -1.0_rk
        if (params%rank == 0) write(*, '(A, A, A, i0, A, i0)') "=== method ", trim(h_tag), "  m=", m, "  slots=", h_slots

        call solve(.true., t_sum)

        ! independent verification with the standard routines on the leaf grid, as in helmholtz_projection
        call verify(v_Linf, v_L2, e_Linf, e_L2, method, m)

        if (params%rank == 0) then
            write(*, '(A, A, A, i0, A, i0, A, es10.3, A, es10.3, A, es10.3, A, es10.3)') "  -> ", trim(h_tag), " m=", m, &
                "  V-cycles=", h_n_vc, "  div Linf=", h_div(1), "  L2=", h_div(2), "  verification Linf=", v_Linf, "  L2=", v_L2
            open(unit=92, file="projection_summary.csv", status="old", position="append")
            write(92, '(A, ",", i0, ",", i0, ",", i0, 18(",", es14.6))') trim(h_tag), m, h_slots, h_n_vc, div0_Linf, div0_L2, b_mean, &
                h_div(1), h_div(2), h_div(2)/b_L2, h_div(3), v_Linf, v_L2, err0_Linf, err0_L2, e_Linf, e_L2, t_sum(1), t_sum(2), t_sum(3), t_sum(4), t_sum(5)
            close(92)
        endif
    end subroutine run_method

    !> poisson_exact_switch (1e30: classic V-cycles only, < 0: none), poisson_exact_m (0: defect correction), sweeps
    subroutine set_scheme(switch, m, sweeps)
        real(kind=rk), intent(in)    :: switch
        integer(kind=ik), intent(in) :: m, sweeps
        params%poisson_exact_switch = switch
        params%poisson_exact_m = m
        params%poisson_exact_sweeps = sweeps
    end subroutine set_scheme

    !> phi = 0, b from S_B, then projection_solve on the leaf grid. With history: the hook records every V-cycle and
    !> t_sum = (total, V-cycles, A, rest, diagnostics), maximum over the ranks
    subroutine solve(with_history, t_sum)
        logical, intent(in)                  :: with_history
        real(kind=rk), intent(out), optional :: t_sum(1:5)
        integer(kind=ik) :: kk, i_free(1:2*max(params%poisson_exact_m, 1))
        real(kind=rk) :: tim(1:3)

        do kk = 1, size(i_free)
            i_free(kk) = S_K + kk
        enddo
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            hvy_tmp(:,:,:,S_XD:nslot,hvy_id) = 0.0_rk
            hvy_tmp(:,:,:,S_K,hvy_id) = hvy_tmp(:,:,:,S_B,hvy_id)
        enddo

        t_hook = 0.0_rk
        t_start = MPI_Wtime()
        if (with_history) then
            call projection_solve(params, hvy_tmp, S_X, S_K, S_W, i_free, tree_ID, params%poisson_order, sync_rhs=.false., &
                init_0=.true., time=0.0_rk, hist=proto_projection_hist, i_b0=S_B, i_hist=S_XD, i_hres=S_RD, timing=tim)
        else
            call projection_solve(params, hvy_tmp, S_X, S_K, S_W, i_free, tree_ID, params%poisson_order, sync_rhs=.false., &
                init_0=.true., time=0.0_rk, timing=tim)
        endif
        if (present(t_sum)) then
            t_sum(1) = MPI_Wtime() - t_start - tim(3)
            t_sum(2:5) = (/ tim(1), tim(2), t_sum(1) - tim(1) - tim(2), tim(3) /)
            call MPI_ALLREDUCE(MPI_IN_PLACE, t_sum, 5, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        endif
    end subroutine solve

    !> u' = u - G phi and div(u') with the standard routines on the pruned leaf grid
    subroutine verify(v_Linf, v_L2, e_Linf, e_L2, method, m)
        real(kind=rk), intent(out)   :: v_Linf, v_L2, e_Linf, e_L2
        character(len=*), intent(in) :: method
        integer(kind=ik), intent(in) :: m
        integer(kind=ik) :: kk

        ! phi is synced with the Poisson predictor by projection_solve
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            hvy_tmp(:,:,:,S_K:S_K+params%dim-1,hvy_id) = hvy_block(:,:,:,1:params%dim,hvy_id)
            call subtract_scalar_gradient(hvy_tmp(:,:,:,S_K:S_K+params%dim-1,hvy_id), dx, Bs, params%g, &
                params%poisson_order, hvy_tmp(:,:,:,S_X,hvy_id))
        enddo
        call sync_ghosts_tree(params, hvy_tmp(:,:,:,S_K:S_K+params%dim-1,:), tree_ID)
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            call compute_divergence(hvy_tmp(:,:,:,S_K:S_K+params%dim-1,hvy_id), dx, Bs, params%g, &
                params%poisson_order, hvy_tmp(:,:,:,S_XD,hvy_id))
        enddo
        call leaf_norms(S_XD, v_Linf, v_L2)

        ! error of u' = u - G phi to the exact divergence-free field (constant offsets are not removed: G phi has none)
        e_Linf = -1.0_rk
        e_L2 = -1.0_rk
        if (has_sol) call vector_error(S_K, e_Linf, e_L2)

        ! projected velocity, overwritten by every method, so the files hold the one of the last method
        if (save_u == 1) then
            do kk = 1, params%dim
                fname = "u"//params%symmetry_vector_component(kk)//"_000000000000.h5"
                call saveHDF5_tree(fname, 0.0_rk, 0, 1, params, hvy_tmp(:,:,:,S_K+kk-1:S_K+kk-1,:), tree_ID)
            enddo
            if (params%rank == 0) write(*, '(A, A)') "  projected velocity saved, method ", trim(method)
        endif

        if (save_fields == 1) then
            if (method == "fgmres") then
                write(fname, '(A, i0, A)') "phi-fgmres", m, "_000000000000.h5"
            else
                write(fname, '(A, A, A)') "phi-", trim(method), "_000000000000.h5"
            endif
            call saveHDF5_tree(fname, 0.0_rk, 0, 1, params, hvy_tmp(:,:,:,S_X:S_X,:), tree_ID)
            if (method == "fgmres") then
                write(fname, '(A, i0, A)') "divu-fgmres", m, "_000000000000.h5"
            else
                write(fname, '(A, A, A)') "divu-", trim(method), "_000000000000.h5"
            endif
            call saveHDF5_tree(fname, 0.0_rk, 0, 1, params, hvy_tmp(:,:,:,S_XD:S_XD,:), tree_ID)
        endif
    end subroutine verify

    !> error of the vector field in slots v..v+dim-1 of hvy_tmp to u_sol in hvy_block(dim+1:2dim):
    !> Linf over all components, L2 of the vector
    subroutine vector_error(v, eLinf, eL2)
        integer(kind=ik), intent(in) :: v
        real(kind=rk), intent(out)   :: eLinf, eL2
        integer(kind=ik) :: kk, ic
        real(kind=rk) :: loc(1:2)
        loc = 0.0_rk
        do kk = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(kk, tree_ID)
            if (.not. block_is_leaf(params, hvy_id)) cycle
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            do ic = 1, params%dim
                loc(1) = max(loc(1), maxval(abs(hvy_tmp(is(1):ie(1),is(2):ie(2),is(3):ie(3),v+ic-1,hvy_id) &
                    - hvy_block(is(1):ie(1),is(2):ie(2),is(3):ie(3),params%dim+ic,hvy_id))))
                loc(2) = loc(2) + product(dx(1:params%dim)) * sum((hvy_tmp(is(1):ie(1),is(2):ie(2),is(3):ie(3),v+ic-1,hvy_id) &
                    - hvy_block(is(1):ie(1),is(2):ie(2),is(3):ie(3),params%dim+ic,hvy_id))**2)
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, loc(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, loc(2), 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        eLinf = loc(1)
        eL2 = sqrt(loc(2))
    end subroutine vector_error

    !> maximum norm and (volume weighted) L2 norm on the leaf blocks
    subroutine leaf_norms(a, n_Linf, n_L2)
        integer(kind=ik), intent(in) :: a
        real(kind=rk), intent(out)   :: n_Linf, n_L2
        real(kind=rk) :: nrm(1:1)
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,a:a,:), tree_ID, "Linfty", nrm, threshold_state_vector=.false.)
        n_Linf = nrm(1)
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,a:a,:), tree_ID, "L2", nrm, threshold_state_vector=.false.)
        n_L2 = nrm(1)
    end subroutine leaf_norms

    !> volume weighted mean on the leaf blocks
    function leaf_mean(a) result(mean)
        integer(kind=ik), intent(in) :: a
        real(kind=rk) :: mean, nrm(1:1)
        call componentWiseNorm_tree(params, hvy_tmp(:,:,:,a:a,:), tree_ID, "Mean", nrm, threshold_state_vector=.false.)
        mean = nrm(1)
    end function leaf_mean

end subroutine proto_projection_krylov
