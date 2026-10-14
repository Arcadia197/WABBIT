subroutine post_wavelet_contributions(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData
    use module_wavelets
    use module_hdf5_wrapper
    use, intrinsic :: ieee_arithmetic

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=chuge)              :: fname_in, fname_out, box_file, table_file
    character(len=chuge), allocatable :: fnames(:)
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgt_ID, tc_length, g, nc, i_arg, i_field
    integer(kind=ik), dimension(3)     :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_ID, ix, iy, iz, even_odd, idim, hvy_tmp_neqn
    integer(kind=ik)                   :: level_me, target_level, Jmax_active, Jmin_active
    integer(kind=ik)                   :: sub0_levels, sub_level, max_sub0_levels
    integer(kind=ik), external         :: count_div2
    integer(kind=ik)                   :: pos_underscore, pos_dot
    logical                            :: norms_only
    ! level and box energies from the wavelet coefficients (Parseval), without reconstruction
    logical                            :: coef_energies
    logical                            :: use_box
    logical                            :: logical_has_extent, logical_has_end
    real(kind=rk), dimension(3)        :: box_origin, box_extent, box_end
    real(kind=rk), dimension(3)        :: box_end_cmd
    ! boxes 1..n_boxes from --box-origin/extent/end (one box) or --box-file; box 0 is the full level contribution u_j
    integer(kind=ik)                   :: n_boxes, i_box
    real(kind=rk), allocatable         :: box_lo(:,:), box_hi(:,:)
    character(len=cshort), allocatable :: box_ids(:)
    ! hvy_tmp holds nc components (one per input field) per slot: slot 1 is the reconstruction helper, the others
    ! start after the offsets idx_dec (stored decomposition, norms-only), idx_coef (coefficients of the current
    ! contribution, boxes) and idx_full (full level contribution u_j, boxes)
    integer(kind=ik)                   :: idx_dec, idx_coef, idx_full
    integer                            :: table_unit
    integer(kind=ik)                   :: lgt_n_file, tc_length_file, dim_file, iteration_file
    integer(kind=ik), dimension(3)     :: Bs_file
    real(kind=rk)                      :: time_file
    real(kind=rk), dimension(3)        :: domain_file
    character(len=chuge)              :: arg

    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename)
    call get_command_argument(2, fname_in)

    ! does the user need help?
    if (fname_in=='--help' .or. fname_in=='--h' .or. fname_in=='-h') then
        if (params%rank==0) then
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " Wabbit postprocessing: compute wavelet contributions per level"
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " This tool decomposes the input field using wavelets, then"
            write(*,'(A)') " reconstructs the field keeping only wavelet coefficients"
            write(*,'(A)') " from each level separately. Outputs one file per level."
            write(*,'(A)') ""
            write(*,'(A)') " ./wabbit-post --wavelet-contributions ux_000.h5 --wavelet=CDF40"
            write(*,'(A)') ""
            write(*,'(A)') " Output: ux-WC0_000.h5, ux-WC1_000.h5, ..., ux-WCJ_000.h5"
            write(*,'(A)') " With sub0 enabled: ux-WC0s1_000.h5, ux-WC0s2_000.h5, ..."
            write(*,'(A)') " Norms of all contributions: table ux-WC_000.t (see --table-file)"
            write(*,'(A)') ""
            write(*,'(A)') " Several fields on the same grid are treated together (one call, one table):"
            write(*,'(A)') " ./wabbit-post --wavelet-contributions ux_000.h5 uy_000.h5 uz_000.h5 --wavelet=CDF40"
            write(*,'(A)') ""
            write(*,'(A)') " To only compute norms per level (no files saved):"
            write(*,'(A)') " ./wabbit-post --wavelet-contributions ux_000.h5 --wavelet=CDF40 --norms-only"
            write(*,'(A)') ""
            write(*,'(A)') "-----------------------------------------------------------"
            write(*,'(A)') " --wavelet=CDF44   - Wavelet type (default: CDF40)"
            write(*,'(A)') " --Jmin=N           - Minimum contribution level to process (default: 0)"
            write(*,'(A)') " --norms-only       - Only compute and print norms, do not save files"
            write(*,'(A)') " --sub0-levels=N    - Additional WC-only periodic sub0 contributions"
            write(*,'(A)') "                     (requires active Jmin=0, max = (times Bs is divisible by 2) - 1)"
            write(*,'(A)') " --box-origin-x/y/z - Box origin for local contributions"
            write(*,'(A)') " --box-extent-x/y/z - Box extent for local contributions"
            write(*,'(A)') " --box-end-x/y/z    - Box end point for local contributions"
            write(*,'(A)') "                     (outside coefficients are set to zero, slice-wise)"
            write(*,'(A)') "                     Per dimension: use either --box-extent OR --box-end"
            write(*,'(A)') "                     For each dimension, specify all 3 (x,y,z) consistently"
            write(*,'(A)') " --box-file=FILE    - Many boxes in one call (requires --norms-only), one box per line:"
            write(*,'(A)') "                       id; x0; y0; z0; x1; y1; z1      (2D: id; x0; y0; x1; y1)"
            write(*,'(A)') "                     with origin (x0,y0,z0) and end (x1,y1,z1). Lines starting with % or #"
            write(*,'(A)') "                     and empty lines are skipped, spaces are ignored (pad it like a table)."
            write(*,'(A)') " --table-file=FILE  - Name of the table (default: first input file with -WC and .t)"
            write(*,'(A)') "                     One row per field, box, level and sub-level, ';'-separated, '%' header."
            write(*,'(A)') "                     Box 0 is the full level contribution u_j. For every box, the table also"
            write(*,'(A)') "                     holds the L2 norm of u_j and the inner product int(u_box*u_j). 1/2*int(u_box*u_j)"
            write(*,'(A)') "                     splits the energy additively among boxes (cross terms are shared equally)"
            write(*,'(A)') " --coefficient-energies - Energies from the wavelet coefficients alone (Parseval), no reconstruction:"
            write(*,'(A)') "                     much less memory and time. Exact for orthogonal wavelets (coiflets), also on"
            write(*,'(A)') "                     adaptive grids; for CDF wavelets the numbers are NOT the energies of u_j."
            write(*,'(A)') "                     A coefficient belongs to a box if its grid point lies in the box. The table holds"
            write(*,'(A)') "                     L2 = sqrt(E), inner_box_full = E (box energies add up exactly), other norms NaN."
            write(*,'(A)') "                     Implies --norms-only, not with --sub0-levels."
            write(*,'(A)') "-----------------------------------------------------------"
        end if
        return
    endif

    ! input fields: all arguments after --wavelet-contributions up to the first option
    nc = 0
    do i_arg = 2, command_argument_count()
        call get_command_argument(i_arg, fname_in)
        if (fname_in(1:2) == "--") exit
        nc = nc + 1
    enddo
    if (nc == 0) call abort(261013_ik, "--wavelet-contributions: no input file given")
    allocate(fnames(1:nc))
    do i_field = 1, nc
        call get_command_argument(i_field+1, fnames(i_field))
    enddo
    ! options are --name=value, so a later argument without a leading - is a misplaced input file
    do i_arg = nc+2, command_argument_count()
        call get_command_argument(i_arg, arg)
        if (arg(1:1) /= "-") call abort(261013_ik, "--wavelet-contributions: input files must come before all options, found "//trim(arg))
    enddo
    fname_in = fnames(1)

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )
    call get_cmd_arg( "--Jmin", params%Jmin, default=0_ik )
    call get_cmd_arg_bool( "--norms-only", norms_only, default=.false. )
    call get_cmd_arg( "--sub0-levels", sub0_levels, default=0_ik )
    call get_cmd_arg( "--box-origin-x", box_origin(1), default=0.0_rk )
    call get_cmd_arg( "--box-origin-y", box_origin(2), default=0.0_rk )
    call get_cmd_arg( "--box-origin-z", box_origin(3), default=0.0_rk )
    call get_cmd_arg( "--box-extent-x", box_extent(1), default=-1.0_rk )
    call get_cmd_arg( "--box-extent-y", box_extent(2), default=-1.0_rk )
    call get_cmd_arg( "--box-extent-z", box_extent(3), default=-1.0_rk )
    call get_cmd_arg( "--box-end-x", box_end_cmd(1), default=-1.0_rk )
    call get_cmd_arg( "--box-end-y", box_end_cmd(2), default=-1.0_rk )
    call get_cmd_arg( "--box-end-z", box_end_cmd(3), default=-1.0_rk )
    call get_cmd_arg( "--box-file", box_file, default="" )
    call get_cmd_arg( "--table-file", table_file, default="" )
    call get_cmd_arg_bool( "--coefficient-energies", coef_energies, default=.false. )
    if (coef_energies) then
        if (sub0_levels /= 0) call abort(261013_ik, "--coefficient-energies does not support --sub0-levels")
        norms_only = .true.
    endif

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    params%forest_size = 1
    params%n_eqn = nc

    call check_file_exists(trim(fname_in))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_in, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! several fields are treated as components of one state vector, so they have to share the grid
    do i_field = 2, nc
        call check_file_exists(trim(fnames(i_field)))
        call read_attributes(fnames(i_field), lgt_n_file, time_file, iteration_file, domain_file, Bs_file, tc_length_file, dim_file)
        if (lgt_n_file /= lgt_n(tree_ID) .or. any(Bs_file /= Bs) .or. tc_length_file /= tc_length .or. dim_file /= params%dim &
            .or. any(abs(domain_file - domain) > 1.0e-12_rk*maxval(domain))) then
            call abort(261013_ik, "--wavelet-contributions: all input files must be on the same grid, but "//trim(fnames(i_field))// &
                " differs from "//trim(fnames(1)))
        endif
    enddo
    if (nc > 1) call check_same_tree()

    ! no cropping is used by default, but the variable is not set automatically
    params%domainSizeCropped = domain

    ! Check box parameters: for each dimension, use either extent or end (not both)
    ! but allow mixing across dimensions
    box_end = 0.0_rk
    do idim = 1, params%dim
        logical_has_extent = box_extent(idim) > 0.0_rk
        logical_has_end = box_end_cmd(idim) > 0.0_rk

        if (logical_has_extent .and. logical_has_end) then
            write(*,'(A,I0,A)') "ERROR: Cannot specify both --box-extent-", idim, &
                " and --box-end-", idim
            call abort(260611_ik, "Cannot use both --box-extent and --box-end for same dimension")
        endif

        ! If end point is specified, validate it's greater than origin
        if (logical_has_end .and. box_end_cmd(idim) <= box_origin(idim)) then
            write(*,'(A,I0,A,F16.8,A,F16.8)') "ERROR: box-end-", idim, " (", box_end_cmd(idim), &
                ") must be > box-origin-", idim, " (", box_origin(idim), ")"
            call abort(260611_ik, "Box end point must be greater than box origin")
        endif

        ! Convert extent to end point if extent was used
        if (logical_has_extent) then
            box_end(idim) = box_origin(idim) + box_extent(idim)
        else if (logical_has_end) then
            box_end(idim) = box_end_cmd(idim)
        endif
    enddo

    ! Check if all dimension has box parameters set
    use_box = all((box_extent(1:params%dim) > 0.0_rk) .or. (box_end_cmd(1:params%dim) > 0.0_rk))

    if (box_file /= "" .and. (any(box_extent > 0.0_rk) .or. any(box_end_cmd > 0.0_rk) .or. any(abs(box_origin) > 0.0_rk))) then
        call abort(261013_ik, "Use either --box-file or --box-origin/--box-extent/--box-end, not both")
    endif
    if (box_file /= "") then
        if (.not. norms_only) call abort(261013_ik, "--box-file requires --norms-only")
        call read_box_file(trim(box_file))
    elseif (use_box) then
        n_boxes = 1
        allocate(box_lo(1:3,1:1), box_hi(1:3,1:1), box_ids(1:1))
        box_lo(:,1) = box_origin
        box_hi(:,1) = box_end
        box_ids(1) = "box"
    else
        n_boxes = 0
        allocate(box_lo(1:3,1:0), box_hi(1:3,1:0), box_ids(1:0))
    endif

    params%Jmax = tc_length
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))

    ! Coarse extension at coarse-fine interfaces: for the lifted CDF wavelets and the coiflets (whose SC are filtered
    ! values too; on adaptive test grids it reduces their level-energy error from up to 57 % to 2 %). Not for coefficient
    ! energies, which are read off the decomposition: coarse extension would delete coefficients at the interfaces
    ! (-2 % at the finest level with, 0.2 % per level without).
    params%useCoarseExtension = (params%isLiftedWavelet .or. params%wavelet(1:7) == "coiflet") .and. .not. coef_energies
    params%useSecurityZone = params%isLiftedWavelet
    if (coef_energies .and. params%wavelet(1:7) /= "coiflet" .and. params%rank == 0) then
        write(*,'(A)') "WARNING: --coefficient-energies with a biorthogonal wavelet: the numbers are not the energies of u_j"
    endif

    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 1

    Bs = params%Bs
    g  = params%g

    ! for full wavelet operations we need 8/7 for 3D or 4/3 for 2D as blocks
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk) * 1.2_rk)+4

    ! slot 1: helper for the reconstruction
    hvy_tmp_neqn = nc
    idx_dec = -1
    idx_coef = -1
    idx_full = -1
    ! for norms-only, we can skip reloading, as we do not need to prune: keep the decomposition
    ! (coefficient energies need neither the copy nor the box slots: hvy_tmp is only used by the decomposition)
    if (norms_only .and. .not. coef_energies) then
        idx_dec = hvy_tmp_neqn
        hvy_tmp_neqn = hvy_tmp_neqn + nc
    endif
    ! with boxes, keep the coefficients of the current contribution (restored for every box) and the
    ! full level contribution u_j (for the inner product)
    if (n_boxes > 0 .and. .not. coef_energies) then
        idx_coef = hvy_tmp_neqn
        idx_full = hvy_tmp_neqn + nc
        hvy_tmp_neqn = hvy_tmp_neqn + 2*nc
    endif

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, neqn_hvy_tmp=hvy_tmp_neqn)

    ! read input data
    call readHDF5vct_tree( fnames, params, hvy_block, tree_ID)

    call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)

    ! perform wavelet decomposition
    if (params%rank == 0) then
        write(*,'(A)') "Performing wavelet decomposition..."
    endif
    call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)
    if (norms_only .and. .not. coef_energies) call copy_block_tmp(idx_dec, to_tmp=.true.)

    ! determine active level range
    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    if (params%Jmin < 0) then
        call abort(260611_ik, "--Jmin must be >= 0")
    endif
    if (params%Jmin > Jmax_active) then
        call abort(260611_ik, "--Jmin exceeds max active level in input")
    endif

    if (params%rank == 0) then
        write(*,'(A,I0,A,I0)') "Active levels: ", params%Jmin, " to ", Jmax_active
        write(*,'(A)') "Computing contributions for each level ..."
    endif

    if (sub0_levels < 0) then
        call abort(260611_ik, "--sub0-levels must be >= 0")
    endif
    if (sub0_levels > 0) then
        if (params%Jmin /= 0) then
            call abort(260611_ik, "--sub0-levels requires --Jmin=0")
        endif
        if (Jmin_active /= 0) then
            call abort(260611_ik, "sub0 requested but active Jmin is not 0")
        endif

        max_sub0_levels = huge(1_ik)
        max_sub0_levels = min(max_sub0_levels, count_div2(Bs(1)))
        max_sub0_levels = min(max_sub0_levels, count_div2(Bs(2)))
        if (params%dim == 3) max_sub0_levels = min(max_sub0_levels, count_div2(Bs(3)))
        ! the SC of level 0 live on Bs/2 points per direction and every sub-level halves this lattice,
        ! which therefore has to stay even: Bs=16 -> 3 sub-levels, Bs=18 -> none
        max_sub0_levels = max_sub0_levels - 1

        if (sub0_levels > max_sub0_levels) then
            write(*,'(A,I0,A,I0)') "Requested --sub0-levels=", sub0_levels, " but max allowed is ", max_sub0_levels
            call abort(260611_ik, "sub0-levels exceeds divide-by-2 limit")
        endif
    endif

    ! even_odd pattern for identifying WC vs SC
    even_odd = mod(params%g + 1, 2)

    ! default table name: first input file with -WC before the last underscore (or before .h5) and .t instead of .h5
    if (table_file == "") then
        pos_dot = index(fname_in, ".h5", back=.true.)
        if (pos_dot == 0) pos_dot = len_trim(fname_in) + 1
        pos_underscore = index(fname_in(1:pos_dot-1), "_", back=.true.)
        ! only an underscore in the file name counts, not one in the directory part
        if (pos_underscore <= index(fname_in(1:pos_dot-1), "/", back=.true.)) pos_underscore = pos_dot
        table_file = fname_in(1:pos_underscore-1)//"-WC"//fname_in(pos_underscore:pos_dot-1)//".t"
    endif
    call open_table()

    if (coef_energies) call output_coefficient_energies()

    ! loop over all active levels (not with coefficient energies, which need no reconstruction)
    do target_level = params%Jmin, Jmax_active
        if (coef_energies) exit
        if (params%rank == 0) then
            write(*,'(A,I0,A,I0)') "Processing level ", target_level, " / ", Jmax_active
        endif

        ! fresh start for each target contribution
        call restore_decomposition(fresh=(target_level == params%Jmin))

        ! zero out wavelet coefficients on all levels except target_level
        ! currently we intentionally keep WC only (SC are discarded for all levels)
        ! future extension could add an additional SC-only output branch
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

            if (level_me /= target_level) then
                ! zero out all values on this block, no matter if SC or WC
                hvy_block(:,:,:,1:nc,hvy_ID) = 0.0_rk
            ! else if (target_level > Jmin_active .and. target_level > params%Jmin) then
            else
                ! on the target level, keep WC but zero SC
                do iz = merge(1, g+1, params%dim==2), merge(1, Bs(3)+g, params%dim==2)
                    do iy = g+1, Bs(2)+g
                        do ix = g+1, Bs(1)+g
                            ! check if this is a scaling coefficient point
                            if (mod(ix,2)==even_odd .and. mod(iy,2)==even_odd .and. &
                                mod(iz,2)==merge(1,even_odd,params%dim==2)) then
                                ! this is a SC point - zero it out
                                hvy_block(ix,iy,iz,1:nc,hvy_ID) = 0.0_rk
                            endif
                        enddo
                    enddo
                enddo
            endif
        enddo

        if (target_level == 0 .and. sub0_levels > 0) then
            ! First, save the regular L0 contribution (WC0), then add sub-level distributions.
            call finalize_and_output(target_level, 0_ik)

            do sub_level = 1, sub0_levels
                ! fresh start for each sub-level contribution
                call restore_decomposition(fresh=.false.)

                ! keep only level-0 blocks before sub0 processing
                do k = 1, hvy_n(tree_ID)
                    hvy_ID = hvy_active(k, tree_ID)
                    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                    level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

                    if (level_me /= 0) then
                        hvy_block(:,:,:,1:nc,hvy_ID) = 0.0_rk
                    endif
                enddo

                ! hierarchical sub-level extraction:
                ! SC0 -> (SC-1,WC-1) -> ... -> (SC-sub_level,WC-sub_level)
                ! while existing coarser WC are kept untouched by sub-level decomposition.
                do k = 1, hvy_n(tree_ID)
                    hvy_ID = hvy_active(k, tree_ID)
                    call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
                    level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

                    if (level_me == 0) then
                        ! input is L0-decomposed, we further decompose it until we get the desired sub-level decomposition
                        call sub0_wavelet_decompose_n(params, hvy_block(:,:,:,1:nc,hvy_ID), 0_ik, sub_level, hvy_tmp(:,:,:,1:nc,hvy_ID))

                        ! keep only WC-sub_level:
                        ! remove WC0..WC(sub_level-1), then remove SC-sub_level
                        do iy = 0, sub_level-1
                            call sub0_delete_coefficients(params, hvy_block(:,:,:,1:nc,hvy_ID), iy, "WC")
                        enddo
                        call sub0_delete_coefficients(params, hvy_block(:,:,:,1:nc,hvy_ID), sub_level, "SC")

                        ! reconstruction from sub-level to level-0 decomposition is deferred to finalize_and_output,
                        ! so it can happen after optional local box filtering
                    endif
                enddo

                call finalize_and_output(0_ik, sub_level)

            enddo

            cycle
        endif

        call finalize_and_output(target_level, 0_ik)
    enddo

    if (params%rank == 0) then
        close(table_unit)
        write(*,'(A,A)') "Table: ", trim(table_file)
        write(*,'(A)') "Done! Created contributions for all levels."
    endif

contains

    !> Norms of the contribution whose coefficients are in hvy_block (sublevel_id 0: regular level contribution):
    !> first the full level contribution u_j (box 0), then every box, one table row per field and box
    subroutine finalize_and_output(level_id, sublevel_id)
        implicit none

        integer(kind=ik), intent(in) :: level_id
        integer(kind=ik), intent(in) :: sublevel_id
        integer(kind=ik) :: n_nonzero_wavelets(1:nc)
        real(kind=rk) :: norms(1:nc,1:7), norm_full(1:nc), inner_box_full(1:nc)

        ! keep the coefficients, they are restored for every box
        if (n_boxes > 0) call copy_block_tmp(idx_coef, to_tmp=.true.)

        ! box 0: full level contribution u_j
        call count_nonzero_wavelets(n_nonzero_wavelets)
        call reconstruct_contribution(sublevel_id)
        call compute_norms(norms)
        norm_full = norms(:,2)
        inner_box_full = norm_full**2
        call write_rows(level_id, sublevel_id, 0_ik, norms, norm_full, inner_box_full, n_nonzero_wavelets)

        if (params%rank == 0) then
            if (sublevel_id == 0) then
                write(*,'(A,I0,A,I0)') "  Level ", level_id, " / ", Jmax_active
            else
                write(*,'(A,I0,A,I0)') "  Level 0 sub-level ", sublevel_id, " / ", sub0_levels
            endif
            if (nc == 1 .and. n_boxes == 0) then
                write(*,'(A,es16.8)') "    L1 norm:     ", norms(1,1)
                write(*,'(A,es16.8)') "    L2 norm:     ", norms(1,2)
                write(*,'(A,es16.8)') "    L3 norm:     ", norms(1,3)
                write(*,'(A,es16.8)') "    L4 norm:     ", norms(1,4)
                write(*,'(A,es16.8)') "    LInfty norm: ", norms(1,5)
                write(*,'(A,es16.8)') "    H1 norm:     ", norms(1,6)
                write(*,'(A,es16.8)') "    Mean:     ", norms(1,7)
                write(*,'(A,I0)')     "    Non-zero wavelets: ", n_nonzero_wavelets(1)
            else
                ! with several fields or boxes, the numbers are in the table only
                do i_field = 1, nc
                    write(*,'(A,I0,A,es16.8,A,I0)') "    field ", i_field, ": full-level L2 norm ", norm_full(i_field), &
                        ", non-zero wavelets ", n_nonzero_wavelets(i_field)
                enddo
                if (n_boxes > 0) write(*,'(A,I0,A)') "    ", n_boxes, " boxes, see table"
            endif
        endif

        if (n_boxes > 0) then
            call copy_block_tmp(idx_full, to_tmp=.true.)

            do i_box = 1, n_boxes
                call copy_block_tmp(idx_coef, to_tmp=.false.)
                call apply_local_box_filter(i_box)

                ! count non-zero wavelet coefficients for the current contribution (for output, not used in reconstruction)
                call count_nonzero_wavelets(n_nonzero_wavelets)
                call reconstruct_contribution(sublevel_id)

                ! inner product of the box contribution u_box with the full level contribution u_j
                call inner_product_full_level(inner_box_full)
                call compute_norms(norms)
                call write_rows(level_id, sublevel_id, i_box, norms, norm_full, inner_box_full, n_nonzero_wavelets)
            enddo
        endif

        ! hvy_block holds u_j, or u_box with a single box from the command line (--box-file requires --norms-only)
        if (.not. norms_only) then
            call prune_fulltree2leafs(params, tree_ID)
            do i_field = 1, nc
                ! append ending before last underscore or if it doesnt exist before last .
                pos_underscore = index(fnames(i_field), "_", back=.true.)
                if (pos_underscore == 0) pos_underscore = index(fnames(i_field), '.', back=.true.)
                if (sublevel_id == 0) then
                    write(fname_out, '(A, A, I0, A)') fnames(i_field)(1:pos_underscore-1), "-WC", level_id, &
                        fnames(i_field)(pos_underscore:LEN_TRIM(fnames(i_field)))
                else
                    write(fname_out, '(A, A, I0, A)') fnames(i_field)(1:pos_underscore-1), "-WC0s", sublevel_id, &
                        fnames(i_field)(pos_underscore:LEN_TRIM(fnames(i_field)))
                endif

                call saveHDF5_tree(fname_out, time, iteration, i_field, params, hvy_block, tree_ID )
                if (params%rank == 0) then
                    write(*,'(A,A)') "  Saved: ", trim(fname_out)
                endif
            enddo
        endif

    end subroutine

    !> Level (and box) energies from the wavelet coefficients alone (Parseval), without any reconstruction.
    !> For an orthogonal wavelet, u_j = sum_k d_k psi_k has ||u_j||^2 = sum_k d_k^2 (orthonormal filters). With the
    !> normalisation used in WABBIT (HD sums to 1, see setup_wavelet), a coefficient carries a factor 2^(-1/2) for every
    !> direction in which it was low-pass filtered, and the volume element is that of its level, so
    !>     ||u_j||^2 = sum_{WC of level j} c^2 * 2^n_even * dx_j^dim,
    !> with n_even the number of directions in which the coefficient sits on a SC position. A coefficient belongs to a box
    !> if its grid point lies in the box (closed box, as in apply_local_box_filter), so box energies add up exactly.
    !> hvy_block holds the full-tree decomposition; every block of level j holds the coefficients of level j.
    subroutine output_coefficient_energies()
        implicit none
        integer(kind=ik) :: k_loc, hvy_id_loc, lgt_id_loc, lvl, i_comp, ix_loc, iy_loc, iz_loc, n_even, mpierr_loc, i_dim
        real(kind=rk) :: x0_loc(3), dx_loc(3), w, xyz(3), c, nan_val
        real(kind=rk) :: norms(1:nc,1:7)
        real(kind=rk), allocatable :: e(:,:,:)
        integer(kind=ik), allocatable :: nnz(:,:,:)
        logical, allocatable :: box_hit(:)

        allocate(e(1:nc, 0:n_boxes, params%Jmin:Jmax_active), nnz(1:nc, 0:n_boxes, params%Jmin:Jmax_active))
        allocate(box_hit(1:n_boxes))
        e = 0.0_rk
        nnz = 0_ik
        xyz = 0.0_rk

        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            call hvy2lgt(lgt_id_loc, hvy_id_loc, params%rank, params%number_blocks)
            lvl = lgt_block(lgt_id_loc, IDX_MESH_LVL)
            if (lvl < params%Jmin) cycle
            call get_block_spacing_origin(params, lgt_id_loc, x0_loc, dx_loc)

            ! boxes that overlap this block at all
            do i_box = 1, n_boxes
                box_hit(i_box) = .true.
                do i_dim = 1, params%dim
                    if (x0_loc(i_dim) + dble(Bs(i_dim)-1)*dx_loc(i_dim) < box_lo(i_dim, i_box) .or. &
                        x0_loc(i_dim) > box_hi(i_dim, i_box)) box_hit(i_box) = .false.
                enddo
            enddo

            do iz_loc = merge(1, g+1, params%dim==2), merge(1, Bs(3)+g, params%dim==2)
                do iy_loc = g+1, Bs(2)+g
                    do ix_loc = g+1, Bs(1)+g
                        n_even = 0
                        if (mod(ix_loc,2) == even_odd) n_even = n_even + 1
                        if (mod(iy_loc,2) == even_odd) n_even = n_even + 1
                        if (params%dim == 3 .and. mod(iz_loc,2) == even_odd) n_even = n_even + 1
                        ! scaling coefficient
                        if (n_even == params%dim) cycle

                        w = 2.0_rk**n_even * product(dx_loc(1:params%dim))
                        xyz(1) = x0_loc(1) + dble(ix_loc-(g+1)) * dx_loc(1)
                        xyz(2) = x0_loc(2) + dble(iy_loc-(g+1)) * dx_loc(2)
                        if (params%dim == 3) xyz(3) = x0_loc(3) + dble(iz_loc-(g+1)) * dx_loc(3)

                        do i_comp = 1, nc
                            c = hvy_block(ix_loc, iy_loc, iz_loc, i_comp, hvy_id_loc)
                            e(i_comp, 0, lvl) = e(i_comp, 0, lvl) + c**2 * w
                            if (abs(c) > 1.0e-12_rk) nnz(i_comp, 0, lvl) = nnz(i_comp, 0, lvl) + 1
                            do i_box = 1, n_boxes
                                if (.not. box_hit(i_box)) cycle
                                if (all(xyz(1:params%dim) >= box_lo(1:params%dim, i_box) .and. &
                                        xyz(1:params%dim) <= box_hi(1:params%dim, i_box))) then
                                    e(i_comp, i_box, lvl) = e(i_comp, i_box, lvl) + c**2 * w
                                    if (abs(c) > 1.0e-12_rk) nnz(i_comp, i_box, lvl) = nnz(i_comp, i_box, lvl) + 1
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, e, size(e), MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr_loc)
        call MPI_ALLREDUCE(MPI_IN_PLACE, nnz, size(nnz), MPI_INTEGER4, MPI_SUM, WABBIT_COMM, mpierr_loc)

        nan_val = ieee_value(nan_val, ieee_quiet_nan)
        if (params%rank == 0) write(*,'(A)') "Level energies from the wavelet coefficients (--coefficient-energies, no reconstruction):"
        do lvl = params%Jmin, Jmax_active
            do i_box = 0, n_boxes
                norms = nan_val
                norms(:,2) = sqrt(e(:, i_box, lvl))
                call write_rows(lvl, 0_ik, i_box, norms, sqrt(e(:, 0, lvl)), e(:, i_box, lvl), nnz(:, i_box, lvl))
            enddo
            if (params%rank == 0) then
                ! same line as in finalize_and_output, log parsers match it exactly
                write(*,'(A,I0,A,I0)') "  Level ", lvl, " / ", Jmax_active
                if (nc == 1 .and. n_boxes == 0) then
                    write(*,'(A,es16.8)') "    L2 norm:     ", sqrt(e(1, 0, lvl))
                    write(*,'(A,I0)')     "    Non-zero wavelets: ", nnz(1, 0, lvl)
                else
                    do i_field = 1, nc
                        write(*,'(A,I0,A,es16.8,A,I0)') "    field ", i_field, ": full-level L2 norm ", sqrt(e(i_field, 0, lvl)), &
                            ", non-zero wavelets ", nnz(i_field, 0, lvl)
                    enddo
                    if (n_boxes > 0) write(*,'(A,I0,A)') "    ", n_boxes, " boxes, see table"
                endif
            endif
        enddo
    end subroutine

    !> Reconstruct the contribution currently held in hvy_block (coefficients) to physical space.
    !> For sub0 outputs, first reconstruct from isolated WC-sub_level to level-0 decomposition.
    subroutine reconstruct_contribution(sublevel_id)
        implicit none
        integer(kind=ik), intent(in) :: sublevel_id

        integer(kind=ik) :: k_loc, hvy_id_loc, lgt_id_loc

        ! ensure reconstruction helper data matches the current coefficients
        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            hvy_tmp(:,:,:,1:nc,hvy_id_loc) = hvy_block(:,:,:,1:nc,hvy_id_loc)
        enddo

        if (sublevel_id > 0_ik) then
            do k_loc = 1, hvy_n(tree_ID)
                hvy_id_loc = hvy_active(k_loc, tree_ID)
                call hvy2lgt(lgt_id_loc, hvy_id_loc, params%rank, params%number_blocks)

                if (lgt_block(lgt_id_loc, IDX_MESH_LVL) == 0) then
                    call sub0_wavelet_reconstruct_n(params, hvy_block(:,:,:,1:nc,hvy_id_loc), sublevel_id, 0_ik, hvy_tmp(:,:,:,1:nc,hvy_id_loc))
                endif
            enddo
        endif

        call wavelet_reconstruct_full_tree(params, hvy_block, hvy_tmp, tree_ID)
    end subroutine

    !> copy hvy_block to the hvy_tmp slot starting after idx (to_tmp) or back
    subroutine copy_block_tmp(idx, to_tmp)
        implicit none
        integer(kind=ik), intent(in) :: idx
        logical, intent(in) :: to_tmp

        integer(kind=ik) :: k_loc, hvy_id_loc

        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            if (to_tmp) then
                hvy_tmp(:,:,:,idx+1:idx+nc,hvy_id_loc) = hvy_block(:,:,:,1:nc,hvy_id_loc)
            else
                hvy_block(:,:,:,1:nc,hvy_id_loc) = hvy_tmp(:,:,:,idx+1:idx+nc,hvy_id_loc)
            endif
        enddo
    end subroutine

    subroutine count_nonzero_wavelets(n_nonzero)
        implicit none
        integer(kind=ik), intent(out) :: n_nonzero(1:nc)

        integer(kind=ik) :: k_loc, hvy_id_loc, i_comp, g_loc(3), mpierr_loc

        g_loc(:) = params%g
        if (params%dim == 2) g_loc(3) = 0

        n_nonzero = 0_ik
        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            do i_comp = 1, nc
                n_nonzero(i_comp) = n_nonzero(i_comp) + count(abs(hvy_block(g_loc(1)+1:g_loc(1)+Bs(1), g_loc(2)+1:g_loc(2)+Bs(2), &
                    g_loc(3)+1:g_loc(3)+Bs(3), i_comp, hvy_id_loc)) > 1.0e-12_rk)
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, n_nonzero, nc, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, mpierr_loc)
    end subroutine

    !> norms of the reconstructed contribution in hvy_block, norms(i_field, :) = L1, L2, L3 (signed), L4, Linfty, H1, mean
    subroutine compute_norms(norms)
        implicit none
        real(kind=rk), intent(out) :: norms(1:nc,1:7)

        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L1",     norms(:,1), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L2",     norms(:,2), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L3",     norms(:,3), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "L4",     norms(:,4), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "Linfty", norms(:,5), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "H1",     norms(:,6), threshold_state_vector=.false., norm_case="leaf")
        call componentWiseNorm_tree(params, hvy_block, tree_ID, "Mean",   norms(:,7), threshold_state_vector=.false., norm_case="leaf")
    end subroutine

    !> int u_box * u_j dV over all leaf blocks, with u_box in hvy_block and u_j in the hvy_tmp slot after idx_full.
    !> Same leaf selection and volume element as componentWiseNorm_tree with norm_case="leaf",
    !> so for a box covering the whole domain the result equals the squared L2 norm.
    subroutine inner_product_full_level(inner)
        implicit none
        real(kind=rk), intent(out) :: inner(1:nc)

        integer(kind=ik) :: k_loc, hvy_id_loc, lgt_id_loc, mpierr_loc, g_loc(3), i_comp
        real(kind=rk) :: x0_loc(3), dx_loc(3)

        g_loc = 0
        g_loc(1:params%dim) = params%g

        inner = 0.0_rk
        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            ! only sum up over leafs
            if (any(hvy_family(hvy_id_loc, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) cycle

            call hvy2lgt(lgt_id_loc, hvy_id_loc, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id_loc, x0_loc, dx_loc)

            do i_comp = 1, nc
                inner(i_comp) = inner(i_comp) + product(dx_loc(1:params%dim)) * sum( &
                    hvy_block(g_loc(1)+1:Bs(1)+g_loc(1), g_loc(2)+1:Bs(2)+g_loc(2), g_loc(3)+1:Bs(3)+g_loc(3), i_comp, hvy_id_loc) * &
                    hvy_tmp(  g_loc(1)+1:Bs(1)+g_loc(1), g_loc(2)+1:Bs(2)+g_loc(2), g_loc(3)+1:Bs(3)+g_loc(3), idx_full+i_comp, hvy_id_loc) )
            enddo
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, inner, nc, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr_loc)
    end subroutine

    !> bring back the full wavelet decomposition into hvy_block: from the stored copy (norms-only), or by reading
    !> and decomposing the files again (the tree was pruned for saving), unless it is still fresh
    subroutine restore_decomposition(fresh)
        implicit none
        logical, intent(in) :: fresh

        if (norms_only) then
            call copy_block_tmp(idx_dec, to_tmp=.false.)
        elseif (.not. fresh) then
            call delete_tree(params, tree_ID)
            call createActiveSortedLists_forest(params)
            call readHDF5vct_tree( fnames, params, hvy_block, tree_ID)
            call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)
            call wavelet_decompose_full_tree(params, hvy_block, tree_ID, hvy_tmp)
        endif
    end subroutine

    !> set all coefficients outside box i_box to zero (slice-wise, the box is closed: origin <= x <= end)
    subroutine apply_local_box_filter(i_box)
        implicit none
        integer(kind=ik), intent(in) :: i_box

        integer(kind=ik) :: k_loc, hvy_id_loc, lgt_id_loc
        integer(kind=ik) :: ix_loc, iy_loc, iz_loc
        real(kind=rk) :: x0_loc(3), dx_loc(3)
        real(kind=rk) :: x, y, z
        real(kind=rk) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(kind=rk) :: lo(3), hi(3)
        logical :: outside_block

        lo = box_lo(:, i_box)
        hi = box_hi(:, i_box)

        do k_loc = 1, hvy_n(tree_ID)
            hvy_id_loc = hvy_active(k_loc, tree_ID)
            call hvy2lgt(lgt_id_loc, hvy_id_loc, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id_loc, x0_loc, dx_loc)

            xmin = x0_loc(1)
            xmax = x0_loc(1) + dble(Bs(1)-1) * dx_loc(1)
            ymin = x0_loc(2)
            ymax = x0_loc(2) + dble(Bs(2)-1) * dx_loc(2)
            zmin = x0_loc(3)
            zmax = x0_loc(3)
            if (params%dim == 3) zmax = x0_loc(3) + dble(Bs(3)-1) * dx_loc(3)

            outside_block = (xmax < lo(1) .or. xmin > hi(1) .or. &
                             ymax < lo(2) .or. ymin > hi(2))
            if (params%dim == 3) outside_block = outside_block .or. (zmax < lo(3) .or. zmin > hi(3))

            if (outside_block) then
                hvy_block(:,:,:, 1:nc, hvy_id_loc) = 0.0_rk
                cycle
            endif

            do iz_loc = merge(1, g+1, params%dim==2), merge(1, Bs(3)+g, params%dim==2)
                z = 0
                if (params%dim == 3) z = dble(iz_loc-(g+1)) * dx_loc(3) + x0_loc(3)

                if (params%dim == 3 .and. (z < lo(3) .or. z > hi(3))) then
                    hvy_block(:, :, iz_loc, 1:nc, hvy_id_loc) = 0.0_rk
                    cycle
                endif

                do iy_loc = g+1, Bs(2)+g
                    y = dble(iy_loc-(g+1)) * dx_loc(2) + x0_loc(2)

                    if (y < lo(2) .or. y > hi(2)) then
                        hvy_block(:, iy_loc, iz_loc, 1:nc, hvy_id_loc) = 0.0_rk
                        cycle
                    endif

                    do ix_loc = g+1, Bs(1)+g
                        x = dble(ix_loc-(g+1)) * dx_loc(1) + x0_loc(1)
                        if (x < lo(1) .or. x > hi(1)) then
                            hvy_block(ix_loc, iy_loc, iz_loc, 1:nc, hvy_id_loc) = 0.0_rk
                        endif
                    enddo
                enddo
            enddo
        enddo
    end subroutine

    !> several fields are read into one tree whose blocks come from the first file, so all files must have the same
    !> blocks in the same order (same treecodes and levels), not only the same number of blocks
    subroutine check_same_tree()
        implicit none
        integer(hid_t) :: file_id
        integer(kind=tsize), allocatable :: tc_first(:), tc_file(:)
        integer(kind=ik), allocatable :: lvl_first(:), lvl_file(:)
        integer(kind=ik) :: i_file, n_blocks_file
        logical :: has_treecodes

        n_blocks_file = lgt_n(tree_ID)
        allocate(tc_first(1:n_blocks_file), tc_file(1:n_blocks_file), lvl_first(1:n_blocks_file), lvl_file(1:n_blocks_file))

        do i_file = 1, nc
            call open_file_hdf5(trim(adjustl(fnames(i_file))), file_id, .false.)
            has_treecodes = dset_exists(file_id, "block_treecode_num")
            if (has_treecodes) then
                call read_dset_mpi_hdf5(file_id, "block_treecode_num", (/0_ik/), (/n_blocks_file-1/), tc_file)
                call read_dset_mpi_hdf5(file_id, "level", (/0_ik/), (/n_blocks_file-1/), lvl_file)
            endif
            call close_file_hdf5(file_id)
            ! old file versions without numerical treecodes: only the block count was checked
            if (.not. has_treecodes) return

            if (i_file == 1) then
                tc_first = tc_file
                lvl_first = lvl_file
            elseif (any(tc_file /= tc_first) .or. any(lvl_file /= lvl_first)) then
                call abort(261013_ik, "--wavelet-contributions: all input files must be on the same grid, but the blocks of "// &
                    trim(fnames(i_file))//" differ from those of "//trim(fnames(1)))
            endif
        enddo
    end subroutine

    !> read the boxes from a ';'-separated file, one box per line: id; x0; y0; [z0;] x1; y1[; z1]
    subroutine read_box_file(fname)
        implicit none
        character(len=*), intent(in) :: fname

        character(len=1024) :: line
        real(kind=rk) :: values(1:6)
        integer :: unit_box, ios, mpierr_loc
        integer(kind=ik) :: i_pass, i_line, n_found

        call check_file_exists(fname)

        n_boxes = 0
        if (params%rank == 0) then
            ! first pass counts the boxes, the second one reads them
            do i_pass = 1, 2
                open(newunit=unit_box, file=fname, status='old', action='read', iostat=ios)
                if (ios /= 0) call abort(261013_ik, "--box-file: cannot open "//fname)
                i_line = 0
                n_found = 0
                do
                    read(unit_box, '(A)', iostat=ios) line
                    if (ios /= 0) exit
                    i_line = i_line + 1
                    if (len_trim(line) == len(line)) call abort(261013_ik, "--box-file: line too long (max 1023 characters)")
                    ! tabs and the CR of Windows line ends count as blanks
                    line = str_replace_text(str_replace_text(line, char(9), " "), char(13), " ")
                    line = adjustl(line)
                    if (len_trim(line) == 0) cycle
                    if (line(1:1) == '%' .or. line(1:1) == '#') cycle
                    n_found = n_found + 1
                    if (i_pass == 2) then
                        call parse_box_line(line, i_line, box_ids(n_found), values)
                        box_lo(1:params%dim, n_found) = values(1:params%dim)
                        box_hi(1:params%dim, n_found) = values(params%dim+1:2*params%dim)
                    endif
                enddo
                close(unit_box)

                if (i_pass == 1) then
                    n_boxes = n_found
                    if (n_boxes == 0) call abort(261013_ik, "--box-file: no box found in "//fname)
                    allocate(box_lo(1:3,1:n_boxes), box_hi(1:3,1:n_boxes), box_ids(1:n_boxes))
                    box_lo = 0.0_rk
                    box_hi = 0.0_rk
                endif
            enddo
        endif

        ! only rank 0 writes, so only the coordinates are needed on all ranks
        call MPI_BCAST(n_boxes, 1, MPI_INTEGER4, 0, WABBIT_COMM, mpierr_loc)
        if (params%rank /= 0) then
            allocate(box_lo(1:3,1:n_boxes), box_hi(1:3,1:n_boxes), box_ids(1:n_boxes))
            box_ids = ""
        endif
        call MPI_BCAST(box_lo, 3*n_boxes, MPI_DOUBLE_PRECISION, 0, WABBIT_COMM, mpierr_loc)
        call MPI_BCAST(box_hi, 3*n_boxes, MPI_DOUBLE_PRECISION, 0, WABBIT_COMM, mpierr_loc)
    end subroutine

    subroutine parse_box_line(line, i_line, id, values)
        implicit none
        character(len=*), intent(in) :: line
        integer(kind=ik), intent(in) :: i_line
        character(len=*), intent(out) :: id
        real(kind=rk), intent(out) :: values(1:6)

        character(len=len(line)) :: rest
        character(len=cshort) :: msg
        integer :: p, ios
        integer(kind=ik) :: i_val, n_val

        n_val = 2*params%dim
        write(msg,'(A,I0,A,I0,A)') "--box-file line ", i_line, ": expected an id and ", n_val, " numbers, "
        values = 0.0_rk

        p = index(line, ';')
        if (p == 0) call abort(261013_ik, trim(msg)//"got: "//trim(line))
        if (len_trim(adjustl(line(1:p-1))) > len(id)) call abort(261013_ik, trim(msg)//"the id is longer than 80 characters: "//trim(line))
        id = adjustl(line(1:p-1))
        if (len_trim(id) == 0) call abort(261013_ik, trim(msg)//"the id is empty: "//trim(line))
        rest = line(p+1:)

        do i_val = 1, n_val
            p = index(rest, ';')
            if (p == 0) p = len_trim(rest) + 1
            if (len_trim(rest(1:p-1)) == 0) call abort(261013_ik, trim(msg)//"got: "//trim(line))
            ! list-directed read would also accept '1 2', '1,5' or '3*2', so allow only the characters of a number
            ios = verify(trim(adjustl(rest(1:p-1))), "0123456789+-.eEdD")
            if (ios == 0) read(rest(1:p-1), *, iostat=ios) values(i_val)
            if (ios /= 0) call abort(261013_ik, trim(msg)//"cannot read '"//trim(adjustl(rest(1:p-1)))//"' in: "//trim(line))
            if (p > len_trim(rest)) then
                rest = ""
            else
                rest = rest(p+1:)
            endif
        enddo
        ! a trailing ';' is fine, more values are not
        if (len_trim(rest) > 0) call abort(261013_ik, trim(msg)//"got more: "//trim(line))

        do i_val = 1, params%dim
            if (values(params%dim+i_val) <= values(i_val)) then
                call abort(261013_ik, "--box-file: box end must be > box origin in every direction: "//trim(line))
            endif
        enddo
    end subroutine

    subroutine open_table()
        implicit none
        integer(kind=ik) :: i_dim
        integer :: ios

        if (params%rank /= 0) return

        open(newunit=table_unit, file=trim(table_file), status='replace', action='write', iostat=ios)
        if (ios /= 0) call abort(261013_ik, "cannot write the table "//trim(table_file)//" (set another one with --table-file)")
        write(table_unit,'(A)') "% wabbit-post --wavelet-contributions: norms of the level contributions u_j (box 0) and of their box parts u_box"
        write(table_unit,'(A,A,A,es16.8,A,I0,A,I0,A,I0)') "% wavelet=", trim(params%wavelet), " time=", time, " Jmin=", params%Jmin, &
            " Jmax_active=", Jmax_active, " sub0_levels=", sub0_levels
        do i_field = 1, nc
            write(table_unit,'(A,I0,A,A)') "% field ", i_field, " = ", trim(fnames(i_field))
        enddo
        write(table_unit,'(A)') "% box 0 = full level contribution u_j (no box)"
        do i_box = 1, n_boxes
            write(table_unit,'(A,I0,A,A,A)', advance='no') "% box ", i_box, " = ", trim(box_ids(i_box)), ": origin"
            do i_dim = 1, params%dim
                write(table_unit,'(1x,es16.8)', advance='no') box_lo(i_dim, i_box)
            enddo
            write(table_unit,'(A)', advance='no') ", end"
            do i_dim = 1, params%dim
                write(table_unit,'(1x,es16.8)', advance='no') box_hi(i_dim, i_box)
            enddo
            write(table_unit,'(A)') ""
        enddo
        write(table_unit,'(A)') "% sublevel 0 = regular level contribution; level 0 with sublevel s = sub0 sub-level s (--sub0-levels)"
        write(table_unit,'(A)') "% L3 is signed ((sum u^3)^(1/3) with the sign kept); L2_full = L2 norm of u_j;"
        write(table_unit,'(A)') "% inner_box_full = int(u_box*u_j) dV (box 0: L2^2); n_nonzero = non-zero wavelet coefficients"
        if (coef_energies) then
            write(table_unit,'(A)') "% --coefficient-energies: E = sum over the wavelet coefficients of the level in the box of c^2*2^n_even*dx_j^dim"
            write(table_unit,'(A)') "% (Parseval, no reconstruction): L2 = sqrt(E), L2_full = sqrt(E of box 0), inner_box_full = E; other norms NaN."
            write(table_unit,'(A)') "% Exact energies of u_j and u_box for orthogonal wavelets only."
        endif
        write(table_unit,'(A7,";",A6,";",A6,";",A9,";",9(A24,";"),A12)') "% field", "box", "level", "sublevel", "L1", "L2", "L3", &
            "L4", "Linf", "H1", "mean", "L2_full", "inner_box_full", "n_nonzero"
    end subroutine

    subroutine write_rows(level_id, sublevel_id, box_id, norms, norm_full, inner, n_nonzero)
        implicit none
        integer(kind=ik), intent(in) :: level_id, sublevel_id, box_id
        real(kind=rk), intent(in) :: norms(1:nc,1:7), norm_full(1:nc), inner(1:nc)
        integer(kind=ik), intent(in) :: n_nonzero(1:nc)

        integer(kind=ik) :: i_comp

        if (params%rank /= 0) return

        do i_comp = 1, nc
            write(table_unit,'(I7,";",I6,";",I6,";",I9,";",9(es24.15e3,";"),I12)') i_comp, box_id, level_id, sublevel_id, &
                norms(i_comp,1:7), norm_full(i_comp), inner(i_comp), n_nonzero(i_comp)
        enddo
    end subroutine

end subroutine


integer(kind=ik) function count_div2(n)
    use module_globals
    implicit none
    integer(kind=ik), intent(in) :: n
    integer(kind=ik) :: m

    count_div2 = 0
    m = n
    do while (m > 0 .and. modulo(m, 2_ik) == 0)
        count_div2 = count_div2 + 1
        m = m / 2
    enddo
end function
