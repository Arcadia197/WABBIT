!> \file
!> \brief Diagnostics of the time-loop cycle (refine -> project -> time step -> adapt), debug output
!> \details Written to cycle_diagnostics.t if [Debug] debug_cycle_diagnostics=1. One row per time step, with the maximum
!! and the RMS of the divergence of the velocity and the kinetic energy after each operation of the cycle:
!!     point 1  after refinement, before the Helmholtz projection
!!     point 2  after the Helmholtz projection, i.e. the field entering the time step
!!     point 3  after the time step (the field that div.t and e_kin.t describe)
!!     point 4  after adaptation; this call writes the row, with the time passed to it
!! The divergence is the one of div.t (order_discretization, ghost nodes from sync_ghosts_RHS_tree) and the kinetic energy
!! is 1/2 int |u|^2 over the active blocks, so the values of point 3 equal div.t and e_kin.t. Energies after adaptation are
!! those of the leaf grid, where the coarse values of lifted wavelets are low-pass filtered: compare energies of the points
!! 1-3 (refined grid) among each other. Divergence in solid bodies is not masked. Columns of points that were not called
!! in a step (e.g. no projection for ACM: point 2 equals point 1) are written as measured; with N_dt_per_grid > 1 point 3
!! holds the last time step on the grid. The calls cost a ghost-node sync and a divergence evaluation each, they are
!! skipped completely without the flag. Nothing is allocated: slot 1 of hvy_tmp, which is free at the four points of the
!! cycle, holds the divergence of a block.

!-------------------------------------------------------------------------------
!> \brief Create cycle_diagnostics.t with its column names
!-------------------------------------------------------------------------------
subroutine init_cycle_diagnostics(params, overwrite)
    implicit none

    type(type_params), intent(in)      :: params
    logical, intent(in)                :: overwrite   !< clear the file if it exists (as for the other .t files)

    character(len=cshort) :: header(1:13)

    if (.not. params%debug_cycle_diagnostics) return

    header = (/ character(len=cshort) :: "time", &
        "divmax_refined", "divrms_refined", "ekin_refined", "divmax_projected", "divrms_projected", "ekin_projected", &
        "divmax_stepped", "divrms_stepped", "ekin_stepped", "divmax_adapted", "divrms_adapted", "ekin_adapted" /)
    call init_t_file('cycle_diagnostics.t', overwrite, header)

end subroutine init_cycle_diagnostics


!-------------------------------------------------------------------------------
!> \brief Measure divergence and kinetic energy at one point of the cycle, write the row at the last point
!-------------------------------------------------------------------------------
subroutine cycle_diagnostics(params, hvy_block, hvy_tmp, tree_ID, point, time)
    implicit none

    type(type_params), intent(in)      :: params
    real(kind=rk), intent(inout)       :: hvy_block(:, :, :, :, :)   !< ghost nodes of the velocity are synced here
    real(kind=rk), intent(inout)       :: hvy_tmp(:, :, :, :, :)     !< slot 1 is overwritten (divergence of the block)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: point   !< 1 refined, 2 projected, 3 stepped, 4 adapted (writes the row)
    real(kind=rk), intent(in)          :: time    !< time written to the row (used at point 4)

    integer(kind=ik) :: k, hvy_id, lgt_id, Bs(1:3), g(1:3), mpicode
    real(kind=rk) :: x0(1:3), dx(1:3), div_max, sums(1:2)   ! sums: integral of div^2, kinetic energy
    real(kind=rk), save :: row(1:12) = -1.0_rk

    if (.not. params%debug_cycle_diagnostics) return

    Bs = params%Bs
    g  = params%g
    if (params%dim == 2) g(3) = 0

    call sync_ghosts_RHS_tree(params, hvy_block(:,:,:,1:params%dim,:), tree_ID)

    div_max = 0.0_rk
    sums    = 0.0_rk
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        call compute_divergence(hvy_block(:,:,:,1:params%dim,hvy_id), dx, Bs, g(1), params%order_discretization, hvy_tmp(:,:,:,1,hvy_id))
        div_max = max(div_max, maxval(abs(hvy_tmp(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), 1, hvy_id))))
        sums(1) = sums(1) + sum(hvy_tmp(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), 1, hvy_id)**2) * product(dx(1:params%dim))
        sums(2) = sums(2) + 0.5_rk * sum(hvy_block(g(1)+1:Bs(1)+g(1), g(2)+1:Bs(2)+g(2), g(3)+1:Bs(3)+g(3), 1:params%dim, hvy_id)**2) &
            * product(dx(1:params%dim))
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE, div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpicode)
    call MPI_ALLREDUCE(MPI_IN_PLACE, sums, 2, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpicode)

    ! domain cropping: the computational domain can be cropped, then its volume is not product(domain)
    row(3*point-2) = div_max
    row(3*point-1) = sqrt(sums(1) / product(params%domainSizeCropped(1:params%dim)))
    row(3*point)   = sums(2)

    if (point == 4) then
        call append_t_file('cycle_diagnostics.t', (/ time, row /))
        row = -1.0_rk
    endif

end subroutine cycle_diagnostics
