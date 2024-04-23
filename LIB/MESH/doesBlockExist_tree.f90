!> \brief given a treecode, look for the block, return .true. and its lgtID if found
! ********************************************************************************************
subroutine doesBlockExist_tree(tcBlock, exists, lgtID, dim, max_level, level, tree_ID)

    implicit none

    !-----------------------------------------------------------------
    !> block treecode we are looking for, tc binary representation
    integer(kind=tsize), intent(in)     :: tcBlock               

    !> true if block with treecode exists
    logical, intent(out)                :: exists                     
    !> light data id of block if found
    integer(kind=ik), intent(out)       :: lgtID                      

    !> level of block we search, defaults to 0
    integer(kind=ik), optional, intent(in)        :: level
    !> index of the tree we are looking at, defaults to 1 (flow)
    integer(kind=ik), optional, intent(in)        :: tree_ID
    !> Dimension of flow, 2 or 3 and defaults to 3
    integer(kind=ik), optional, intent(in)        :: dim
    !> Max level, should be params%Jmax, defaults to maxdigits
    integer(kind=ik), optional, intent(in)        :: max_level
    !-----------------------------------------------------------------

    ! variables for optionals and defaults
    integer(kind=ik)                    :: n_level, n_tree_ID, n_dim, max_tclevel
    integer(kind=ik)                    :: k, i1, i2, imid            ! loop variables
    integer(kind=tsize)                 :: num_treecode               ! unique tc-id
    integer(kind=ik)                    :: id_a_now(3)                ! array for comparisons

    ! Set default for Tree_ID
    n_tree_ID = 1; if (present(tree_ID)) n_tree_ID = tree_ID
  
    ! Set default for level and handle negative cases
    n_dim = 3; if (present(dim)) n_dim = dim
    max_tclevel = maxdigits; if (present(max_level)) max_tclevel = max_level
    n_level = 0; if (present(level)) n_level = level
    if (n_level < 0) n_level = max_tclevel + n_level + 1

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    exists = .false.
    lgtID = -1

    ! construct array of what we want to compare
    ! tree_ID and level are combined into one number for performance purposes
    ! as max_level is 31 for 2D, we shift tree_ID by 100
    ! this might be confusing but can be entangled easily
    id_a_now = (/-1, -1, n_level + n_tree_ID * 100/)
    call set_tc(id_a_now, tcBlock)

    ! ! construct array of what we want to compare
    ! id_a_now = (/-1, -1, n_level, n_tree_ID/)
    ! call set_tc(id_a_now, tcBlock)

    !> 2nd: binary search. start with the entire interval, then choose either right or left half
    i1 = 1
    i2 = lgt_n(n_tree_ID)

    ! if the interval is only two entries, check both and return the winner
    do while ( abs(i2-i1) >= 2 )
        ! cut interval in two parts
        imid = (i1+i2) / 2
        if (tc_id_lower(lgt_sortednumlist(2:4, imid,n_tree_ID), id_a_now, .true.)) then
            i1 = imid
        else
            i2 = imid
        end if
    end do

    if ( all(lgt_sortednumlist(2:4, i1,n_tree_ID) == id_a_now) .and. lgt_sortednumlist(1, i1,n_tree_ID) > 0) then
        ! found the block we're looking for
        exists = .true.
        lgtID = int( lgt_sortednumlist(1,i1,n_tree_ID), kind=ik)
        return
    end if

    if ( all(lgt_sortednumlist(2:4, i2,n_tree_ID) == id_a_now) .and. lgt_sortednumlist(1, i2,n_tree_ID) > 0) then
        ! found the block we're looking for
        exists = .true.
        lgtID = int( lgt_sortednumlist(1,i2,n_tree_ID), kind=ik)
        return
    end if

    ! write(*, '("treecode - id=", i0)') lgtID

end subroutine doesBlockExist_tree