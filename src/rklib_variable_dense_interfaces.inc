    ! variable step dense output interfaces

    pure module function rkdp54_dense_coefficients(me) result(b)
        implicit none
        class(rkdp54_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable :: b !! dense output coefficients
    end function rkdp54_dense_coefficients

    pure module function rkdp85_dense_coefficients(me) result(b)
        implicit none
        class(rkdp85_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable :: b !! dense output coefficients
    end function rkdp85_dense_coefficients

    module subroutine rkdp85_dense_extra_stages(me)
        implicit none
        class(rkdp85_class),intent(inout) :: me
    end subroutine rkdp85_dense_extra_stages

