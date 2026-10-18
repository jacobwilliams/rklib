    ! variable step dense output interfaces

    pure module subroutine rkdp54_dense_coefficients(me,b,form)
        implicit none
        class(rkdp54_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkdp54_dense_coefficients

    pure module subroutine rkdp85_dense_coefficients(me,b,form)
        implicit none
        class(rkdp85_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkdp85_dense_coefficients

    module subroutine rkdp85_dense_extra_stages(me)
        implicit none
        class(rkdp85_class),intent(inout) :: me
    end subroutine rkdp85_dense_extra_stages

