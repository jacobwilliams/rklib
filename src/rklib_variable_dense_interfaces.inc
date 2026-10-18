    ! variable step dense output interfaces

    pure module subroutine rkdp54_dense_coefficients(me,b,form)
        implicit none
        class(rkdp54_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkdp54_dense_coefficients

    pure module subroutine rkv65e_dense_coefficients(me,b,form)
        implicit none
        class(rkv65e_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv65e_dense_coefficients

    module subroutine rkv65e_dense_extra_stages(me)
        implicit none
        class(rkv65e_class),intent(inout) :: me
    end subroutine rkv65e_dense_extra_stages

    pure module subroutine rkv65r_dense_coefficients(me,b,form)
        implicit none
        class(rkv65r_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv65r_dense_coefficients

    module subroutine rkv65r_dense_extra_stages(me)
        implicit none
        class(rkv65r_class),intent(inout) :: me
    end subroutine rkv65r_dense_extra_stages

    pure module subroutine rkv76e_dense_coefficients(me,b,form)
        implicit none
        class(rkv76e_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv76e_dense_coefficients

    module subroutine rkv76e_dense_extra_stages(me)
        implicit none
        class(rkv76e_class),intent(inout) :: me
    end subroutine rkv76e_dense_extra_stages

    pure module subroutine rkv76r_dense_coefficients(me,b,form)
        implicit none
        class(rkv76r_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv76r_dense_coefficients

    module subroutine rkv76r_dense_extra_stages(me)
        implicit none
        class(rkv76r_class),intent(inout) :: me
    end subroutine rkv76r_dense_extra_stages

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

    pure module subroutine rkv87e_dense_coefficients(me,b,form)
        implicit none
        class(rkv87e_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv87e_dense_coefficients

    module subroutine rkv87e_dense_extra_stages(me)
        implicit none
        class(rkv87e_class),intent(inout) :: me
    end subroutine rkv87e_dense_extra_stages

    pure module subroutine rkv87r_dense_coefficients(me,b,form)
        implicit none
        class(rkv87r_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv87r_dense_coefficients

    module subroutine rkv87r_dense_extra_stages(me)
        implicit none
        class(rkv87r_class),intent(inout) :: me
    end subroutine rkv87r_dense_extra_stages

    pure module subroutine rkv98e_dense_coefficients(me,b,form)
        implicit none
        class(rkv98e_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv98e_dense_coefficients

    module subroutine rkv98e_dense_extra_stages(me)
        implicit none
        class(rkv98e_class),intent(inout) :: me
    end subroutine rkv98e_dense_extra_stages

    pure module subroutine rkv98r_dense_coefficients(me,b,form)
        implicit none
        class(rkv98r_class),intent(in) :: me
        real(wp),dimension(:,:),allocatable,intent(out) :: b !! dense output coefficients
        integer,intent(out) :: form !! the form of `b`
    end subroutine rkv98r_dense_coefficients

    module subroutine rkv98r_dense_extra_stages(me)
        implicit none
        class(rkv98r_class),intent(inout) :: me
    end subroutine rkv98r_dense_extra_stages

