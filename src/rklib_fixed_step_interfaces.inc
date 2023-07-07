        ! fixed step interfaces

        module subroutine euler(me,t,x,h,xf)
            implicit none
            class(euler_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine euler
        module subroutine midpoint(me,t,x,h,xf)
            implicit none
            class(midpoint_class),intent(inout)  :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine midpoint
        module subroutine heun(me,t,x,h,xf)
            implicit none
            class(heun_class),intent(inout)      :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine heun
        module subroutine rk3(me,t,x,h,xf)
            implicit none
            class(rk3_class),intent(inout)       :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk3
        module subroutine rk4(me,t,x,h,xf)
            implicit none
            class(rk4_class),intent(inout)       :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk4
        module subroutine rks4(me,t,x,h,xf)
            implicit none
            class(rks4_class),intent(inout)      :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rks4
        module subroutine rks5(me,t,x,h,xf)
            implicit none
            class(rks5_class),intent(inout)      :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rks5
        module subroutine rkb6(me,t,x,h,xf)
            implicit none
            class(rkb6_class),intent(inout)      :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rkb6
        module subroutine rk7(me,t,x,h,xf)
            implicit none
            class(rk7_class),intent(inout)       :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk7
        module subroutine rk8_10(me,t,x,h,xf)
            implicit none
            class(rk8_10_class),intent(inout)    :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk8_10
        module subroutine rk8_12(me,t,x,h,xf)
            implicit none
            class(rk8_12_class),intent(inout)    :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rk8_12
        module subroutine rkcv8(me,t,x,h,xf)
            implicit none
            class(rkcv8_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rkcv8
        module subroutine rkz10(me,t,x,h,xf)
            implicit none
            class(rkz10_class),intent(inout)     :: me
            real(wp),intent(in)                  :: t   !! initial time
            real(wp),dimension(me%n),intent(in)  :: x   !! initial state
            real(wp),intent(in)                  :: h   !! time step
            real(wp),dimension(me%n),intent(out) :: xf  !! state at time `t+h`
        end subroutine rkz10