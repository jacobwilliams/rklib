        pure module function rkbs32_order(me) result(p)
            implicit none
            class(rkbs32_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkbs32_order
        pure module function rkf45_order(me) result(p)
            implicit none
            class(rkf45_class),intent(in)  :: me
            integer                        :: p    !! order of the method
        end function rkf45_order
        pure module function rkck54_order(me) result(p)
            implicit none
            class(rkck54_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkck54_order
        pure module function rkdp54_order(me) result(p)
            implicit none
            class(rkdp54_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkdp54_order
        pure module function rkt54_order(me) result(p)
            implicit none
            class(rkt54_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkt54_order
        pure module function rks54_order(me) result(p)
            implicit none
            class(rks54_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rks54_order
        pure module function rkdp65_order(me) result(p)
            implicit none
            class(rkdp65_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkdp65_order
        pure module function rkc65_order(me) result(p)
            implicit none
            class(rkc65_class),intent(in)  :: me
            integer                        :: p    !! order of the method
        end function rkc65_order
        pure module function rktp64_order(me) result(p)
            implicit none
            class(rktp64_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rktp64_order
        pure module function rkv65e_order(me) result(p)
            implicit none
            class(rkv65e_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv65e_order
        pure module function rktf65_order(me) result(p)
            implicit none
            class(rktf65_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rktf65_order
        pure module function rkv65r_order(me) result(p)
            implicit none
            class(rkv65r_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv65r_order
        pure module function rkv76e_order(me) result(p)
            implicit none
            class(rkv76e_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv76e_order
        pure module function rkv76r_order(me) result(p)
            implicit none
            class(rkv76r_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv76r_order
        pure module function rkf78_order(me) result(p)
            implicit none
            class(rkf78_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf78_order
        pure module function rkdp87_order(me) result(p)
            implicit none
            class(rkdp87_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkdp87_order
        pure module function rkv87e_order(me) result(p)
            implicit none
            class(rkv87e_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv87e_order
        pure module function rkv87r_order(me) result(p)
            implicit none
            class(rkv87r_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv87r_order
        pure module function rkk87_order(me) result(p)
            implicit none
            class(rkk87_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkk87_order
        pure module function rkv78_order(me) result(p)
            implicit none
            class(rkv78_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkv78_order
        pure module function rktp75_order(me) result(p)
            implicit none
            class(rktp75_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rktp75_order
        pure module function rktmy7_order(me) result(p)
            implicit none
            class(rktmy7_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rktmy7_order
        pure module function rktp86_order(me) result(p)
            implicit none
            class(rktp86_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rktp86_order
        pure module function rkf89_order(me) result(p)
            implicit none
            class(rkf89_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkf89_order
        pure module function rkv89_order(me) result(p)
            implicit none
            class(rkv89_class),intent(in) :: me
            integer                       :: p    !! order of the method
        end function rkv89_order
        pure module function rkt98a_order(me) result(p)
            implicit none
            class(rkt98a_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkt98a_order
        pure module function rkv98e_order(me) result(p)
            implicit none
            class(rkv98e_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv98e_order
        pure module function rkv98r_order(me) result(p)
            implicit none
            class(rkv98r_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkv98r_order
        pure module function rkf108_order(me) result(p)
            implicit none
            class(rkf108_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkf108_order
        pure module function rkc108_order(me) result(p)
            implicit none
            class(rkc108_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkc108_order
        pure module function rkb109_order(me) result(p)
            implicit none
            class(rkb109_class),intent(in) :: me
            integer                        :: p    !! order of the method
        end function rkb109_order
        pure module function rks1110a_order(me) result(p)
            implicit none
            class(rks1110a_class),intent(in) :: me
            integer                         :: p    !! order of the method
        end function rks1110a_order
        pure module function rkf1210_order(me) result(p)
            implicit none
            class(rkf1210_class),intent(in) :: me
            integer                         :: p    !! order of the method
        end function rkf1210_order
        pure module function rko129_order(me) result(p)
            implicit none
            class(rko129_class),intent(in)  :: me
            integer                         :: p    !! order of the method
        end function rko129_order
        pure module function rkf1412_order(me) result(p)
            implicit none
            class(rkf1412_class),intent(in) :: me
            integer                         :: p    !! order of the method
        end function rkf1412_order
