pure module function euler_properties(me) result(p)
    implicit none
    class(euler_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function euler_properties

pure module function midpoint_properties(me) result(p)
    implicit none
    class(midpoint_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function midpoint_properties

pure module function heun_properties(me) result(p)
    implicit none
    class(heun_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function heun_properties

pure module function rkssp22_properties(me) result(p)
    implicit none
    class(rkssp22_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkssp22_properties

pure module function rk3_properties(me) result(p)
    implicit none
    class(rk3_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk3_properties

pure module function rkssp33_properties(me) result(p)
    implicit none
    class(rkssp33_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkssp33_properties

pure module function rkssp53_properties(me) result(p)
    implicit none
    class(rkssp53_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkssp53_properties

pure module function rk4_properties(me) result(p)
    implicit none
    class(rk4_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk4_properties

pure module function rks4_properties(me) result(p)
    implicit none
    class(rks4_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rks4_properties

pure module function rkr4_properties(me) result(p)
    implicit none
    class(rkr4_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkr4_properties

pure module function rkls44_properties(me) result(p)
    implicit none
    class(rkls44_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkls44_properties

pure module function rkls54_properties(me) result(p)
    implicit none
    class(rkls54_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkls54_properties

pure module function rkssp54_properties(me) result(p)
    implicit none
    class(rkssp54_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkssp54_properties

pure module function rks5_properties(me) result(p)
    implicit none
    class(rks5_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rks5_properties

pure module function rk5_properties(me) result(p)
    implicit none
    class(rk5_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk5_properties

pure module function rkc5_properties(me) result(p)
    implicit none
    class(rkc5_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkc5_properties

pure module function rkl5_properties(me) result(p)
    implicit none
    class(rkl5_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkl5_properties

pure module function rklk5a_properties(me) result(p)
    implicit none
    class(rklk5a_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rklk5a_properties

pure module function rklk5b_properties(me) result(p)
    implicit none
    class(rklk5b_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rklk5b_properties

pure module function rkb6_properties(me) result(p)
    implicit none
    class(rkb6_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkb6_properties

pure module function rk7_properties(me) result(p)
    implicit none
    class(rk7_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk7_properties

pure module function rk8_10_properties(me) result(p)
    implicit none
    class(rk8_10_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk8_10_properties

pure module function rkcv8_properties(me) result(p)
    implicit none
    class(rkcv8_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkcv8_properties

pure module function rk8_12_properties(me) result(p)
    implicit none
    class(rk8_12_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rk8_12_properties

pure module function rks10_properties(me) result(p)
    implicit none
    class(rks10_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rks10_properties

pure module function rkz10_properties(me) result(p)
    implicit none
    class(rkz10_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkz10_properties

pure module function rko10_properties(me) result(p)
    implicit none
    class(rko10_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rko10_properties

pure module function rkh10_properties(me) result(p)
    implicit none
    class(rkh10_class),intent(in) :: me
    type(rklib_properties) :: p !! properties of the method
end function rkh10_properties

