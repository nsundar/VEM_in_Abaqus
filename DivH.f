!====================================================================
        DivH(6,6) = (4.d0/dia**2.d0)*area
        DivH(6,7) = (2.d0/dia**2.d0)*s13
        DivH(6,8) = (2.d0/dia**2.d0)*s12
        DivH(6,9) = (4.d0/dia**2.d0)*s12
        DivH(6,12) = (4.d0/dia**2.d0)*s13
        DivH(7,7) = (1.d0/dia**2.d0)*s17
        DivH(7,8) = (1.d0/dia**2.d0)*s15
        DivH(7,9) = (2.d0/dia**2.d0)*s15
        DivH(7,12) = (2.d0/dia**2.d0)*s17
        DivH(8,8) = (1.d0/dia**2.d0)*s16
        DivH(8,9) = (2.d0/dia**2.d0)*s16
        DivH(8,12) = (2.d0/dia**2.d0)*s15
        DivH(9,9) = (4.d0/dia**2.d0)*s16
        DivH(9,12) = (4.d0/dia**2.d0)*s15
        DivH(12,12) = (4.d0/dia**2.d0)*s17
        
        do i = 6,12
        do j = 6,12
        DivH(j,i) = DivH(i,j)
        end do !j = 6,12
        end do ! i=6,12
      
