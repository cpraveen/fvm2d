      subroutine check_string (str1, str2)
      implicit none
      character str1*(*), str2*(*)

      if(str1.ne.str2)then
         print*,'Expecting input ', str1, ' but found ', str2
         stop
      endif

      end
