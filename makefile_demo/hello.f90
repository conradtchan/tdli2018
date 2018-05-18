! This is how we define a module
module hello
contains
   ! This module contains the subroutine 'say_hello'
   subroutine say_hello
      ! print the text 'hello'
      print*, 'hello'
   end subroutine say_hello

end module hello
