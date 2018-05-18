program makefile_demo

   ! We usually don't want ALL of our code in one file,
   ! so we break it up into different modules.
   ! Here, we're importing some subroutines from the modules called
   ! hello and world
   use hello, only:say_hello
   use world, only:say_world

   ! Always use 'implicit none'
   ! To learn why, see:
   !    http://www.personal.psu.edu/jhm/f90/statements/implicit.html
   implicit none

   ! Call the subroutines that we've imported from other modules
   ! The first subroutine prints 'hello',
   ! and the second subroutine prints 'world'
   call say_hello
   call say_world

end program makefile_demo
