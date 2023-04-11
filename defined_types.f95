module defined_types
  type :: mat2
     real(8), allocatable, dimension(:,:) :: mat2
  end type mat2
  type :: vec
     real(8), allocatable, dimension(:) :: vec
  end type vec
  type :: mat3
     real(8), allocatable, dimension(:,:,:) :: mat3
  end type mat3
  type :: vecmat2
     type(mat2), allocatable, dimension(:) :: vecmat2
  end type vecmat2
  type :: vecvec
     type(vec), allocatable, dimension(:) :: vecvec
  end type vecvec
  type :: vecvecmat2
     type(vecmat2), allocatable, dimension(:) :: vecvecmat2
  end type vecvecmat2
  type :: vecvecvec
     type(vecvec), allocatable, dimension(:) :: vecvecvec
  end type vecvecvec
  type :: vecmat3
     type(mat3), allocatable, dimension(:) :: vecmat3
  end type vecmat3

end module defined_types
