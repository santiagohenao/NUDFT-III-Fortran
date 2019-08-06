program NUDFT
implicit none

    integer, parameter :: general_kind=4

    real(kind=general_kind), parameter :: pi=atan(1.)*4.0


	real(kind=general_kind), parameter :: T_low=0.1,T_hi=10.0,dT=0.0001
    real(kind=general_kind), allocatable :: dat(:,:),x(:),y(:),T(:)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call import_OGLE_data(5,dat)

    allocate(x(size(dat(:,1))),y(size(dat(:,1))))
	x=dat(:,1)
	y=dat(:,2)
    !normalize101
    y=(y-mean(y))/std(y)

    print*,NFT(5.612)


contains

    real(kind=general_kind) function mean(arr)
        implicit none
        real(kind=general_kind), intent(in), dimension(:) :: arr

        mean=sum(arr)/size(arr)
    end function mean

    real(kind=general_kind) function std(arr)
        implicit none
        real(kind=general_kind), intent(in), dimension(:) :: arr

        std=sqrt(sum((arr-mean(arr))**2)/(size(arr)-1))
    end function std

    subroutine print_array(arr,arr_size)
        implicit none
        real(kind=general_kind), intent(in) :: arr(*)
        integer, intent(in) :: arr_size
        integer :: i
        do i = 1,arr_size
            print*,arr(i)
        end do
    end subroutine print_array

    integer function count_file_lines(filename)
        implicit none
        character (len=21), intent(in) :: filename

        count_file_lines=0
        open(50,file=filename)
        do
            read(50,*,end=10)
            count_file_lines=count_file_lines+1
        end do
        10 close(50)
    end function count_file_lines

    character(len=4) function str4(k)
        implicit none
        ! left-pad an integer with 4 zeros
        integer, intent(in) :: k
        write (str4, '(I4.4)') k
    end function str4

    character(len=21) function OGLE_file_string(k)
        implicit none
        ! ogle filename generator
        integer, intent(in) :: k
        OGLE_file_string="OGLE-LMC-CEP-"//str4(k)//".dat"
    end function OGLE_file_string

    subroutine import_OGLE_data(k,M)
        implicit none
        integer, intent(in) :: k
		character (len=21) :: filename
		real(kind=general_kind), allocatable :: M(:,:)
		integer i,row,col

        filename=OGLE_file_string(k)

		row=count_file_lines(filename)
		col=3 ! specific for OGLE files
		allocate(M(row,col))

		open(50,file=filename)
		do i=1,row
			read(50,*) M(i,:)
		end do
		close(50)
	end subroutine import_OGLE_data

    complex(kind=general_kind) function NFT(T_)
        implicit none
        ! complex nonuniform fourier transform
        real(kind=general_kind), intent(in) :: T_
        NFT=sum(y*exp(cmplx(0.0,-2*pi)*x/T_))
    end function NFT




end program NUDFT
