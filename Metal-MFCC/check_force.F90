module check_comparm
implicit none
integer,allocatable::check_flag(:)
integer(kind=8),allocatable::check_atom(:)
integer(kind=8)::check_num,check_unitnum
real(kind=8),allocatable::checked_Force(:,:)
real(kind=8),allocatable::ref_Force(:,:)
integer(kind=8),allocatable::check_unit(:)
end module check_comparm

subroutine checkforce(job_type)
use comparm
use check_comparm
use mpi
integer(kind=8)::i,j,k
character(len=50)::filename
character(len=1)::job_type
real(kind=8)::tmp
if(job_type.eq.'w')then
    if(.not.allocated(checked_Force)) allocate(checked_Force(3,natom))
    if(.not.allocated(ref_Force)) allocate(ref_Force(3,natom))
    if(.not.allocated(check_flag)) allocate(check_flag(unit_num))
    open(201,file='mfcc.force')
        read(201,*)
        read(201,*)
        read(201,*)
        read(201,*)
        read(201,*)
        do i=1,natom
            read(201,'(10x,3f16.8)')(checked_Force(j,i),j=1,3)
        end do
    close(201)
end if
do i=1,nres
    check_flag=0
    do j=resstart(i),resend(i)
        do k=1,natom
            if(atomdis(j,k).le.8.0)then
                check_flag(unit_index(k))=1
            end if
        end do
    end do
    do j=2,unit_num
        if((check_flag(j-1).eq.1).and.(check_flag(j+1).eq.1))then
            if(check_flag(j).eq.0)then
            check_flag(j)=2
            end if
        end if
    end do
    gauss_tag=0
    do j=1,natom
        if(check_flag(unit_index(j)).ne.0)then
            gauss_tag(j)=1
        end if
    end do
    write(filename,'(A9,i3.3)')'CHECK_Res',i
    if(job_type.eq.'w')then
    call g09_job(filename)
    else if(job_type.eq.'r')then
    call read_data(filename)
    do j=resstart(i),resend(i)
        do k=1,3
            ref_Force(k,j)=tmpforce(k,j)
        end do
    end do
    end if
end do
if(job_type.eq.'r')then
    do i=1,natom
        tmp=0
        do j=1,3
            tmp=(checked_Force(j,i)-ref_Force(j,i))**2+tmp
        end do
        tmp=sqrt(tmp)*627.51/0.529
        if(tmp.gt.15)then
            check_num=check_num+1
            check_atom(check_num)=i
        end if
    end do
    open(201,file='CHECK_Force')
        write(201,'(i10)')check_num
        do i=1,check_num
            write(201,'(i10,2x,a4,2x,a4,2x,i5,2x,6f12.4)')&
                              check_atom(i),atomname(check_atom(i)),&
                              resname(check_atom(i)),resnum(check_atom(i)),&
                                    (checked_Force(j,check_atom(i)),j=1,3),&
                                        (ref_Force(j,check_atom(i)),j=1,3)
        end do
    close(201)
end if
end subroutine

subroutine get_check_point
use comparm
use mpi
use check_comparm
integer(kind=8)::i
character(len=80)::pline
integer::ist
open(201,file='input')
do while(.true.)
    read(201,'(A80)',iostat=ist)pline
    if(ist .ne. 0) exit
    if(pline(1:12).eq.'CHECK PONIT:')then
        read(pline,'(15x,i5)')check_num
        if(.not.allocated(check_atom)) allocate(check_atom(check_num))
        do i=1,check_num
            read(201,'(i5)')check_atom(i)
        end do
    end if
end do
close(201)
if(.not.allocated(check_unit)) allocate(check_unit(unit_num))
if(.not.allocated(check_flag)) allocate(check_flag(unit_num))
end subroutine

subroutine check_convergence(check_point,radius,job_type)
use comparm
use mpi
use check_comparm
integer(kind=8)::i,check_point
real(kind=8)::dis,radius
character(len=1)::job_type
character(len=50)::filename
check_flag=0
check_unitnum=0
check_unit=0
do i=1,natom
    dis=atomdis(i,check_point)
    if(dis.le.radius)then
        check_flag(unit_index(i))=1
    end if
end do
do i=2,unit_num-1
    if(check_flag(i-1).eq.1.and.check_flag(i+1).eq.1)then
        if(check_flag(i).eq.0)then
            check_flag(i)=2
        end if
    end if
end do
do i=1,unit_num
    if(check_flag(i).ne.0)then
        check_unitnum=check_unitnum+1
        check_unit(check_unitnum)=i
    end if
end do
gauss_tag=0
do i=1,natom
    if(check_flag(unit_index(i)).ne.0)then
        gauss_tag(i)=1
    end if
end do
write(filename,'(A3,i4.4,A4,f3.1)')'BP_',check_point,'_Rd_',radius
if(job_type.eq.'w')then
    call g09_job(filename)
else if(job_type.eq.'r')then
    call read_data(filename)
    write(*,'(A10)')'CHECK unit'
    write(*,'(15i5)')(check_unit(i),i=1,check_unitnum)
    write(*,'(F5.2,3f14.8)')radius,(tmpforce(i,check_point),i=1,3)
end if
end subroutine

subroutine creat_sp_gjf(job_type)
use comparm
use mpi
integer(kind=8)::i,j,k,cycle_num
integer(kind=8)::test_num
integer(kind=8),allocatable::test_unit(:)
character(len=80)::pline
character(len=50)::filename
character(len=1)::job_type

open(777,file='input')
do while(.true.)
    read(777,'(A)',iostat=ist)pline
    if(ist.ne.0)then
        exit
    end if
    if(pline(1:8).eq.'TEST Gjf')then
        read(777,*)cycle_num
        do i=1,cycle_num
            read(777,*)test_num
            if(.not.allocated(test_unit)) allocate(test_unit(test_num))
            read(777,*)test_unit
            gauss_tag=0
            do j=1,test_num
                do k=unit_pt(1,test_unit(j)),unit_pt(2,test_unit(j))
                    gauss_tag(k)=1
                end do
            end do
            write(filename,'(A5,i2.2)')'TEST_',i
            if(job_type.eq.'w')then
                call g09_job(filename)
            else if(job_type.eq.'r')then
                call read_data(filename)
                write(*,'(A11,10i5)')'TEST unit: ',test_unit
                write(*,'(i5,3f14.8)')i,(tmpforce(j,bad_point)*627.51/0.529,j=1,3)
            end if
            deallocate(test_unit)
        end do
    else if(pline(1:8).eq.'TEST End')then
        exit
    end if
end do
close(777)
end subroutine
