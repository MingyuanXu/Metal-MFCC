program Metal_MB
use comparm
use mpi
use check_comparm

implicit none
integer(kind=4)::i
integer(kind=8)::m,n
real(kind=8)::radius
character(len=50)::result_file
character(len=80)::sort_cmd

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

if(myid.eq.0)then
    read(*,nml=qm_ctrl)
    read(*,nml=job_ctrl)
    result_file=trim(adjustl(basename))//'.result'
!    open(101,file=result_file)
    open(101,file='mfcc.out')
    open(103,file='Structure_correction')
    open(102,file=trim(adjustl(basename))//'.cmd')
    call init_system
    call cut_unit
    call QM_center
    call cut_fragment
    call adjust_Fragment
    if(trim(adjustl(check_conver)).eq.'Single')then
        call get_check_point
        write(*,*)check_num
        do m=1,check_num
            write(*,*)'Check the Force of atom :',check_atom(m)
            do n=4,9
                radius=n*1.0d0
                call check_convergence(check_atom(m),radius,'w')
            end do
        end do
!        do i=4,9
!            radius=i*1.0d0
!            call check_convergence(bad_point,radius,'w')
!        end do
!        call creat_sp_gjf('w')
    else if(trim(adjustl(check_conver)).eq.'All')then
        call checkforce('w')
    else 
        if(if_fullQM.eq.'yes')then
            call Full_QM_calculation('w')
        else 
            if(trim(adjustl(Frag_method)).eq.'GMFCC_S')then
            call unitment_job2('w')
            else
            call unitment_job('w')
            end if
!!!  !       call MFCC_B2_correct('w')
            call Body2_interaction('w')
            if(b2_cutoff.eq.0)then
                call H_Bond_correction('w')
            end if
            if(center_atom.ne.0)then
                call QM_correct('w')
            end if
            call all_sys
            close(102)
        end if
    end if
    if(trim(adjustl(software)).eq.'g09'.or.trim(adjustl(software)).eq.'g16')then
        sort_cmd='sort -n -k 3 -r '//trim(adjustl(basename))//&
            ".cmd > sort.cmd"
!        sort_cmd='cat '//trim(adjustl(basename))//&
!            ".cmd|awk '{print $1,$2}'>sort.cmd"
    else if(trim(adjustl(software)).eq.'orca')then
        sort_cmd='sort -n -k 3 -r '//trim(adjustl(basename))//&
            ".cmd  > sort.cmd"
    else if(trim(adjustl(software)).eq.'dftb+')then
        sort_cmd='sort -n -k 5 -r '//trim(adjustl(basename))//&
            ".cmd  > sort.cmd"
    end if
    call system(sort_cmd)
    call arrange_task
    open(102,file='GMFCC_pairs')
        write(102,'(40i2)')GMFCC_pairs
    close(102)
end if
if(trim(adjustl(if_fullQM)).ne.'yes')then
    call mpi_barrier(mpi_comm_world,ierr)
    call mpi_bcast(ncmd,1,mpi_integer4,0,mpi_comm_world,ierr)
    call mpi_bcast(task_id,ncmd,mpi_integer4,0,mpi_comm_world,ierr)
    call mpi_bcast(task_cmd,ncmd*256,mpi_character,0,mpi_comm_world,ierr)
    call mpi_bcast(ifrun, 30, mpi_character, 0, mpi_comm_world, ierr) 
    if(trim(adjustl(ifrun)).eq.'yes')then
        if(stable_mode.eq.0)then
!            do i=myid+1,ncmd-1,nproc
!               call system(cmd(i))
!               write(*,*)myid,cmd(i)
!           end do
!        else if(stable_mode.eq.1)then
!            do i=myid+1,ncmd-1,nproc
!                call system(cmd(2*i+1))
!                call system(cmd(2*i+2))
!            end do
            do i=1,ncmd
                if(task_id(i).eq.myid)then
                    write(*,*)task_cmd(i)
                    call system(task_cmd(i))
                end if
            end do
        end if
    end if
else 
    if(myid.eq.0)then
        if(trim(adjustl(ifrun)).eq.'yes')then
            if(trim(adjustl(software)).eq.'g09'.or.trim(adjustl(software)).eq.'g16')then
                call system(trim(adjustl(software))//' Full_QM.gjf')
            else if(trim(adjustl(software)).eq.'dftb+')then
                call system('cd Full_QM && dftb+ >Full_QM.log')
            end if
        end if
    end if
end if

call mpi_barrier(mpi_comm_world,ierr)
!=====================================
!read_data
!=====================================
if(myid.eq.0)then
    if(trim(adjustl(check_conver)).eq.'Single')then
!        write(*,'(A11,i5,A20)')'Bad point: ',bad_point,' convergence test'
!        do i=4,9
!            radius=i*1.0d0
!            call check_convergence(bad_point,radius,'r')
!        end do
!        write(*,'(A)')'Additional Test:'
!        call creat_sp_gjf('r')
        do m=1,check_num
            do n=4,9
                radius=n*1.0d0
                call check_convergence(check_atom(m),radius,'r')
            end do
        end do
    else if(trim(adjustl(check_conver)).eq.'All')then
        call checkforce('r')
    else
        if(trim(adjustl(if_fullQM)).eq.'yes')then
            call Full_QM_calculation('r')
        else 
            if(trim(adjustl(Frag_method)).eq.'GMFCC_S')then
                call unitment_job2('r')
            else
                call unitment_job('r')
                write(*,*)'hello'
            end if
            if(center_atom.ne.0)then
                call QM_correct('r')
            end if
            call Body2_interaction('r')
            if(b2_cutoff.eq.0)then
                call H_Bond_correction('r')
            end if
            call Energy_Force_cal
        end if
        call print_result
        close(101)
        close(103)
    end if
end if

call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)

end program

subroutine print_result
use comparm
use mpi
implicit none
integer(kind=8)::i,j
!write(101,'(a)')'Energy:'
!write(101,'(f20.10)')energy
write(101,'(a)')'Force:'
do i=1,natom
    write(101,'(i10,3f16.8)')i,(force(j,i),j=1,3)
end do
end subroutine

subroutine Full_QM_calculation(job_type)
use comparm
use mpi
implicit none
character(len=1)::job_type
character(len=50)::filename
integer(kind=8)::i,j
gauss_tag=1
filename='Full_QM'

if(job_type.eq.'w')then
call G09_job(filename)
else if(job_type.eq.'r')then
if(.not.allocated(Force)) allocate(Force(3,natom))
if(.not.allocated(bg_force)) allocate(bg_force(3,natom))
Energy=0
Force=0
call read_data(filename)
Energy=tmpenergy
Force=tmpforce
if(trim(adjustl(add_bg)).eq.'yes')then
    bg_force=tmp_bgforce
    if(add_bg.eq.'yes')then
        open(701,file='BG_Force.dat')
            do i=1,bg_num
                write(701,'(i10,3f14.8)')i,(bg_force(j,i),j=1,3)
            end do
        close(701)
    end if
end if
write(101,'(A16)')'Total energy is:'
write(101,'(F16.8)')Energy
end if
end subroutine
