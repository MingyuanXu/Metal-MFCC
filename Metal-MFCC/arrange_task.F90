subroutine arrange_task
    use comparm
    use mpi
    implicit none
    integer(kind=8),allocatable::cost(:)
!    integer(kind=8),allocatable::task_id(:)
    integer(kind=8),allocatable::proc_job(:)
    character(len=80)::pline
    integer::ierrs
    integer::max_proc,min_proc,move_job_id
    integer::i,j
    open(102,file='sort.cmd')
        ncmd=0
        do while(.true.)
            read(102,'(A)',iostat=ierrs) pline
            if(ierrs .ne. 0 ) exit
            ncmd=ncmd+1
        end do
    close(102)
!    if(.not.allocated(task_cmd)) allocate(task_cmd(ncmd))
    if(.not.allocated(cost)) allocate(cost(ncmd))
!    if(.not.allocated(task_id)) allocate(task_id(ncmd)) 
    if(.not.allocated(proc_job)) allocate(proc_job(nproc))
    open(102,file='sort.cmd')
        do i=1,ncmd
            if(trim(adjustl(software)).eq.'g09'.or.trim(adjustl(software)).eq.'g16')then
                read(102,'(A50,i20)')task_cmd(i),cost(i)
            else if(trim(adjustl(software)).eq.'dftb+')then
                read(102,'(A100,i20)')task_cmd(i),cost(i)
            end if
        end do
    close(102)
    j=0
    proc_job=0
    do i=1,ncmd
        task_id(i)=j 
        proc_job(j+1)=proc_job(j+1)+cost(i)
        j=j+1
        if(j.eq.nproc)then
            j=0
        end if
    end do
    do while(.true.) 
        max_proc=1
        min_proc=1
        do i=2,nproc
            if(proc_job(i).gt.proc_job(max_proc))then
                max_proc=i
            end if 
            if(proc_job(i).lt.proc_job(min_proc))then
                min_proc=i
            end if
        end do 
        do i=1,ncmd
            if(task_id(i).eq.(max_proc-1))then
                move_job_id=i
            end if
        end do
        if ((proc_job(max_proc)-proc_job(min_proc)).gt.cost(move_job_id))then
            task_id(move_job_id)=min_proc-1
            proc_job=0
            do i=1,ncmd
                proc_job(task_id(i)+1)=proc_job(task_id(i)+1)+cost(i)
            end do
        else if((proc_job(max_proc)-proc_job(min_proc))&
                .le.cost(move_job_id))then
            exit
        end if
    end do    
    open(102,file='Job_arrange')
    do i=1,ncmd
        write(102,'(A100,2i10)')task_cmd(i),cost(i),task_id(i)
    end do
    close(102)
    write(*,*)max_proc-1,proc_job(max_proc),min_proc-1,proc_job(min_proc),cost(move_job_id)
    do i=1,nproc
        write(*,*)i,proc_job(i)
    end do
end subroutine arrange_task

