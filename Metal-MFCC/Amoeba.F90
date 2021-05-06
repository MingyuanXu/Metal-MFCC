module amoeba
use mpi
use comparm

integer(kind=8),allocatable::Amber2Amoeba(:),Amoeba2Amber(:),Amoeba_list(:,:)
real(kind=8),allocatable::   Amoeba_Force(:,:)
real(kind=8)::Amoeba_Energy
character(len=80)::Fullsys,subsys,Amoebaprm,Binpath
namelist /amoeba_ctrl/ Fullsys,subsys,Amoebaprm,Binpath

end module

subroutine Translate_atomindex_To_Amoeba
use comparm
use amoeba 
use mpi
implicit none

integer(kind=8)::i,j
character(len=80)::filename
character(len=200)::cmd
logical::fexist
real(kind=8),allocatable::tmp(:,:)

read(*,nml=amoeba_ctrl)
cmd='rm '//trim(adjustl(Fullsys))//'.xyz* '//trim(adjustl(Fullsys))//'.seq*'
call system(cmd)
cmd='rm '//trim(adjustl(subsys))//'.xyz* '//trim(adjustl(subsys))//'.seq*' 
call system(cmd)
!cmd='ambpdb -p '//trim(adjustl(Fullsys))//'.parm7 < '//&
!    trim(adjustl(Fullsys))//'.crd > '//trim(adjustl(Fullsys))//'.pdb '
!call system(cmd)
cmd='ambpdb -p '//trim(adjustl(subsys))//'.parm7 < '//&
    trim(adjustl(subsys))//'.crd > '//trim(adjustl(subsys))//'.pdb '
call system(cmd)

!inquire(file=trim(adjustl(Fullsys))//'.pdb',exist=fexist)
!if(.not.fexist)then
!    write(*,'(A)')trim(adjustl(Fullsys))//".pdb doesn't exist!"
!end if

inquire(file=trim(adjustl(subsys))//'.pdb',exist=fexist)
if(.not.fexist)then
    write(*,'(A)')trim(adjustl(subsys))//".pdb doesn't exist!"
end if

filename=trim(adjustl(Fullsys))//'.xyz'
cmd=trim(adjustl(Binpath))//&
'/pdbxyz '//trim(adjustl(Fullsys))//'.pdb '//trim(adjustl(Amoebaprm))
write(*,*)cmd
call system(cmd)

filename=trim(adjustl(subsys))//'.xyz'
cmd=trim(adjustl(Binpath))//&
'/pdbxyz '//trim(adjustl(subsys))//'.pdb '//trim(adjustl(Amoebaprm))
write(*,*)cmd
call system(cmd)

filename=trim(adjustl(basename))//'.trans2Amoeba'
inquire(file=trim(adjustl(filename)),exist=fexist)
if(fexist)then
    if(.not.allocated(Amber2Amoeba)) allocate(Amber2Amoeba(natom))
    if(.not.allocated(Amoeba2Amber)) allocate(Amoeba2Amber(natom))
    open(501,file=trim(adjustl(filename))) 
        do i=1,natom 
        read(501,'(2i10)')Amber2Amoeba(i),Amoeba2Amber(i) 
        end do
    close(501)
else 
    if(.not.allocated(Amber2Amoeba)) allocate(Amber2Amoeba(natom))
    if(.not.allocated(Amoeba2Amber)) allocate(Amoeba2Amber(natom))
    filename=trim(adjustl(subsys))//'.xyz'
   ! write(*,*)filename
    if(.not.allocated(tmp)) allocate(tmp(3,natom))
    open(501,file=filename)
        read(501,*)
        do i=1,natom
            read(501,'(11x,3f12.6)')tmp(1,i),tmp(2,i),tmp(3,i)
        end do
    close(501)
    write(*,*)tmp(1,1),tmp(2,1),tmp(3,1)
    write(*,*)atomcrd(1,2),atomcrd(2,2),atomcrd(3,2)
    do i=1,natom
        do j=1,natom
            if((abs(atomcrd(1,i)-tmp(1,j)).le.0.001).and.(abs(atomcrd(2,i)-tmp(2,j)).le.0.001)&
            .and.(abs(atomcrd(3,i)-tmp(3,j)).le.0.001))then
           ! write(*,*)i,j
                Amber2Amoeba(i)=j
                Amoeba2Amber(j)=i
            end if
        end do
    end do
    filename=trim(adjustl(basename))//'.trans2Amoeba'
    open(501,file=trim(adjustl(filename)))
        do i=1,natom
            write(501,'(2i10)')Amber2Amoeba(i),Amoeba2Amber(i)
        end do
    close(501)
end if
end subroutine

subroutine print_GMFCC_pairs
use comparm
use amoeba
use mpi
integer(kind=8)::i,j

if(.not.allocated(Amoeba_list)) allocate(Amoeba_list(natom,natom))
Amoeba_list=0
open(501,file='GMFCC_pairs')
    do i=1,natom
        do j=1,natom
            if(GMFCC_pairs(i,j).eq.1)then
                Amoeba_list(Amber2Amoeba(i),Amber2Amoeba(j))=1
                Amoeba_list(Amber2Amoeba(j),Amber2Amoeba(i))=1
            end if
        end do
    end do
    write(501,'(40i2)')Amoeba_list
close(501)

end subroutine

subroutine cal_Amoeba_interaction
use comparm
use amoeba
use mpi
integer(kind=8)::i,j
character(len=80)::filename
integer(kind=4)::cmdstatus
character(len=200)::cmd
logical::fexist

if(.not.allocated(Amoeba_Force)) allocate(Amoeba_Force(3,natom))

filename=trim(adjustl(subsys))//'.key'
inquire(file=filename,exist=fexist)
if(.not.fexist)then
    open(501,file=filename)
        write(501,'(A)')' '
        write(501,'(A16,A40)')'parameters      ',trim(adjustl(Amoebaprm))
        write(501,'(A)')'verbose'
        write(501,'(A)')'randomseed       123456789'
        write(501,'(A)')'neighbor-list'
        write(501,'(A)')'polar-eps        0.001'
    close(501)
end if
filename='GMFCC_ctrl'
inquire(file=filename,exist=fexist)
if(.not.fexist)then
    open(501,file='GMFCC_ctrl')
        write(501,'(A)')'yes'
        write(501,'(i10)')natom
        write(501,'(A)')subsys
    close(501)
end if

filename=trim(adjustl(subsys))//'.run'
inquire(file=filename,exist=fexist)
if(.not.fexist)then
    open(501,file=trim(adjustl(filename)))
        write(501,'(A)')trim(adjustl(Binpath))//'/dynamic '//trim(adjustl(subsys))//&
        ' 1 1.0 1.0 2 300'
    close(501)
end if

cmd='chmod +x '//trim(adjustl(filename))
cmdstatus=system(cmd)
if (cmdstatus .ne. 0 )then
    write(*,*)'Failed to give '//trim(adjustl(filename))//' power!'
    stop 
end if
cmd='./'//trim(adjustl(filename))
cmdstatus= system(cmd)
if (cmdstatus .ne. 0 )then
    write(*,*)'Failed to run'//trim(adjustl(filename))//'!'
    stop 
end if

Amoeba_Force=0

filename=trim(adjustl(subsys))//'.Mpole'
inquire(file=filename,exist=fexist)
if(fexist)then
    open(501,file=trim(adjustl(filename)))
        read(501,'(13x,f14.8)')tmpenergy
        read(501,*)
        do i=1,natom
            read(501,'(10x,3f14.8)') (tmpforce(j,i),j=1,3)       
        end do
    close(501)
    Amoeba_Force=Amoeba_Force+tmpforce
    Amoeba_Energy=Amoeba_Energy+tmpenergy
else 
    write(*,*)'Not find '//trim(adjustl(filename))//'!'
    stop
end if

filename=trim(adjustl(subsys))//'.Polar'
inquire(file=filename,exist=fexist)
if(fexist)then
    open(501,file=trim(adjustl(filename)))
        read(501,'(13x,f14.8)')tmpenergy
        read(501,*)
        do i=1,natom
            read(501,'(10x,3f14.8)')(tmpforce(j,i),j=1,3)
        end do
    close(501)
    Amoeba_Force=Amoeba_Force+tmpforce
    Amoeba_Energy=Amoeba_Energy+tmpenergy
else 
    write(*,*)'Not find '//trim(adjustl(filename))//'!'
    stop
end if

filename=trim(adjustl(subsys))//'.vdw14-7'
inquire(file=filename,exist=fexist)
if(fexist)then
    open(501,file=trim(adjustl(filename)))
        read(501,'(11x,f14.8)')tmpenergy
        read(501,*)
        do i=1,natom
            read(501,'(10x,3f14.8)')(tmpforce(j,i),j=1,3)
        end do
    close(501)
    Amoeba_Force=Amoeba_Force+tmpforce
    Amoeba_Energy=Amoeba_Energy+tmpenergy
    Amoeba_Force=Amoeba_Force*0.529/627.51
    Amoeba_Energy=Amoeba_Energy/627.51
    do i=1,natom
        do j=1,3
            Force(j,Amoeba2Amber(i))=Force(j,Amoeba2Amber(i))-Amoeba_Force(j,i)
        end do
    end do
    Energy=Energy+Amoeba_Energy
else 
    write(*,*)'Not find '//trim(adjustl(filename))//'!'
    stop
end if
end subroutine
