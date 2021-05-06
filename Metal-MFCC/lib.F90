subroutine distance(i,j,dis)

use comparm
use mpi
implicit none

integer(kind=8)::i,j,k
real(kind=8)::dis
dis=0
do k=1,3
dis=dis+(atomcrd(k,i)-atomcrd(k,j))**2
end do
dis=sqrt(dis)
end subroutine

subroutine all_sys
use comparm
use mpi
implicit none
character(len=50)::filename
gauss_tag=1
filename='ALL_sys'
call g09_job(filename)

end subroutine

subroutine adjust_charge
use comparm
use mpi

implicit none

integer(kind=8)::i,listnum,atom1,atom2,ist
real(kind=8)::delta_charge
character(len=80)::pline

if(.not.allocated(crg)) allocate(crg(natom))
crg=atomcrg
open(201,file='input')
do while(.true.)
    read(201,'(A80)',iostat=ist)pline
    if(ist.ne.0)then
        exit
    end if
    if(pline(1:25).eq.'Charge Transfer')then
    read(201,'(i5)')listnum
        do i=1,listnum
            read(201,*)atom1,atom2,delta_charge
            crg(atom1)=crg(atom1)-delta_charge
            crg(atom2)=crg(atom2)+delta_charge
        end do
    end if
end do
close(201)

end subroutine
subroutine cal_Angle(i,j,k,theta)
    use comparm
    use mpi
    implicit none
    integer(kind=8)::i,j,k,m
    real(kind=8)::theta,cos_theta,a(3),b(3),ra,rb,ab
    real(kind=8),parameter::pt999=0.999999d0
    real(kind=8),parameter::pi=3.1415926
    do m=1,3
        a(m)=atomcrd(m,i)-atomcrd(m,j)
        b(m)=atomcrd(m,k)-atomcrd(m,j)
    end do
    ra=sqrt(a(1)**2+a(2)**2+a(3)**2)
    rb=sqrt(b(1)**2+b(2)**2+b(3)**2)
    ab=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
    cos_theta=ab/(ra*rb)
    if(cos_theta.gt.pt999)then
        cos_theta=pt999
    else if(cos_theta.lt.-1*pt999)then
        cos_theta=-1*pt999
    end if
    theta = acos(cos_theta)/pi*180
    return
end subroutine

subroutine judge_H_bond
use comparm
use mpi
implicit none
integer(kind=8)::i,j,ist1,ist2
character(len=80)::pline

open(201,file=trim(adjustl(basename))//'.hbn')
do while (.true.)
    read(201,'(A80)',iostat=ist1)pline
    if(ist1.ne.0)then
        exit
    end if
    if(pline(2:25).eq.'Structure file analyzed:')then
        do while(.true.)
            read(201,'(i5,16x,i5)',iostat=ist2)i,j
            if(ist2.ne.0)then
                exit
            end if
            bond_num=bond_num+1
            H_bond(1,bond_num)=i
            H_bond(2,bond_num)=j
        end do
    end if
end do
close(201)
open(201,file='H_bond.info')
write(201,'(A)')'H Bond list'
write(201,'(i5)')bond_num
do i=1,bond_num
    write(201,'(2i5,a5,i5,a5,i5,a5,i5,a5,2i5)')i,&
        H_bond(1,i),atomname(H_bond(1,i)),H_bond(2,i),atomname(H_bond(2,i)),&
        resnum(H_bond(1,i)),resname(H_bond(1,i)),resnum(H_bond(2,i)),resname(H_bond(2,i)),&
        unit_index(H_bond(1,i)),unit_index(H_bond(2,i))
end do
close(201)
end subroutine

subroutine cal_vdw(i,j)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,k,vdw_index
real(kind=8)::F_vdw(3),vdw_energy,dis

dis=atomdis(i,j)
F_vdw=0
vdw_index=ico(ntypes*(iac(i)-1)+iac(j))
if(vdw_index.ge.0)then
    vdw_energy=(vdwa(vdw_index)/dis**12-vdwb(vdw_index)/dis**6)/627.51
    do k=1,3
    F_vdw(k)=(12*vdwa(vdw_index)/dis**13-6*vdwb(vdw_index)/dis**7)*(atomcrd(k,i)-atomcrd(k,j))/dis
    F_vdw(k)=F_vdw(k)*rfactor/627.51
!    F_vdw(k)=F_vdw(k)
    tmpforce(k,i)=tmpforce(k,i)+F_vdw(k)
    tmpforce(k,j)=tmpforce(k,j)-F_vdw(k)
    end do
else if(vdw_index.lt.0)then
    vdw_index=abs(vdw_index)
    vdw_energy=(sola(vdw_index)/dis**12-solb(vdw_index)/dis**10)/627.51
    do k=1,3
    F_vdw(k)=(12*sola(vdw_index)/dis**13-10*solb(vdw_index)/dis**11)*(atomcrd(k,i)-atomcrd(k,j))/dis
    F_vdw(k)=F_vdw(k)*rfactor/627.51
!    F_vdw(k)=F_vdw(k)
    tmpforce(k,i)=tmpforce(k,i)+F_vdw(k)
    tmpforce(k,j)=tmpforce(k,j)-F_vdw(k)
    end do
end if
!tmpenergy=tmpenergy+vdw_energy
end subroutine


