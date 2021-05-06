subroutine cut_unit
use comparm
use mpi
implicit none
integer(kind=8)::i,j,Xstart,Xend,tmp
real(kind=8)::dis
logical::test
if(.not. allocated(restype)) allocate(restype(nres))
if(.not. allocated(unit_pt)) allocate(unit_pt(2,3*nres))
if(.not. allocated(unit_index)) allocate(unit_index(natom))
if(.not. allocated(unit_charge)) allocate(unit_charge(3*nres))
if(.not. allocated(unit_type))  allocate(unit_type(3*nres))

unit_num=0
restype=' '
unit_pt=0
unit_index=0
unit_charge=0
Xstart=0
Xend=0
unit_type=' '

do i=1,chain_num
    do j=chain_start(i),chain_end(i)
        restype(j)='P'
    end do
    if(residue(chain_start(i)).eq.'ACE')then
        unit_num=unit_num+1
        unit_pt(1,unit_num)=resstart(chain_start(i))
        if(residue(chain_start(i)+1).eq.'PRO')then
            unit_pt(2,unit_num)=selectC(chain_start(i)+1)-1
            Xstart=unit_pt(2,unit_num)+1
            unit_charge(unit_num)=0
            unit_type(unit_num)='P'
        else
            unit_pt(2,unit_num)=selectCA(chain_start(i)+1)-1
            Xstart=unit_pt(2,unit_num)+1
            unit_charge(unit_num)=0
            unit_type(unit_num)='H'
        end if
        Xend=selectC(chain_end(i)-1)-1
        do j=Xstart,Xend
            if((atomname(j).eq.'CA').and.(residue(resnum(j)).ne.'PRO'))then
                unit_num=unit_num+1
                unit_pt(1,unit_num)=selectCA(resnum(j))
                unit_pt(2,unit_num)=selectC(resnum(j))-1
                unit_charge(unit_num)=int_rescrg(resnum(j))
                unit_type(unit_num)='R'
            else if(atomname(j).eq.'C')then
                if(residue(resnum(j)+1).eq.'PRO')then
                    unit_num=unit_num+1
                    unit_pt(1,unit_num)=j
                    unit_pt(2,unit_num)=selectC(resnum(j)+1)-1
                    unit_charge(unit_num)=0
                    unit_type(unit_num)='S'
                else if(residue(resnum(j)+1).ne.'PRO')then
                    unit_num=unit_num+1
                    unit_pt(1,unit_num)=j
                    unit_pt(2,unit_num)=selectCA(resnum(j)+1)-1
                    unit_charge(unit_num)=0
                    unit_type(unit_num)='B'
                end if
            end if
        end do
        unit_num=unit_num+1
        unit_pt(1,unit_num)=Xend+1
        unit_pt(2,unit_num)=resend(chain_end(i))
        unit_charge(unit_num)=0
        unit_type(unit_num)='E'
    else 
        unit_num=unit_num+1
        unit_pt(1,unit_num)=resstart(chain_start(i))
        unit_pt(2,unit_num)=selectC(chain_start(i))-1
        unit_charge(unit_num)=int_rescrg(chain_start(i))
        unit_type(unit_num)='H'
        Xstart=unit_pt(2,unit_num)+1
        if(residue(chain_end(i)).eq.'PRO')then
            Xend=selectC(chain_end(i)-1)
        else 
            Xend=selectCA(chain_end(i))-1
        end if
        do j=Xstart,Xend
            if((atomname(j).eq.'CA').and.(residue(resnum(j)).ne.'PRO'))then
                unit_num=unit_num+1
                unit_pt(1,unit_num)=selectCA(resnum(j))
                unit_pt(2,unit_num)=selectC(resnum(j))-1
                unit_charge(unit_num)=int_rescrg(resnum(j))
                unit_type(unit_num)='R'
            else if(atomname(j).eq.'C')then
                if(residue(resnum(j)+1).eq.'PRO')then
                    unit_num=unit_num+1
                    unit_pt(1,unit_num)=j
                    unit_pt(2,unit_num)=selectC(resnum(j)+1)-1
                    unit_charge(unit_num)=0
                    unit_type(unit_num)='S'
                else if(residue(resnum(j)+1).ne.'PRO')then
                    unit_num=unit_num+1
                    unit_pt(1,unit_num)=j
                    unit_pt(2,unit_num)=selectCA(resnum(j)+1)-1
                    unit_charge(unit_num)=0
                    unit_type(unit_num)='B'
                end if
            end if
        end do 
        unit_num=unit_num+1
        unit_pt(1,unit_num)=Xend+1
        unit_pt(2,unit_num)=resend(chain_end(i))
        unit_charge(unit_num)=int_rescrg(chain_end(i))
        unit_type(unit_num)='E'
    end if
end do

do i=1,nres
    if(restype(i).ne.'P')then
        restype(i)='L'
        unit_num=unit_num+1
        unit_pt(1,unit_num)=resstart(i)
        unit_pt(2,unit_num)=resend(i)
        unit_charge(unit_num)=int_rescrg(i)
        unit_type(unit_num)='L'
        do j=1,special_resnum
            if(i.eq.special_res(j))then
                unit_type(unit_num)='T'
            end if
        end do
    end if
end do
tmp=unit_num
combine_num=0
combine_unit=0
test=.true.
if(test)then
do i=1,tmp-2
    do j=i+2,tmp
        if(unit_type(i).eq.'R'.and.unit_type(j).eq.'R')then
            dis=0
            call cal_unitdis(i,j,dis)
            if(dis.le.1.5)then
                if(unit_charge(i).ne.0.and.unit_charge(j).ne.0)then
                    combine_num=combine_num+1
                    unit_num=unit_num+1
                    unit_charge(i)=0
                    unit_pt(1,i)=unit_pt(1,i)
                    unit_pt(2,i)=selectCG(resnum(unit_pt(1,i)))-1
                    unit_type(i)='M'
                    !write(*,*)i,unit_pt(1,i),unit_pt(2,i),unit_charge(i),unit_type(i)
                    unit_charge(unit_num)=int_rescrg(resnum(unit_pt(1,i)))
                    unit_pt(1,unit_num)=selectCG(resnum(unit_pt(1,i)))
                    unit_pt(2,unit_num)=selectC(resnum(unit_pt(1,i)))-1
                    unit_type(unit_num)='N'
                    !write(*,*)unit_num,unit_pt(1,unit_num),unit_pt(2,unit_num),&
                    !unit_charge(unit_num),unit_type(unit_num)
                    combine_unit(1,combine_num)=unit_num
                    unit_num=unit_num+1
                    unit_charge(j)=0
                    unit_pt(1,j)=unit_pt(1,j)
                    unit_pt(2,j)=selectCG(resnum(unit_pt(1,j)))-1
                    unit_type(j)='M'
                    !write(*,*)j,unit_pt(1,j),unit_pt(2,j),unit_charge(j),unit_type(j)
                    unit_charge(unit_num)=int_rescrg(resnum(unit_pt(1,j)))
                    unit_pt(1,unit_num)=selectCG(resnum(unit_pt(1,j)))
                    unit_pt(2,unit_num)=selectC(resnum(unit_pt(1,j)))-1
                    unit_type(unit_num)='N'
                    !write(*,*)unit_num,unit_pt(1,unit_num),unit_pt(2,unit_num),unit_charge(unit_num),unit_type(unit_num)
                    combine_unit(2,combine_num)=unit_num
                    !write(*,*)combine_unit(1,combine_num),combine_unit(2,combine_num)
                end if
            end if
        end if
    end do
end do
end if

open(501,file='unit_info')
    write(501,'(A)')'FRAG_INFORMATION'
    do i=1,unit_num
        write(501,'(i10,2X,2A4,4x,2I8,4x,2A4,i3,a2)')&
              i,resname(unit_pt(1,i)),resname(unit_pt(2,i)),&
              unit_pt(1,i),unit_pt(2,i),atomname(unit_pt(1,i)),atomname(unit_pt(2,i)),&
              unit_charge(i),unit_type(i)
    end do
close(501)
do i=1,unit_num
    do j=unit_pt(1,i),unit_pt(2,i)
        unit_index(j)=i
    end do
end do
if(.not. allocated(gauss_unit))              allocate(gauss_unit(unit_num))
if(.not.allocated(unit_connect))            allocate(unit_connect(unit_num,unit_num))
end subroutine 

subroutine QM_center
use comparm
use mpi
implicit none

integer(kind=8)::i,j,k,tmp
real(kind=8)::dis
logical::tmp_tag(unit_num)
character(len=80)::pline

if(.not.allocated(QM_tag)) allocate(QM_tag(unit_num))
if(.not.allocated(coordinate_atom)) allocate(coordinate_atom(natom))
coordinate_atom=0
QM_tag=0
tmp=0
tmp_tag=.false.
if(center_atom.ne.0)then
do i=1,unit_num
    do j=unit_pt(1,i),unit_pt(2,i)
        dis=atomdis(j,center_atom)
        if(dis.le.2.6d0)then
            tmp_tag(i)=.true.
        end if 
        if(dis.le.2.8d0)then
            coordinate_atom(j)=1
        end if
    end do
end do

do i=1,natom
    if(coordinate_atom(i).eq.1)then
    do j=1,unit_num
        do k=unit_pt(1,j),unit_pt(2,j)
            dis=atomdis(i,k)
            if(dis.le.2.4)then
                tmp_tag(j)=.true.
            end if
        end do
    end do
    end if
end do

do i=1,natom
    if(coordinate_atom(i).eq.1)then
        dis=atomdis(i,center_atom)
        if(dis.gt.2.0)then
           coordinate_atom(i)=0
        end if
    end if
end do
!tmp_tag(12)=.true.
!tmp_tag(16)=.false.
if(tmp_tag(1).eqv..true.) then
    tmp=tmp+1
    QM_tag(i)=tmp
end if

do i=2,unit_num
    if(restype(resnum(unit_pt(1,i))).eq.'P')then
        if(tmp_tag(i).eqv..true.)then
            if(tmp_tag(i-1).eqv..true.)then
                QM_tag(i)=tmp
            else 
                tmp=tmp+1
                QM_tag(i)=tmp
            end if
        end if
    else 
        if(tmp_tag(i).eqv..true.)then
            tmp=tmp+1
            QM_tag(i)=tmp
        end if
    end if
end do

open (501,file='input')
    do while(.true.)
        read(501,'(A80)',iostat=ierr)pline
        if(ierr.ne.0)then
            exit
        end if
        if(pline(1:12).eq.'Combine_unit')then
            read(501,*)i,j
!            write(*,*)i,j,unit_index(i),unit_index(j),QM_tag(unit_index(i)),QM_tag(unit_index(j))
            if(i.le.j)then
                QM_tag(unit_index(j))=QM_tag(unit_index(i))
                do k=unit_index(j),unit_num
                    if(QM_tag(k).ne.0)then
                        QM_tag(k)=QM_tag(k)-1
                        tmp=tmp-1
                    end if
                end do
            else 
                QM_tag(unit_index(i))=QM_tag(unit_index(j))
                do k=unit_index(i)+1,unit_num
                    if(QM_tag(k).ne.0)then
                        QM_tag(k)=QM_tag(k)-1
                    end if
                end do
                tmp=tmp-1
            end if
        end if
    end do
close(501)
end if
QM_ligand=tmp
open(501,file='QM_center')
do i=1,QM_ligand
    write(501,'(A9,i10)')'QM_Fragment',i
    do j=1,unit_num
        if(QM_tag(j).eq.i)then
            write(501,'(i10)')j
        end if
    end do
end do
close(501)
end subroutine
!===============================================================
subroutine cut_fragment 
use comparm
use mpi
integer(kind=8)::i,j,k,targ,ori,sp_tmp,fragment_mask(unit_num)
real(kind=8)::dis
if(.not.allocated(fragment)) allocate(fragment(unit_num))
if(.not.allocated(frag_type)) allocate(frag_type(unit_num))

fragment=0
frag_num=0
if(trim(adjustl(frag_mode)).eq.'up')then
    frag_num=0
    do i=1,chain_num
        frag_num=frag_num+1
        fragment(unit_index(resstart(chain_start(i))))=frag_num
        do j=unit_index(resstart(chain_start(i)))+1,unit_index(resend(chain_end(i)))
            if(unit_type(j).eq.'B'.or.unit_type(j).eq.'S'.or.unit_type(j).eq.'E')then
                frag_num=frag_num+1
            end if
            fragment(j)=frag_num
        end do
    end do
else if(trim(adjustl(frag_mode)).eq.'down')then
    frag_num=0
    do i=1,chain_num
        frag_num=frag_num+1
        fragment(unit_index(resstart(chain_start(i))))=frag_num
        do j=unit_index(resstart(chain_start(i)))+1,unit_index(resend(chain_end(i)))
            if(unit_type(j).eq.'R'.or.unit_type(j).eq.'S'.or.unit_type(j).eq.'E'.or.unit_type(j).eq.'M')then
                frag_num=frag_num+1
            end if
            fragment(j)=frag_num
        end do
    end do
end if
do i=1,unit_num
    if(unit_type(i).eq.'T')then
        targ=i
        ori=i
100     if(fragment(targ).eq.0)then
            if(ori.eq.targ)then
                frag_num=frag_num+1
            end if
            fragment_mask=0
            fragment(targ)=frag_num
            do j=1,unit_num
                call cal_unitdis(targ,j,dis)
                if(dis.le.3.0.and.fragment(j).eq.0.and.j.ne.targ)then
                    fragment(j)=fragment(targ)
                    fragment_mask(j)=1
                end if
            end do
            do j=1,unit_num
                if(fragment(j).eq.frag_num.and.fragment_mask(j).eq.1)then
                    do k=1,unit_num
                        if(fragment(k).eq.0)then
                            call cal_unitdis(j,k,dis)
                            if(dis.le.3.0 .and.unit_type(j).eq.'T')then
                                fragment(k)=fragment(targ)
                            end if
                            if(dis.le.3.0 .and. unit_type(k).eq.'T')then
                                targ=k
                                goto 100
                            end if
                        end if
                    end do
                end if
            end do
        end if
    end if
end do
sp_tmp=frag_num

do i=1,unit_num
    if(unit_type(i).eq.'L'.and.fragment(i).eq.0)then
        do j=1,unit_num
            if(unit_type(j).eq.'L'.and.fragment(j).eq.0.and.i.ne.j)then
                call cal_unitdis(i,j,dis)
                if (dis.le.3.0)then
                    frag_num=frag_num+1
                    fragment(i)=frag_num
                    fragment(j)=frag_num
                    write(*,*)frag_num,i,j
                    goto 102
                end if
            end if
        end do
102    end if
end do
do i=1,unit_num
    if(unit_type(i).eq.'L'.and.fragment(i).eq.0)then
        do j=sp_tmp+1,frag_num
            do k=1,unit_num
                if(fragment(k).eq.j)then
                    call cal_unitdis(i,k,dis)
                    if(dis.le.3.0)then
                        fragment(i)=j
                        goto 101
                    end if
               end if
           end do
        end do
101 end if
end do
do i=1,unit_num
    if(unit_type(i).eq.'L'.and.fragment(i).eq.0)then
        frag_num=frag_num+1
        fragment(i)=frag_num
    end if
end do
do i=1,combine_num
    frag_num=frag_num+1
    fragment(combine_unit(1,i))=frag_num
    fragment(combine_unit(2,i))=frag_num
end do

end subroutine
!===============================================================
subroutine adjust_Fragment
use comparm
use mpi
integer(kind=8)::i,j!,tmp1,tmp2,k
if(.not.allocated(frag_type)) allocate(frag_type(frag_num))
frag_type=0
do i=1,unit_num
    if(QM_tag(i).ne.0)then
        frag_type(fragment(i))=1
        if(unit_pt(1,i).eq.center_atom)then
            frag_type(fragment(i))=2
        end if
    else if(unit_type(i).eq.'L'.or.unit_type(i).eq.'T')then
        frag_type(fragment(i))=-1
    end if
end do

open(199,file='Fragment_arrange')
do i=1,frag_num
    write(199,'(A)')'+++++++++++++'
    write(199,'(A15,i5,i2)')'Frag_index:',i,frag_type(i)
    do j=1,unit_num
        if(fragment(j).eq.i)then
            write(199,'(i10,2X,2A4,4x,2I8,4x,2A4,i3,a2)')j,resname(unit_pt(1,j)),resname(unit_pt(2,j)),&
            unit_pt(1,j),unit_pt(2,j),atomname(unit_pt(1,j)),atomname(unit_pt(2,j)),unit_charge(j),unit_type(j)            
        end if
    end do
end do
close(199)
end subroutine
!===============================================================
subroutine identify_unitment(qmstart,qmend)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,qmstart,qmend
do i=qmstart,qmend
    do j=qmstart,qmend
        if(gauss_tag(i).eq.1.and.gauss_tag(j).eq.1)then
            if(unit_index(i).ne.0.and.unit_index(j).ne.0)then
                line_connect(unit_index(i),unit_index(j))=1
            end if
        end if
    end do
end do
end subroutine
!===============================================================
subroutine G09_job(filename)
use comparm
use mpi
implicit none
character(len=50)::filename
if(trim(adjustl(software)).eq.'g09'.or.trim(adjustl(software)).eq.'g16')then
    if(stable_mode.eq.0)then
        call G09_submit(filename,0)
    else if(stable_mode.eq.1)then
        call G09_submit(filename,1)
        call G09_submit(filename,2)
    end if
else if(trim(adjustl(software)).eq.'orca')then
    call Orca_submit(filename)
else if(trim(adjustl(software)).eq.'dftb+')then
    call dftb_submit(filename)
end if
end subroutine
subroutine read_data(filename)
use comparm
use mpi
use dftb3
implicit none
character(len=50)::filename
if(trim(adjustl(software)).eq.'g09'.or.trim(adjustl(software)).eq.'g16')then
    call gau_read_data(filename)
else if(trim(adjustl(software)).eq.'dftb+')then
    call dftb3_read_data(filename)
end if
end subroutine
!===============================================================
