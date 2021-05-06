subroutine G09_submit(filename,filetype)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,charge,spin,H_num,electron_num,QM_atom,QM_atom_old,ist
integer(kind=4)::filetype
logical::fexist,exist_flag
character(len=50)::filename,tmpname
character(len=50)::addkw,ee_gmfcc
character(len=80)::pline
electron_num=0
QM_atom=0
exist_flag=.false.
do i=1,natom
    if(gauss_tag(i).eq.1)then
        QM_atom=QM_atom+1
    end if
end do
QM_atom_old=0
if(filetype.ne.1)then
    tmpname=trim(adjustl(filename))//'.gjf'
    inquire(file=trim(adjustl(tmpname)),exist=fexist)
    if(fexist)then
        open(202,file=trim(adjustl(tmpname)))
        do while(.true.)
            read(202,'(A80)',iostat=ist)pline
            if(ist.ne.0)then
                exit
            end if
            if(pline(1:10).eq.'QM atom:  ')then
                read(pline,'(10x,i5)')QM_atom_old
                if(QM_atom_old.eq.QM_atom)then
                    exist_flag=.true.
                    exit
                end if
            end if
        end do
        close(202)
    end if
    open(202,file=trim(adjustl(filename))//'.gjf')
else if(filetype.eq.1)then
    tmpname=trim(adjustl(filename))//'_stable.gjf'
    inquire(file=trim(adjustl(tmpname)),exist=fexist)
    if(fexist)then
        open(202,file=trim(adjustl(tmpname)))
        do while(.true.)
            read(202,'(A80)',iostat=ist)pline
            if(ist.ne.0)then
                exit
            end if
            if(pline(1:10).eq.'QM atom:  ')then
                read(pline,'(10x,i5)')QM_atom_old
                if(QM_atom_old.eq.QM_atom)then
                    exist_flag=.true.
                    exit
                end if
            end if
        end do
        close(202)
    end if
    open(202,file=trim(adjustl(filename))//'_stable.gjf')
end if
open(203,file=trim(adjustl(filename))//'.pdb')
write(202,'(A)')'%nproc='//trim(adjustl(qmproc))
write(202,'(A)')'%mem='//trim(adjustl(qmmem))
write(202,'(A)')'%chk='//trim(adjustl(filename))//'.chk'
write(202,*)
if(filetype.eq.0)then
    addkw='force '
else if(filetype.eq.1)then
    addkw=' stable=opt '
else if(filetype.eq.2)then
    addkw='force guess=read '
end if
ee_gmfcc=' charge Prop=(Field,read) '
if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
    addkw=trim(adjustl(addkw))//' '//trim(adjustl(ee_gmfcc))
end if
if(trim(adjustl(filename)).eq.'Full_QM')then
    addkw=' force '
    if(add_bg.eq.'yes')then
        addkw=trim(adjustl(addkw))//' '//trim(adjustl(ee_gmfcc))
    end if
end if
inquire(file=trim(adjustl(filename))//'.chk',exist=fexist)
if(fexist.and.exist_flag)then
    addkw=trim(adjustl(addkw))//' guess=read'
end if
if(trim(adjustl(qmmethod)).eq.'PM7'.or.trim(adjustl(qmmethod)).eq.'AM1'&
    .or.trim(adjustl(qmmethod)).eq.'PM3'.or.trim(adjustl(qmmethod)).eq.'PM6')then
write(202,'(A)')'#p '//trim(adjustl(qmmethod))//&
' nosymm '//trim(adjustl(addkw))//' '//trim(adjustl(add_word))
else 
write(202,'(A)')'#p '//trim(adjustl(qmmethod))//'/'//trim(adjustl(qmbasis))//&
' nosymm '//trim(adjustl(addkw))//' '//trim(adjustl(add_word))
end if
write(202,*)
write(202,'(A10,i5)')'QM atom:  ',QM_atom
write(202,'(A)')
if(trim(adjustl(filename)).ne.'ALL_sys')then
    do i=1,natom
        if(gauss_tag(i).eq.1)then
            do j=i,natom
                if(gauss_tag(j).eq.1)then
                    GMFCC_pairs(i,j)=1
                    GMFCC_pairs(j,i)=1
                end if
            end do
        end if
    end do
end if
charge=0
gauss_unit=0
do i=1,natom
    if(gauss_tag(i).eq.1.and.unit_index(i).ne.0)then
        gauss_unit(unit_index(i))=1
    end if
end do
do i=1,unit_num
    if(gauss_unit(i).eq.1)then
        charge=charge+unit_charge(i)
    end if
end do
spin=1
do i=1,natom
   if(gauss_tag(i).eq.1.and.elementnum(i).eq.29.and.resname(i).eq.'CU2')then
       spin=2
   else if(gauss_tag(i).eq.1.and.elementnum(i).eq.26.and.resname(i).eq.'FE3')then
       spin=2
   end if
end do
write(202,'(2i5)')charge,spin
do i=1,natom
    if(gauss_tag(i).eq.1)then
       write(202,'(A2,5x,3f14.8)')element(i),(atomcrd(j,i),j=1,3)
       write(203,299) 'ATOM  ',i, atomname(i), resname(i)(1:3), &
        resnum(i), (atomcrd(j,i),j=1,3)
        electron_num=electron_num+elementnum(i)
    end if
end do
H_num=0
do i=1,nbona
    if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).ne.1)then
        call add_H(ib(i),jb(i))
        H_num=H_num+1
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,jb(i))=1
                GMFCC_pairs(jb(i),j)=1
            end if
        end do
    else if(gauss_tag(ib(i)).ne.1.and.gauss_tag(jb(i)).eq.1)then
        call add_H(jb(i),ib(i))
        H_num=H_num+1
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,ib(i))=1
                GMFCC_pairs(ib(i),j)=1
            end if
        end do
    end if
end do
do i=1,nbonh
    if(gauss_tag(ibh(i)).eq.1.and.gauss_tag(jbh(i)).ne.1)then
        call add_H(ibh(i),jbh(i))
        H_num=H_num+1
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,jbh(i))=1
                GMFCC_pairs(jbh(i),j)=1
            end if
        end do
    else if(gauss_tag(ibh(i)).ne.1.and.gauss_tag(jbh(i)).eq.1)then
        call add_H(jbh(i),ibh(i))
        H_num=H_num+1
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,ibh(i))=1
                GMFCC_pairs(ibh(i),j)=1
            end if
        end do
    end if
end do
write(202,'(A)')
if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
    do i=1,natom
        if(gauss_tag(i).eq.0)then
            write(202,'(3f14.8,5x,f14.8)')(atomcrd(j,i),j=1,3),crg(i)
        end if
    end do
    if(add_bg.eq.'yes')then
        do i=1,bg_num
            write(202,'(3f14.8,5x,f14.8)')(bg_crd(j,i),j=1,3),bg_crg(i)
        end do
    end if
    write(202,'(A)')
    do i=1,natom
        if(gauss_tag(i).eq.0)then
            write(202,'(3f14.8)')(atomcrd(j,i),j=1,3)
        end if
    end do
    if(add_bg.eq.'yes')then
        do i=1,bg_num
            write(202,'(3f14.8)')(bg_crd(j,i),j=1,3)
        end do
    end if
    write(202,'(A)')
end if
close (202)
if(trim(adjustl(filename)).ne.'ALL_sys')then
    if(filetype.ne.1)then
        write(102,'(A50,i20)')trim(adjustl(software))//' '//trim(adjustl(filename))//'.gjf',electron_num**3
    else 
        write(102,'(A50,i20)')trim(adjustl(software))//' '//trim(adjustl(filename))//'_stable.gjf',electron_num**3
    end if
end if
299 format(a6,1x,I4,1x,a4,1x,a3,2x,I4,4x,3f8.3)
close(203)
end subroutine
!=============================================================================
subroutine add_H(i,j)
use comparm
use mpi
implicit none
integer(kind=8)::i,j
real(kind=8)::dis1,bond_length,x,y,z

dis1=atomdis(i,j)

if(elementnum(i).eq.7)then
    bond_length = 1.00000
else if (elementnum(i).eq.8)then
    bond_length = 0.99619
else if (elementnum(i).eq.6)then
    bond_length = 1.09000
else if (elementnum(i).eq.16)then
    bond_length = 1.31000
end if

x=bond_length*(atomcrd(1,j)-atomcrd(1,i))/dis1 +atomcrd(1,i)
y=bond_length*(atomcrd(2,j)-atomcrd(2,i))/dis1 +atomcrd(2,i)
z=bond_length*(atomcrd(3,j)-atomcrd(3,i))/dis1 +atomcrd(3,i)
if(trim(adjustl(software)).eq.'sqm')then
write(202,'(i8,A2,5x,3F14.8)')1,'H ',x,y,z
else 
write(202,'(A2,5x,3F14.8)')'H ',x,y,z
end if
end subroutine
!===============================================================================
subroutine gau_read_data(filename)
use comparm
use mpi
implicit none

integer(kind=8)::i,j,H_num
character(len=50)::filename
logical::   fexist,file_tag
character(len=80)::pline
real(kind=8)::HForce(3),Mass_sum,tmp(3),tmpcrg
integer(kind=4)::ist
tmpenergy=0
tmpforce=0
tmpcharge=0
if(add_bg.eq.'yes')then
    tmp_bgfield=0
    tmp_bgforce=0
end if

inquire(file=trim(adjustl(filename))//'.log',exist=fexist)

if(fexist)then
    call check_file(filename,file_tag)
    if(file_tag.eqv..false.)then
        open(202,file=trim(adjustl(filename))//'.log')
        open(201,file=trim(adjustl(filename))//'.data')
        H_num=0
        do i=1,nbona
            if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).ne.1)then
                H_num=H_num+1
!                gauss_tag(jb(i))=-1
            else if(gauss_tag(ib(i)).ne.1.and.gauss_tag(jb(i)).eq.1)then
                H_num=H_num+1
!                gauss_tag(ib(i))=-1
            end if
        end do
        do i=1,nbonh
            if(gauss_tag(ibh(i)).eq.1.and.gauss_tag(jbh(i)).ne.1)then
                H_num=H_num+1
!                gauss_tag(jbh(i))=-1
            else if(gauss_tag(ibh(i)).ne.1.and.gauss_tag(jbh(i)).eq.1)then
                H_num=H_num+1
!                gauss_tag(ibh(i))=-1
            end if
        end do

        Mass_sum=0
        do i=1,natom
            if(gauss_tag(i).eq.1)then
                Mass_sum=Mass_sum+mass(i)
            end if
        end do

        do while(.true.)
            read(202,'(A80)',iostat=ist)pline
            if(ist.ne.0)then
                exit
            end if
            if(pline(2:10).eq.'SCF Done:')then
                if(trim(adjustl(qmmethod)).eq.'HF')then
                    read(pline,'(20x,F16.8)')tmpenergy
                else
                    read(pline,'(23x,F16.8)')tmpenergy
                end if
                write(201,'(A12,4X,F20.10)')'QM_energy=  ',tmpenergy
            end if
            if(pline(2:27).eq.'Self energy of the charges')then
                read(pline,'(30x,F20.10)')tmpcenergy
                write(201,'(A16,F20.10)')'charge energy = ',tmpcenergy
            end if
            if(pline(2:19).eq.'Mulliken charges:')then
                read(202,*)
                write(201,*)'Mulliken charge:'
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        read(202,'(10x,f12.6)')tmpcharge(i)
                        write(201,'(i10,F12.6)')i,tmpcharge(i)
                    end if
                end do
            end if
            if(pline(38:59).eq.'Forces (Hartrees/Bohr)')then
                read(202,*)
                read(202,*)
                write(201,'(A)')'Force'
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        read(202,'(23x,3F15.9)')(tmpforce(j,i),j=1,3)
                    end if
                end do
                !=======================Keep Energy conservation======================
                !do i=1,H_num
                !    read(202,'(23x,3f15.9)')(HForce(j),j=1,3)
                !    do k=1,natom
                !        if(gauss_tag(k).eq.1)then
                !            do j=1,3
                !                tmpforce(j,k)=tmpforce(j,k)+HForce(j)*mass(k)/Mass_sum
                !            end do
                !        end if
                !    end do
                !end do
                !=====================================================================
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        write(201,'(i10,3F15.9)')i,(tmpforce(j,i),j=1,3)
                    end if
                end do
                do i=1,nbona
                    if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).eq.0)then
                        read(202,'(23x,3f15.9)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,jb(i))=tmpforce(j,jb(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')jb(i),(tmpforce(j,jb(i)),j=1,3)
                    else if(gauss_tag(ib(i)).eq.0.and.gauss_tag(jb(i)).eq.1)then
                        read(202,'(23x,3f15.9)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ib(i))=tmpforce(j,ib(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')ib(i),(tmpforce(j,ib(i)),j=1,3)
                    end if
                end do
                do i=1,nbonh
                    if(gauss_tag(ibh(i)).eq.1.and.gauss_tag(jbh(i)).eq.0)then
                        read(202,'(23x,3f15.9)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ibh(i))=tmpforce(j,ibh(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')jbh(i),(tmpforce(j,jbh(i)),j=1,3)
                    else if(gauss_tag(ibh(i)).eq.0.and.gauss_tag(jbh(i)).eq.1)then
                        read(202,'(23x,3f15.9)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ibh(i))=tmpforce(j,ibh(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')ibh(i),(tmpforce(j,ibh(i)),j=1,3)
                    end if
                end do
                write(201,*)
            end if
            if(pline(42:55).eq.'Electric Field'.and.trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                read(202,*)
                read(202,*)
                write(201,'(A)') 'Field'  
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        read(202,*)
                    end if
                end do 
                do i=1,H_num
                    read(202,*)
                end do 
                if(trim(adjustl(if_field)).eq.'yes')then
                    do i=1,natom
                        if(gauss_tag(i).eq.0)then
                            read(202,'(24x,3F14.6)')(tmpfield(j,i),j=1,3)
                            write(201,'(I10,3F15.9)')i,(tmpfield(j,i),j=1,3)
                            do j=1,3
                                tmpforce(j,i)=tmpfield(j,i)*crg(i)
                            end do
                        end if
                    end do
                    if(add_bg.eq.'yes')then
                        write(201,'(A)')'bg_field'
                        do i=1,bg_num
                            read(202,'(24x,3f14.6)')(tmp_bgfield(j,i),j=1,3)
                            write(201,'(i10,3f15.9)')i,(tmp_bgfield(j,i),j=1,3)
                            do j=1,3
                                tmp_bgforce(j,i)=tmp_bgfield(j,i)*bg_crg(i)
                            end do
                        end do
                    end if
                end if
            end if
        end do
        close(202)
        close(201)
        !=============================================================
        !do i=1,natom
        !     if(gauss_tag(i).eq.0.and.coordinate_atom(i).eq.1)then
        !         do j=1,3
        !             tmpforce(j,i)=0
        !         end do
        !     end if
        !end do
        !=============================================================
        open(200,file=trim(adjustl(filename))//'.tmp')
            do i=1,natom
                write(200,'(3f14.8)')(tmpforce(j,i),j=1,3)
            end do
        close(200)
        tmp=0
        tmpcrg=0
        do i=1,natom
            do j=1,3
                tmp(j)=tmp(j)+tmpforce(j,i)
            end do
        end do
        if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
        if(add_bg.eq.'yes')then
            do i=1,bg_num
                do j=1,3
                    tmp(j)=tmp(j)+tmp_bgforce(j,i)
                end do
            end do
        end if
        end if
        if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
            do i=1,natom
                if(gauss_tag(i).eq.0)then
                    tmpcrg=abs(crg(i))+tmpcrg
                end if
            end do
            if(add_bg.eq.'yes')then
                do i=1,bg_num
                    tmpcrg=abs(bg_crg(i))+tmpcrg
                end do
            end if
            do i=1,natom
                if(gauss_tag(i).eq.0)then
                    do j=1,3
                        tmpforce(j,i)=tmpforce(j,i)-tmp(j)*abs(crg(i))/tmpcrg
                    end do
                end if
            end do
            if(add_bg.eq.'yes')then
            do i=1,bg_num
                do j=1,3
                    tmp_bgforce(j,i)=tmp_bgforce(j,i)-tmp(j)*abs(bg_crg(i))/tmpcrg
                end do
            end do
            end if
        end if
        !write(*,'(A20,3f6.3)')trim(adjustl(filename)),(tmp(i)*627.51/0.529,i=1,3)
    else
        write(*,*)'Error:'
        write(*,*)trim(adjustl(filename))//'.log end with error'
        write(101,*)'End with error'
        stop
    end if
else 
    write(*,*)'Error:'
    write(*,*)trim(adjustl(filename))//'.log does not exist'
    write(101,*)'End with error'
    stop
end if
return
end subroutine

subroutine check_file(file_name,file_flag)
implicit none
integer(kind=4)::ist
character(len=50)::file_name
character(len=80)::pline
logical::file_flag
file_flag=.false.
open (1001,file=trim(adjustl(file_name))//'.log')
do while(.true.)
    read(1001,'(A80)',iostat=ist)pline
    if(ist.ne.0)exit
    if(pline(2:18).eq.'Error termination')then
        write(101,'(A)')"File "//trim(adjustl(file_name))//&
                       'terminated with error'
        file_flag=.true.
    end if
end do
close(1001)
return
end subroutine 

