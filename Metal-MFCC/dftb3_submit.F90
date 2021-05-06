subroutine set_dftb_parameter
use comparm
use dftb3
use mpi
implicit none
integer(kind=8)::i
if(.not.allocated(MaxAngularMomentum)) allocate(MaxAngularMomentum(natom))
if(.not.allocated(HubbardDerivs))      allocate(HubbardDerivs(natom))
do i=1,natom
    if(elementnum(i).eq.1)then
        MaxAngularMomentum(i)='s'
        HubbardDerivs(i)=-0.1857
    else if(elementnum(i).eq.6)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.1492
    else if(elementnum(i).eq.8)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.1575
    else if(elementnum(i).eq.7)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.1535
    else if(elementnum(i).eq.9)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.1623
    else if(elementnum(i).eq.11)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.0454
    else if(elementnum(i).eq.12)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.02
    else if(elementnum(i).eq.15)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.14
    else if(elementnum(i).eq.16)then
        MaxAngularMomentum(i)='d'
        HubbardDerivs(i)=-0.11
    else if(elementnum(i).eq.17)then
        MaxAngularMomentum(i)='d'
        HubbardDerivs(i)=-0.0697
    else if(elementnum(i).eq.19)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.0339
    else if(elementnum(i).eq.20)then
        MaxAngularMomentum(i)='p'
        HubbardDerivs(i)=-0.0340
    else if(elementnum(i).eq.35)then
        MaxAngularMomentum(i)='d'
        HubbardDerivs(i)=-0.0573
    else if(elementnum(i).eq.30)then
        MaxAngularMomentum(i)='d'
        HubbardDerivs(i)=-0.03
    else if(elementnum(i).eq.53)then
        MaxAngularMomentum(i)='d'
        HubbardDerivs(i)=-0.0433
    end if
end do
end subroutine

subroutine dftb_submit(filename)
use comparm
use mpi
use dftb3
implicit none
integer(kind=8)::i,j,k,charge,H_num,electron_num,QM_atom,MM_atom,H_index
logical::exist_flag
character(len=50)::filename,cn,clist
character(len=256)::cmd
logical::judge_flag
if(.not.allocated(element_tag)) allocate(element_tag(natom))
element_Tag=0
element_list=''
element_MaxAngM=''
element_Hubbard=0
element_typenum=0
QM_atom=0
MM_atom=0
electron_num=0
exist_flag=.false.
charge=0
gauss_unit=0
!spin=1
call system('mkdir '//trim(adjustl(filename)))
open(202,file=trim(adjustl(filename))//'.coord')
open(203,file=trim(adjustl(filename))//'.field')
open(204,file=trim(adjustl(filename))//'.hsd')
open(206,file=trim(adjustl(filename))//'.pdb')

write(202,'(A6,2X,A5)')'XNATOM','C'
write(202,'(A5)')'XLIST'

do i=1,natom
    if(gauss_tag(i).eq.1)then
        judge_flag=.false.
        QM_atom=QM_atom+1
        do j=1,element_typenum 
            if(element(i).eq.element_list(j))then
                judge_flag=.true.
                element_tag(i)=j
                write(202,'(i8,i8,5x,3f14.8)')QM_atom,element_tag(i),(atomcrd(k,i),k=1,3)
                write(206,299)'ATOM  ',i, atomname(i), resname(i)(1:3), &
                               resnum(i), (atomcrd(k,i),k=1,3)
                electron_num=electron_num+elementnum(i)
            end if
        end do
        if(judge_flag.eqv..false.)then
            element_typenum=element_typenum+1
            if(trim(adjustl(element(i))).eq.'H')then
                H_index=element_typenum
            end if
            element_list(element_typenum)=element(i)
            element_tag(i)=element_typenum
            write(202,'(i8,i8,5x,3f14.8)')QM_atom,element_tag(i),(atomcrd(k,i),k=1,3)
            write(206,299)'ATOM  ',i, atomname(i), resname(i)(1:3), &
                           resnum(i), (atomcrd(k,i),k=1,3)
            electron_num=electron_num+elementnum(i)
            element_MaxAngM(element_typenum)=MaxAngularMomentum(i)
            element_Hubbard(element_typenum)=HubbardDerivs(i)
        end if
        if(trim(adjustl(filename)).ne.'ALL_sys')then
            do j=i,natom
                if(gauss_tag(j).eq.1)then
                    GMFCC_pairs(i,j)=1
                    GMFCC_pairs(j,i)=1
                end if
            end do
        end if
        if(unit_index(i).ne.0)then
            gauss_unit(unit_index(i))=1
        end if
    else if(gauss_tag(i).eq.0)then
        if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
            MM_atom=MM_atom+1
            write(203,'(3f14.8,5x,f14.8)')(atomcrd(j,i),j=1,3),crg(i)
        end if
    end if
end do
if(trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
    if(add_bg.eq.'yes')then
        do i=1,bg_num
            MM_atom=MM_atom+1
            write(203,'(3f14.8,5x,f14.8)')(bg_crd(j,i),j=1,3),bg_crg(i)
        end do
    end if
end if
close(203)

do i=1,unit_num
    if(gauss_unit(i).eq.1)then
        charge=charge+unit_charge(i)
    end if
end do

H_num=0

do i=1,nbona
    if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).ne.1)then
        H_num=H_num+1
        call dftb_add_H(ib(i),jb(i),QM_atom+H_num,H_index)
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,jb(i))=1
                GMFCC_pairs(jb(i),j)=1
            end if
        end do
    else if(gauss_tag(ib(i)).ne.1.and.gauss_tag(jb(i)).eq.1)then
        H_num=H_num+1
        call dftb_add_H(jb(i),ib(i),QM_atom+H_num,H_index)
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
        H_num=H_num+1
        call dftb_add_H(ibh(i),jbh(i),QM_atom+H_num,H_index)
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,jbh(i))=1
                GMFCC_pairs(jbh(i),j)=1
            end if
        end do
    else if(gauss_tag(ibh(i)).ne.1.and.gauss_tag(jbh(i)).eq.1)then
        H_num=H_num+1
        call dftb_add_H(jbh(i),ibh(i),QM_atom+H_num,H_index)
        electron_num=electron_num+1
        do j=1,natom
            if(gauss_tag(j).eq.1)then
                GMFCC_pairs(j,ibh(i))=1
                GMFCC_pairs(ibh(i),j)=1
            end if
        end do
    end if
end do
close (202)
close(206)
write(cn,'(i5)')QM_atom+H_num
open(206,file='info')
write(206,'(i10)')H_num
close(206)
cmd="sed -i 's/XNATOM/"//trim(adjustl(cn))//"/g' "//trim(adjustl(filename))//'.coord'
call system(trim(adjustl(cmd)))
write(clist,'(20A2)')element_list
cmd="sed -i 's/XLIST/"//trim(adjustl(clist))//"/g' "//trim(adjustl(filename))//'.coord'
call system(trim(adjustl(cmd)))
write(204,'(A)')'Geometry=GenFormat{'
write(204,'(A)')'<<<"'//trim(adjustl(filename))//'.coord"'
write(204,'(A)')'}'
write(204,'(A)')''
write(204,'(A)')'Hamiltonian=DFTB{'
write(204,'(A)')'    SCC=Yes'
write(204,'(A)')'    SCCTolerance=1.0e-6'
write(204,'(A)')'    MaxSCCIterations=10000'
write(204,'(A)')'    Mixer=Broyden{'
write(204,'(A)')'        MixingParameter=0.1'
write(204,'(A)')'#        CacheIterations=-1'
write(204,'(A)')'        InverseJacobiweight=0.01'
write(204,'(A)')'        MinimalWeight=1'
write(204,'(A)')'        MaximalWeight=100000.'
write(204,'(A)')'        WeightFactor=0.01'
write(204,'(A)')'    }'
write(204,'(A)')''
write(204,'(A)')'    SlaterKosterFiles=Type2FileNames{'
write(204,'(A)')"        Prefix='"//trim(adjustl(workpath))//trim(adjustl(dftb_prmpath))//"'"
write(204,'(A)')'        Separator="-"'
write(204,'(A)')'        Suffix=".skf"'
write(204,'(A)')'        LowerCaseTypeName=No'
write(204,'(A)')'    }'
write(204,'(A)')''
write(204,'(A)')'    MaxAngularMomentum={'
do i=1,element_typenum
    write(204,'(A)')'        '//element_list(i)//'="'//element_MaxAngM(i)//'"'
end do
write(204,'(A)')'    }'
write(204,'(A12,i5)')'    charge=',charge
write(204,'(A)')'    SpinPolarisation={}'
write(204,'(A)')'    Filling=Fermi{'
write(204,'(A)')'        Temperature[k]=300'
write(204,'(A)')'    }'
if(Frag_method.eq.'EE-GMFCC')then
    if((trim(adjustl(filename)).ne.'Full_QM' .and. trim(adjustl(filename)).ne.'ALL_sys')&
    .or.(add_bg.eq.'yes'))then
        write(204,'(A)')'    ElectricField={'
        write(204,'(A)')'        PointCharges={'
        write(204,'(A)')'            CoordsAndCharges [Angstrom]=DirectRead{'
        write(204,'(A25,i10)')'                Records=',MM_atom
        write(204,'(A)')'                File="'//trim(adjustl(filename))//'.field"'
        write(204,'(A)')'            }'
        write(204,'(A)')'        }'
        write(204,'(A)')'    }'
    end if
end if
write(204,'(A)')'    OrbitalResolvedSCC=No'
write(204,'(A)')'    ReadInitialCharges=No'
write(204,'(A)')'    Eigensolver=DivideAndConquer{}'
write(204,'(A)')'    OldSKInterpolation=No'
write(204,'(A)')'    ThirdOrderFull=Yes'
write(204,'(A)')'    DampXH=Yes'
write(204,'(A)')'    DampXHExponent=4.00'
write(204,'(A)')'    HubbardDerivs={'
do i=1,element_typenum
    write(204,'(8x,A3,f8.4)')trim(adjustl(element_list(i)))//'=',element_Hubbard(i)
end do
write(204,'(A)')'    }'
write(204,'(A)')'    Dispersion = DftD3{'
write(204,'(A)')'        Damping = BeckeJohnson{'
write(204,'(A)')'            a1 = 0.746'
write(204,'(A)')'            a2 = 4.191'
write(204,'(A)')'        }'
write(204,'(A)')'        s8 = 3.209'
write(204,'(A)')'    }'
write(204,'(A)')'}'
write(204,'(A)')''
write(204,'(A)')'Options={'
write(204,'(A)')'    WriteAutotestTag=Yes'
write(204,'(A)')'    WriteDetailedXML=Yes'
write(204,'(A)')'    writeResultsTag=Yes'
write(204,'(A)')'    RandomSeed=0'
write(204,'(A)')'}'
write(204,'(A)')''
write(204,'(A)')'Analysis={'
write(204,'(A)')'    CalculateForces=Yes'
write(204,'(A)')'    WriteEigenvectors=No'
write(204,'(A)')'    AtomResolvedEnergies=No'
write(204,'(A)')'    writeBandOut=Yes'
!if(Frag_method.eq.'EE-GMFCC')then
!    if((trim(adjustl(filename)).ne.'Full_QM' .and. trim(adjustl(filename)).ne.'ALL_sys')&
!    .or.(add_bg.eq.'yes'))then
!        write(204,'(A)')'    ElectrostaticPotential={'
!        write(204,'(A)')'       OutputFile="'//trim(adjustl(filename))//'.Efield"'
!        write(204,'(A)')'       Softening=1E-6'
!        write(204,'(A)')'       Points={'
!        do i=1,natom
!            if(gauss_tag(i).eq.0)then
!                write(204,'(3f14.8)')(atomcrd(j,i),j=1,3)
!            end if
!        end do
!        if(add_bg.eq.'yes')then
!            do i=1,bg_num
!                write(203,'(3f14.8)')(bg_crd(j,i),j=1,3)
!            end do
!        end if
!        write(204,'(A)')'       }'
!        write(204,'(A)')'    }'
!    end if
!end if
write(204,'(A)')'}'
write(204,'(A)')''
write(204,'(A)')'ParserOptions={'
write(204,'(A)')'    ParserVersion=5'
write(204,'(A)')'}'
write(204,'(A)')'  ' 
close(204)
call system('mv info '//trim(adjustl(filename))//'.* '//trim(adjustl(filename)))
call system('cd '//trim(adjustl(workpath))//trim(adjustl(filename))//&
            '&& cp '//trim(adjustl(filename))//'.hsd dftb_in.hsd')
if(trim(adjustl(filename)).ne.'ALL_sys')then
    write(102,'(A100,i20)')'cd '//trim(adjustl(workpath))//trim(adjustl(filename))//'&& dftb+ >'&
    //trim(adjustl(filename))//'.out',electron_num**3
end if
299 format(a6,1x,I4,1x,a4,1x,a3,2x,I4,4x,3f8.3)
deallocate (element_tag)
end subroutine
!=============================================================================
subroutine dftb_add_H(i,j,A_index,H_index)
use comparm
use mpi
use dftb3
implicit none
integer(kind=8)::i,j,A_index,H_index
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

write(202,'(i10,i10,5x,3F14.8)')A_index,H_index,x,y,z

end subroutine
!===============================================================================
subroutine dftb3_read_data(filename)
use comparm
use mpi
implicit none

integer(kind=8)::i,j,H_num
character(len=50)::filename
logical::   fexist,file_tag
character(len=80)::pline
real(kind=8)::HForce(3),Mass_sum,tmp(3),tmpcrg
integer(kind=4)::ist
real(kind=8)::Edp,EH,E_scc,Erp,TE,T_ele_E
tmpenergy=0
tmpforce=0
tmpcharge=0
if(add_bg.eq.'yes')then
    tmp_bgfield=0
    tmp_bgforce=0
end if
tmpD_pairs=0
inquire(file=trim(adjustl(filename))//'/detailed.out',exist=fexist)

if(fexist)then
    call dftb3_check_file(filename,file_tag)
    if(file_tag.eqv..false.)then
        H_num=0
        open(202,file=trim(adjustl(filename))//'/detailed.out')
        open(201,file=trim(adjustl(filename))//'/calinfo.data')
        open(203,file=trim(adjustl(filename))//'/info')
            read(203,'(i10)')H_num
        close(203)

        Mass_sum=0
        do i=1,natom
            if(gauss_tag(i).eq.1)then
                Mass_sum=Mass_sum+mass(i)
            end if
        end do
        do i=1,natom
            do j=i,natom
                if(gauss_tag(i).ne.0.or.gauss_tag(j).ne.0)then
                    tmpD_pairs(i,j)=tmpD_pairs(i,j)+1
                    tmpD_pairs(j,i)=tmpD_pairs(j,i)+1
                end if
            end do
        end do
        do while(.true.)
            read(202,'(A80)',iostat=ist)pline
            if(ist.ne.0)then
                exit
            end if
            if(pline(1:14).eq.'Total energy:')then
                read(pline,'(30x,F20.10)')TE
            end if
            if(pline(1:9).eq.'Energy H0')then
                read(pline,'(30x,F20.10)')EH
            end if
            if(pline(1:10).eq.'Energy SCC')then
                read(pline,'(30x,F20.10)')E_scc
            end if
            if(pline(1:23).eq.'Total Electronic energy')then
                read(pline,'(30x,f20.10)')T_ele_E
            end if
            if(pline(1:16).eq.'Repulsive energy')then
                read(pline,'(30x,f20.10)')Erp
                tmpenergy=T_ele_E+Erp 
                write(201,'(A12,4X,F20.10)')'QM_energy=  ',tmpenergy
            end if
            if(pline(1:17).eq.'Dispersion energy')then
                read(pline,'(30x,f20.10)')Edp
            end if
            tmpcenergy=0
            if(pline(1:21).eq.' Atomic gross charges')then
                read(202,*)
                write(201,*)'Atomic gross charge:'
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        read(202,'(6x,f16.8)')tmpcharge(i)
                        write(201,'(i10,F12.6)')i,tmpcharge(i)
                    end if
                end do
            end if
            
            if(pline(1:12).eq.'Total Forces')then
                write(201,'(A)')'Force'
                do i=1,natom
                    if(gauss_tag(i).eq.1)then
                        read(202,'(3F20.12)')(tmpforce(j,i),j=1,3)
                        write(201,'(i10,3F15.9)')i,(tmpforce(j,i),j=1,3)
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
                do i=1,nbona
                    if(gauss_tag(ib(i)).eq.1.and.gauss_tag(jb(i)).eq.0)then
                        read(202,'(3f20.12)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,jb(i))=tmpforce(j,jb(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')jb(i),(tmpforce(j,jb(i)),j=1,3)
                    else if(gauss_tag(ib(i)).eq.0.and.gauss_tag(jb(i)).eq.1)then
                        read(202,'(3f20.12)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ib(i))=tmpforce(j,ib(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')ib(i),(tmpforce(j,ib(i)),j=1,3)
                    end if
                end do
                do i=1,nbonh
                    if(gauss_tag(ibh(i)).eq.1.and.gauss_tag(jbh(i)).eq.0)then
                        read(202,'(3f20.12)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ibh(i))=tmpforce(j,ibh(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')jbh(i),(tmpforce(j,jbh(i)),j=1,3)
                    else if(gauss_tag(ibh(i)).eq.0.and.gauss_tag(jbh(i)).eq.1)then
                        read(202,'(3f20.12)')(HForce(j),j=1,3)
                        do j=1,3
                            tmpforce(j,ibh(i))=tmpforce(j,ibh(i))+HForce(j)
                        end do
                        write(201,'(i10,3F15.9)')ibh(i),(tmpforce(j,ibh(i)),j=1,3)
                    end if
                end do
                write(201,*)
            end if
            if(pline(1:26).eq.'Forces on external charges'.and.&
                trim(adjustl(Frag_method)).eq.'EE-GMFCC')then
                write(201,'(A)') 'Force on external charges'  
                do i=1,natom
                    if(gauss_tag(i).eq.0)then
                        read(202,'(3F20.12)')(tmp(j),j=1,3)
                        do j=1,3
                            tmpforce(j,i)=tmpforce(j,i)+tmp(j)
                        end do
                        write(201,'(I10,3F15.9)')i,(tmp(j),j=1,3)
                    end if
                end do
                if(add_bg.eq.'yes')then
                    write(201,'(A)')'bg_force'
                    do i=1,bg_num
                        read(202,'(3f20.12)')(tmp_bgforce(j,i),j=1,3)
                        write(201,'(i10,3f15.9)')i,(tmp_bgforce(j,i),j=1,3)
                    end do
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
        open(200,file=trim(adjustl(filename))//'/calinfo.forcetmp')
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
        write(*,*)filename,tmp*627.51/0.529
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
        write(*,*)filename,tmp*627.51/0.529
        !write(*,'(A20,3f6.3)')trim(adjustl(filename)),(tmp(i)*627.51/0.529,i=1,3)
    else
        write(*,*)'Error:'
        write(*,*)trim(adjustl(filename))//'/detailed.out end with error'
        write(101,*)'End with error'
        stop
    end if
else 
    write(*,*)'Error:'
    write(*,*)trim(adjustl(filename))//'/detailed.out does not exist'
    write(101,*)'Program Stop!!!'
    stop
end if
return
end subroutine

subroutine dftb3_check_file(file_name,file_flag)
implicit none
integer(kind=4)::ist
character(len=50)::file_name
character(len=80)::pline
logical::file_flag
file_flag=.false.
open (1001,file=trim(adjustl(file_name))//'/detailed.out')
do while(.true.)
    read(1001,'(A80)',iostat=ist)pline
    if(ist.ne.0)exit
    if(pline(1:24).eq.' -> SCC is NOT converged')then
        write(101,'(A)')"Job "//trim(adjustl(file_name))//&
                       'terminated with error'
        file_flag=.true.
    end if
end do
close(1001)
return
end subroutine 


