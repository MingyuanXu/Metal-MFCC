subroutine init_system
use comparm
use mpi
use amoeba
use dftb3
implicit none

integer(kind=8)::pts(20)
integer(kind=8),allocatable::respt(:)
character(len=60)::prmfile,crdfile
character(len=80)::pline
integer(kind=8)::i,j,k,ferr,iatom,jatom,katom,latom
integer(kind=4)::ierrs
real(kind=8)::dis
dis=0

prmfile=trim(adjustl(basename))//'.parm7'
crdfile=trim(adjustl(basename))//'.crd'

open(201,file=prmfile)
do i=1,6
    read(201,*)
end do

read(201,'(10I8)')pts
natom  = pts(1)
ntypes = pts(2)
nico   = ntypes*ntypes
nvdwp  = (ntypes +1)*ntypes/2
nbonh  = pts(3)
nbona  = pts(4)
ntheth = pts(5)
ntheta = pts(6)
nphih  = pts(7)
nphia  = pts(8)
nnb    = pts(11)
nres   = pts(12)
nphb   = pts(20)
numbnd = pts(16)

if(.not. allocated(iac))                    allocate(iac(natom))
if(.not. allocated(ico))                    allocate(ico(nico))
if(.not. allocated(vdwa))                   allocate(vdwa(nvdwp))
if(.not. allocated(vdwb))                   allocate(vdwb(nvdwp))
if(.not. allocated(sola))                   allocate(sola(nphb))
if(.not. allocated(solb))                   allocate(solb(nphb))
if(.not. allocated(elementnum))             allocate(elementnum(natom))
if(.not. allocated(element))                allocate(element(natom))
if(.not. allocated(rescrg))                 allocate(rescrg(nres))
if(.not. allocated(int_rescrg))             allocate(int_rescrg(nres))
if(.not. allocated(atomcrg))                allocate(atomcrg(natom))
if(.not. allocated(atomcrd))                allocate(atomcrd(3,natom))
if(.not. allocated(atomname))               allocate(atomname(natom))
if(.not. allocated(resname))                allocate(resname(natom))
if(.not. allocated(resnum))                 allocate(resnum(natom))
if(.not. allocated(residue))                allocate(residue(nres))
if(.not. allocated(respt))                  allocate(respt(nres))
if(.not. allocated(resstart))               allocate(resstart(nres))
if(.not. allocated(resend))                 allocate(resend(nres))
if(.not. allocated(selectC))                allocate(selectC(nres))
if(.not. allocated(selectCA))               allocate(selectCA(nres))
if(.not. allocated(selectCB))               allocate(selectCB(nres)) 
if(.not. allocated(selectN))                allocate(selectN(nres))
if(.not. allocated(selectCG))               allocate(selectCG(nres))
if(.not. allocated(ibh))                    allocate(ibh(nbonh))
if(.not. allocated(jbh))                    allocate(jbh(nbonh))
if(.not. allocated(icbh))                   allocate(icbh(nbonh))
if(.not. allocated(ib))                     allocate(ib(nbona))
if(.not. allocated(jb))                     allocate(jb(nbona))
if(.not. allocated(icb))                    allocate(icb(nbona))
if(.not. allocated(tmpforce))               allocate(tmpforce(3,natom))
if(.not. allocated(tmpcharge))              allocate(tmpcharge(natom))
if(.not. allocated(tmpfield))               allocate(tmpfield(3,natom))
if(.not. allocated(gauss_tag))              allocate(gauss_tag(natom))
if(.not. allocated(GMFCC_pairs))            allocate(GMFCC_pairs(natom,natom))
if(.not. allocated(iph))                    allocate(iph(nphih))
if(.not. allocated(jph))                    allocate(jph(nphih))
if(.not. allocated(kph))                    allocate(kph(nphih))
if(.not. allocated(lph))                    allocate(lph(nphih))
if(.not. allocated(icph))                   allocate(icph(nphih))
if(.not. allocated(ip))                     allocate(ip(nphia))
if(.not. allocated(jp))                     allocate(jp(nphia))
if(.not. allocated(kp))                     allocate(kp(nphia))
if(.not. allocated(lp))                     allocate(lp(nphia))
if(.not. allocated(icp))                    allocate(icp(nphia))
if(.not. allocated(numex))                  allocate(numex(natom))
if(.not. allocated(natex))                  allocate(natex(nnb))
if(.not. allocated(numex))                  allocate(numex(natom))
if(.not. allocated(connect_map))            allocate(connect_map(natom,natom))
if(.not. allocated(atomdis))                allocate(atomdis(natom,natom))
if(.not. allocated(mass))                   allocate(mass(natom))
if(.not. allocated(Bond_EQ))                allocate(Bond_EQ(numbnd))
if(.not. allocated(atomv))                  allocate(atomv(3,natom))
if(.not. allocated(restrain_F))             allocate(restrain_F(3,natom))
if(.not.allocated(QM_charge))               allocate(QM_charge(natom))
if(.not.allocated(D_pairs))                 allocate(D_pairs(natom,natom))
if(.not.allocated(tmpD_pairs))              allocate(tmpD_pairs(natom,natom))
D_pairs=0
GMFCC_pairs=0
atomv=0
restrain_F=0
QM_charge=0

do while(.true.) 
    read(201,'(A80)',iostat=ferr)pline
    if(ferr /= 0) exit
    if(pline(2:15) .eq. 'FLAG ATOM_NAME') then
        read(201,*)
        read(201,'(20A4)') atomname
    else if (pline(2:12) .eq.'FLAG CHARGE') then
        read(201,*)
        read(201,'(5E16.8)') atomcrg
        atomcrg = atomcrg / 18.2223
    else if (pline(2:10).eq.'FLAG MASS')then
        read(201,*)
        read(201,'(5E16.8)') mass
    else if (pline(2:19) .eq. 'FLAG ATOMIC_NUMBER')then
        read(201,*)
        read(201,'(10I8)') elementnum
    else if (pline(2:21) .eq. 'FLAG ATOM_TYPE_INDEX') then
        read(201,*)
        read(201,'(10I8)') iac
    else if (pline(2:26) .eq. 'FLAG NONBONDED_PARM_INDEX')then
        read(201,*)
        read(201,'(10I8)') ico
    else if(pline(2:17).eq.'FLAG HBOND_ACOEF')then
        read(201,*)
        read(201,'(5E16.8)') sola
    else if(pline(2:17).eq.'FLAG HBOND_BCOEF')then
        read(201,*)
        read(201,'(5E16.8)') solb
    else if (pline(2:25).eq. 'FLAG LENNARD_JONES_ACOEF')then
        read(201,*)
        read(201,'(5E16.8)')vdwa
    else if (pline(2:25).eq.'FLAG LENNARD_JONES_BCOEF') then
        read(201,*)
        read(201,'(5E16.8)')vdwb
    else if (pline(2:26) .eq. 'FLAG RESIDUE_POINTER') then
        read(201,*)
        read(201,'(10I8)')respt
        rescrg=0.0d0
        do i=1,nres-1
            do j=respt(i),respt(i+1)-1
                resname(j) = residue(i)
                resnum(j)=i
                rescrg(i)=rescrg(i)+atomcrg(j)
            end do
            resstart(i) =respt(i)
            resend(i)=respt(i+1)-1
        end do
        do j =respt(nres),natom
           resname(j)=residue(nres)
           resnum(j)=i
           rescrg(nres) = rescrg(nres) + atomcrg(j)
        end do
        resstart(nres) = respt(nres)
        resend(nres)   = natom
    else if (pline(2:26).eq. 'FLAG RESIDUE_LABEL') then
        read(201,*)
        read(201,'(20A4)')residue
    else if (pline(2:24).eq.'FLAG BONDS_INC_HYDROGEN') then
        read(201,*)
        read(201,'(10i8)')(ibh(i),jbh(i),icbh(i),i=1,nbonh)
        do i=1,nbonh
            ibh(i)=abs(ibh(i))/3+1
            jbh(i)=abs(jbh(i))/3+1
        end do
    else if (pline(2:28).eq.'FLAG BONDS_WITHOUT_HYDROGEN') then
        read(201,*)
        read(201,'(10i8)')(ib(i),jb(i),icb(i),i=1,nbona)
        do i=1,nbona
            ib(i)=abs(ib(i))/3+1
            jb(i)=abs(jb(i))/3+1
        end do
    else if (pline(2:22).eq.'FLAG BOND_EQUIL_VALUE')then
        read(201,*)
        read(201,'(5E16.8)')Bond_EQ
    else if (pline(2:28).eq.'FLAG DIHEDRALS_INC_HYDROGEN') then
        read (201,*)
        read (201,'(10i8)')(iph(i),jph(i),kph(i),lph(i),icph(i),i=1,nphih)
    else if (pline(2:32).eq.'FLAG DIHEDRALS_WITHOUT_HYDROGEN')then
        read (201,*)
        read (201,'(10i8)') (ip(i),jp(i),kp(i),lp(i),icp(i),i=1,nphia)
    else if (pline(2:27).eq.'FLAG NUMBER_EXCLUDED_ATOMS')then
        read(201,*)
        read(201,'(10I8)') numex
    else if (pline(2:25).eq.'FLAG EXCLUDED_ATOMS_LIST')then
        read (201,*)
        read (201,'(10I8)') natex
    end if
end do  
close (201)

if(.not.allocated(ex_atom_start))      allocate(ex_atom_start(natom))
if(.not.allocated(ex_atom_end))        allocate(ex_atom_end(natom))
ex_atom_start(1)=1
ex_atom_end(1)=numex(1)
do i=2,natom
    ex_atom_start(i)=ex_atom_end(i-1)+1
    ex_atom_end(i)  =ex_atom_end(i-1)+numex(i)
end do

connect_map=0
do i=1,natom
    connect_map(i,i)=-1
end do
do i=1,natom-1
    do j=i+1,natom
        do k=ex_atom_start(i),ex_atom_end(i)
            if(j.eq.natex(k))then
                connect_map(i,j)=1
                connect_map(j,i)=1
            end if
        end do
    end do
end do
do i=1,nphih
    iatom=abs(iph(i))/3+1
    jatom=abs(jph(i))/3+1
    katom=abs(kph(i))/3+1
    latom=abs(lph(i))/3+1
    if(kph(i).ge.0)then
        connect_map(iatom,latom)=2
        connect_map(latom,iatom)=2
    end if
end do
do i=1,nphia
    iatom=abs(ip(i))/3+1
    jatom=abs(jp(i))/3+1
    katom=abs(kp(i))/3+1
    latom=abs(lp(i))/3+1
    if(kp(i).ge.0)then
        connect_map(iatom,latom)=2
        connect_map(latom,iatom)=2
    end if
end do

!===========================================================================

open (201,file=crdfile)
read(201,*)
read(201,*)
read(201,'(6F12.7)')atomcrd
close(201)
if(add_bg.eq.'yes')then
    open(201,file='bg_crg.dat')
        read(201,'(i10)')bg_num
        if(.not.allocated(bg_crg))                  allocate(bg_crg(bg_num))
        if(.not.allocated(bg_crd))                  allocate(bg_crd(3,bg_num))
        if(.not.allocated(bg_force))                allocate(bg_force(3,bg_num))
        if(.not.allocated(tmp_bgfield))             allocate(tmp_bgfield(3,bg_num))
        if(.not.allocated(tmp_bgforce))             allocate(tmp_bgforce(3,bg_num))
        do i=1,bg_num
            read(201,'(4f14.8)')(bg_crd(j,i),j=1,3),bg_crg(i)
        end do
        bg_force=0
    close(201)
end if
do i=1,nbonh
    call distance(ibh(i),jbh(i),dis)
    if((dis.gt.Bond_EQ(icbh(i))*1.2)&
        .or.(dis.le.Bond_EQ(icbh(i))*0.8))then
        if(elementnum(ibh(i)).eq.1)then
    !        call restrain_H(jbh(i),ibh(i),Bond_EQ(icbh(i)))
        else if(elementnum(jbh(i)).eq.1)then
    !        call restrain_H(ibh(i),jbh(i),Bond_EQ(icbh(i)))
        end if
        write(*,'(A33,i5,A5,i5,A10)')'Warning!!!maybe the Bond between ',ibh(i),&
                                 ' and ',jbh(i),' is Wrong!'
        write(*,'(A)')'Please do the convergence test to check the force of this bond'
        write(*,'(A11,2x,f14.6)')'Restrain_F:',sqrt(restrain_F(1,ibh(i))**2+&
                    restrain_F(2,ibh(i))**2+restrain_F(3,ibh(i))**2)
        write(103,'(A33,i5,A5,i5,A10)')'Warning!!!maybe the Bond between ',ibh(i),&
                                 ' and ',jbh(i),' is Wrong!'
        write(103,'(A)')'Please do the convergence test to check the force of this bond'
        write(103,'(A11,2x,f14.6)')'Restrain_F:',sqrt(restrain_F(1,ibh(i))**2+&
                    restrain_F(2,ibh(i))**2+restrain_F(3,ibh(i))**2)
    end if
end do
open(201,file='MFCC.crd')
    write(201,'(A)')'default'
    write(201,'(i6)') natom
    write(201,'(6f12.7)')atomcrd
close(201)
do i=1,natom-1
    do j=i,natom
        call distance(i,j,dis)
        atomdis(i,j)=dis
        atomdis(j,i)=dis
    end do
end do
!===========================================================================

open(201,file='input')
special_resnum=0
chain_num=0
do while(.true.)
read(201,'(a)',IOSTAT=ierrs)pline
    IF(ierrs.NE.0)THEN
        EXIT
    END IF
    if(pline(1:10).eq.'chain_info')then
        read(201,'(i10)')chain_num
        if(.not.allocated(chain_start)) allocate(chain_start(chain_num))
        if(.not.allocated(chain_end))   allocate(chain_end(chain_num))
        if(.not.allocated(chain_tag))   allocate(chain_tag(nres))
        chain_tag=0
        write(*,'(A)')'The information of protein:'
        do i=1,chain_num
            read(201,'(2i5)')chain_start(i),chain_end(i)
            write(*,*)'the',i,'th chain: ',chain_start(i),chain_end(i)
            do j=chain_start(i),chain_end(i)
                chain_tag(j)=i
            end do
        end do
    else if(pline(1:11).eq.'special_res')then
        read(201,'(i10)')special_resnum
        if(.not.allocated(special_res)) allocate(special_res(special_resnum))
        read(201,'(5i5)')special_res
    end if
end do
close(201)

!===========================================================================
do i=1,natom
    if (atomname(i).eq.'ZN')then
        atomname(i)='Zn'
    else if (atomname(i).eq.'CU')then
        atomname(i)='Cu'
    else if(atomname(i).eq.'FE')then
        atomname(i)='Fe'
    end if
    if(elementnum(i).eq.1)then
        element(i)='H '
    else if (elementnum(i).eq.7)then
        element(i)='N '
    else if (elementnum(i).eq.6)then
        element(i)='C '
    else if (elementnum(i).eq.8)then
        element(i)='O '
    else if (elementnum(i).eq.16)then
        element(i)='S '
    else if (elementnum(i).eq.11)then
        element(i)='Na'
    else if (elementnum(i).eq.17)then
        element(i)='Cl'
    else if (elementnum(i).eq.30)then
        element(i)='Zn'
    else if (elementnum(i).eq.20)then
        element(i)='Ca'
    else if (elementnum(i).eq.12)then
        element(i)='Mg'
    else if (elementnum(i).eq.29)then
        element(i)='Cu'
    else if(elementnum(i).eq.26)then
        element(i)='Fe'
    else if (elementnum(i).eq. -1)then
        element(i)=trim(adjustl(atomname(i)))
    end if 
end do
do i=1,natom
     if (atomname(i) .eq. 'C   ')then
        selectC(resnum(i))=i
     else if (atomname(i).eq.'N   ')then
        selectN(resnum(i))=i
     else if (atomname(i).eq.'CA  ')then
        selectCA(resnum(i)) =i
     else if(atomname(i).eq.'CB  ') then
        selectCB(resnum(i))=i
     else if(atomname(i).eq.'CG')then
        selectCG(resnum(i))=i
     end if
end do
if(trim(adjustl(software)).eq.'dftb+')then
    if(trim(adjustl(software)).eq.'dftb+')then
        read(*,nml=dftb3_ctrl)
        call set_dftb_parameter
    end if
end if
do i=1,nres
    int_rescrg(i)=NINT(rescrg(i))
!    write(*,*)residue(i),int_rescrg(i)
end do
!=====================================================================
!interface to correct background_charge
call adjust_charge
write(*,'(2A10)')trim(adjustl(Frag_method)),trim(adjustl(charge_type))
if(trim(adjustl(charge_type)).eq.'EPB')then
!    open(201,file='EPB_charge')
!        do i=1,natom
!            read(201,'(f16.8)')crg(i)
!        end do
!    close(201)
call cal_EPB_charge
else if(trim(adjustl(charge_type)).eq.'Amber')then
write(*,*)'here'
!if Amber,do nothing     
else if(trim(adjustl(charge_type)).eq.'Amoeba')then
    call Translate_atomindex_To_Amoeba
end if
!======================================================================
end subroutine
subroutine restrain_H(i,j,length)
use comparm
use mpi
implicit none
integer(kind=8)::i,j,m
real(kind=8)::vector(3),length,dis
call distance(i,j,dis)
if(dis.gt.length*1.2)then
    do m=1,3
        vector(m)=(atomcrd(m,j)-atomcrd(m,i))/dis
        atomcrd(m,j)=atomcrd(m,i)+vector(m)*length*1.2
!        restrain_F(m,j)=-100*(dis-length)*vector(m)
!        restrain_F(m,i)=100*(dis-length)*vector(m)
    end do
else if(dis.lt.length*0.8)then
    do m=1,3
        vector(m)=(atomcrd(m,j)-atomcrd(m,i))/dis
        atomcrd(m,j)=atomcrd(m,i)+vector(m)*length*0.8
!        restrain_F(m,j)=100*(dis-length)*vector(m)
!        restrain_F(m,i)=-100*(dis-length)*vector(m)
    end do
end if

end subroutine
subroutine cal_unitdis(i,j,dis)
    use comparm
    use mpi
    implicit none
    integer(kind=8)::i,j,m,n
    real(kind=8)::dis
    dis=999.0
    do m=unit_pt(1,i),unit_pt(2,i)
        do n=unit_pt(1,j),unit_pt(2,j)
            if(atomdis(m,n).le.dis)then
                dis=atomdis(m,n)
            end if
        end do
    end do
end subroutine
subroutine cal_EPB_charge
use comparm
use mpi
implicit none
integer(kind=8)::             Polar_bond_num,subtract_num
integer(kind=8)::             i,j
real(kind=8)::                tmp1,tmp2,dis,dELE_potential,tmp3,dq
integer(kind=8),allocatable:: Polar_bond(:,:),subtract_interaction(:,:),&
                              Polar_bondfactor(:)
real(kind=8),allocatable::    EPB_charge(:),ELE_potential(:),subtract_factor(:),&
                              Polar_factor(:),charge_back(:),equil_bond(:),gas_crg(:)
real(kind=8),parameter::      scale_f=0.83333333333d0
logical::                     iteration_flag,fexist
integer(kind=8)::             iteration_num
if(.not.allocated(EPB_charge))             allocate(EPB_charge(natom))
if(.not.allocated(charge_back))            allocate(charge_back(natom))
if(.not.allocated(gas_crg))                allocate(gas_crg(natom))
if(.not.allocated(ELE_potential))          allocate(ELE_potential(natom))
EPB_charge=crg
charge_back=crg
gas_crg=crg
ELE_potential=0
inquire(file='jcgparam.inp',exist=fexist)
if(fexist)then
open(501,file='jcgparam.inp') 
    read(501,*)Polar_bond_num
    if(.not.allocated(Polar_bond))         allocate(Polar_bond(Polar_bond_num,2))
    if(.not.allocated(equil_bond))         allocate(equil_bond(Polar_bond_num))
    if(.not.allocated(Polar_bondfactor))   allocate(Polar_bondfactor(Polar_bond_num))
    if(.not.allocated(Polar_factor))       allocate(Polar_factor(Polar_bond_num))
    do i=1,Polar_bond_num
        read(501,'(3i10,4f10.5)')Polar_bond(i,1),Polar_bond(i,2),&
                                 Polar_bondfactor(i),equil_bond(i),&
                                 Polar_factor(i),tmp1,tmp2
        gas_crg(Polar_bond(i,1))=tmp1
        gas_crg(Polar_bond(i,2))=tmp2
    end do
close(501)
EPB_charge=gas_crg
else 
    write(*,*)'Error:jcgparam.inp does not exist!'
end if
inquire(file='EPB_CHARGE',exist=fexist)
if(fexist)then
    open(501,file='EPB_CHARGE')
        read(501,'(5F15.8)')EPB_charge
    close(501)
end if    
inquire(file='sceeparm.txt',exist=fexist)
if(fexist)then
    open(501,file='sceeparm.txt')
        read(501,*)subtract_num
        if(.not.allocated(subtract_interaction))&
            allocate(subtract_interaction(subtract_num,2))
        if(.not.allocated(subtract_factor))&
            allocate(subtract_factor(subtract_num))
        do i=1,subtract_num
            read(501,*)subtract_interaction(i,1),subtract_interaction(i,2),&
                        subtract_factor(i)
        end do
    close(501)
else 
    write(*,*)'Error:sceeparm.txt does not exist!'
end if
do i=1,subtract_num
    connect_map(subtract_interaction(i,1),subtract_interaction(i,2))=3
    connect_map(subtract_interaction(i,2),subtract_interaction(i,1))=3
end do
iteration_num=0
open(501,file='Iteration_EPB_charge')
write(501,'(A15,i10)')'Iteration_num:',iteration_num
write(501,'(5f16.8)')EPB_charge
do while(.true.)
    iteration_num=iteration_num+1
    charge_back=EPB_charge 
    ELE_potential=0
    do i=1,natom
        do j=1,natom
            if(connect_map(i,j).eq.0)then
                dis=atomdis(i,j)
                ELE_potential(i)=ELE_potential(i)+charge_back(j)/dis*18.2223
            else if(connect_map(i,j).eq.2)then
                dis=atomdis(i,j)
                ELE_potential(i)=ELE_potential(i)+charge_back(j)/dis*scale_f*18.2223
            end if    
        end do
    end do
    open(502,file='ELE_Potential')
        !read(502,*)
        write(502,'(f16.8)') ELE_potential
    close(502)
!    open(502,file='Check')
    do i=1,Polar_bond_num
        dELE_potential=ELE_potential(Polar_bond(i,1))-ELE_potential(Polar_bond(i,2))
!       write(502,'(2i10,3f16.8)')Polar_bond(i,1),Polar_bond(i,2),ELE_potential(Polar_bond(i,1)),&
!            ELE_potential(Polar_bond(i,2)),dELE_potential
        dq=dELE_potential/(2*Polar_factor(i)*(4.8*equil_bond(i))**2)
        dq=dq/2*Polar_bondfactor(i)*18.2223
!       write(502,'(f16.8)')dq
        EPB_charge(Polar_bond(i,1))=gas_crg(Polar_bond(i,1))+dq
        EPB_charge(Polar_bond(i,2))=gas_crg(Polar_bond(i,2))-dq
    end do
!    close(502)
    iteration_flag=.true.
    do i=1,natom
        tmp3=abs(charge_back(i)-EPB_charge(i))
        if(tmp3.gt.0.001)then
            iteration_flag=.false.
        end if
    end do
    write(501,'(A15,i10)')'Iteration_num:',iteration_num
    write(501,'(5f16.8)')EPB_charge
    if(iteration_flag.eqv..true.)then
        exit
    end if
end do
crg=EPB_charge
close(501)
open(501,file='EPB_CHARGE')
    write(501,'(5F15.8)')EPB_charge
close(501)
end subroutine
