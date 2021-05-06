module comparm
use mpi
integer(kind=8)::               natom,nres,ntypes,nico,nvdwp,nbona,&
                                nbonh,nphb,ntheth,ntheta,nnb,nphih,nphia,numbnd
integer(kind=8),allocatable::   iac(:),ico(:),elementnum(:),numex(:),natex(:)
real(kind=8),allocatable::      atomcrg(:),rescrg(:),crg(:)
integer(kind=8),allocatable::   int_rescrg(:)
real(kind=8),allocatable::      atomcrd(:,:),atomv(:,:),restrain_F(:,:),bgcrd(:,:),mass(:)
real(kind=8),allocatable::      vdwa(:),vdwb(:),sola(:),solb(:)
character(len=4),allocatable::  atomname(:),resname(:),residue(:)
integer(kind=8),allocatable::   resstart(:),resend(:),resnum(:)
character(len=2),allocatable::  element(:)
integer(kind=8),allocatable::   selectC(:),selectN(:),selectCA(:),selectCB(:),selectCG(:)
integer(kind=8),allocatable::   ex_atom_start(:),ex_atom_end(:)
integer(kind=8),allocatable::   ibh(:),jbh(:),icbh(:),ib(:),jb(:),icb(:)
integer(kind=8),allocatable::   ith(:),jth(:),kth(:),icth(:)
integer(kind=8),allocatable::   it(:),jt(:),kt(:),ict(:)
integer(kind=8),allocatable::   iph(:),jph(:),kph(:),lph(:),icph(:)
integer(kind=8),allocatable::   ip(:),jp(:),kp(:),lp(:),icp(:)
real(kind=8),allocatable::      Bond_EQ(:)
integer(kind=8),allocatable::   connect_map(:,:)
real(kind=8),allocatable::      atomdis(:,:)
real(kind=8),parameter::        rfactor=0.529177210920d0
real(kind=8),parameter::        ae=627.509d0
!===============================================
integer(kind=8)::               unit_num,HB_num
integer(kind=8),allocatable::   unit_pt(:,:),unit_index(:),unit_charge(:)
integer(kind=8),parameter::     combine_max=10
integer(kind=8)::               combine_num
integer(kind=8)::               combine_unit(2,combine_max)
integer(kind=8)::               frag_num
integer(kind=8),allocatable::   fragment(:),frag_type(:)
character(len=1),allocatable::  restype(:),unit_type(:)
!===============================================
integer(kind=8)::chain_num
integer(kind=8),allocatable::   chain_start(:),chain_end(:)
integer(kind=8),allocatable::   chain_tag(:)
integer(kind=8),allocatable::   GMFCC_pairs(:,:)
integer(kind=8)::special_resnum
integer(kind=8),allocatable::   special_res(:)
!==============================================
character(len=30)::             qmproc,qmmem,qmmethod,qmbasis
real(kind=8)::                  b2_cutoff
integer(kind=8)::               center_atom,stable_mode,strategy,bad_point
namelist /qm_ctrl/  qmproc,qmmem,qmmethod,qmbasis,center_atom,stable_mode
character(len=30):: basename,ifrun,ifread,if_field,Frag_method,if_b3,if_cal_charge
character(len=60)::             add_word
character(len=3)::              add_bg,if_fullQM
character(len=6)::              check_conver
character(len=4)::              frag_mode
character(len=6)::              charge_type
character(len=10)::             software
namelist /job_ctrl/ basename,ifrun,ifread,if_field,add_bg,add_word,&
                    strategy,bad_point,Frag_method,b2_cutoff,frag_mode,if_fullQM,charge_type,&
                    check_conver,software,if_cal_charge
!==============================================
integer(kind=8),allocatable::   gauss_tag(:),gauss_unit(:),B1_tag(:),B1_utag(:)
integer(kind=8),allocatable::   QM_tag(:)
integer(kind=8)::               B2_num,QM_2b_num,QM_ligand,B2_unum,B2_uf_num
real(kind=8),allocatable::      B2_E(:),B2_F(:,:,:),B1_E(:),B1_F(:,:,:),&
                                Body2_Force(:,:),B2_uE(:),B2_uF(:,:,:),B1_uF(:,:,:),B1_uE(:),&
                                B2_uf_E(:),B2_uf_F(:,:,:)
real(kind=8)::                  Body2_Energy,QM_correct_E,QM_ALL_E,Body3_Energy
real(kind=8),allocatable::      QM_1B_E(:),QM_2B_E(:),&
                                QM_charge(:),QM_1B_charge(:),QM_2B_charge(:)
real(kind=8),allocatable::      QM_1B_F(:,:,:),QM_2B_F(:,:,:),QM_correct_F(:,:)
real(kind=8),allocatable::      QM_ALL_F(:,:)
real(kind=8),allocatable::      HB_energy(:),HB_Force(:,:,:),HB_charge(:,:),&
                                HB1_Force(:,:,:),HB1_energy(:),HB1_charge(:,:)
integer(kind=8),allocatable::   HB_connect(:,:)   
logical,allocatable::           HB_tag(:)
integer(kind=8),allocatable::   B2_connect(:,:),B2_uconnect(:,:),unit_connect(:,:),B2_ufconnect(:,:)
!==============================================
integer(kind=8)::               Re_cal_num
integer(kind=8),allocatable::   Re_connect(:,:)
real(kind=8),allocatable::      Re_1B_E(:),Re_1B_F(:,:,:),Re_2B_E(:),Re_2B_F(:,:,:)
integer(kind=8)::               add_num
integer(kind=8),allocatable::   add_connect(:,:)
real(kind=8),allocatable::      add_1B_E(:),add_1B_F(:,:,:),add_2B_E(:),add_2B_F(:,:,:)
!===============================================
integer(kind=8)::               bg_num
real(kind=8),allocatable::      bg_crg(:),bg_crd(:,:),bg_force(:,:)
integer(kind=8),allocatable::   bg_elementnum(:),bg_element(:)
!===============================================
real(kind=8),allocatable::      fforce(:,:),cforce(:,:),Force(:,:),Dcount_F(:,:)
integer(kind=8),allocatable::   line_connect(:,:)
real(kind=8)::                  fenergy,cenergy,Dcount_E,Energy
real(kind=8),allocatable::      tmpforce(:,:),tmpfield(:,:),&
                                tmp_bgforce(:,:),tmp_bgfield(:,:),tmpcharge(:)
real(kind=8)::                  tmpenergy,tmpcenergy
integer(kind=8),allocatable::   connect_tag(:,:),connect_judge(:,:)
integer(kind=8),allocatable::   coordinate_atom(:)
!===============================================
integer(kind=8),parameter::     Hbond_max=1000
integer(kind=8)::               H_bond(2,Hbond_max),bond_num
!integer(kind=8),allocatable::   vdw_judge(:,:)
integer(kind=8),allocatable::   D_pairs(:,:)
integer(kind=8),allocatable::   tmpD_pairs(:,:)
integer(kind=8),allocatable::   B2D_pairs(:,:)
integer(kind=8),allocatable::   B2FD_pairs(:,:,:)
integer(kind=8),allocatable::   B2UD_pairs(:,:,:)
integer(kind=8),allocatable::   B2UFD_pairs(:,:,:)
integer(kind=8),allocatable::   B1FD_pairs(:,:,:)
integer(kind=8),allocatable::   B1UD_pairs(:,:,:)
real(kind=8),allocatable::      B2_BGF(:,:)
real(kind=8),allocatable::      B2F_BGF(:,:,:),B2U_BGF(:,:,:),B2UF_BGF(:,:,:)
real(kind=8),allocatable::      B1F_BGF(:,:,:),B1U_BGF(:,:,:)

!===============================================
integer(kind=4)::               ierr,nproc,myid,istatus(MPI_STATUS_SIZE)
!===============================================
integer::ncmd
integer::task_id(1000)
character(len=250)::task_cmd(1000)
end module
!==============================================
module dftb3
use mpi
integer(kind=8),allocatable::               element_tag(:)
character(len=2)::              element_list(20) 
character(len=1),allocatable::              MaxAngularMomentum(:)
real(kind=8),allocatable::                  HubbardDerivs(:)
character(len=1)::                          element_MaxAngM(20)
real(kind=8)::                              element_Hubbard(20)
character(len=200)::                         dftb_prmpath
character(len=200)::                         workpath
integer(kind=8)::                           element_typenum
namelist /dftb3_ctrl/ dftb_prmpath,workpath
end module

