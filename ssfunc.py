import numpy as np
from sympy import *
import math

def pf_calc(a,b):
    c = a + b

    return c


def lsm_sys_droop_gfl(sysData):
    xdot = Matrix([0]*33)
    
    delta_1 = symbols('delta_1')
    Po_1 = symbols('Po_1')
    Qo_1 = symbols('Qo_1')
    phid_1 = symbols('phid_1')
    phiq_1 = symbols('phiq_1')
    gammad_1 = symbols('gammad_1')
    gammaq_1 = symbols('gammaq_1')
    iid_1 = symbols('iid_1')
    iiq_1 = symbols('iiq_1')
    vcd_1 = symbols('vcd_1')
    vcq_1 = symbols('vcq_1')
    iod_1 = symbols('iod_1')
    ioq_1 = symbols('ioq_1')
    delta_2 = symbols('delta_2')
    epsilonL_2 = symbols('epsilonL_2')
    wf_2 = symbols('wf_2')
    Po_2 = symbols('Po_2')
    Qo_2 = symbols('Qo_2')
    phid_2 = symbols('phid_2')
    phiq_2 = symbols('phiq_2')
    gammad_2 = symbols('gammad_2')
    gammaq_2 = symbols('gammaq_2')
    iid_2 = symbols('iid_2')
    iiq_2 = symbols('iiq_2')
    vcd_2 = symbols('vcd_2')
    vcq_2 = symbols('vcq_2')
    iod_2 = symbols('iod_2')
    ioq_2 = symbols('ioq_2')
    ibranchD_1 = symbols('ibranchD_1')
    ibranchQ_1 = symbols('ibranchQ_1')
    ibranchD_2 = symbols('ibranchD_2')
    ibranchQ_2 = symbols('ibranchQ_2')
    iloadD = symbols('iloadD')
    iloadQ = symbols('iloadQ')
    
    wbase = sysData.wbase
    Rx = sysData.Rx
    
    mp_1 = sysData.generator[0].mp
    mq_1 = sysData.generator[0].mq
    P0_1 = sysData.generator[0].P0
    Q0_1 = sysData.generator[0].Q0
    w0_1 = sysData.generator[0].w0
    V0_1 = sysData.generator[0].V0
    Rv_1 = sysData.generator[0].Rv
    Lv_1 = sysData.generator[0].Lv
    KpV_1 = sysData.generator[0].KpV
    KiV_1 = sysData.generator[0].KiV
    KpC_1 = sysData.generator[0].KpC
    KiC_1 = sysData.generator[0].KiC
    Rt_1 = sysData.generator[0].Rt
    Lt_1 = sysData.generator[0].Lt
    Rd_1 = sysData.generator[0].Rd
    Cf_1 = sysData.generator[0].Cf
    Rg_1 = sysData.generator[0].Rg
    Lg_1 = sysData.generator[0].Lg
    wc_1 = sysData.generator[0].wc
    
    mp_2 = sysData.generator[1].mp
    mq_2 = sysData.generator[1].mq
    P0_2 = sysData.generator[1].P0
    Q0_2 = sysData.generator[1].Q0
    w0_2 = sysData.generator[1].w0
    V0_2 = sysData.generator[1].V0
    KpL_2 = sysData.generator[1].KpL
    KiL_2 = sysData.generator[1].KiL
    KpS_2 = sysData.generator[1].KpS
    KiS_2 = sysData.generator[1].KiS
    KpC_2 = sysData.generator[1].KpC
    KiC_2 = sysData.generator[1].KiC
    Rt_2 = sysData.generator[1].Rt
    Lt_2 = sysData.generator[1].Lt
    Rd_2 = sysData.generator[1].Rd
    Cf_2 = sysData.generator[1].Cf
    Rg_2 = sysData.generator[1].Rg
    Lg_2 = sysData.generator[1].Lg
    wc_2 = sysData.generator[1].wc
    wcf_2 = sysData.generator[1].wcf

    Rbranch_1 = sysData.branch[0].R
    Lbranch_1 = sysData.branch[0].L
    Rbranch_2 = sysData.branch[1].R
    Lbranch_2 = sysData.branch[1].L
    Rload = sysData.load[0].R
    Lload = sysData.load[0].L
    
    winv_1 = w0_1 - mp_1*(Po_1 - P0_1)
    vod_1 = vcd_1 + Rd_1*(iid_1 - iod_1)
    voq_1 = vcq_1 + Rd_1*(iiq_1 - ioq_1)
    vdVirtual_1 = Rv_1*iod_1 - winv_1*Lv_1*ioq_1
    vqVirtual_1 = Rv_1*ioq_1 + winv_1*Lv_1*iod_1
    vodRef_1 = V0_1 - mq_1*(Qo_1 - Q0_1) - vdVirtual_1
    voqRef_1 = 0 - vqVirtual_1
    iidRef_1 = KpV_1*(vodRef_1 - vod_1) + KiV_1*phid_1
    iiqRef_1 = KpV_1*(voqRef_1 - voq_1) + KiV_1*phiq_1
    vidRef_1 = KpC_1*(iidRef_1 - iid_1) + KiC_1*gammad_1 - w0_1*Lt_1*iiq_1
    viqRef_1 = KpC_1*(iiqRef_1 - iiq_1) + KiC_1*gammaq_1 + w0_1*Lt_1*iid_1
    ioD_1 = (iod_1*cos(0) - ioq_1*sin(0))
    ioQ_1 = (iod_1*sin(0) + ioq_1*cos(0))

    vod_2 = vcd_2 + Rd_2*(iid_2 - iod_2)
    voq_2 = vcq_2 + Rd_2*(iiq_2 - ioq_2)
    winv_2 = KpL_2*voq_2 + KiL_2*epsilonL_2 + w0_2
    Pref_2 = (w0_2 - wf_2)/mp_2 + P0_2
    Qref_2 = (V0_2 - sqrt(vod_2*vod_2+voq_2*voq_2))/mq_2 + Q0_2
    iidRef_2 = KpS_2*(Pref_2 - Po_2) + KiS_2*phid_2
    iiqRef_2 = KpS_2*(Qo_2 - Qref_2) + KiS_2*phiq_2
    vidRef_2 = KpC_2*(iidRef_2 - iid_2) + KiC_2*gammad_2 - w0_2*Lt_2*iiq_2
    viqRef_2 = KpC_2*(iiqRef_2 - iiq_2) + KiC_2*gammaq_2 + w0_2*Lt_2*iid_2
    ioD_2 = (iod_2*cos(delta_2) - ioq_2*sin(delta_2))
    ioQ_2 = (iod_2*sin(delta_2) + ioq_2*cos(delta_2))  
    
    wcom = winv_1

    vbD_1 = Rx*(ioD_1 - ibranchD_1)
    vbQ_1 = Rx*(ioQ_1 - ibranchQ_1)
    vbD_2 = Rx*(ioD_2 - ibranchD_2)
    vbQ_2 = Rx*(ioQ_2 - ibranchQ_2)
    vbD_3 = Rx*(ibranchD_1 + ibranchD_2 - iloadD)
    vbQ_3 = Rx*(ibranchQ_1 + ibranchQ_2 - iloadQ)   
    vbd_1 = vbD_1*cos(0) + vbQ_1*sin(0)
    vbq_1 = -vbD_1*sin(0) + vbQ_1*cos(0)
    vbd_2 = vbD_2*cos(delta_2) + vbQ_2*sin(delta_2)
    vbq_2 = -vbD_2*sin(delta_2) + vbQ_2*cos(delta_2)

    x = Matrix([Po_1,Qo_1,phid_1,phiq_1,gammad_1,gammaq_1,iid_1,iiq_1,vcd_1,vcq_1,iod_1,ioq_1,
                delta_2,epsilonL_2,wf_2,Po_2,Qo_2,phid_2,phiq_2,gammad_2,gammaq_2,iid_2,iiq_2,vcd_2,vcq_2,iod_2,ioq_2,
                ibranchD_1,ibranchQ_1,ibranchD_2,ibranchQ_2,iloadD,iloadQ])
    
    #xdot[0] = wbase*(winv_1 - wcom)
    xdot[0] = -wc_1*Po_1 + wc_1*(vod_1*iod_1 + voq_1*ioq_1)
    xdot[1] = -wc_1*Qo_1 + wc_1*(voq_1*iod_1 - vod_1*ioq_1)
    xdot[2] = vodRef_1 - vod_1
    xdot[3] = voqRef_1 - voq_1
    xdot[4] = iidRef_1 - iid_1
    xdot[5] = iiqRef_1 - iiq_1
    xdot[6] = wbase*(vidRef_1 - vod_1 - Rt_1*iid_1 + winv_1*Lt_1*iiq_1)/Lt_1
    xdot[7] = wbase*(viqRef_1 - voq_1 - Rt_1*iiq_1 - winv_1*Lt_1*iid_1)/Lt_1
    xdot[8] = wbase*(iid_1 - iod_1 + winv_1*Cf_1*vcq_1)/Cf_1
    xdot[9] = wbase*(iiq_1 - ioq_1 - winv_1*Cf_1*vcd_1)/Cf_1
    xdot[10] = wbase*(vod_1 - vbd_1 - Rg_1*iod_1 + winv_1*Lg_1*ioq_1)/Lg_1
    xdot[11] = wbase*(voq_1 - vbq_1 - Rg_1*ioq_1 - winv_1*Lg_1*iod_1)/Lg_1
    xdot[12] = wbase*(winv_2 - wcom)
    xdot[13] = voq_2
    xdot[14] = -wcf_2*wf_2 + wcf_2*winv_2
    xdot[15] = -wc_2*Po_2 + wc_2*(vod_2*iod_2 + voq_2*ioq_2)
    xdot[16] = -wc_2*Qo_2 + wc_2*(voq_2*iod_2 - vod_2*ioq_2)
    xdot[17] = Pref_2 - Po_2
    xdot[18] = Qo_2 - Qref_2
    xdot[19] = iidRef_2 - iid_2
    xdot[20] = iiqRef_2 - iiq_2
    xdot[21] = wbase*(vidRef_2 - vod_2 - Rt_2*iid_2 + winv_2*Lt_2*iiq_2)/Lt_2
    xdot[22] = wbase*(viqRef_2 - voq_2 - Rt_2*iiq_2 - winv_2*Lt_2*iid_2)/Lt_2
    xdot[23] = wbase*(iid_2 - iod_2 + winv_2*Cf_2*vcq_2)/Cf_2
    xdot[24] = wbase*(iiq_2 - ioq_2 - winv_2*Cf_2*vcd_2)/Cf_2
    xdot[25] = wbase*(vod_2 - vbd_2 - Rg_2*iod_2 + winv_2*Lg_2*ioq_2)/Lg_2
    xdot[26] = wbase*(voq_2 - vbq_2 - Rg_2*ioq_2 - winv_2*Lg_2*iod_2)/Lg_2    
    xdot[27] = wbase*(vbD_1 - vbD_3 - Rbranch_1*ibranchD_1 + wcom*Lbranch_1*ibranchQ_1)/Lbranch_1;
    xdot[28] = wbase*(vbQ_1 - vbQ_3 - Rbranch_1*ibranchQ_1 - wcom*Lbranch_1*ibranchD_1)/Lbranch_1;
    xdot[29] = wbase*(vbD_2 - vbD_3 - Rbranch_2*ibranchD_2 + wcom*Lbranch_2*ibranchQ_2)/Lbranch_2;
    xdot[30] = wbase*(vbQ_2 - vbQ_3 - Rbranch_2*ibranchQ_2 - wcom*Lbranch_2*ibranchD_2)/Lbranch_2;
    xdot[31] = wbase*(vbD_3 - Rload*iloadD + wcom*Lload*iloadQ)/Lload;
    xdot[32] = wbase*(vbQ_3 - Rload*iloadQ - wcom*Lload*iloadD)/Lload;
    return x, xdot

def lsm_sys_2droop(sysData):
    xdot = Matrix([0]*31)
    
    delta_1 = symbols('delta_1')
    Po_1 = symbols('Po_1')
    Qo_1 = symbols('Qo_1')
    phid_1 = symbols('phid_1')
    phiq_1 = symbols('phiq_1')
    gammad_1 = symbols('gammad_1')
    gammaq_1 = symbols('gammaq_1')
    iid_1 = symbols('iid_1')
    iiq_1 = symbols('iiq_1')
    vcd_1 = symbols('vcd_1')
    vcq_1 = symbols('vcq_1')
    iod_1 = symbols('iod_1')
    ioq_1 = symbols('ioq_1')
    delta_2 = symbols('delta_2')
    Po_2 = symbols('Po_2')
    Qo_2 = symbols('Qo_2')
    phid_2 = symbols('phid_2')
    phiq_2 = symbols('phiq_2')
    gammad_2 = symbols('gammad_2')
    gammaq_2 = symbols('gammaq_2')
    iid_2 = symbols('iid_2')
    iiq_2 = symbols('iiq_2')
    vcd_2 = symbols('vcd_2')
    vcq_2 = symbols('vcq_2')
    iod_2 = symbols('iod_2')
    ioq_2 = symbols('ioq_2')
    ibranchD_1 = symbols('ibranchD_1')
    ibranchQ_1 = symbols('ibranchQ_1')
    ibranchD_2 = symbols('ibranchD_2')
    ibranchQ_2 = symbols('ibranchQ_2')
    iloadD = symbols('iloadD')
    iloadQ = symbols('iloadQ')
    
    wbase = sysData.wbase
    Rx = sysData.Rx
    
    mp_1 = sysData.generator[0].mp
    mq_1 = sysData.generator[0].mq
    P0_1 = sysData.generator[0].P0
    Q0_1 = sysData.generator[0].Q0
    w0_1 = sysData.generator[0].w0
    V0_1 = sysData.generator[0].V0
    Rv_1 = sysData.generator[0].Rv
    Lv_1 = sysData.generator[0].Lv
    KpV_1 = sysData.generator[0].KpV
    KiV_1 = sysData.generator[0].KiV
    KpC_1 = sysData.generator[0].KpC
    KiC_1 = sysData.generator[0].KiC
    Rt_1 = sysData.generator[0].Rt
    Lt_1 = sysData.generator[0].Lt
    Rd_1 = sysData.generator[0].Rd
    Cf_1 = sysData.generator[0].Cf
    Rg_1 = sysData.generator[0].Rg
    Lg_1 = sysData.generator[0].Lg
    wc_1 = sysData.generator[0].wc
    
    mp_2 = sysData.generator[0].mp
    mq_2 = sysData.generator[0].mq
    P0_2 = sysData.generator[0].P0
    Q0_2 = sysData.generator[0].Q0
    w0_2 = sysData.generator[0].w0
    V0_2 = sysData.generator[0].V0
    Rv_2 = sysData.generator[0].Rv
    Lv_2 = sysData.generator[0].Lv
    KpV_2 = sysData.generator[0].KpV
    KiV_2 = sysData.generator[0].KiV
    KpC_2 = sysData.generator[0].KpC
    KiC_2 = sysData.generator[0].KiC
    Rt_2 = sysData.generator[0].Rt
    Lt_2 = sysData.generator[0].Lt
    Rd_2 = sysData.generator[0].Rd
    Cf_2 = sysData.generator[0].Cf
    Rg_2 = sysData.generator[0].Rg
    Lg_2 = sysData.generator[0].Lg
    wc_2 = sysData.generator[0].wc

    Rbranch_1 = sysData.branch[0].R
    Lbranch_1 = sysData.branch[0].L
    Rbranch_2 = sysData.branch[1].R
    Lbranch_2 = sysData.branch[1].L
    Rload = sysData.load[0].R
    Lload = sysData.load[0].L
    
    winv_1 = w0_1 - mp_1*(Po_1 - P0_1)
    vod_1 = vcd_1 + Rd_1*(iid_1 - iod_1)
    voq_1 = vcq_1 + Rd_1*(iiq_1 - ioq_1)
    vdVirtual_1 = Rv_1*iod_1 - winv_1*Lv_1*ioq_1
    vqVirtual_1 = Rv_1*ioq_1 + winv_1*Lv_1*iod_1
    vodRef_1 = V0_1 - mq_1*(Qo_1 - Q0_1) - vdVirtual_1
    voqRef_1 = 0 - vqVirtual_1
    iidRef_1 = KpV_1*(vodRef_1 - vod_1) + KiV_1*phid_1
    iiqRef_1 = KpV_1*(voqRef_1 - voq_1) + KiV_1*phiq_1
    vidRef_1 = KpC_1*(iidRef_1 - iid_1) + KiC_1*gammad_1 - w0_1*Lt_1*iiq_1
    viqRef_1 = KpC_1*(iiqRef_1 - iiq_1) + KiC_1*gammaq_1 + w0_1*Lt_1*iid_1
    ioD_1 = (iod_1*cos(0) - ioq_1*sin(0))
    ioQ_1 = (iod_1*sin(0) + ioq_1*cos(0))

    winv_2 = w0_2 - mp_2*(Po_2 - P0_2)
    vod_2 = vcd_2 + Rd_2*(iid_2 - iod_2)
    voq_2 = vcq_2 + Rd_2*(iiq_2 - ioq_2)
    vdVirtual_2 = Rv_2*iod_2 - winv_2*Lv_2*ioq_2
    vqVirtual_2 = Rv_2*ioq_2 + winv_2*Lv_2*iod_2
    vodRef_2 = V0_2 - mq_2*(Qo_2 - Q0_2) - vdVirtual_2
    voqRef_2 = 0 - vqVirtual_2
    iidRef_2 = KpV_2*(vodRef_2 - vod_2) + KiV_2*phid_2
    iiqRef_2 = KpV_2*(voqRef_2 - voq_2) + KiV_2*phiq_2
    vidRef_2 = KpC_2*(iidRef_2 - iid_2) + KiC_2*gammad_2 - w0_2*Lt_2*iiq_2
    viqRef_2 = KpC_2*(iiqRef_2 - iiq_2) + KiC_2*gammaq_2 + w0_2*Lt_2*iid_2
    ioD_2 = (iod_2*cos(delta_2) - ioq_2*sin(delta_2))
    ioQ_2 = (iod_2*sin(delta_2) + ioq_2*cos(delta_2))
    
    wcom = winv_1

    vbD_1 = Rx*(ioD_1 - ibranchD_1)
    vbQ_1 = Rx*(ioQ_1 - ibranchQ_1)
    vbD_2 = Rx*(ioD_2 - ibranchD_2)
    vbQ_2 = Rx*(ioQ_2 - ibranchQ_2)
    vbD_3 = Rx*(ibranchD_1 + ibranchD_2 - iloadD)
    vbQ_3 = Rx*(ibranchQ_1 + ibranchQ_2 - iloadQ)   
    vbd_1 = vbD_1*cos(0) + vbQ_1*sin(0)
    vbq_1 = -vbD_1*sin(0) + vbQ_1*cos(0)
    vbd_2 = vbD_2*cos(delta_2) + vbQ_2*sin(delta_2)
    vbq_2 = -vbD_2*sin(delta_2) + vbQ_2*cos(delta_2)

    x = Matrix([Po_1,Qo_1,phid_1,phiq_1,gammad_1,gammaq_1,iid_1,iiq_1,vcd_1,vcq_1,iod_1,ioq_1,
                delta_2,Po_2,Qo_2,phid_2,phiq_2,gammad_2,gammaq_2,iid_2,iiq_2,vcd_2,vcq_2,iod_2,ioq_2,
                ibranchD_1,ibranchQ_1,ibranchD_2,ibranchQ_2,iloadD,iloadQ])
    
    #xdot[0] = wbase*(winv_1 - wcom)
    xdot[0] = -wc_1*Po_1 + wc_1*(vod_1*iod_1 + voq_1*ioq_1)
    xdot[1] = -wc_1*Qo_1 + wc_1*(voq_1*iod_1 - vod_1*ioq_1)
    xdot[2] = vodRef_1 - vod_1
    xdot[3] = voqRef_1 - voq_1
    xdot[4] = iidRef_1 - iid_1
    xdot[5] = iiqRef_1 - iiq_1
    xdot[6] = wbase*(vidRef_1 - vod_1 - Rt_1*iid_1 + winv_1*Lt_1*iiq_1)/Lt_1
    xdot[7] = wbase*(viqRef_1 - voq_1 - Rt_1*iiq_1 - winv_1*Lt_1*iid_1)/Lt_1
    xdot[8] = wbase*(iid_1 - iod_1 + winv_1*Cf_1*vcq_1)/Cf_1
    xdot[9] = wbase*(iiq_1 - ioq_1 - winv_1*Cf_1*vcd_1)/Cf_1
    xdot[10] = wbase*(vod_1 - vbd_1 - Rg_1*iod_1 + winv_1*Lg_1*ioq_1)/Lg_1
    xdot[11] = wbase*(voq_1 - vbq_1 - Rg_1*ioq_1 - winv_1*Lg_1*iod_1)/Lg_1
    xdot[12] = wbase*(winv_2 - wcom)
    xdot[13] = -wc_2*Po_2 + wc_2*(vod_2*iod_2 + voq_2*ioq_2)
    xdot[14] = -wc_2*Qo_2 + wc_2*(voq_2*iod_2 - vod_2*ioq_2)
    xdot[15] = vodRef_2 - vod_2
    xdot[16] = voqRef_2 - voq_2
    xdot[17] = iidRef_2 - iid_2
    xdot[18] = iiqRef_2 - iiq_2
    xdot[19] = wbase*(vidRef_2 - vod_2 - Rt_2*iid_2 + winv_2*Lt_2*iiq_2)/Lt_2
    xdot[20] = wbase*(viqRef_2 - voq_2 - Rt_2*iiq_2 - winv_2*Lt_2*iid_2)/Lt_2
    xdot[21] = wbase*(iid_2 - iod_2 + winv_2*Cf_2*vcq_2)/Cf_2
    xdot[22] = wbase*(iiq_2 - ioq_2 - winv_2*Cf_2*vcd_2)/Cf_2
    xdot[23] = wbase*(vod_2 - vbd_2 - Rg_2*iod_2 + winv_2*Lg_2*ioq_2)/Lg_2
    xdot[24] = wbase*(voq_2 - vbq_2 - Rg_2*ioq_2 - winv_2*Lg_2*iod_2)/Lg_2  
    xdot[25] = wbase*(vbD_1 - vbD_3 - Rbranch_1*ibranchD_1 + wcom*Lbranch_1*ibranchQ_1)/Lbranch_1;
    xdot[26] = wbase*(vbQ_1 - vbQ_3 - Rbranch_1*ibranchQ_1 - wcom*Lbranch_1*ibranchD_1)/Lbranch_1;
    xdot[27] = wbase*(vbD_2 - vbD_3 - Rbranch_2*ibranchD_2 + wcom*Lbranch_2*ibranchQ_2)/Lbranch_2;
    xdot[28] = wbase*(vbQ_2 - vbQ_3 - Rbranch_2*ibranchQ_2 - wcom*Lbranch_2*ibranchD_2)/Lbranch_2;
    xdot[29] = wbase*(vbD_3 - Rload*iloadD + wcom*Lload*iloadQ)/Lload;
    xdot[30] = wbase*(vbQ_3 - Rload*iloadQ - wcom*Lload*iloadD)/Lload;
    return x, xdot

