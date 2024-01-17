import numpy as np
import pandas as pd
import math
class Data:
    def __init__(self,baseVA,wbase,isFrequencyKnown,Rx,bus,branch,load,generator):
        self.baseVA = baseVA
        self.wbase = wbase
        self.isFrequencyKnown = isFrequencyKnown
        self.Rx = Rx
        self.bus = bus
        self.branch = branch
        self.load = load
        self.generator = generator

class Branch:
    def __init__(self,FromBus,ToBus,R,L):
        self.FromBus = FromBus
        self.ToBus = ToBus
        self.R = R
        self.L = L

class Load:
    def __init__(self,Bus,R,L):
        self.Bus = Bus
        self.R = R
        self.L = L

class GFM_IBR:
    def __init__(self,mp,mq,P0,Q0,w0,V0,Rv,Lv,KpV,KiV,KpC,KiC,Rt,Lt,Rd,Cf,Rg,Lg,wc):
        self.mp = mp
        self.mq = mq
        self.P0 = P0
        self.Q0 = Q0
        self.w0 = w0
        self.V0 = V0        
        self.Rv = Rv
        self.Lv = Lv
        self.KpV = KpV
        self.KiV = KiV
        self.KpC = KpC
        self.KiC = KiC
        self.Rt = Rt
        self.Lt = Lt
        self.Rd = Rd
        self.Cf = Cf
        self.Rg = Rg
        self.Lg = Lg
        self.wc = wc

class GFL_IBR:
    def __init__(self,mp,mq,P0,Q0,w0,V0,KpL,KiL,KpS,KiS,KpC,KiC,Rt,Lt,Rd,Cf,Rg,Lg,wc,wcf):
        self.mp = mp
        self.mq = mq
        self.P0 = P0
        self.Q0 = Q0
        self.w0 = w0
        self.V0 = V0        
        self.KpL = KpL
        self.KiL = KiL
        self.KpS = KpS
        self.KiS = KiS
        self.KpC = KpC
        self.KiC = KiC
        self.Rt = Rt
        self.Lt = Lt
        self.Rd = Rd
        self.Cf = Cf
        self.Rg = Rg
        self.Lg = Lg
        self.wc = wc
        self.wcf = wcf

def case_3bus_droop_gfl():
    baseVA = 100e6
    wbase = 2*math.pi*60
    isFrequencyKnown = False
    Rx = 2e3
    bus = pd.DataFrame([[1,'Bus1','PQ-Ref',13.8e3],
                        [2,'Bus2','PQ',13.8e3],
                        [3,'Bus3','PQ',13.8e3]])
    branch = [Branch(1,3,0.01,0.05),
              Branch(2,3,0.01,0.05)]
    load = [Load(3,0.90,0.20)]
    generator = [GFM_IBR(0.05,0.05,0.5,0.2,1.0,1.0,0.00,0.00,0.6,3.6,0.6,12.0,0.02,0.10,100,0.05,0.01,0.05,31.4159),
                 GFL_IBR(0.05,0.05,0.5,0.2,1.0,1.0,1.20,32.0,0.6,3.6,0.6,06.0,0.02,0.10,100,0.05,0.01,0.05,31.4159,100)]
    sysData = Data(baseVA,wbase,isFrequencyKnown,Rx,bus,branch,load,generator)    
    return sysData

def case_3bus_2droop():
    baseVA = 100e6
    wbase = 2*math.pi*60
    isFrequencyKnown = False
    Rx = 2e3
    bus = pd.DataFrame([[1,'Bus1','PQ-Ref',13.8e3],
                        [2,'Bus2','PQ',13.8e3],
                        [3,'Bus3','PQ',13.8e3]])
    branch = [Branch(1,3,0.01,0.05),
              Branch(2,3,0.01,0.05)]
    load = [Load(3,0.90,0.20)]
    generator = [GFM_IBR(0.05,0.05,0.2,0.2,1.0,1.0,0.00,0.00,0.4,9.6,0.3,6.0,0.02,0.10,100,0.05,0.01,0.05,31.4159),
                 GFM_IBR(0.05,0.05,0.8,0.2,1.0,1.0,0.00,0.00,0.4,9.6,0.3,6.0,0.02,0.10,100,0.05,0.01,0.05,31.4159)]
    sysData = Data(baseVA,wbase,isFrequencyKnown,Rx,bus,branch,load,generator)    
    return sysData
    
    
