# -*- coding: cp1251 -*-
import math, os
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='Arial') # ��� �������� ������
import PUmodel, maplepy

# ������ �������� SD1_s_rel0. ��������� ���� ����� ������ ����������� ���������� ��������� (������� S1_T0>0).

class Maplesim_model():
    """������ ���� � Maplesim
    ������� ����� ����� PUmodel.msim"""
    pu=None # ��'��� ����� PUmodel.PU
    ms=None # ��'��� ����� maplepy.MapleInterface4Maplesim
    
    def prepareParams(self):
        """ϳ��������� �������� Maplesim �����"""
        pu=self.pu # ����
        rs=pu.suckerRodString.items # ������ ������
        n=len(rs) # ������� ������
        
        M=pu.suckerRodString.massList() # ���� ������
        M2=pu.suckerRodString.fluidMassList() # ���� ����� ��� ������
        C2=pu.suckerRodString.hydrodynamicRodResistanceList()
        W=pu.suckerRodString.weightForceList()
        C=[r.stiffness() for r in rs]
        D=[r.materialRodResistance() for r in rs]
        v=pu.pumpingUnit.velAvg()
        
        AN=pu.suckerRodString.wellAngleList() # ������ ���� alfa (�������� n)
        ANS=pu.suckerRodString.secAngleList() # ���� ��������� ������ �� ��������
        AN1=pu.suckerRodString.angleTop()
        AN2=pu.suckerRodString.angleBottom()
                        
#         for p in 'M M2 C2 W C D v AN ANS AN1 AN2'.split():
#             print p, eval(p) # �����������
        
        # ���������� ��������� ��� ����� Maplesim
        modelParams={'S1_amplitude':pu.pumpingUnit.A, # ��������� ��������-��������
                     'S1_freqHz':pu.pumpingUnit.ns,
                     'S1_phase':pu.pumpingUnit.fi,
                     'HC1_A':pu.pump.pistonArea(), # ��������� ��������� �������
                     'CV1_Ropen':pu.pump.valveLossCoefficient(v),
                     'CV2_Ropen':pu.pump.valveLossCoefficient(v),
                     'CP1_D':pu.well.diametr, # ?
                     'CP1_L':pu.suckerRodString.length(),
                     #������ VerticalPipe1 ����� ��������������� FP2_P
                     #'VerticalPipe1_z':pu.suckerRodString.height(),
                     'AP1_P':0.0, # ���� �� ����
                     'AP2_P':pu.outPres, # ���� �� ����
                     # ������������� ����, ���� 䳺 �� ������� �� ��� ���� �����
                     # ����������� ����� ��������. ���� � ����. �������, ���� 䳺 � ����� ������� ��������
                     'FP2_P':pu.fluidWeight()/pu.pump.pistonArea()-pu.outTubeFluidWeight()/pu.pump.pistonArea()-pu.outTubePres,
                     'rhoFluid':pu.well.density,
                     'nuFluid':pu.well.viscosity,
                     # ������� ����� ������ �����
                     'FI1_A':(sum(M2)/pu.well.density)/pu.suckerRodString.length(),
                     # ������ FI1_L=0, ���� �� ����������� ������� �����
                     'FI1_L':pu.suckerRodString.length(),
                     'TF1_fs':1.25*pu.pump.pistonFrictionForce(),
                     'TF1_fc':pu.pump.pistonFrictionForce()
                     }
        #W[0]+=15000.0 # !!! ���� ����� � ���������        
        for i in range(n): # ��� ����� ������ ni (i=0, 1, 2...)
            modelParams['n'+str(i)+'_M1_m']=M[i]
            # !!! ����� ��� � ����� ����������� ����������� ��������
            modelParams['n'+str(i)+'_G2_k']=10*C2[i] # �������� � TF1_d!!! ������ ������� �� ��� ���� ����
            modelParams['n'+str(i)+'_CF1_f_constant']=-W[i] # ����� ����
            modelParams['n'+str(i)+'_TF1_fs']=0.2*W[i]*math.tan(AN[i]) # ���� ����� �� ���� max
            modelParams['n'+str(i)+'_TF1_fc']=0.16*W[i]*math.tan(AN[i]) # ���� ����� �� ����
            modelParams['n'+str(i)+'_TF1_d']=10*C2[i] # ����. �'������ ����� (�������� � G2_k!!!) ������ �������
            modelParams['n'+str(i)+'_G1_k']=0.16*math.sin(AN1[i]) # ������� ���� ����� �� ��������� ������
            modelParams['n'+str(i)+'_G4_k']=0.16*math.sin(AN2[i]) # ������� ���� ����� �� �������� ������
            modelParams['n'+str(i)+'_SD1_c']=C[i]
            modelParams['n'+str(i)+'_SD1_d']=D[i] # ������� ���� ������, ������� ���� �� ������������� �����������
        
#         for p in sorted(modelParams.keys()):
#             print p, modelParams[p]
            
        self.modelParams=modelParams
            
    def createMaplesimModel(self, execute=False, resFileName=None):
        """������� Maplesim ������. ������ ��, ���� execute=True"""       
        ms=maplepy.MapleInterface4Maplesim()
        ms.path=os.getcwd().replace('\\', '/')
        ms.filenameMaplesim="PUmodel.msim"
        ms.paramDictMaplesim=self.modelParams # {} ���� �� ���������� �����
        #print ms.getCode()
        if execute: ms.execute=True
        if resFileName: ms.resultCSVfile=resFileName
        ms.runMaple()
        self.ms=ms
            
    def drawDynamometerCard(self, S0=1.0):
        """������� ����-������� � ������������
        ������ ������ � PUmodel.msim �� ���� ����� Length � Force
        �� ���� ����������� (������ ������)
        S0 - ����� ����������� ������ ������� ������ (���� S0=1, �������� ����)"""    
        results=self.ms.readCSVfile(5) # ������ ����������
        
        # ���� ����� ������� �����
        T=1.0/self.modelParams['S1_freqHz'] # �����
        tend=results[-1][0] # ������� ���
        tstart=tend-T # ������� ���������� ������
        
        f=[] # ���� ��� ����������
        s=[] # ����������
        f1=[]
        s1=[]
        for res in results: # ��� ������� ����
            if res[0]>=tstart: # �������� � ���� tstart
                f.append(res[1]/S0) # !!!����������
                s.append(res[2]) # ����������
                f1.append(res[3])
                s1.append(res[4])
        plt.plot(s, f, 'k-')
        plt.plot(s1, f1, 'b-')
        plt.grid(True)
        plt.xlabel(u"S, �") #������ �� x
        if S0==1.0:
            plt.ylabel(u"F, �") #������ �� y
        else:
            plt.ylabel("s, Pa") #������ �� y
        plt.savefig("MaplesimDynCard.png")
        return min(f),max(f),min(s1),max(s1),min(f1),max(f1)
        
    def drawDynamometerCardMulti(self, filenames=[],S0L=[],freqL=[],styleL=[]):
        """������� ����-������� � ������������ ��� �������� �����-����������
        ���.������������ drawDynamometerCard"""    
        ms=maplepy.MapleInterface4Maplesim()
        ms.path=os.getcwd().replace('\\', '/')
        for name,S0,freq,style in zip(filenames,S0L,freqL,styleL):
            ms.resultCSVfile=name
            results=ms.readCSVfile(3) # ������ ����������
            
            # ���� ����� ������� �����
            T=1.0/freq # �����
            tend=results[-1][0] # ������� ���
            tstart=tend-T # ������� ���������� ������
            
            f=[] # ���� ��� ����������
            s=[] # ����������
            
            delta=0
            #if name=='13.csv': delta=0 # !!!�������� ��������� �����������, ���� ������� 
            
            for res in results: # ��� ������� ����
                if res[0]>=tstart: # �������� � ���� tstart
                    f.append(res[1]/S0-delta) # !!!����������
                    s.append(res[2]) # ����������
            plt.plot(s, f, style)
        plt.grid(True)
        plt.xlabel(u"S, �") #������ �� x
        plt.ylabel(u"F, �") #������ �� y
        plt.savefig("MaplesimDynCard.png")
                    
    def optimize(self):
        """��������� ������ ��� ����� ������� �� ���������"""
        parList=[i/10.0 for i in [1,2,3,4,5,6,7,8,9,10]] # ������ ������� ���������
        fl=open("ampl-freq.csv", "w")
        self.pu=PUmodel.PU()
        for par in parList: # ��� ������� ��������
            self.pu.pumpingUnit.ns=par # ������ �������� ���������
            self.prepareParams()
            #������� ����� ���������� x_MaplesimResults.csv
            #self.createMaplesimModel(True, str(par)+'_MaplesimResults.csv')
            self.createMaplesimModel(True)
            fl.write("%f;%f;%f;%f;%f;%f\n"%self.drawDynamometerCard())
            fl.flush()
        fl.close() 

class Maplesim_model2(Maplesim_model): # ��������� Maplesim_model
    """�������� ������ ���� � Maplesim
    ������� ����� ����� PUmodel.msim"""
    
    def prepareParams(self):
        """ϳ��������� ��������� Maplesim �����"""
        pu=self.pu # ����
        
        M=pu.suckerRodString.mass() # ���� ������
        W=pu.suckerRodString.weight() # ���� ������
        WF=pu.fluidWeight() # ���� �����
        C=pu.suckerRodString.stiffness() # ��������� ������
        D=pu.suckerRodString.damping()# ������������ ���������� �����
        #D=50*D # �����, �������� � 50! (��� �������� psi=0.5)
        #D=psi*math.sqrt(C*M)/(2*math.pi) # ?
        
        # ���������� ��������� ��� ����� Maplesim
        modelParams={'S1_amplitude':pu.pumpingUnit.A, # ��������� ��������-��������
                     'S1_freqHz':pu.pumpingUnit.ns,
                     'S1_phase':pu.pumpingUnit.fi,
                     'SD1_c':C, # ��������� ��������� �������
                     'SD1_d':D,
                     'M1_m':M,
                     'CF1_f_constant':-W,
                     'G1_k':-WF
                     }
        
        for p in sorted(modelParams.keys()):
            print p, modelParams[p]
            
        self.modelParams=modelParams

class Maplesim_model3(Maplesim_model): # ��������� Maplesim_model
    """������ ��������-�������� � Maplesim
    ������� ����� ����� PUmodel.msim"""
    
    def prepareParams(self):
        """ϳ��������� �������� Maplesim �����"""
        pu=self.pu # ����
        # ���������� ��������� ��� ����� Maplesim
        modelParams={'p1_L':pu.pumpingUnit.p1['L'], # ��������� �����
                     'p1_m':pu.pumpingUnit.p1['m'],
                     'p1_Iz':pu.pumpingUnit.p1['Iz'],
                     'p2_L':pu.pumpingUnit.p2['L'],
                     'p2_m':pu.pumpingUnit.p2['m'],
                     'p2_Iz':pu.pumpingUnit.p2['Iz'],
                     'p3_L':pu.pumpingUnit.p3['L'],
                     'p3_m':pu.pumpingUnit.p3['m'],
                     'p3_Iz':pu.pumpingUnit.p3['Iz'],
                     'p4_L':pu.pumpingUnit.p4['L'],
                     'p4_m':pu.pumpingUnit.p4['m'],
                     'p4_Iz':pu.pumpingUnit.p4['Iz'],
                     'p5_L':pu.pumpingUnit.p5['L'], # ����� ���������
                     'p5_m':pu.pumpingUnit.p5['m'], # ���� ���������
                     'p5_Iz':pu.pumpingUnit.p5['Iz'],
                     'dx':pu.pumpingUnit.dx, # ������ �� �����
                     'dy':pu.pumpingUnit.dy,
                     'r_k':pu.pumpingUnit.r, # ����� ���������
                     'const1_k':-pu.pumpingUnit.omega # ������ �������� ���������, ���/�
                     # ���� "-" ��� ������� ��������� (������� �.15)
                     # ��� ��� ����� ������
                     # �� ������ FF1_teta=pi/2
                     }
        
        for p in sorted(modelParams.keys()):
            print p, modelParams[p]
            
        self.modelParams=modelParams
                
if __name__ == '__main__':
    model=Maplesim_model() # !!! ������ ����� ����
    #model.optimize()
    
    model.pu=PUmodel.PU()
    print "Fmin,Fmax=", model.pu.polishedRodForce()
    #print "Rp=", model.pu.balancing()
    model.prepareParams()
    
    model.createMaplesimModel(execute=1)
    #model.drawDynamometerCard(S0=math.pi*0.019**2/4) #-math.pi*0.0085**2/4)
    #model.drawDynamometerCard()
    # ����� �������� �������� ���������� � ���������� maplesim plots � .csv
    #model.drawDynamometerCardMulti(['1.csv','2.csv','3.csv'],[326.0e-6,326.0e-6,326.0e-6],[0.166666,0.166666,0.166666],['k-','k--','k:'])
    
    #f=model.modelParams['S1_freqHz']
    model.drawDynamometerCardMulti(['testAPI.csv'],[1.0],[0.10833],['k-'])
    #model.drawDynamometerCardMulti(['2_.csv','1_.csv','2.csv','1.csv'],[1.0,1.0,1.0,1.0],[0.10833,0.10833,f,f],['k--','k--','k-','k-'])
    