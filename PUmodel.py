# -*- coding: cp1251 -*-
#������� �.�. ���������� ������ �����: ������� ������� ��� �����. � �: ���� ���-�� ������ � ��� ��� ����� � ���� ��. �.�. �������, 2003. � 816 �.
import math

# ��� ����������� �������� �������� ������ �� ��������
# ������������ ��������� � ������ �:
# ������� ��������, ��������, �����
# ��� ��� � �� ���� ����� ������ ���� � �������
def eqYoungsModulus(E,D,d):
    """������������ ������ ���� �������� ������ ������� D
    � ����� ���������� �������� k=E*A/L=Eeq*Aeq/L
    �������� (�������� D) � �������� ������ (����. �������� D,
    �������� d � ������� ���� E)"""
    A=math.pi*(D**2.0-d**2.0)/4.0 # ����� ������ �����
    Aeq=math.pi*(D**2.0)/4.0 # ����� ������ �������� ������
    return E*A/Aeq

def eqDensity(Ro,D,d,Roo=0.0):
    """����������� ������� �������� ������ ������� D
    � ����� ���������� ���� m=Ro*A*L+Roo*Ao*L=Roeq*Aeq*L
    �������� (�������� D) � �������� ������ (����. �������� D,
    �������� d, �������� ��� Ro � ��� ������ Roo)
    ���� Roo=0 - ������� ������ ��������� �������"""
    A=math.pi*(D**2.0-d**2.0)/4.0 # ����� ������ �����
    Ao=math.pi*(d**2.0)/4.0 # ����� ������ ������
    Aeq=math.pi*(D**2.0)/4.0 # ����� ������ �������� ������
    return (Ro*A+Roo*Ao)/Aeq
    

class SuckerRod(object):
    """���� �������� ����� ��� ������ ���������� �����. ������� Ѳ"""
    diametrs=[0.016,0.019,0.022,0.025] # ������� ���
    length=8.0 # ������� (�� �������������)
    materials={"�����":{"������ ��������":2.1e11,
                        "�������":8000.0}, # ��� ������� ���� �������� ��� ���������� ��� ���� � �������
               "�����������":{"������ ��������":0.5e11,
                              "�������":2100.0},
               "�������22���":{"������ ��������":eqYoungsModulus(2.1e11,0.022,0.0085),
                              "�������":eqDensity(8000.0,0.022,0.0085,1000.0)},
               "�������������22���":{"������ ��������":eqYoungsModulus(0.5e11,0.022,0.0085),
                              "�������":eqDensity(2100.0,0.022,0.0085,1000.0)}}
    #!!!���� Roo=0, �� ������� ������ � ������������ � ��������� ������� (��� �� ����������)
    alfa=0.0 # ��� ��������� ����������� ����� ������ �� �������� (���.)
    
    def __init__(self, PU, diametr, length=8.0, material="�����", alfa=0.0):
        """�����������"""
        self.PU=PU # ������������ �������� ������� ���������
        self.diametr=diametr # ������ ��� ������
        self.length=length # ������� ������ ��� ������
        self.material=material # �������
        self.alfa=alfa
    
    def area(self):
        """����� ����������� ������ S=pi*r^2"""
        return math.pi*(self.diametr/2.0)**2
    
    def volume(self):
        """��'��, ���� ����� ������ � ���� V=L*pi*r^2"""
        return self.length*math.pi*(self.diametr/2.0)**2
    
    def mass(self):
        """���� m=V*ro"""
        return self.volume()*self.materials[self.material]["�������"]
    
    def fluidVolume(self):
        """��'�� ����� � ��� ������� ������ (������)"""
        A=self.PU.well.sectionalArea()-self.area() # �����! ��� ����� ���, � �� ��������
        return A*self.length
    
    def fluidMass(self):
        """���� ����� � ��� ������� ������ (������)"""
        return self.PU.well.density*self.fluidVolume()
    
    def weight(self):
        """���� � ����"""
        fluidDensity=self.PU.well.density # ������� �����
        Farh=fluidDensity*9.81*self.volume() # ���� ��������
        w=self.mass()*9.81 # ���� � �����
        return w-Farh
    
    def stiffness(self):
        """��������� k=E*S/L"""
        E=self.materials[self.material]["������ ��������"]
        return E*self.area()/self.length
    
    def hydrodynamicRodResistance(self):
        """ó������������ ��� '�' ������ (������) �� ϳ���������
        ����������� �� ��� ���� ����. ���� ����� F=c*v"""
        msh=self.PU.well.diametr/self.diametr
        Msh=1.0/((math.log(msh)*(msh**2+1)/(msh**2-1))-1) # ���������� ���
        # v=2nS - �������� �����
        R=math.pi**2*self.PU.well.viscosity*self.PU.well.density*Msh
        return R*self.length
    
    def hydrodynamicRodResistance2(self):
        """ó������������ ��� '�' ������ (������) �� ��������
        ����������� �� ��� ���� ����. ���� ����� F=c*v"""
        m=self.diametr/self.PU.well.diametr
        R=16.9*self.PU.well.viscosity*self.PU.well.density*m**5.49
        R=R*1000 # ������� ?
        return R*self.length
    
    def materialRodResistance(self,sigmaA=None):
        """������������ ���������� ����� ��������, �*�/�
        sigmaA - �������� ��������� (��)
        sigmaA=None - �� �������� �� ����. ������."""
        k=self.stiffness() #��������� �������
        #��������� ���������� ���������� ��������
        if self.material.startswith("�����"):
            if sigmaA!=None: #���� �������� �� ����. ������.
                psi=6e-5*sigmaA/1000000.0
            else: #���� �� �������� �� ����. ������.
                #psi<0.13 (������ � �������� �������)[Damping of materials and members in structures]
                #psi=0.25...0.5 (������� ���������)
                #psi=0.38...0.88(������� ��������� � ���������)
                psi=0.01
        if self.material.startswith("�����������"):
            if sigmaA!=None: #���� �������� �� ����. ������.
                psi=4e-4*sigmaA/1000000.0 #���������� "�������� � �������"
            else: #���� �� �������� �� ����. ������.
                psi=0.08
        omega=self.PU.pumpingUnit.omega #������ ������� �������� ����
        #m=self.mass() # ���� �������
        #omega=math.sqrt(k/m) # ������ ������� �������
        #[�������� � ������� �6 �130 �13,19]
        c=psi*k/(2*math.pi*omega) #k=E*S/L
        
#         #!�������� ���������� ������ ����� ������
#         nz=int(self.length/8) # ������ �������� ������ �� ������
#         kz=1000000.0 # ��������� ��������� ������� ?
#         cz=0.2*kz/(2*math.pi*omega) #��� ���� ����� ��������� ������� (psi=0,006...0,4)
#         czs=1/(nz*(1/cz)) # ��� ���� ����� ��� ������ �� ������
#         return 1/(1/c+1/czs)# ��� ���� ����� ������ (�������+�������)
        return 20.0*c # ����� �������� � 20.0 (��� �����. ����� 50)

###############################################################   
class SuckerRodString(object):
    """���� ������ �������� �����"""
    items=[] # ������ ����� (������) [������1, ������2,...](��������� � �����)
    
    def __init__(self, PU, sections):
        """�����������.
        PU - ����
        sections=[(������, �������, �������, ��� ��������� �� ��������(���.),...]"""
        self.PU=PU
        for diametr, length, material, alfa in sections: # ��� ����� ������
            # �������� ������ (������ ���������� �����)
            self.items.append(SuckerRod(PU, diametr, length, material, alfa))
    
    def length(self,start=0,end=None):
        """�������. start, end - ������� ����"""
        return sum([r.length for r in self.items[start:end]])
    
    def xList(self):
        """������ ��������� ����� ����� ������"""
        X=[]
        x=0.0 # ���������� 
        for r in self.items:
            x=x+r.length # �������� ������� ���� ������
            X.append(x)
        return X
    
    def height(self):
        """������ �������� ������ ����� �� ���������"""
        return sum([r.length*math.cos(r.alfa) for r in self.items])
    
    def volume(self,start=0,end=None):
        return sum([r.volume() for r in self.items[start:end]])
    
    def stiffness(self):
        """�������� ���� ��������� ������ 1/k=1/k1+1/k2+..."""
        return 1/sum([1/r.stiffness() for r in self.items])
    
    def damping(self):
        """�������� ���� ����������� ������ 1/d=1/d1+1/d2+...
        [�������� � ������� �6 �177]"""
        return 1/sum([1/r.materialRodResistance() for r in self.items])
    
    def mass(self,start=0,end=None):
        return sum([r.mass() for r in self.items[start:end]])
    
    def massList(self):
        """������ ��� ������"""
        return [r.mass() for r in self.items]
    
    def fluidMassList(self):
        """������ ��� ����� ��� �����"""
        return [r.fluidMass() for r in self.items]
    
    def weight(self,start=0,end=None):
        return sum([r.weight() for r in self.items[start:end]])
    
    def wellAngleList(self):
        """������ ���� ����������� ����� ������"""
        return [x.alfa for x in self.items] 
    
#     def secAngle(self,L1,a1,L2,a2):
#         """������� ��� ��������� �� �������� ������ �� �������� ����������� �
#         ��������� L1,L2 �� ���������� �� �������� a1,a2 (���)
#         ʳ��� ������ ���������� ������� ������
#         """
#         x=L1*math.sin(a1)/2+L2*math.sin(a2)/2 # ������
#         y=L1*math.cos(a1)/2+L2*math.cos(a2)/2 # ������
#         a=math.atan(x/y) # ��� ��������� �� �������� ������
#         return a
        
    def secAngleList(self):
        """������� ������ ���� ������ an_i
        �� ������ ����������� ����� ������ alfa_i.
             an0          an1          an2
        ------------O------------O------------O
                  alfa0        alfa1        alfa2
        ��������� ��������� �� ����� ������ ������ ������.
        """
        angles=[self.items[0].alfa] # ��� ������� ������=���� ������� ������ �����������  
        aprev=self.items[0].alfa # ��� ����� ������ ����������� (������). �� ������� 0.
        for r in self.items[1:]: # ��� ����� ������ ��� �����
            a=(aprev+r.alfa)/2 # ����������� ��������� ���
            angles.append(a) # �������� � ������
            aprev=r.alfa # �����'����� �� ��������� ������
        return angles
      
    def angleTop(self):
        """���� �� �������� ������� � ������������"""
        AN=self.wellAngleList()
        ANS=self.secAngleList()
        return [ans-an for an,ans in zip(AN,ANS)]
        
    def angleBottom(self):
        """���� �� ������������ � ������� �������
        ������ �������� 0"""
        AN=self.wellAngleList()
        ANS=self.secAngleList()     
        return [an-ans for an,ans in zip(AN[:-1],ANS[1:])]+[0.0]  
    
    def hydrodynamicRodResistanceList(self):
        """������ ������������� ����� ������"""
        return [r.hydrodynamicRodResistance() for r in self.items]
    
    def weightForceList(self):
        """C����� ������� ��� �� ���� ������
        (� ����������� ���� �����������)"""
        return [r.weight()*math.cos(r.alfa) for r in self.items]

        
class Well(object):
    """���� ����������� � ���"""
    diameters=[60.3-2*5.0, 73.0-2*5.5, 73.0-2*7.0, 88.9-2*6.5, 114.3-2*7.0]
    diametr=diameters[4]/1000.0 #4 ������ ��� (��������� �� ����������)
    density=1000.0 #1000.0 ������� ����� [�����]
    viscosity=0.000002 #0.000001 # ���������� �'������ (����� �� �����!!!)
    bulkModulus=2.2e9 # ��'����� ������ ��������
    
    def __init__(self,PU):
        self.PU=PU
        self.length=PU.suckerRodString.length()
    def sectionalArea(self):
        """����� ����������� ������ ������ ���"""
        return math.pi*(self.diametr/2.0)**2
    def area(self):
        """����� ����������� ������ S=pi*r^2"""
        D=114.3/1000
        d=(114.3-2*7.0)/1000
        return math.pi*(D/2.0)**2-math.pi*(d/2.0)**2
    def stiffness(self):
        """��������� k=E*S/L"""
        E=2.1e11 # ������ ��������
        return E*self.area()/self.length
    def volume(self):
        return self.length*self.area()
    def mass(self):
        return self.volume()*8000.0
    def weight(self):
        """���� � ����"""
        Farh=self.density*9.81*self.volume() # ���� ��������
        w=self.mass()*9.81 # ���� � �����
        return w-Farh
    def materialRodResistance(self):
        """������������ ���������� ����� ��������, �*�/�"""
        k=self.stiffness() # ���������
        psi=0.2 # ���������� ���������� ��������
        omega=self.PU.pumpingUnit.omega #������ ������� �������� ����
        return psi*k/(2*math.pi*omega)
    
class Pump(object):
    """���� �������������� ������"""
    diameters=[27, 32, 38, 44, 50, 57, 57, 63, 70, 95]
    diametr=38/1000.0 #38 ������ ��������
    valveDiametr=0.025 # ������ ������ ���� �������
    openValveResistance=0.00001 # ��� ��������� ������� ???
    pistonGap=0.1e-3 #0.1e-3 ����� �� ��������� � ��������
    # ��� v - �������� ��������, vv - �������� � ��� �������  !!!
    
    def __init__(self,PU):
        self.PU=PU # ����
        
    def volumeFlowRate(self,v):
        """���������� ������ ������ (��'���� �������) Q=A*v"""
        return self.pistonArea()*v #self.PU.pumpingUnit.velAvg()
    
    def pistonArea(self):
        """����� ����������� ������ ��������"""
        return math.pi*(self.diametr/2.0)**2 # !!! ���������
    
    def valveDischargeCoefficient(self,Re):
        """���������� ������� �������. �������� �� ����� ����������. ó��������"""
        if Re<=225:
            mu=0.0846*Re**0.2872
        elif Re<=30000:
            mu=0.4
        elif Re<=300000:
            mu=0.0085*Re**0.3764
        else: mu=1.0
        return mu
    
    def valveLossCoefficient(self,Re):
        """���������� ������ ������� (�������� �����)"""
        return 1.0/self.valveDischargeCoefficient(Re)**2.0
    
    def valveV(self,v):
        """�������� ������ � ���, v - �������� ��������"""
        return self.volumeFlowRate(v)/(math.pi*(self.valveDiametr/2.0)**2)
        
    def valveRe(self,vv):
        """����� ����������, vv - �������� ������ � ���"""
        nu=self.PU.well.viscosity
        Re=abs(vv)*self.valveDiametr/nu
        return Re
    
    def dynamicPressure(self, vv):
        """��������� ����, vv - �������� ������ � ���"""
        ro=self.PU.well.density
        return vv*vv*ro/2
    
    def dPvalve(self,vv):
        """������� ����� �� ������. ó��������
        vv - �������� ������ � ���"""
        Re=self.valveRe(vv)
        dP=self.valveLossCoefficient(Re)*self.dynamicPressure(vv)
        return dP
    
    def valveHydroResistForce(self,v):
        """���� �������������� ����� ������� ������.
        ĳ� �� ��� ���� ���� � �����. ���������� ����� ����
        v - �������� ��������"""
        vv=self.valveV(v)
        return self.dPvalve(vv)*self.pistonArea() # ������� ������ ������
    
    def pistonFrictionForce(self):
        """���� ����� �������� �� ����� �������, �.
        ������� �������. �������� 1000 �
        ĳ� �� ��� ���� ���� � �����. ���������� ����� ����"""
        return 1.84*self.diametr/self.pistonGap-137 # �� ��� ���������� �����
 
class PumpingUnit():
    """���� ��������-��������"""
    stroke=3.0 #2.0 # ������� ���� ����� ������
    n=6.5 #10.0 # ������� ������� (�������� ����) �� �������
    fi=math.pi # ���������� ��� (������ � ������� ����� ��� � ��������, ���� ������ �����������)
    # ����� �������� (�������, ����, ������ ������� ������� �� z, ��� ��������� ����� ����� ���)
    # ���: ���8-3-4000 (������������)
    p1=dict(L=1.29,m=1982.0,Iz=635.0) # ��� ���������
    p2=dict(L=3.0,m=499.0,Iz=728.0) # ��� ������ � ��������
    p3=dict(L=2.0,m=643.0,Iz=443.0) # ���� ����� ���������
    p4=dict(L=2.29,m=911.0,Iz=1106.0) # ������ ����� ���������
    p5=dict(L=0.789,m=4*750.0,Iz=0.0) # ��������� (L - ����� ������ ���) ��������, �� ���� ����������� � ����� (Iz=0.0)
    # �����! �������� �� ������������� p5 �������� �������� balancing
    dx=1.345 # ������������� ������� �� ����� ����� ��������� � ����������� ���� ���������
    dy=3.012 # 3.0, 3.035 ? # ����������� ������� �� ����� (�������������)
    r=0.84 # ����� ��������� (�������� �� ����)
    # ������� �� ����. �� �.42 ������� ?
    
    def __init__(self):
        self.A=0.5*self.stroke # ��������
        self.ns=self.n/60.0 # �������, ��
        self.omega=2*math.pi*self.ns # ������ �������, ���/�
    
    def velAvg(self):
        """������� ��������, �/�"""
        return 2*self.stroke*self.n/60.0
    
    def motionX(self, t):
        """����� ����� ���� ����� ������.
        ������� ���������� X ����� ������ � ������ ���� t"""
        #!!!!!��������
        #v=self.A*self.omega*math.cos(self.omega*t+self.fi) # ��������
        x=self.A*math.sin(self.omega*t+self.fi) # ����������
        return x
     
class PU():
    """���� ������������ �������� ������� ���������"""
    suckerRodString=None
    well=None
    pump=None
    outPres=100000.0 # 100000.0 ���� �� ���� � ���
    outTubePres=500000.0 # 500000.0 ���� � ���������� �������
    heightDynamic = None # ������� ���������� ���� (���� None, �� �� �������������)
    # ��������� heightDynamic ������ ����� Fmax (���. ������� �-�� 9.90)
    timeEnd=20.0 #����� ����
    timeStep=0.1 # 0.01 ���� ����
    
    def __init__(self):
        """�����������"""
        h=1500.0 # ������� ������
        sections=[(0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0),
                  (0.019, h/8.0, "�����", 0.0)] # 0.2
        self.suckerRodString=SuckerRodString(PU=self,sections=sections) # ������ �����
        self.well=Well(PU=self) # �����������
        self.pump=Pump(PU=self) # �����
        self.pumpingUnit=PumpingUnit() # �������-��������
        if self.heightDynamic==None:
            self.heightDynamic=self.suckerRodString.height()-100.0 # -100 ��������
        r,m=self.balancing() # ����� � ���� ��������
        self.pumpingUnit.p5['L']=r
        self.pumpingUnit.p5['m']=m
              
    def fluidVolume(self):
        """��'�� ������ ����� ��� ��������� V=H*S"""
        return self.suckerRodString.height()*self.pump.pistonArea()
    
    def fluidWeight(self): # �� ����������� ������� �����
        """���� ����� ��� ��������� w=V*ro*g.
        ĳ� �� ������� �� ��� ���� �����. ���������� ����"""
        return self.fluidVolume()*self.well.density*9.81
    
    def outPresForce(self):
        """������������ �� ����� �� ����.
        ĳ� �� ������� �� ��� ���� �����. ���������� ����"""
        return self.pump.pistonArea()*self.outPres
    
    def outTubePresForce(self):
        """������������ �� ����� � ���������� �������.
        ĳ� �� ������� �� ��� ���� �����. ���������� �����"""
        return self.pump.pistonArea()*self.outTubePres
    
    def outTubeFluidWeight(self):
        """���� ����� � ���������� �������.
        ĳ� �� ������� �� ��� ���� �����. ���������� �����"""
        v=self.pump.pistonArea()*(self.suckerRodString.height()-self.heightDynamic) # ���
        return v*self.well.density*9.81
    
    def tubePressureFrictLossForse(self,v):
        """������������ �� ������ ����� �� ����� ������ �� ������ ���
        ĳ� �� ������� �� ��� ���� �����. ���������� ����. �� ��� ����� Maplesim"""
        v=self.pump.volumeFlowRate(v)/self.well.sectionalArea() # ������� �������� � ���
        D=self.well.diametr
        L=self.suckerRodString.length()
        nu=self.well.viscosity
        ro=self.well.density
        Re=abs(v)*D/nu # ����� ����������
        fL=64/Re # ��� ����������
        fT=0.316/math.pow(Re, 1.0/4) # ��� ������������� (������� �������)
        if Re<=2000: # ���� ��������� �����
            f=fL
        elif Re>=4000: # ���� ������������ �����
            f=fT
        else: # � ������ �������
            f=fL+(fT-fL)*(Re-2000)/2000
        dP=(f*L/D)*v*v*ro/2 # ������ �����
        return self.pump.pistonArea()*dP # ������ �������
        
    def polishedRodForce(self):
        """������� ������ (Fmin,Fmax) ������������ �� ����������� ���������
        ֳ ������� ���, ���� ��������� �����=������ ������
        ����� ������������� ������� ����������� �������"""
        F1=self.suckerRodString.weight() # ���� ����� � ����
        F2=self.fluidWeight() # ���� �����
        k=(self.pumpingUnit.stroke*self.pumpingUnit.n**2)/1790
        Fmin=F1*(1-k) # ������� ������
        Fmax=(F1+F2)*(1+k) # ������� �������
        return (Fmin,Fmax)
    
    def sigmaPr(self,d):
        """��������� ���������� (������)"""
        A=math.pi*(d/2.0)**2 # ����� ������ ������� ������
        Fmin,Fmax=self.polishedRodForce()
        return (Fmax-0.5625*Fmin)/A
    
    def polishedRodForce2(self):
        """������� ������ (Fmin,Fmax) ������������ �� ���� ����������� ���������"""
        F1=self.suckerRodString.weight() # ���� ����� � ����
        F2=self.fluidWeight() # ���� �����
        Fmin=F1
        Fmax=F1+F2
        return (Fmin,Fmax)
    
    def balancing(self,Pmin=None,Pmax=None):
        """��������� ���������� �������������� [�������]. ��� �������� ���� �� � ���
        ������� ����� � ������� ���� ��������"""
        if Pmin==None and Pmax==None: # ���� �� ������
            Pmin,Pmax=self.polishedRodForce() # ����������� �� ����������� ���������
        Mp=750.0 # ���� 1 ���������
        A0=52000.0 # ������ ��� ���� �������������� ������ �������
        S=self.pumpingUnit.stroke
        K_list=[2,4,6,8] # ������ ������� ��������
        A=Pmin*S/2+Pmax*S/2 # ������ ���, �� ��������������
        R_list=[(A-A0)/(2*9.81*Mp*k) for k in K_list] # ������ ������ ��������
        #print "Rp=", R_list
        Rmax=self.pumpingUnit.p1['L']-0.2 # ������������ ����� ��������� ?��������
        for i,r in enumerate(R_list):
            if r<Rmax: break # ��������� ����� � ��������� ������� ��������
        return r, K_list[i]*Mp
    
    def info(self): # ����� ����������
        s=""
        s+="���� ����� � ���� %d �\n"%self.suckerRodString.weight()
        s+="���� ����� %d �\n"%self.fluidWeight()
        s+="���� ��������� Fmin=%d Fmax=%d\n"%self.polishedRodForce2()
        s+="��������� Fmin=%d Fmax=%d\n"%self.polishedRodForce()
        s+="��������� ���������� %d\n"%self.sigmaPr(d=self.suckerRodString.items[0].diametr)
        s+="������ �������� %f\n"%self.pump.diametr
        s+="������� ������ %d\n"%self.suckerRodString.length()
        s+="��������� ������ %d �/�\n"%self.suckerRodString.stiffness()
        s+="���� ������ %d ��\n"%self.suckerRodString.mass()
        n0=math.sqrt(self.suckerRodString.stiffness()/self.suckerRodString.mass()) # ������ ������� �������
        s+="������ ������� ������ %f ���/�\n"%n0
        s+="�������� ����� %d �\n"%self.heightDynamic
        s+="���� ����� � ���� ������� %f �\n"%self.outTubeFluidWeight()
        s+="������������ ���������� ����� �������� %f \n"%self.suckerRodString.damping()
        s+="������������ ���������� ����� �������� ������ 1 %f \n"%self.suckerRodString.items[0].materialRodResistance()
        # � ���������� ������� ������ ����� �� 0 (�� � ���������)
        s+="��������� ������ 1 %f \n"%self.suckerRodString.items[0].stiffness()
        s+="���� ������ 1 %f \n"%self.suckerRodString.items[0].mass()
        
        s+="���� ��� %f \n"%self.well.mass()
        s+="��������� ��� %f \n"%self.well.stiffness()
        s+="���. ����. ��. �����. ��� %f \n"%self.well.materialRodResistance()
        print s 

if __name__ == '__main__':
    pu=PU()
    pu.info()
#     A=math.pi*(0.019/2.0)**2 # ����� ������ ������� ������
#     Fmin,Fmax=7000,30000
#     print (Fmax-0.5625*Fmin)/A