# -*- coding: CP1251 -*-

# -��������� ��������� ������� ����� �������� �� ����������
# -���������� ������ ������ K, D
# -��������� ���������� ���������
# -��������� ���������� ��������� ���������� ����������� ����� ����������
# -���������� ������� �� ������ a.dt=0.01. ����� a.dt=timeStep

import math, os
import PUmodel
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='Arial') # ��� �������� ������

class Automaton(object):
    """����������� ������� ��� ����������� ���� ����� �������� �������"""
    def __init__(self):
        """����������� ��������"""
        self.x=None # ����������
        self.x0=None # ��������� ����������
        self.prevx=None # ��������� ����������
        self.bc=None # �������� ����� (0 - ���������)
        self.f=0.0 # ������� ����
        self.fs=0.0 # ���� ������ ������
        self.left=None # ���� ���� (������)
        self.right=None # ������ ���� (�����)
        self.delta=None # ��������� ������� �� ������� ����� (rSection.length)
        self.k=None # ���������� �������� ������ ������ (rSection.stiffness())
        self.c=0.0 # ���������� ����� �������� ������ ������
        self.m=1.0 # ���� ������ ������
        self.a=0.0 # �����������
        self.v=0.0 # ��������
        self.prevv=0.0 # ��������� ��������
        self.eps=0.0 # ���'����
        self.dx=1.0 # ���� ����������
        self.area=None # ����� ����������� ������ ������� (�������) ������ (rSection.area())
        self.domain=None # ��'��� ����� CA_model
        self.rSection=None # ������ ������
        self.alfa=0.0 # ��� ��������� ����������� �� �������� ��� ��������

    def rule(self):
        """������� �������� ��������"""
        
        if self.bc!=None: # ���� ������ �������� �����
            self.x=self.bc # ������ x
        
        if self.domain.t!=0.0: # ���� ��� �� 0
            self.v=self.velocity() # ����������� ��������
            self.a=self.acceleration() # ����������� �����������
        else: self.v=self.a=0.0
        
        if self.bc!=None: # ���� ������ �������� �����
            return # ����� (�� ������� ��������� x)
        
        Fk=self.rightForce()-self.leftForce() # ���� ��������
        Fv=self.rightForce2()-self.leftForce2() # ���� �����������
        Fa=self.dynamicForce() # ���� �������
        # ������ ����
        F=self.fs+self.f+self.pistonForce()+self.frictionForce()+self.hydrodynamicResistanceForce() 
        
        # �������� ��� ��� ��� ��������: m*a+c*v+k*u=F
        # ���� ����>0 �� ���� ���������� ����
        # �������� �� x:
        # [0]---[1]---[2]--->+x,+v,+a,+f,������,����� 
        # ��������� ���'����
        self.eps=Fa+Fv+Fk-F
        #print "eps=",self.eps
        
        # ���� ���'���� ����������� �� 0
        if self.eps>0:
            self.x=self.x-self.dx # �������� x �� ���� ��������
        if self.eps<0:
            self.x=self.x+self.dx # �������� x �� ���� ��������
        
    def velocity(self):
        return (self.x-self.prevx)/self.domain.dt # ��������
        
    def acceleration(self):
        return (self.v-self.prevv)/self.domain.dt # �����������
               
    def leftForce(self):
        """������� �������� ���� ������� ���� (������)"""
        if self.left: # ���� � ���� ����
            dLeftBegin=self.left.delta # ��������� ������� ������� ����
            dLeft=self.x-self.left.x # ������� ������� ����
            fLeft=(dLeftBegin-dLeft)*self.left.k # ���� ������� ����
        else:
            fLeft=0.0
        return fLeft
        
    def leftForce2(self):
        """������� �������� ���� �������� ���� (������)"""
        if self.left: # ���� � ���� ����
            vrel=self.v-self.left.v # ������� ��������
            fLeft=-self.left.rSection.materialRodResistance(None)*vrel
        else:
            fLeft=0.0
        return fLeft
    
    def rightForce(self):
        """������� �������� ���� ������� ������ (�����)"""
        if self.right: # ���� � ������ ����
            dRightBegin=self.delta # ��������� ������� ������� ������
            dRight=self.right.x-self.x # ������� ������� ������ 
            fRight=(dRightBegin-dRight)*self.k # ���� ������� ������
        else:
            fRight=0.0
        return fRight
           
    def rightForce2(self):
        """������� �������� ���� �������� ������ (�����)"""
        if self.right: # ���� � ������ ����
            vrel=self.right.v-self.v # ������� ��������
            fRight=-self.rSection.materialRodResistance(None)*vrel
        else:
            fRight=0.0
        return fRight
        
    def dynamicForce(self):
        """���� ������� ����� F=m*a"""
        return self.m*self.a
    
    def dynamicFluidForce(self):
        """���� ������� �����. ĳ� �� �������"""
        k=self.domain.PU.pump.pistonArea()/self.domain.PU.well.sectionalArea()
        a=self.a*k # ����������� ����� ����� ����������� �������� � k ���
        return -a*sum(self.domain.PU.suckerRodString.fluidMassList())
            
    def normalForce(self):
        """���� ���������� ��� � ���� �� ����, ����� � ������� ������"""
        F1=self.f*math.tan(self.alfa) # ��������� ���� �� ���� (f - ������ ���� �� ����)
        F2=F3=0.0
        if self.left:
            F2=self.leftForce()*math.sin(self.an1) # ��������� ���� �� ����� ������
        if self.right:
            F3=self.rightForce()*math.sin(self.an2) # ��������� ���� �� ������� ������
        return F1+abs(F2+F3) # ���� ��������� ������� �� ������ ��� � ������
        
    def frictionForce(self):
        """���� �����"""
        # ���� ����� �� ��������
        #Ft=C*(Psh*math.sin(a1)-F12*math.sin(a12-a1)) # ???
        #Ft=C*(Psh*math.sin(a1)+2*F12*math.tan(a2/2-a1/2)) # a12-a1=a2/2-a1/2
        #Ft=C*(Psh*math.sin(a1)+F12*(a2-a1)) # ������� �������
        
        v=self.v #��������
        Fn=self.normalForce() #��������� ���� # ��� �������� �������� ������ ��� 1000.0
        fc=0.16 #0.25 # ���������� ����� �������� (�����-����� � �����������)
        #F=-fc*abs(Fn)*math.tanh(v/0.01) # �������� ������
        
        Fc=abs(Fn)*fc # Coulomb friction force
        fmax=0.2 # ���������� ����� �������� (�����-����� � ����� �����������)
        Fmax=abs(Fn)*fmax # Maximum sliding friction force
        vs=0.1 # Stribeck sliding velocity coefficient
        n=1.0 # Decay exponent
        # ������� ���� �����
        F=-math.tanh(v/0.01)*(Fc+(Fmax-Fc)*math.exp(-(abs(v)/vs)**n))
        return F
    
    def hydrodynamicResistanceForce(self):
        """C��� �������������� ����� ������"""
        v=self.v # ��������
        if v<=0: # ���� ��� �����
            k=self.domain.PU.pump.pistonArea()/self.domain.PU.well.sectionalArea()
            vu=v-v*k # �������� ����� ������� �����
            F=-self.c*vu*10.0 #�����!!! ��������
            # F=0.0 # ���� �������� ����� � ����� ���
        else: # ���� ��� ����
            F=-self.c*v*10.0  #�����!!! ��������
        return F*math.tanh(abs(v)/0.01) # ��������� ��� 0
    
    def pistonForce(self):
        """���� ������� ��� �� ������� ������"""
        v=self.v # ��������
        F=0.0 # ��� ��� ����� ��� ���������� 
        if self.right!=None: return F # ���� �� ������� ����� (�����)
        
        if v<0: # ���� ��� �������� �����
            # ��������� ����:
            F+=self.domain.PU.fluidWeight()
            F+=self.domain.PU.outPresForce()
            F+=abs(self.domain.PU.tubePressureFrictLossForse(v)) # ��� abs ��� ���� �������
            F+=self.dynamicFluidForce()
            # ��������� �����:
            F+=-self.domain.PU.outTubePresForce()
            F+=-self.domain.PU.outTubeFluidWeight()
        if v!=0:
            # ������� ���� ���� �� ��� ���� (!) ����� � ����
            # � ��������� ����� ���� ��������
            F+=-math.copysign(1, v)*abs(self.domain.PU.pump.valveHydroResistForce(v))       
            F+=-math.copysign(1, v)*self.domain.PU.pump.pistonFrictionForce() # ������ ����� ��������
        
        F=F*math.tanh(abs(v)/0.01) # �-��� ������������ ��� 0 (F=max, ���� v=0.01) 
        #��� F=F*(1.0-math.exp(v/0.01))
        return F
         
    def rightStress(self):
        """������� �������� ���������� � ������ ������ (�����)"""
        if self.right:
            return self.rightForce()/self.area
        else:
            return 0.0
        
################################################################
class CA_model(object):
    """������� �������� (������� ��������)"""
    def __init__(self,PU):
        """������� ������� �������� �� ��'����� ����� PU (����)
        ����� ������� ��������:
                ������1                      ������2
        O--------------------O-----------------------------------O
        m0=0      d0         m1                d1                m2 
        c0=0      k0         c1                k1                c2
        f0        a0         f1                a1                f2
       alfa0=0   ans0       alfa1             ans1              alfa2                                       
        ������ ������� (0) ���������.
        """
        
        self.PU=PU # ����
        rs=PU.suckerRodString.items # ������ ������
        n=len(rs)+1 # ������� ����� (��������)
        # ������ ������� (0) ���������.
        X=[0.0]+PU.suckerRodString.xList() # c����� ��������� �����
        M=[0.0]+PU.suckerRodString.massList() # ���� �����
        Fs=[0.0]+PU.suckerRodString.weightForceList() # ������ ��� �� ���� ������
        C=[0.0]+PU.suckerRodString.hydrodynamicRodResistanceList() # ����������� ����� ������

        AN=[0.0]+PU.suckerRodString.wellAngleList()
        #ANS=PU.suckerRodString.secAngleList()
        AN1=[None]+pu.suckerRodString.angleTop()
        AN2=[0.0]+pu.suckerRodString.angleBottom()
        
        self.dt=self.PU.timeStep # ���� ����
        self.T=self.timeList(PU.timeEnd, PU.timeStep) # ������ ������� ���� [t0,t1,t2,...]
        self.t=self.T[0] # ������� �������� ���� (0.0)
        
        BCdict={} #������� ��������� ���� (���. self.BCList)
        # ����� ���� ����� ������
        BCdict[0]=[PU.pumpingUnit.motionX(t) for t in self.T]
        self.BC=CA_model.BCList(BCdict, n, self.T, default=None) # ������ ����� ��������� ���� [[bc0,bc1,bc2,...],...]
                
        F=[0.0]*n # ������ ����� ������� ��� (���� �������)
        Fdict=dict(zip(range(n),F)) # ������� ��� (���. self.BCList)
        self.F=CA_model.BCList(Fdict, n, self.T, default=0.0) # ������ ����� ������� ��� [[f0,f1,f2,...],...]

        self.p=[Automaton() for i in range(n)] # ������ ��������
        
        # ��� ������� �������� ��� ���������� ������ ����������
        for i in range(n-1):
            self.p[i].rSection=rs[i]
            self.p[i].k=rs[i].stiffness()
            self.p[i].delta=rs[i].length
            self.p[i].area=rs[i].area()
        
        # ������ ���� ���������� ��� ��������
        for i in range(n):
            self.p[i].bc=self.BC[0][i]
            self.p[i].f=self.F[0][i]
            self.p[i].fs=Fs[i]
            self.p[i].m=M[i]
            self.p[i].c=C[i]
            self.p[i].alfa=AN[i]
            self.p[i].x=X[i]
            self.p[i].x0=X[i]
            self.p[i].domain=self
            self.p[i].an1=AN1[i]
            self.p[i].an2=AN2[i]
            
        for i in range(n): # ������ �����
            if i!=0: # ���� �� ������
                self.p[i].left=self.p[i-1] # ������ ����� �����
            if i!=len(self.p)-1: # ���� �� �������
                self.p[i].right=self.p[i+1] # ������ ������� �����
        
        # ������� ������ ����������
        # ���������, ������ ����� ������� x [[x0,x1,x2,...],...]
        self.history={'T':self.T,'X':[],'V':[],'A':[],'Fl':[],'Fr':[],'Sr':[]}
                      
    def run(self):
        """�������� ������� �������� � ���������� ������ ����"""
        eps=0.01*len(self.p) # ������� ����������
        e=eps+1 # ���� ���'����
        eprev=eps+1 # ��������� ���� ���'����
        # ���������� ���� ���������� (�������� �� delta)
        dx=min([a.delta for a in self.p if a.delta!=None])/10.0 # !!! �������� � ����� 10-100
        niter=0 # ������� �������� (��� ���������� ���������)
        while e>eps: # ���� ���� ���'���� > eps
            for a in self.p: # ��� ������� ��������
                a.dx=dx # ������ ���� ����������
                a.rule() # ����������� ������� �������� ��������
            e=sum([abs(a.eps) for a in self.p]) # ���� ���'����
            #print e
            if e>=eprev: # ���� ���� ���'���� >= ����������
                dx=dx*0.935 # �������� ���� !!! �������� � ����� 0.1-0.99
                #print "dx=",dx
            eprev=e # �����'����� ���� ���'����
            niter=niter+1
        #print niter
        
#         a=self.p[1]
#         print 'lf=',a.leftForce()
#         print 'rf=',a.rightForce()
#         print 'lf2=',a.leftForce2()
#         print 'lvrel=',a.v-a.left.v
#         print 'rf2=',a.rightForce2()
#         print 'rvrel=',a.right.v-a.v
            
    def runDynamic(self):
        """�������� �������� ������� ��������"""
        # �������� �������� ���������, ���������� � ����������
        for a in self.p: # ��� ������� ��������
            a.prevx=a.x # ��������� ����������
            a.prevv=0.0 # ��������� ��������
        for i,t in enumerate(self.T): # ��� ������� �������� ����
            self.t=t # ������� �������� ����
            print "solve t=",t
            for a,f,bc in zip(self.p, self.F[i], self.BC[i]): # ��� ������� ��������
                a.f=f # ������ ���� � ������ ���� t
                a.bc=bc # ������ �������� ����� � ������ ���� t
            self.run() # ������ ������������ ����
            for a in self.p: # ��� ������� ��������
                a.prevx=a.x # ���'������� ������ x �� v
                a.prevv=a.v
            self.appendHistory()# �����'����� �����
        self.writeHistory()
            
    def appendHistory(self):
        """���� � ������ ����� ������ ������� ��� ��������� ����"""
        self.history['X'].append([a.x for a in self.p])
        self.history['V'].append([a.v for a in self.p])
        self.history['A'].append([a.a for a in self.p])
        self.history['Fl'].append([a.leftForce() for a in self.p])
        self.history['Fr'].append([a.rightForce() for a in self.p])
        self.history['Sr'].append([a.rightStress() for a in self.p])
        
    def draw(self):
        """���� ��������� �������� �������������� (� ���������� �����)"""
        s=""
        prev=0
        for a in self.p:
            s+="-"*int(round(a.x-prev)-1)+"O"
            prev=a.x # ���������� ������������ ��������
        print s
        #time.sleep(1)
        
    def drawPlot(self, itemIndexes, key, x0scale=1, toFile=True):
        """���� ��������� �������� x[i] �� t, �� i - ������ ��������
        itemIndexes - ������ ������� ��������
        key - ���� ����� x
        x0scale - ������� �������� �� ������� ������� ������ �� �� X (1...0):
                  1 - ��������, 0 - ������ ��� 0.
                  ����� ��� 0-�� ����� �� ��������.
                  ��������������� ��� ���������� ���������� ��������� �����,
                  ���� ���� ��� ������� �������� �� �������"""
        plt.figure()
        plt.gca().invert_yaxis() # �������� ��� y
        T=self.history['T'] # ������ ������� ����
        x00=self.history[key][0][0] # x 0-�� ����� ��� t=0
        for i in itemIndexes: # ��� ������� �������
            x0i=self.history[key][0][i] # x i-�� ����� ��� t=0
            x0=(1-x0scale)*(x0i-x00) # �������� ������� �� �� X
            X=[x[i]-x0 for x in self.history[key]] # ������ ����� ����� "i"
            plt.plot(T, X, 'k-') # (T, X, 'r-', T, X,'ro')
        plt.grid(True)
        plt.xlabel(u't, c') #������ �� x
        plt.ylabel(key) #������ �� y
        if toFile:
            fileName='CAmodelPlot{0}{1}.png'.format(itemIndexes,key)
            plt.savefig(fileName)
        else:
            plt.show()
        
    def drawDynamometerCard(self, itemIndexes, keyX, keyY, signY=-1, toFile=True, tstart=None, tend=None):    
        """���� ������������ - ��������� Y(X), ���������, ����(����������)"""
        plt.figure()
        plt.gca().invert_xaxis() # �������� ��� x
        plt.grid(True)
        plt.xlabel(u"S, �") #������ �� x
        plt.ylabel(u"F, �") #������ �� y        
        for i in itemIndexes: # ��� ������� �����
            # ������� ����������� � �������� ���� � �����
            tstarti=self.history['T'].index(tstart) if tstart in self.T else None 
            tendi=self.history['T'].index(tend)+1 if tend in self.T else None
            X=[x[i] for x in self.history[keyX][tstarti:tendi]] # ������ ����� X ����� "i"
            x0=self.history[keyX][0][i] # ��������� ��������
            X=[x-x0 for x in X] # ������� ��������� ����������
            Y=[signY*y[i] for y in self.history[keyY][tstarti:tendi]] # ������ ����� Y ����� "i"
            plt.plot(X, Y, 'k-')
        if toFile:
            fileName='CAmodelDynCard{0}{1}.png'.format(keyX,keyY)
            plt.savefig(fileName)
        else:
            plt.show()
        
    def writeHistory(self):
        """������ ������ � ������� ����"""
        import pickle
        f = open('CAmodelHistory.pkl', 'wb')
        pickle.dump(self.history, f) # �������� ������
        f.close()
        
    def readHistory(self):
        """���� ������ � �������� �����"""
        import pickle
        fileName='CAmodelHistory.pkl'
        if not os.path.isfile(fileName): # ���� ����� ����
            print "No file "+fileName
            return False # �����, ���������� False
        f = open(fileName, 'rb')
        self.history=pickle.load(f) # ����������� ������
        f.close()
        print "Warning! Old results"
        return True
       
    def writeCSV(self,itemIndexes,key):
        """������ �������� ������ � ���� CSV.
        itemIndexes - ������ ������� ��������
        key - ���� �����"""
        import csv
        csv_file=open("RodHistory_"+key+".csv", "wb")
        writer = csv.writer(csv_file,delimiter = ';')
        for time,data in zip(self.T,self.history[key]):
            X=[data[i] for i in itemIndexes] # ���� ��� �� ��������� ���������
            writer.writerow([time]+X) # �������� �����
        csv_file.close()

    def info(self):
        """�������� �������� ��� �������� ��������"""
        attrs=["rSection","k","delta","area","bc","x","prevx","v","prevv",\
               "a","m","alfa","an1","an2","c","eps","f"]
        for attr in attrs:
            print attr+"=", [a.__getattribute__(attr) for a in self.p]
        print "leftForce=", [a.leftForce() for a in self.p]
        print "rightForce=", [a.rightForce() for a in self.p]
         
                    
    @staticmethod
    def timeList(timeEnd,timeStep):
        """������� ������ ������� ����"""
        nTimeSteps=int(timeEnd/timeStep)
        return [x*timeStep for x in range(nTimeSteps+1)]  

    @staticmethod
    def BCList(BC, n, timeList, default=None):
        """
        ������� ������ ��������� ���� (��� ���) ��� �������� 1,2,... � ������:
        [[bc0,bc1,...],...],
        �� ������ ������ ������� ������� �������� ����.
        
        BC - ������� ��������� ����, �� ����� - ������� ��, � ���������� ������ ����:
        1.ĳ��� (������� ��������): 0.0
        2.�������� ��� ���-��������: {0.0:1.0,0.05:2.0}
        3.������ ����� �������: [1.0,2.0,3.0,4.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]
             ���: [math.sin(x) for x in range(nt)]
        n - ������� ��������
        timeList - ������ ������� ����
        default - �������� ������ �� �������������"""
        
        nt=len(timeList) # ������� ������� ����
        res=[] # ���������
        for i in range(nt):
            res.append([default]*n) #  ���������� res ������(!) ��������
                            
        for k in BC: # ��� ������� �������� �����
            if type(BC[k]).__name__=='dict': # ���� ������ ��� ���������� ��������� ����
                i=0
                prev=None # ���������
                for t in res: # ��� ������� �������� ����
                    if BC[k].has_key(timeList[i]):             
                        t[k]=BC[k][timeList[i]]
                        prev=t[k]
                    else:
                        t[k]=prev
                    i=i+1
            elif type(BC[k]).__name__=='list': # ���� ������ ��� ������� ��������� ����
                i=0
                for t in res: # ��� ������� �������� ����
                    t[k]=BC[k][i]
                    i=i+1  
            else: # ���� ���������
                for t in res: # ��� ������� �������� ����
                    t[k]=BC[k]
        #print res            
        return res 
    
    def resultsCAmodel(self):
        """�������� ���������� ��������� �����"""
        if self.history['X']==[]: # ���� ���� ����������
            oldHistExist=self.readHistory() # ��������� ���� ����������
            if not oldHistExist: return
            
        tend=self.T[-1]
        T=1.0/self.PU.pumpingUnit.ns # �����
        tstart=tend-T # ����� ������� �����
        # ������ ��������� �������� �� tstart � ������ self.T
        dt=[abs(t-tstart) for t in self.T]
        tstart=self.T[dt.index(min(dt))]
        self.drawDynamometerCard([0,1,2,3,4,5,6,7], "X", "Fr", tstart=tstart, tend=tend)
        self.drawPlot([0,1,2,3,4,5,6,7,8],'X', x0scale=0.001) 
        self.drawPlot([0,1,2,3,4,5,6,7],'Sr')
        
if __name__ == '__main__':
    pu=PUmodel.PU()
    model=CA_model(pu)  
    #model.info()
    #model.run() # �������
    #model.info()
    
    model.runDynamic() # �������
    #model.info()
    model.resultsCAmodel()
    