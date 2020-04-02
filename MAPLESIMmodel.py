# -*- coding: cp1251 -*-
import math, os
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='Arial') # для підтримки Юнікоду
import PUmodel, maplepy

# Бажано задавати SD1_s_rel0. Визначити його можна шляхом моделювання непрацюючої установки (задавши S1_T0>0).

class Maplesim_model():
    """Модель СШНУ у Maplesim
    Потребує файлу моделі PUmodel.msim"""
    pu=None # об'єкт класу PUmodel.PU
    ms=None # об'єкт класу maplepy.MapleInterface4Maplesim
    
    def prepareParams(self):
        """Підготовлює парметри Maplesim моделі"""
        pu=self.pu # СШНУ
        rs=pu.suckerRodString.items # список секцій
        n=len(rs) # кількість секцій
        
        M=pu.suckerRodString.massList() # маса секцій
        M2=pu.suckerRodString.fluidMassList() # маса рідини біля секцій
        C2=pu.suckerRodString.hydrodynamicRodResistanceList()
        W=pu.suckerRodString.weightForceList()
        C=[r.stiffness() for r in rs]
        D=[r.materialRodResistance() for r in rs]
        v=pu.pumpingUnit.velAvg()
        
        AN=pu.suckerRodString.wellAngleList() # список кутів alfa (довжиною n)
        ANS=pu.suckerRodString.secAngleList() # кути відхилення секцій від вертикалі
        AN1=pu.suckerRodString.angleTop()
        AN2=pu.suckerRodString.angleBottom()
                        
#         for p in 'M M2 C2 W C D v AN ANS AN1 AN2'.split():
#             print p, eval(p) # надрукувати
        
        # підготувати параметри для моделі Maplesim
        modelParams={'S1_amplitude':pu.pumpingUnit.A, # параметри верстата-гойдалки
                     'S1_freqHz':pu.pumpingUnit.ns,
                     'S1_phase':pu.pumpingUnit.fi,
                     'HC1_A':pu.pump.pistonArea(), # параметри гідравлічної частини
                     'CV1_Ropen':pu.pump.valveLossCoefficient(v),
                     'CV2_Ropen':pu.pump.valveLossCoefficient(v),
                     'CP1_D':pu.well.diametr, # ?
                     'CP1_L':pu.suckerRodString.length(),
                     #Замість VerticalPipe1 можна використовувати FP2_P
                     #'VerticalPipe1_z':pu.suckerRodString.height(),
                     'AP1_P':0.0, # тиск на вході
                     'AP2_P':pu.outPres, # тиск на гирлі
                     # гідростатичний тиск, який діє на плунжер під час ходу вверх
                     # враховується також гідростат. тиск в затр. просторі, який діє з другої сторони плунжера
                     'FP2_P':pu.fluidWeight()/pu.pump.pistonArea()-pu.outTubeFluidWeight()/pu.pump.pistonArea()-pu.outTubePres,
                     'rhoFluid':pu.well.density,
                     'nuFluid':pu.well.viscosity,
                     # середня площа січення труби
                     'FI1_A':(sum(M2)/pu.well.density)/pu.suckerRodString.length(),
                     # ввести FI1_L=0, якщо не моделюється інерція рідини
                     'FI1_L':pu.suckerRodString.length(),
                     'TF1_fs':1.25*pu.pump.pistonFrictionForce(),
                     'TF1_fc':pu.pump.pistonFrictionForce()
                     }
        #W[0]+=15000.0 # !!! сила тертя в ущільненні        
        for i in range(n): # для кожної секції ni (i=0, 1, 2...)
            modelParams['n'+str(i)+'_M1_m']=M[i]
            # !!! увага тут і нижче коефіцієнти демпфування збільшено
            modelParams['n'+str(i)+'_G2_k']=10*C2[i] # узгодити з TF1_d!!! гасить вібрації під час ходу вниз
            modelParams['n'+str(i)+'_CF1_f_constant']=-W[i] # увага знак
            modelParams['n'+str(i)+'_TF1_fs']=0.2*W[i]*math.tan(AN[i]) # сила тертя від ваги max
            modelParams['n'+str(i)+'_TF1_fc']=0.16*W[i]*math.tan(AN[i]) # сила тертя від ваги
            modelParams['n'+str(i)+'_TF1_d']=10*C2[i] # коеф. в'язкого тертя (узгодити з G2_k!!!) гасить вібрації
            modelParams['n'+str(i)+'_G1_k']=0.16*math.sin(AN1[i]) # множник сили тертя від верхнього натягу
            modelParams['n'+str(i)+'_G4_k']=0.16*math.sin(AN2[i]) # множник сили тертя від нижнього натягу
            modelParams['n'+str(i)+'_SD1_c']=C[i]
            modelParams['n'+str(i)+'_SD1_d']=D[i] # повинен бути більшим, оскільки існує ще конструкційне демпфування
        
#         for p in sorted(modelParams.keys()):
#             print p, modelParams[p]
            
        self.modelParams=modelParams
            
    def createMaplesimModel(self, execute=False, resFileName=None):
        """Створює Maplesim модель. Виконує її, якщо execute=True"""       
        ms=maplepy.MapleInterface4Maplesim()
        ms.path=os.getcwd().replace('\\', '/')
        ms.filenameMaplesim="PUmodel.msim"
        ms.paramDictMaplesim=self.modelParams # {} якщо не передавати нічого
        #print ms.getCode()
        if execute: ms.execute=True
        if resFileName: ms.resultCSVfile=resFileName
        ms.runMaple()
        self.ms=ms
            
    def drawDynamometerCard(self, S0=1.0):
        """Створює файл-рисунок з динамограмою
        Першою пробою в PUmodel.msim має бути проба Length і Force
        на місці динамометра (верхня штанга)
        S0 - площа поперечного січення верхньої штанги (якщо S0=1, виводить силу)"""    
        results=self.ms.readCSVfile(5) # список результатів
        
        # рисує тільки останній період
        T=1.0/self.modelParams['S1_freqHz'] # період
        tend=results[-1][0] # кінцевий час
        tstart=tend-T # початок останнього періоду
        
        f=[] # сила або напруження
        s=[] # переміщення
        f1=[]
        s1=[]
        for res in results: # для кожного часу
            if res[0]>=tstart: # починати з часу tstart
                f.append(res[1]/S0) # !!!напруження
                s.append(res[2]) # переміщення
                f1.append(res[3])
                s1.append(res[4])
        plt.plot(s, f, 'k-')
        plt.plot(s1, f1, 'b-')
        plt.grid(True)
        plt.xlabel(u"S, м") #надпис осі x
        if S0==1.0:
            plt.ylabel(u"F, Н") #надпис осі y
        else:
            plt.ylabel("s, Pa") #надпис осі y
        plt.savefig("MaplesimDynCard.png")
        return min(f),max(f),min(s1),max(s1),min(f1),max(f1)
        
    def drawDynamometerCardMulti(self, filenames=[],S0L=[],freqL=[],styleL=[]):
        """Створює файл-рисунок з динамограмою для багатьох файлів-результатів
        див.документацію drawDynamometerCard"""    
        ms=maplepy.MapleInterface4Maplesim()
        ms.path=os.getcwd().replace('\\', '/')
        for name,S0,freq,style in zip(filenames,S0L,freqL,styleL):
            ms.resultCSVfile=name
            results=ms.readCSVfile(3) # список результатів
            
            # рисує тільки останній період
            T=1.0/freq # період
            tend=results[-1][0] # кінцевий час
            tstart=tend-T # початок останнього періоду
            
            f=[] # сила або напруження
            s=[] # переміщення
            
            delta=0
            #if name=='13.csv': delta=0 # !!!поправка практичної динамограми, якщо потрібно 
            
            for res in results: # для кожного часу
                if res[0]>=tstart: # починати з часу tstart
                    f.append(res[1]/S0-delta) # !!!напруження
                    s.append(res[2]) # переміщення
            plt.plot(s, f, style)
        plt.grid(True)
        plt.xlabel(u"S, м") #надпис осі x
        plt.ylabel(u"F, Н") #надпис осі y
        plt.savefig("MaplesimDynCard.png")
                    
    def optimize(self):
        """Розраховує модель для різних значень її параметрів"""
        parList=[i/10.0 for i in [1,2,3,4,5,6,7,8,9,10]] # список значень параметра
        fl=open("ampl-freq.csv", "w")
        self.pu=PUmodel.PU()
        for par in parList: # для кожного значення
            self.pu.pumpingUnit.ns=par # змінити значення параметра
            self.prepareParams()
            #Створює файли результатів x_MaplesimResults.csv
            #self.createMaplesimModel(True, str(par)+'_MaplesimResults.csv')
            self.createMaplesimModel(True)
            fl.write("%f;%f;%f;%f;%f;%f\n"%self.drawDynamometerCard())
            fl.flush()
        fl.close() 

class Maplesim_model2(Maplesim_model): # успадковує Maplesim_model
    """Спрощена модель СШНУ у Maplesim
    Потребує файлу моделі PUmodel.msim"""
    
    def prepareParams(self):
        """Підготовлює параметри Maplesim моделі"""
        pu=self.pu # СШНУ
        
        M=pu.suckerRodString.mass() # маса колони
        W=pu.suckerRodString.weight() # вага колони
        WF=pu.fluidWeight() # вага рідини
        C=pu.suckerRodString.stiffness() # жорсткість колони
        D=pu.suckerRodString.damping()# еквівалентний коефіцієнт опору
        #D=50*D # увага, збільшено в 50! (або прийняти psi=0.5)
        #D=psi*math.sqrt(C*M)/(2*math.pi) # ?
        
        # підготувати параметри для моделі Maplesim
        modelParams={'S1_amplitude':pu.pumpingUnit.A, # параметри верстата-гойдалки
                     'S1_freqHz':pu.pumpingUnit.ns,
                     'S1_phase':pu.pumpingUnit.fi,
                     'SD1_c':C, # параметри гідравлічної частини
                     'SD1_d':D,
                     'M1_m':M,
                     'CF1_f_constant':-W,
                     'G1_k':-WF
                     }
        
        for p in sorted(modelParams.keys()):
            print p, modelParams[p]
            
        self.modelParams=modelParams

class Maplesim_model3(Maplesim_model): # успадковує Maplesim_model
    """Модель верстата-гойдалки у Maplesim
    Потребує файлу моделі PUmodel.msim"""
    
    def prepareParams(self):
        """Підготовлює парметри Maplesim моделі"""
        pu=self.pu # СШНУ
        # підготувати параметри для моделі Maplesim
        modelParams={'p1_L':pu.pumpingUnit.p1['L'], # параметри ланок
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
                     'p5_L':pu.pumpingUnit.p5['L'], # радіус противаги
                     'p5_m':pu.pumpingUnit.p5['m'], # маса противаги
                     'p5_Iz':pu.pumpingUnit.p5['Iz'],
                     'dx':pu.pumpingUnit.dx, # відстані між осями
                     'dy':pu.pumpingUnit.dy,
                     'r_k':pu.pumpingUnit.r, # радіус кривошипа
                     'const1_k':-pu.pumpingUnit.omega # кутова швидкість кривошипів, рад/с
                     # знак "-" для відємного дизаксіалу (Архіпов с.15)
                     # тоді хід вверх довший
                     # ще бажано FF1_teta=pi/2
                     }
        
        for p in sorted(modelParams.keys()):
            print p, modelParams[p]
            
        self.modelParams=modelParams
                
if __name__ == '__main__':
    model=Maplesim_model() # !!! задати вірний клас
    #model.optimize()
    
    model.pu=PUmodel.PU()
    print "Fmin,Fmax=", model.pu.polishedRodForce()
    #print "Rp=", model.pu.balancing()
    model.prepareParams()
    
    model.createMaplesimModel(execute=1)
    #model.drawDynamometerCard(S0=math.pi*0.019**2/4) #-math.pi*0.0085**2/4)
    #model.drawDynamometerCard()
    # іншим способом побудови динамограм є збереження maplesim plots у .csv
    #model.drawDynamometerCardMulti(['1.csv','2.csv','3.csv'],[326.0e-6,326.0e-6,326.0e-6],[0.166666,0.166666,0.166666],['k-','k--','k:'])
    
    #f=model.modelParams['S1_freqHz']
    model.drawDynamometerCardMulti(['testAPI.csv'],[1.0],[0.10833],['k-'])
    #model.drawDynamometerCardMulti(['2_.csv','1_.csv','2.csv','1.csv'],[1.0,1.0,1.0,1.0],[0.10833,0.10833,f,f],['k--','k--','k-','k-'])
    