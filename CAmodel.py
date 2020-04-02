# -*- coding: CP1251 -*-

# -Добавлено можливість задання різних відстаней між автоматами
# -Виправлено розміри списків K, D
# -Добавлено розрахунок напружень
# -Необхідно задаватись величиною подрібнення початкового кроку наближення
# -Виправлена помилка де завжди a.dt=0.01. Тепер a.dt=timeStep

import math, os
import PUmodel
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', family='Arial') # для підтримки Юнікоду

class Automaton(object):
    """Абстрактний автомат для моделювання руху вузла пружного стержня"""
    def __init__(self):
        """Конструктор автомата"""
        self.x=None # координата
        self.x0=None # початкова координата
        self.prevx=None # попередня координата
        self.bc=None # гранична умова (0 - нерухомий)
        self.f=0.0 # зовнішня сила
        self.fs=0.0 # вага секції справа
        self.left=None # лівий сусід (верхній)
        self.right=None # правий сусід (нижній)
        self.delta=None # початкова відстань до правого сусіда (rSection.length)
        self.k=None # коефіцієнт жорсткісті секції справа (rSection.stiffness())
        self.c=0.0 # коефіцієнт опору матеріалу секції справа
        self.m=1.0 # маса секції справа
        self.a=0.0 # прискорення
        self.v=0.0 # швидкість
        self.prevv=0.0 # попередня швидкість
        self.eps=0.0 # нев'язка
        self.dx=1.0 # крок наближення
        self.area=None # площа поперечного січення пружини (стержня) справа (rSection.area())
        self.domain=None # об'єкт класу CA_model
        self.rSection=None # секція справа
        self.alfa=0.0 # кут відхилення свердловини від вертикалі біля автомата

    def rule(self):
        """Правило поведінки автомата"""
        
        if self.bc!=None: # якщо задана гранична умова
            self.x=self.bc # задати x
        
        if self.domain.t!=0.0: # якщо час не 0
            self.v=self.velocity() # розрахувати швидкість
            self.a=self.acceleration() # розрахувати прискорення
        else: self.v=self.a=0.0
        
        if self.bc!=None: # якщо задана гранична умова
            return # вийти (не потрібно знаходити x)
        
        Fk=self.rightForce()-self.leftForce() # сили пружності
        Fv=self.rightForce2()-self.leftForce2() # сили демпфування
        Fa=self.dynamicForce() # сила інерції
        # зовнішні сили
        F=self.fs+self.f+self.pistonForce()+self.frictionForce()+self.hydrodynamicResistanceForce() 
        
        # рівновага усіх сил біля автомата: m*a+c*v+k*u=F
        # якщо сила>0 то вона направлена вниз
        # Напрямок осі x:
        # [0]---[1]---[2]--->+x,+v,+a,+f,правий,нижній 
        # обчислити нев'язку
        self.eps=Fa+Fv+Fk-F
        #print "eps=",self.eps
        
        # якщо нев'язка відрізняється від 0
        if self.eps>0:
            self.x=self.x-self.dx # збільшити x на малу величину
        if self.eps<0:
            self.x=self.x+self.dx # зменшити x на малу величину
        
    def velocity(self):
        return (self.x-self.prevx)/self.domain.dt # швидкість
        
    def acceleration(self):
        return (self.v-self.prevv)/self.domain.dt # прискорення
               
    def leftForce(self):
        """Повертає значення сили пружини зліва (вверху)"""
        if self.left: # якщо є лівий сусід
            dLeftBegin=self.left.delta # початкова довжина пружини зліва
            dLeft=self.x-self.left.x # довжина пружини зліва
            fLeft=(dLeftBegin-dLeft)*self.left.k # сила пружини зліва
        else:
            fLeft=0.0
        return fLeft
        
    def leftForce2(self):
        """Повертає значення сили демпфера зліва (вверху)"""
        if self.left: # якщо є лівий сусід
            vrel=self.v-self.left.v # відносна швидкість
            fLeft=-self.left.rSection.materialRodResistance(None)*vrel
        else:
            fLeft=0.0
        return fLeft
    
    def rightForce(self):
        """Повертає значення сили пружини справа (внизу)"""
        if self.right: # якщо є правий сусід
            dRightBegin=self.delta # початкова довжина пружини справа
            dRight=self.right.x-self.x # довжина пружини справа 
            fRight=(dRightBegin-dRight)*self.k # сила пружини справа
        else:
            fRight=0.0
        return fRight
           
    def rightForce2(self):
        """Повертає значення сили демпфера справа (внизу)"""
        if self.right: # якщо є правий сусід
            vrel=self.right.v-self.v # відносна швидкість
            fRight=-self.rSection.materialRodResistance(None)*vrel
        else:
            fRight=0.0
        return fRight
        
    def dynamicForce(self):
        """Сила інерції штанг F=m*a"""
        return self.m*self.a
    
    def dynamicFluidForce(self):
        """Сила інерції рідини. Діє на плунжер"""
        k=self.domain.PU.pump.pistonArea()/self.domain.PU.well.sectionalArea()
        a=self.a*k # прискорення рідини менше прискорення плунжера у k раз
        return -a*sum(self.domain.PU.suckerRodString.fluidMassList())
            
    def normalForce(self):
        """Сума нормальних сил у вузлі від ваги, лівого і правого натягу"""
        F1=self.f*math.tan(self.alfa) # нормальна сила від ваги (f - осьова сила від ваги)
        F2=F3=0.0
        if self.left:
            F2=self.leftForce()*math.sin(self.an1) # нормальна сила від лівого натягу
        if self.right:
            F3=self.rightForce()*math.sin(self.an2) # нормальна сила від правого натягу
        return F1+abs(F2+F3) # якщо нехтувати зазором між стінкою НКТ і вузлом
        
    def frictionForce(self):
        """Сила тертя"""
        # сила тертя на інтервалі
        #Ft=C*(Psh*math.sin(a1)-F12*math.sin(a12-a1)) # ???
        #Ft=C*(Psh*math.sin(a1)+2*F12*math.tan(a2/2-a1/2)) # a12-a1=a2/2-a1/2
        #Ft=C*(Psh*math.sin(a1)+F12*(a2-a1)) # формула Песляка
        
        v=self.v #швидкість
        Fn=self.normalForce() #нормальна сила # для перевірки напрямку введіть тут 1000.0
        fc=0.16 #0.25 # коефіцієнт тертя ковзання (сталь-сталь зі змащуванням)
        #F=-fc*abs(Fn)*math.tanh(v/0.01) # спрощена модель
        
        Fc=abs(Fn)*fc # Coulomb friction force
        fmax=0.2 # коефіцієнт тертя ковзання (сталь-сталь з малим змащуванням)
        Fmax=abs(Fn)*fmax # Maximum sliding friction force
        vs=0.1 # Stribeck sliding velocity coefficient
        n=1.0 # Decay exponent
        # сумарна сила тертя
        F=-math.tanh(v/0.01)*(Fc+(Fmax-Fc)*math.exp(-(abs(v)/vs)**n))
        return F
    
    def hydrodynamicResistanceForce(self):
        """Cила гідродинамічного опору секції"""
        v=self.v # швидкість
        if v<=0: # якщо хід вверх
            k=self.domain.PU.pump.pistonArea()/self.domain.PU.well.sectionalArea()
            vu=v-v*k # швидкість штанг відносно рідини
            F=-self.c*vu*10.0 #увага!!! збільшено
            # F=0.0 # якщо швидкості штанг і рідини рівні
        else: # якщо хід вниз
            F=-self.c*v*10.0  #увага!!! збільшено
        return F*math.tanh(abs(v)/0.01) # згладжуємо біля 0
    
    def pistonForce(self):
        """Сума зовнішніх сил на плунжері насоса"""
        v=self.v # швидкість
        F=0.0 # для усіх вузлів крім останнього 
        if self.right!=None: return F # якщо не останній вузол (насос)
        
        if v<0: # якщо рух плунжера вгору
            # направлені вниз:
            F+=self.domain.PU.fluidWeight()
            F+=self.domain.PU.outPresForce()
            F+=abs(self.domain.PU.tubePressureFrictLossForse(v)) # тут abs про всяк випадок
            F+=self.dynamicFluidForce()
            # направлені вверх:
            F+=-self.domain.PU.outTubePresForce()
            F+=-self.domain.PU.outTubeFluidWeight()
        if v!=0:
            # наступні сили діють під час РУХУ (!) вгору і вниз
            # і направлені проти руху плунжера
            F+=-math.copysign(1, v)*abs(self.domain.PU.pump.valveHydroResistForce(v))       
            F+=-math.copysign(1, v)*self.domain.PU.pump.pistonFrictionForce() # модель тертя плунжера
        
        F=F*math.tanh(abs(v)/0.01) # ф-ція згладжування біля 0 (F=max, якщо v=0.01) 
        #або F=F*(1.0-math.exp(v/0.01))
        return F
         
    def rightStress(self):
        """Повертає значення напруження в пружині справа (внизу)"""
        if self.right:
            return self.rightForce()/self.area
        else:
            return 0.0
        
################################################################
class CA_model(object):
    """Система автоматів (пружний стержень)"""
    def __init__(self,PU):
        """Створює систему автоматів за об'єктом класу PU (СШНУ)
        Схема системи автоматів:
                секція1                      секція2
        O--------------------O-----------------------------------O
        m0=0      d0         m1                d1                m2 
        c0=0      k0         c1                k1                c2
        f0        a0         f1                a1                f2
       alfa0=0   ans0       alfa1             ans1              alfa2                                       
        Перший автомат (0) допоміжний.
        """
        
        self.PU=PU # СШНУ
        rs=PU.suckerRodString.items # список секцій
        n=len(rs)+1 # кількість вузлів (автоматів)
        # Перший автомат (0) допоміжний.
        X=[0.0]+PU.suckerRodString.xList() # cписок координат вузлів
        M=[0.0]+PU.suckerRodString.massList() # маса вузлів
        Fs=[0.0]+PU.suckerRodString.weightForceList() # список сил від ваги секцій
        C=[0.0]+PU.suckerRodString.hydrodynamicRodResistanceList() # гідродинамічні опори секцій

        AN=[0.0]+PU.suckerRodString.wellAngleList()
        #ANS=PU.suckerRodString.secAngleList()
        AN1=[None]+pu.suckerRodString.angleTop()
        AN2=[0.0]+pu.suckerRodString.angleBottom()
        
        self.dt=self.PU.timeStep # крок часу
        self.T=self.timeList(PU.timeEnd, PU.timeStep) # список значень часу [t0,t1,t2,...]
        self.t=self.T[0] # поточне значення часу (0.0)
        
        BCdict={} #словник граничних умов (див. self.BCList)
        # закон руху точки підвіски
        BCdict[0]=[PU.pumpingUnit.motionX(t) for t in self.T]
        self.BC=CA_model.BCList(BCdict, n, self.T, default=None) # список історії граничних умов [[bc0,bc1,bc2,...],...]
                
        F=[0.0]*n # список інших зовнішніх сил (якщо потрібно)
        Fdict=dict(zip(range(n),F)) # словник сил (див. self.BCList)
        self.F=CA_model.BCList(Fdict, n, self.T, default=0.0) # список історії зовнішніх сил [[f0,f1,f2,...],...]

        self.p=[Automaton() for i in range(n)] # список автоматів
        
        # для кожного автомата крім останнього задати властивості
        for i in range(n-1):
            self.p[i].rSection=rs[i]
            self.p[i].k=rs[i].stiffness()
            self.p[i].delta=rs[i].length
            self.p[i].area=rs[i].area()
        
        # задаємо інші властивості усіх автоматів
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
            
        for i in range(n): # задаємо сусідів
            if i!=0: # якщо не перший
                self.p[i].left=self.p[i-1] # задати лівого сусіда
            if i!=len(self.p)-1: # якщо не останній
                self.p[i].right=self.p[i+1] # задати правого сусіда
        
        # словник історій результатів
        # наприклад, список історії значень x [[x0,x1,x2,...],...]
        self.history={'T':self.T,'X':[],'V':[],'A':[],'Fl':[],'Fr':[],'Sr':[]}
                      
    def run(self):
        """Еволюція системи автоматів в конкретний момент часу"""
        eps=0.01*len(self.p) # похибка розрахунку
        e=eps+1 # сума нев'язок
        eprev=eps+1 # попередня сума нев'язок
        # початковий крок наближення (залежить від delta)
        dx=min([a.delta for a in self.p if a.delta!=None])/10.0 # !!! Змінювати в межах 10-100
        niter=0 # кількість ітерацій (для оптимізації алгоритму)
        while e>eps: # поки сума нев'язок > eps
            for a in self.p: # для кожного автомата
                a.dx=dx # змінити крок наближення
                a.rule() # застосувати правило поведінки автомата
            e=sum([abs(a.eps) for a in self.p]) # сума нев'язок
            #print e
            if e>=eprev: # якщо сума нев'язок >= попередньої
                dx=dx*0.935 # зменшити крок !!! Змінювати в межах 0.1-0.99
                #print "dx=",dx
            eprev=e # запам'ятати суму нев'язок
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
        """Динамічна еволюція системи автоматів"""
        # початкові значення координат, швидкостей і прискорень
        for a in self.p: # для кожного автомата
            a.prevx=a.x # попередня координата
            a.prevv=0.0 # попередня швидкість
        for i,t in enumerate(self.T): # для кожного значення часу
            self.t=t # поточне значення часу
            print "solve t=",t
            for a,f,bc in zip(self.p, self.F[i], self.BC[i]): # для кожного автомата
                a.f=f # задати силу в момент часу t
                a.bc=bc # задати граничну умову в момент часу t
            self.run() # знайти урівноважений стан
            for a in self.p: # для кожного автомата
                a.prevx=a.x # зап'ямятати поточні x та v
                a.prevv=a.v
            self.appendHistory()# запам'ятати історії
        self.writeHistory()
            
    def appendHistory(self):
        """Додає в списки історії список значень для поточного часу"""
        self.history['X'].append([a.x for a in self.p])
        self.history['V'].append([a.v for a in self.p])
        self.history['A'].append([a.a for a in self.p])
        self.history['Fl'].append([a.leftForce() for a in self.p])
        self.history['Fr'].append([a.rightForce() for a in self.p])
        self.history['Sr'].append([a.rightStress() for a in self.p])
        
    def draw(self):
        """Рисує положення автоматів псевдографікою (в текстовому режимі)"""
        s=""
        prev=0
        for a in self.p:
            s+="-"*int(round(a.x-prev)-1)+"O"
            prev=a.x # координата попереднього автомата
        print s
        #time.sleep(1)
        
    def drawPlot(self, itemIndexes, key, x0scale=1, toFile=True):
        """Рисує залежність величини x[i] від t, де i - індекс автомата
        itemIndexes - список індексів автоматів
        key - ключ історії x
        x0scale - масштаб відстаней між першими точками кривих по осі X (1...0):
                  1 - реальний, 0 - відстані рівні 0.
                  Крива для 0-го вузла не зміщується.
                  Використовується для покращення візуалізації переміщень вузлів,
                  якщо вони малі відносно відстаней між вузлами"""
        plt.figure()
        plt.gca().invert_yaxis() # обернути вісь y
        T=self.history['T'] # список значень часу
        x00=self.history[key][0][0] # x 0-го вузла для t=0
        for i in itemIndexes: # для кожного індексу
            x0i=self.history[key][0][i] # x i-го вузла для t=0
            x0=(1-x0scale)*(x0i-x00) # величина зміщення по осі X
            X=[x[i]-x0 for x in self.history[key]] # список історії вузла "i"
            plt.plot(T, X, 'k-') # (T, X, 'r-', T, X,'ro')
        plt.grid(True)
        plt.xlabel(u't, c') #надпис осі x
        plt.ylabel(key) #надпис осі y
        if toFile:
            fileName='CAmodelPlot{0}{1}.png'.format(itemIndexes,key)
            plt.savefig(fileName)
        else:
            plt.show()
        
    def drawDynamometerCard(self, itemIndexes, keyX, keyY, signY=-1, toFile=True, tstart=None, tend=None):    
        """Рисує динамограмму - залежність Y(X), наприклад, сила(переміщення)"""
        plt.figure()
        plt.gca().invert_xaxis() # обернути вісь x
        plt.grid(True)
        plt.xlabel(u"S, м") #надпис осі x
        plt.ylabel(u"F, Н") #надпис осі y        
        for i in itemIndexes: # для кожного вузла
            # індекси початкового і кінцевого часу в історії
            tstarti=self.history['T'].index(tstart) if tstart in self.T else None 
            tendi=self.history['T'].index(tend)+1 if tend in self.T else None
            X=[x[i] for x in self.history[keyX][tstarti:tendi]] # список історії X вузла "i"
            x0=self.history[keyX][0][i] # початкове значення
            X=[x-x0 for x in X] # відносно початкової координати
            Y=[signY*y[i] for y in self.history[keyY][tstarti:tendi]] # список історії Y вузла "i"
            plt.plot(X, Y, 'k-')
        if toFile:
            fileName='CAmodelDynCard{0}{1}.png'.format(keyX,keyY)
            plt.savefig(fileName)
        else:
            plt.show()
        
    def writeHistory(self):
        """Зберігає історію в бінарний файл"""
        import pickle
        f = open('CAmodelHistory.pkl', 'wb')
        pickle.dump(self.history, f) # зберегти історію
        f.close()
        
    def readHistory(self):
        """Читає історію з бінарного файлу"""
        import pickle
        fileName='CAmodelHistory.pkl'
        if not os.path.isfile(fileName): # якщо файлу немає
            print "No file "+fileName
            return False # вийти, повернувши False
        f = open(fileName, 'rb')
        self.history=pickle.load(f) # завантажити історію
        f.close()
        print "Warning! Old results"
        return True
       
    def writeCSV(self,itemIndexes,key):
        """Записує значення історій у файл CSV.
        itemIndexes - список індексів автоматів
        key - ключ історії"""
        import csv
        csv_file=open("RodHistory_"+key+".csv", "wb")
        writer = csv.writer(csv_file,delimiter = ';')
        for time,data in zip(self.T,self.history[key]):
            X=[data[i] for i in itemIndexes] # лише дані за вказаними індексами
            writer.writerow([time]+X) # записати рядок
        csv_file.close()

    def info(self):
        """Виводить значення усіх атрибутів автоматів"""
        attrs=["rSection","k","delta","area","bc","x","prevx","v","prevv",\
               "a","m","alfa","an1","an2","c","eps","f"]
        for attr in attrs:
            print attr+"=", [a.__getattribute__(attr) for a in self.p]
        print "leftForce=", [a.leftForce() for a in self.p]
        print "rightForce=", [a.rightForce() for a in self.p]
         
                    
    @staticmethod
    def timeList(timeEnd,timeStep):
        """Повертає список значень часу"""
        nTimeSteps=int(timeEnd/timeStep)
        return [x*timeStep for x in range(nTimeSteps+1)]  

    @staticmethod
    def BCList(BC, n, timeList, default=None):
        """
        Повертає список граничних умов (або сил) для автоматів 1,2,... у вигляді:
        [[bc0,bc1,...],...],
        де перший список відповідає першому значенню часу.
        
        BC - словник граничних умов, де ключі - індекси КА, а значеннями мужуть бути:
        1.Дійсні (постійне значення): 0.0
        2.Словники пар час-значення: {0.0:1.0,0.05:2.0}
        3.Списки історії значень: [1.0,2.0,3.0,4.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]
             або: [math.sin(x) for x in range(nt)]
        n - кількість автоматів
        timeList - список значень часу
        default - значення списків за замовчуванням"""
        
        nt=len(timeList) # кількість значень часу
        res=[] # результат
        for i in range(nt):
            res.append([default]*n) #  заповнюємо res різними(!) списками
                            
        for k in BC: # для кожного заданого вузла
            if type(BC[k]).__name__=='dict': # якщо задана для конкретних інтервалів часу
                i=0
                prev=None # попередній
                for t in res: # для кожного значення часу
                    if BC[k].has_key(timeList[i]):             
                        t[k]=BC[k][timeList[i]]
                        prev=t[k]
                    else:
                        t[k]=prev
                    i=i+1
            elif type(BC[k]).__name__=='list': # якщо задана для кожного інтервалу часу
                i=0
                for t in res: # для кожного значення часу
                    t[k]=BC[k][i]
                    i=i+1  
            else: # якщо константа
                for t in res: # для кожного значення часу
                    t[k]=BC[k]
        #print res            
        return res 
    
    def resultsCAmodel(self):
        """Виводить результати автоматної моделі"""
        if self.history['X']==[]: # якщо немає результатів
            oldHistExist=self.readHistory() # прочитати старі результати
            if not oldHistExist: return
            
        tend=self.T[-1]
        T=1.0/self.PU.pumpingUnit.ns # період
        tstart=tend-T # тільки останній період
        # знайти найближче значення до tstart у списку self.T
        dt=[abs(t-tstart) for t in self.T]
        tstart=self.T[dt.index(min(dt))]
        self.drawDynamometerCard([0,1,2,3,4,5,6,7], "X", "Fr", tstart=tstart, tend=tend)
        self.drawPlot([0,1,2,3,4,5,6,7,8],'X', x0scale=0.001) 
        self.drawPlot([0,1,2,3,4,5,6,7],'Sr')
        
if __name__ == '__main__':
    pu=PUmodel.PU()
    model=CA_model(pu)  
    #model.info()
    #model.run() # статика
    #model.info()
    
    model.runDynamic() # динаміка
    #model.info()
    model.resultsCAmodel()
    