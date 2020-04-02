# -*- coding: cp1251 -*-
#Мищенко И.Т. Скважинная добыча нефти: Учебное пособие для вузов. — М: ФГУП Изд-во «Нефть и газ» РГУ нефти и газа им. И.М. Губкина, 2003. — 816 с.
import math

# Для моделювання пустотілої трубчатої штанги її замінюємо
# еквівалентною суцільною з такими ж:
# зовнішнім діаметром, пружністю, масою
# але тоді у неї буде інший модуль Юнга і густина
def eqYoungsModulus(E,D,d):
    """Еквівалентний модуль Юнга суцільної штанги діамером D
    з умови однаковості пружності k=E*A/L=Eeq*Aeq/L
    суцільної (діаметром D) і трубчатої штанги (зовн. діаметром D,
    внутрішнім d і модулем Юнга E)"""
    A=math.pi*(D**2.0-d**2.0)/4.0 # площа січення труби
    Aeq=math.pi*(D**2.0)/4.0 # площа січення суцільної штанги
    return E*A/Aeq

def eqDensity(Ro,D,d,Roo=0.0):
    """Еквівалентна густина суцільної штанги діамером D
    з умови однаковості маси m=Ro*A*L+Roo*Ao*L=Roeq*Aeq*L
    суцільної (діаметром D) і трубчатої штанги (зовн. діаметром D,
    внутрішнім d, густиною тіла Ro і тіла отвору Roo)
    Якщо Roo=0 - пустотілі штанги наповенені повітрям"""
    A=math.pi*(D**2.0-d**2.0)/4.0 # площа січення труби
    Ao=math.pi*(d**2.0)/4.0 # площа січення отвору
    Aeq=math.pi*(D**2.0)/4.0 # площа січення суцільної штанги
    return (Ro*A+Roo*Ao)/Aeq
    

class SuckerRod(object):
    """Клас насосних штанг або секції однотипних штанг. Одиниці СІ"""
    diametrs=[0.016,0.019,0.022,0.025] # діаметри тіла
    length=8.0 # довжина (за замовчуванням)
    materials={"сталь":{"модуль пружності":2.1e11,
                        "густина":8000.0}, # тут густина дещо завищена для врахування мас муфт і головок
               "склопластик":{"модуль пружності":0.5e11,
                              "густина":2100.0},
               "стальТр22Екв":{"модуль пружності":eqYoungsModulus(2.1e11,0.022,0.0085),
                              "густина":eqDensity(8000.0,0.022,0.0085,1000.0)},
               "склопластикТр22Екв":{"модуль пружності":eqYoungsModulus(0.5e11,0.022,0.0085),
                              "густина":eqDensity(2100.0,0.022,0.0085,1000.0)}}
    #!!!якщо Roo=0, то пустотілі штанги є герметичними і наповенені повітрям (але це недоцільно)
    alfa=0.0 # кут відхилення свердловини внизу секції від вертикалі (рад.)
    
    def __init__(self, PU, diametr, length=8.0, material="сталь", alfa=0.0):
        """Конструктор"""
        self.PU=PU # свердловинна штангова насосна установка
        self.diametr=diametr # діаметр тіла штанги
        self.length=length # довжина штанги або секції
        self.material=material # матеріал
        self.alfa=alfa
    
    def area(self):
        """Площа поперечного січення S=pi*r^2"""
        return math.pi*(self.diametr/2.0)**2
    
    def volume(self):
        """Об'єм, який займає штанга в рідині V=L*pi*r^2"""
        return self.length*math.pi*(self.diametr/2.0)**2
    
    def mass(self):
        """Маса m=V*ro"""
        return self.volume()*self.materials[self.material]["густина"]
    
    def fluidVolume(self):
        """Об'єм рідини в НКТ навколо штанги (секції)"""
        A=self.PU.well.sectionalArea()-self.area() # увага! тут площа НКТ, а не плунжера
        return A*self.length
    
    def fluidMass(self):
        """Маса рідини в НКТ навколо штанги (секції)"""
        return self.PU.well.density*self.fluidVolume()
    
    def weight(self):
        """Вага в рідині"""
        fluidDensity=self.PU.well.density # густина рідини
        Farh=fluidDensity*9.81*self.volume() # сила Архімеда
        w=self.mass()*9.81 # вага в повітрі
        return w-Farh
    
    def stiffness(self):
        """Жорсткість k=E*S/L"""
        E=self.materials[self.material]["модуль пружності"]
        return E*self.area()/self.length
    
    def hydrodynamicRodResistance(self):
        """Гідродинамічний опір 'с' штанги (секції) за Пірвердяном
        Враховується під час ходу вниз. Сила опору F=c*v"""
        msh=self.PU.well.diametr/self.diametr
        Msh=1.0/((math.log(msh)*(msh**2+1)/(msh**2-1))-1) # коефіцієнт Мшт
        # v=2nS - швидкість штанг
        R=math.pi**2*self.PU.well.viscosity*self.PU.well.density*Msh
        return R*self.length
    
    def hydrodynamicRodResistance2(self):
        """Гідродинамічний опір 'с' штанги (секції) за Андреєвим
        Враховується під час ходу вниз. Сила опору F=c*v"""
        m=self.diametr/self.PU.well.diametr
        R=16.9*self.PU.well.viscosity*self.PU.well.density*m**5.49
        R=R*1000 # одиниці ?
        return R*self.length
    
    def materialRodResistance(self,sigmaA=None):
        """Еквівалентний коефіцієнт опору матеріалу, Н*с/м
        sigmaA - амплітуда напружень (Па)
        sigmaA=None - не залежить від ампл. напруж."""
        k=self.stiffness() #жорсткість стержня
        #визначаємо коефіцієнт поглинання матеріалу
        if self.material.startswith("сталь"):
            if sigmaA!=None: #якщо залежить від ампл. напруж.
                psi=6e-5*sigmaA/1000000.0
            else: #якщо не залежить від ампл. напруж.
                #psi<0.13 (метали в пружному діапазоні)[Damping of materials and members in structures]
                #psi=0.25...0.5 (металічні структури)
                #psi=0.38...0.88(металічні структури зі зєднаннями)
                psi=0.01
        if self.material.startswith("склопластик"):
            if sigmaA!=None: #якщо залежить від ампл. напруж.
                psi=4e-4*sigmaA/1000000.0 #справочник "Вибрации в технике"
            else: #якщо не залежить від ампл. напруж.
                psi=0.08
        omega=self.PU.pumpingUnit.omega #кутова частота зовнішньої сили
        #m=self.mass() # маса стержня
        #omega=math.sqrt(k/m) # власна частота стержня
        #[Вибрации в технике т6 с130 ф13,19]
        c=psi*k/(2*math.pi*omega) #k=E*S/L
        
#         #!уточнити врахування окремо опору зєднань
#         nz=int(self.length/8) # кількіст муфтових зєднань на секції
#         kz=1000000.0 # жорсткість муфтового зєднання ?
#         cz=0.2*kz/(2*math.pi*omega) #екв коеф опору муфтового зєднання (psi=0,006...0,4)
#         czs=1/(nz*(1/cz)) # екв коеф опору усіх зєднань на секції
#         return 1/(1/c+1/czs)# екв коеф опору секції (матеріал+зєднання)
        return 20.0*c # увага збільшено у 20.0 (для спрощ. моделі 50)

###############################################################   
class SuckerRodString(object):
    """Клас колони насосних штанг"""
    items=[] # список штанг (секцій) [штанга1, штанга2,...](починаючи з гирла)
    
    def __init__(self, PU, sections):
        """Конструктор.
        PU - СШНУ
        sections=[(діаметр, довжина, матеріал, кут відхилення від вертикалі(рад.),...]"""
        self.PU=PU
        for diametr, length, material, alfa in sections: # для кожної секції
            # добавити штангу (секцію однотипних штанг)
            self.items.append(SuckerRod(PU, diametr, length, material, alfa))
    
    def length(self,start=0,end=None):
        """Довжина. start, end - індекси зрізу"""
        return sum([r.length for r in self.items[start:end]])
    
    def xList(self):
        """Список координат нижніх кінців секцій"""
        X=[]
        x=0.0 # координата 
        for r in self.items:
            x=x+r.length # добавити довжину даної секції
            X.append(x)
        return X
    
    def height(self):
        """Висота проекції колони штанг на вертикаль"""
        return sum([r.length*math.cos(r.alfa) for r in self.items])
    
    def volume(self,start=0,end=None):
        return sum([r.volume() for r in self.items[start:end]])
    
    def stiffness(self):
        """Сумарний коеф жорсткості колони 1/k=1/k1+1/k2+..."""
        return 1/sum([1/r.stiffness() for r in self.items])
    
    def damping(self):
        """Сумарний коеф демпфування колони 1/d=1/d1+1/d2+...
        [Вибрации в технике т6 с177]"""
        return 1/sum([1/r.materialRodResistance() for r in self.items])
    
    def mass(self,start=0,end=None):
        return sum([r.mass() for r in self.items[start:end]])
    
    def massList(self):
        """Список мас секцій"""
        return [r.mass() for r in self.items]
    
    def fluidMassList(self):
        """Список мас рідини біля вузлів"""
        return [r.fluidMass() for r in self.items]
    
    def weight(self,start=0,end=None):
        return sum([r.weight() for r in self.items[start:end]])
    
    def wellAngleList(self):
        """Список кутів свердловини внизу секцій"""
        return [x.alfa for x in self.items] 
    
#     def secAngle(self,L1,a1,L2,a2):
#         """Повертає кут відхилення від вертикалі секції між ділянками свердловини з
#         довжинами L1,L2 та відхиленням від вертикалі a1,a2 (рад)
#         Кінці секції торкаються середин ділянок
#         """
#         x=L1*math.sin(a1)/2+L2*math.sin(a2)/2 # ширина
#         y=L1*math.cos(a1)/2+L2*math.cos(a2)/2 # висота
#         a=math.atan(x/y) # кут відхилення від вертикалі секції
#         return a
        
    def secAngleList(self):
        """Повертає список кутів секцій an_i
        за кутами свердловини внизу секції alfa_i.
             an0          an1          an2
        ------------O------------O------------O
                  alfa0        alfa1        alfa2
        Розраховує приблизно за умови рівності довжин секцій.
        """
        angles=[self.items[0].alfa] # кут верхньої секції=куту верхньої ділянки свердловини  
        aprev=self.items[0].alfa # кут першої ділянки свердловини (вверху). Як правило 0.
        for r in self.items[1:]: # для кожної ділянки крім першої
            a=(aprev+r.alfa)/2 # розрахувати наближено кут
            angles.append(a) # добавити в список
            aprev=r.alfa # запам'ятати як попередню секцію
        return angles
      
    def angleTop(self):
        """Кути між верхньою штангою і свердловиною"""
        AN=self.wellAngleList()
        ANS=self.secAngleList()
        return [ans-an for an,ans in zip(AN,ANS)]
        
    def angleBottom(self):
        """Кути між свердловиною і нижньою штангою
        Останнє значення 0"""
        AN=self.wellAngleList()
        ANS=self.secAngleList()     
        return [an-ans for an,ans in zip(AN[:-1],ANS[1:])]+[0.0]  
    
    def hydrodynamicRodResistanceList(self):
        """Список гідродинамічних опорів секцій"""
        return [r.hydrodynamicRodResistance() for r in self.items]
    
    def weightForceList(self):
        """Cписок осьових сил від ваги секцій
        (з врахуванням кута свердловини)"""
        return [r.weight()*math.cos(r.alfa) for r in self.items]

        
class Well(object):
    """Клас свердловини і НКТ"""
    diameters=[60.3-2*5.0, 73.0-2*5.5, 73.0-2*7.0, 88.9-2*6.5, 114.3-2*7.0]
    diametr=diameters[4]/1000.0 #4 діаметр НКТ (однаковий на інтервалах)
    density=1000.0 #1000.0 густина суміші [Белов]
    viscosity=0.000002 #0.000001 # кінематична в'язкість (змінна по висоті!!!)
    bulkModulus=2.2e9 # об'ємний модуль пружності
    
    def __init__(self,PU):
        self.PU=PU
        self.length=PU.suckerRodString.length()
    def sectionalArea(self):
        """Площа поперечного січення отвору НКТ"""
        return math.pi*(self.diametr/2.0)**2
    def area(self):
        """Площа поперечного січення S=pi*r^2"""
        D=114.3/1000
        d=(114.3-2*7.0)/1000
        return math.pi*(D/2.0)**2-math.pi*(d/2.0)**2
    def stiffness(self):
        """Жорсткість k=E*S/L"""
        E=2.1e11 # модуль пружності
        return E*self.area()/self.length
    def volume(self):
        return self.length*self.area()
    def mass(self):
        return self.volume()*8000.0
    def weight(self):
        """Вага в рідині"""
        Farh=self.density*9.81*self.volume() # сила Архімеда
        w=self.mass()*9.81 # вага в повітрі
        return w-Farh
    def materialRodResistance(self):
        """Еквівалентний коефіцієнт опору матеріалу, Н*с/м"""
        k=self.stiffness() # жорсткість
        psi=0.2 # коефіцієнт поглинання матеріалу
        omega=self.PU.pumpingUnit.omega #кутова частота зовнішньої сили
        return psi*k/(2*math.pi*omega)
    
class Pump(object):
    """Клас свердловинного насоса"""
    diameters=[27, 32, 38, 44, 50, 57, 57, 63, 70, 95]
    diametr=38/1000.0 #38 діаметр плунжера
    valveDiametr=0.025 # діаметр отвору сідла клапана
    openValveResistance=0.00001 # опір відкритого клапана ???
    pistonGap=0.1e-3 #0.1e-3 зазор між плунжером і циліндром
    # далі v - швидкість плунжера, vv - швидкість в сідлі клапана  !!!
    
    def __init__(self,PU):
        self.PU=PU # СШНУ
        
    def volumeFlowRate(self,v):
        """Теоретична подача насоса (об'ємна витрата) Q=A*v"""
        return self.pistonArea()*v #self.PU.pumpingUnit.velAvg()
    
    def pistonArea(self):
        """Площа поперечного січення плунжера"""
        return math.pi*(self.diametr/2.0)**2 # !!! приблизно
    
    def valveDischargeCoefficient(self,Re):
        """Коефіцієнт витрати клапана. Залежить від числа Рейнольдса. Гіматудінов"""
        if Re<=225:
            mu=0.0846*Re**0.2872
        elif Re<=30000:
            mu=0.4
        elif Re<=300000:
            mu=0.0085*Re**0.3764
        else: mu=1.0
        return mu
    
    def valveLossCoefficient(self,Re):
        """Коефіцієнт витрат клапана (місцевого опору)"""
        return 1.0/self.valveDischargeCoefficient(Re)**2.0
    
    def valveV(self,v):
        """Швидкість потоку в сідлі, v - швидкість плунжера"""
        return self.volumeFlowRate(v)/(math.pi*(self.valveDiametr/2.0)**2)
        
    def valveRe(self,vv):
        """Число Рейнольдса, vv - швидкість потоку в сідлі"""
        nu=self.PU.well.viscosity
        Re=abs(vv)*self.valveDiametr/nu
        return Re
    
    def dynamicPressure(self, vv):
        """Динамічний тиск, vv - швидкість потоку в сідлі"""
        ro=self.PU.well.density
        return vv*vv*ro/2
    
    def dPvalve(self,vv):
        """Перепад тисків на клапані. Гіматудінов
        vv - швидкість потоку в сідлі"""
        Re=self.valveRe(vv)
        dP=self.valveLossCoefficient(Re)*self.dynamicPressure(vv)
        return dP
    
    def valveHydroResistForce(self,v):
        """Сила гідродинамічного опору клапана насоса.
        Діє під час ходу вниз і вверх. Направлена проти руху
        v - швидкість плунжера"""
        vv=self.valveV(v)
        return self.dPvalve(vv)*self.pistonArea() # повертає завжди додатнє
    
    def pistonFrictionForce(self):
        """Сила тертя плунжера об стінки циліндра, Н.
        Формула Сердюка. Орієнтовно 1000 Н
        Діє під час ходу вниз і вверх. Направлена проти руху"""
        return 1.84*self.diametr/self.pistonGap-137 # під час змащування водою
 
class PumpingUnit():
    """Клас верстата-гойдалки"""
    stroke=3.0 #2.0 # довжина ходу точки підвіски
    n=6.5 #10.0 # кількість гойдань (подвійних ходів) за хвилину
    fi=math.pi # початковий кут (почати з верхньої точки або з середини, якщо колона стабілізована)
    # ланки механізму (довжина, маса, момент інерції відносно осі z, яка проходить через центр мас)
    # Тип: СКД8-3-4000 (дезаксіальний)
    p1=dict(L=1.29,m=1982.0,Iz=635.0) # два кривошипи
    p2=dict(L=3.0,m=499.0,Iz=728.0) # два шатуна і траверса
    p3=dict(L=2.0,m=643.0,Iz=443.0) # заднє плече балансира
    p4=dict(L=2.29,m=911.0,Iz=1106.0) # переднє плече балансира
    p5=dict(L=0.789,m=4*750.0,Iz=0.0) # противаги (L - радіус центра мас) Приймаємо, що маса зосереджена в точці (Iz=0.0)
    # Увага! Значення за замовчуванням p5 уточнити функцією balancing
    dx=1.345 # горизонтальна відстань між осями опори балансира і тихохідного валу редуктора
    dy=3.012 # 3.0, 3.035 ? # вертикальна відстань між осями (розраховується)
    r=0.84 # радіус кривошипа (залежить від ходу)
    # примітка до табл. на с.42 Архіпов ?
    
    def __init__(self):
        self.A=0.5*self.stroke # амплітуда
        self.ns=self.n/60.0 # частота, Гц
        self.omega=2*math.pi*self.ns # кутова частота, рад/с
    
    def velAvg(self):
        """Середня швидкість, м/с"""
        return 2*self.stroke*self.n/60.0
    
    def motionX(self, t):
        """Описує закон руху точки підвіски.
        Повертає координату X точки підвіски в момент часу t"""
        #!!!!!уточнити
        #v=self.A*self.omega*math.cos(self.omega*t+self.fi) # швидкість
        x=self.A*math.sin(self.omega*t+self.fi) # переміщення
        return x
     
class PU():
    """Клас свердловинної штангової насосної установки"""
    suckerRodString=None
    well=None
    pump=None
    outPres=100000.0 # 100000.0 тиск на гирлі в НКТ
    outTubePres=500000.0 # 500000.0 тиск в затрубному просторі
    heightDynamic = None # глибина динамічного рівня (якщо None, то за замовчуванням)
    # зменшення heightDynamic зменшує тільки Fmax (див. Мищенко ф-ла 9.90)
    timeEnd=20.0 #кінець часу
    timeStep=0.1 # 0.01 крок часу
    
    def __init__(self):
        """Конструктор"""
        h=1500.0 # довжина колони
        sections=[(0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0),
                  (0.019, h/8.0, "сталь", 0.0)] # 0.2
        self.suckerRodString=SuckerRodString(PU=self,sections=sections) # колона штанг
        self.well=Well(PU=self) # свердловина
        self.pump=Pump(PU=self) # насос
        self.pumpingUnit=PumpingUnit() # верстат-гойдалка
        if self.heightDynamic==None:
            self.heightDynamic=self.suckerRodString.height()-100.0 # -100 орієнтовно
        r,m=self.balancing() # радіус і маса противаг
        self.pumpingUnit.p5['L']=r
        self.pumpingUnit.p5['m']=m
              
    def fluidVolume(self):
        """Об'єм стовба рідини над плунжером V=H*S"""
        return self.suckerRodString.height()*self.pump.pistonArea()
    
    def fluidWeight(self): # не враховується інерція рідини
        """Вага рідини над плунжером w=V*ro*g.
        Діє на плунжер під час ходу вверх. Направлена вниз"""
        return self.fluidVolume()*self.well.density*9.81
    
    def outPresForce(self):
        """Навантаження від тиску на гирлі.
        Діє на плунжер під час ходу вверх. Направлене вниз"""
        return self.pump.pistonArea()*self.outPres
    
    def outTubePresForce(self):
        """Навантаження від тиску в затрубному просторі.
        Діє на плунжер під час ходу вверх. Направлене вверх"""
        return self.pump.pistonArea()*self.outTubePres
    
    def outTubeFluidWeight(self):
        """Вага рідини в затрубному просторі.
        Діє на плунжер під час ходу вверх. Направлена вверх"""
        v=self.pump.pistonArea()*(self.suckerRodString.height()-self.heightDynamic) # обєм
        return v*self.well.density*9.81
    
    def tubePressureFrictLossForse(self,v):
        """Навантаження від втрати тиску від тертя потоку по довжині НКТ
        Діє на плунжер під час ходу вверх. Направлене вниз. Не для моделі Maplesim"""
        v=self.pump.volumeFlowRate(v)/self.well.sectionalArea() # середня швидкість в НКТ
        D=self.well.diametr
        L=self.suckerRodString.length()
        nu=self.well.viscosity
        ro=self.well.density
        Re=abs(v)*D/nu # число Рейнольдса
        fL=64/Re # для ламінарного
        fT=0.316/math.pow(Re, 1.0/4) # для турбулентного (формула Блазіуса)
        if Re<=2000: # якщо ламінарний режим
            f=fL
        elif Re>=4000: # якщо турбулентний режим
            f=fT
        else: # в іншому випадку
            f=fL+(fT-fL)*(Re-2000)/2000
        dP=(f*L/D)*v*v*ro/2 # втрати тиску
        return self.pump.pistonArea()*dP # завжди додатня
        
    def polishedRodForce(self):
        """Повертає кортеж (Fmin,Fmax) розрахований за наближеними формулами
        Ці формули вірні, коли динамічний рівень=довжині спуску
        Тобто розраховується найбільш небезпечний випадок"""
        F1=self.suckerRodString.weight() # вага штанг в рідині
        F2=self.fluidWeight() # вага рідини
        k=(self.pumpingUnit.stroke*self.pumpingUnit.n**2)/1790
        Fmin=F1*(1-k) # формула Миллса
        Fmax=(F1+F2)*(1+k) # формула Кемлера
        return (Fmin,Fmax)
    
    def sigmaPr(self,d):
        """Приведене напруження (Круман)"""
        A=math.pi*(d/2.0)**2 # площа січення верхньої штанги
        Fmin,Fmax=self.polishedRodForce()
        return (Fmax-0.5625*Fmin)/A
    
    def polishedRodForce2(self):
        """Повертає кортеж (Fmin,Fmax) розрахований за дуже наближеними формулами"""
        F1=self.suckerRodString.weight() # вага штанг в рідині
        F2=self.fluidWeight() # вага рідини
        Fmin=F1
        Fmax=F1+F2
        return (Fmin,Fmax)
    
    def balancing(self,Pmin=None,Pmax=None):
        """Параметри орієнтовного урівноважування [Архіпов]. Для верстатів типу СК і СКД
        Повертає радіус і сумарну масу противаг"""
        if Pmin==None and Pmax==None: # якщо не задано
            Pmin,Pmax=self.polishedRodForce() # розрахувати за наближеними формулами
        Mp=750.0 # маса 1 противаги
        A0=52000.0 # робота сил ваги неурівноважених частин вестата
        S=self.pumpingUnit.stroke
        K_list=[2,4,6,8] # список кількості противаг
        A=Pmin*S/2+Pmax*S/2 # робота сил, що урівноважуються
        R_list=[(A-A0)/(2*9.81*Mp*k) for k in K_list] # список радіусів противаг
        #print "Rp=", R_list
        Rmax=self.pumpingUnit.p1['L']-0.2 # максимальний радіус противаги ?уточнити
        for i,r in enumerate(R_list):
            if r<Rmax: break # вибирають радіус з мінімальною кількістю противаг
        return r, K_list[i]*Mp
    
    def info(self): # друкує інформацію
        s=""
        s+="вага штанг в рідині %d Н\n"%self.suckerRodString.weight()
        s+="вага рідини %d Н\n"%self.fluidWeight()
        s+="дуже наближено Fmin=%d Fmax=%d\n"%self.polishedRodForce2()
        s+="наближено Fmin=%d Fmax=%d\n"%self.polishedRodForce()
        s+="приведене напруження %d\n"%self.sigmaPr(d=self.suckerRodString.items[0].diametr)
        s+="діаметр плунжера %f\n"%self.pump.diametr
        s+="довжина колони %d\n"%self.suckerRodString.length()
        s+="жорсткість колони %d Н/м\n"%self.suckerRodString.stiffness()
        s+="маса колони %d кг\n"%self.suckerRodString.mass()
        n0=math.sqrt(self.suckerRodString.stiffness()/self.suckerRodString.mass()) # власна частота стержня
        s+="власна частота колони %f рад/с\n"%n0
        s+="диамічний рівень %d м\n"%self.heightDynamic
        s+="вага рідини в затр просторі %f Н\n"%self.outTubeFluidWeight()
        s+="еквівалентний коефіцієнт опору матеріалу %f \n"%self.suckerRodString.damping()
        s+="еквівалентний коефіцієнт опору матеріалу секції 1 %f \n"%self.suckerRodString.items[0].materialRodResistance()
        # зі збільшенням довжини секції прямує до 0 (як і жорсткість)
        s+="жорсткість секції 1 %f \n"%self.suckerRodString.items[0].stiffness()
        s+="маса секції 1 %f \n"%self.suckerRodString.items[0].mass()
        
        s+="маса НКТ %f \n"%self.well.mass()
        s+="жорсткість НКТ %f \n"%self.well.stiffness()
        s+="екв. коеф. оп. матер. НКТ %f \n"%self.well.materialRodResistance()
        print s 

if __name__ == '__main__':
    pu=PU()
    pu.info()
#     A=math.pi*(0.019/2.0)**2 # площа січення верхньої штанги
#     Fmin,Fmax=7000,30000
#     print (Fmax-0.5625*Fmin)/A