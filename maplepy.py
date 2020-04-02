# -*- coding: CP1251 -*-

class MapleInterface(object):
    """Інтерфейс до Maple"""
    mapleExePath=r"c:\Program Files\Maple 18\bin.win\cmaple.exe"
    codeTemplate="" # шаблон коду Maple
    def getCode(self):
        """Повертає код Maple"""
        return ""
    def runMaple(self):
        """Виконує код Maple"""
        import subprocess,tempfile,os
        tf=tempfile.NamedTemporaryFile(delete=False, suffix='.mpl') # тимчасовий файл .mpl
        tf.write(self.getCode()) # записати код
        tf.close()
        subprocess.Popen(self.mapleExePath+' '+tf.name).communicate() # виконати код Maple
        os.unlink(tf.name)# видалити тимчасовий файл    

class MapleInterface4Maplesim(MapleInterface):
    """Інтерфейс до Maplesim"""
    path="c:/9" # шлях до файлів моделі та результатів (розділювач тільки /)
    filenameMaplesim="testAPI.msim" # файл моделі Maplesim
    # приклад доступу до параметру SD1_c підсистеми n: n_SD1_c
    paramDictMaplesim={"SD1_c":1, "SD1_d":0.1} # словник параметрів моделі Maplesim
    execute=False # чи виконувати відразу симуляцію
    resultCSVfile="testAPI.csv" # файл результатів CSV
    # шаблони коду Maple
    codeTemplate1=r"""A :=MapleSim:-LinkModel('filename'="{filenameMaplesim}"):
A:-SetParameters([{setParams}]):
"""
    codeTemplate2=r"""simData:=A:-Simulate(output = datapoint):
ExportMatrix("{resultCSVfile}", simData, target=delimited, delimiter=";"):
"""
    def paramString(self,paramDictMaplesim):
        """Повертає рядок параметрів MapleSim"""
        paramList=[k+" = "+str(v) for k,v in paramDictMaplesim.iteritems()] # перетворити в список
        return ", ".join(paramList) # об'єднати в рядок
    
    def getCode(self):
        """Повертає код Maple"""
        setParams=self.paramString(self.paramDictMaplesim)
        filenameMaplesim=self.path+'/'+self.filenameMaplesim
        if self.execute: # якщо відразу виконувати симуляцію
            allCodeTemplate=self.codeTemplate1+self.codeTemplate2
            resultCSVfile=self.path+'/'+self.resultCSVfile
            code=allCodeTemplate.format(filenameMaplesim=filenameMaplesim, setParams=setParams, resultCSVfile=resultCSVfile)
        else:
            code=self.codeTemplate1.format(filenameMaplesim=filenameMaplesim, setParams=setParams)
        return code
        
    def readCSVfile(self,nc=2):
        """Читає файл результатів CSV та повертає список з результатами
        у вигляді [[t1,x1,y1,...],[t2,x2,y2,...],...]
        nc - кількість змінних (мінімум 2)"""
        import csv
        resultCSVfile=self.path+'/'+self.resultCSVfile
        csv_file=open(resultCSVfile, "rb")
        reader=csv.reader(csv_file,delimiter = ';')
        resultList=[]
        for row in reader: # для кожного рядка
            oneResult=[] # список результатів з одного рядка
            for i in range(nc): # для кожної змінної
                oneResult.append(float(row[i])) # добавити значення змінної
            resultList.append(oneResult) # добавити список результатів з одного рядка
        csv_file.close()
        return resultList

# class MapleInterface4Maplesim2(MapleInterface4Maplesim):
#     """Інтерфейс до Maplesim. Для моделей з підсистемами"""
#     paramDictMaplesim={'SD1_c':2} # параметри моделі
#     paramDictMaplesim2={'n':{'SD1_d':0.2}} # параметри підсистем моделі
#     CodeTemplate3="""A:-SetSubsystemName("{subSystemName}"):
# A:-SetParameters([{setParams}]):
# """
# 
#     def getSubCode(self):
#         """Повертає код Maple для вводу параметрів підсистем"""
#         codeList=[]
#         for sname in self.paramDictMaplesim2:
#             setParams=self.paramString(self.paramDictMaplesim2[sname])
#             subCode=self.CodeTemplate3.format(subSystemName=sname,setParams=setParams)
#             codeList.append(subCode)
#         return '\n'.join(codeList)
#     
#     def getCode(self):
#         """Повертає код Maple"""
#         setParams=self.paramString(self.paramDictMaplesim)
#         filenameMaplesim=self.path+'/'+self.filenameMaplesim
#         if self.execute: # якщо відразу виконувати симуляцію
#             allCodeTemplate=self.codeTemplate1+self.getSubCode()+self.codeTemplate2
#             resultCSVfile=self.path+'/'+self.resultCSVfile
#             code=allCodeTemplate.format(filenameMaplesim=filenameMaplesim, setParams=setParams, resultCSVfile=resultCSVfile)
#         else:
#             allCodeTemplate=self.codeTemplate1+self.getSubCode()
#             code=allCodeTemplate.format(filenameMaplesim=filenameMaplesim, setParams=setParams)
#         return code

if __name__ == '__main__':        
    ms=MapleInterface4Maplesim()
    print ms.getCode()
    #ms.runMaple()
    #print ms.readCSVfile(3)