# -*- coding: CP1251 -*-

class MapleInterface(object):
    """��������� �� Maple"""
    mapleExePath=r"c:\Program Files\Maple 18\bin.win\cmaple.exe"
    codeTemplate="" # ������ ���� Maple
    def getCode(self):
        """������� ��� Maple"""
        return ""
    def runMaple(self):
        """������ ��� Maple"""
        import subprocess,tempfile,os
        tf=tempfile.NamedTemporaryFile(delete=False, suffix='.mpl') # ���������� ���� .mpl
        tf.write(self.getCode()) # �������� ���
        tf.close()
        subprocess.Popen(self.mapleExePath+' '+tf.name).communicate() # �������� ��� Maple
        os.unlink(tf.name)# �������� ���������� ����    

class MapleInterface4Maplesim(MapleInterface):
    """��������� �� Maplesim"""
    path="c:/9" # ���� �� ����� ����� �� ���������� (��������� ����� /)
    filenameMaplesim="testAPI.msim" # ���� ����� Maplesim
    # ������� ������� �� ��������� SD1_c ��������� n: n_SD1_c
    paramDictMaplesim={"SD1_c":1, "SD1_d":0.1} # ������� ��������� ����� Maplesim
    execute=False # �� ���������� ������ ���������
    resultCSVfile="testAPI.csv" # ���� ���������� CSV
    # ������� ���� Maple
    codeTemplate1=r"""A :=MapleSim:-LinkModel('filename'="{filenameMaplesim}"):
A:-SetParameters([{setParams}]):
"""
    codeTemplate2=r"""simData:=A:-Simulate(output = datapoint):
ExportMatrix("{resultCSVfile}", simData, target=delimited, delimiter=";"):
"""
    def paramString(self,paramDictMaplesim):
        """������� ����� ��������� MapleSim"""
        paramList=[k+" = "+str(v) for k,v in paramDictMaplesim.iteritems()] # ����������� � ������
        return ", ".join(paramList) # ��'������ � �����
    
    def getCode(self):
        """������� ��� Maple"""
        setParams=self.paramString(self.paramDictMaplesim)
        filenameMaplesim=self.path+'/'+self.filenameMaplesim
        if self.execute: # ���� ������ ���������� ���������
            allCodeTemplate=self.codeTemplate1+self.codeTemplate2
            resultCSVfile=self.path+'/'+self.resultCSVfile
            code=allCodeTemplate.format(filenameMaplesim=filenameMaplesim, setParams=setParams, resultCSVfile=resultCSVfile)
        else:
            code=self.codeTemplate1.format(filenameMaplesim=filenameMaplesim, setParams=setParams)
        return code
        
    def readCSVfile(self,nc=2):
        """���� ���� ���������� CSV �� ������� ������ � ������������
        � ������ [[t1,x1,y1,...],[t2,x2,y2,...],...]
        nc - ������� ������ (����� 2)"""
        import csv
        resultCSVfile=self.path+'/'+self.resultCSVfile
        csv_file=open(resultCSVfile, "rb")
        reader=csv.reader(csv_file,delimiter = ';')
        resultList=[]
        for row in reader: # ��� ������� �����
            oneResult=[] # ������ ���������� � ������ �����
            for i in range(nc): # ��� ����� �����
                oneResult.append(float(row[i])) # �������� �������� �����
            resultList.append(oneResult) # �������� ������ ���������� � ������ �����
        csv_file.close()
        return resultList

# class MapleInterface4Maplesim2(MapleInterface4Maplesim):
#     """��������� �� Maplesim. ��� ������� � �����������"""
#     paramDictMaplesim={'SD1_c':2} # ��������� �����
#     paramDictMaplesim2={'n':{'SD1_d':0.2}} # ��������� �������� �����
#     CodeTemplate3="""A:-SetSubsystemName("{subSystemName}"):
# A:-SetParameters([{setParams}]):
# """
# 
#     def getSubCode(self):
#         """������� ��� Maple ��� ����� ��������� ��������"""
#         codeList=[]
#         for sname in self.paramDictMaplesim2:
#             setParams=self.paramString(self.paramDictMaplesim2[sname])
#             subCode=self.CodeTemplate3.format(subSystemName=sname,setParams=setParams)
#             codeList.append(subCode)
#         return '\n'.join(codeList)
#     
#     def getCode(self):
#         """������� ��� Maple"""
#         setParams=self.paramString(self.paramDictMaplesim)
#         filenameMaplesim=self.path+'/'+self.filenameMaplesim
#         if self.execute: # ���� ������ ���������� ���������
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