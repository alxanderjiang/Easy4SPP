#文件名: brdc_proc.py
#广播星历解析与卫星位置运算库
#注: 本文件当前仅支持: RINEX-3版本GPS, BDS, GAL星历
from satpos import *
from sppp import *
#from sppp_multiGNSS import reconstruct_obs_mat,check_obs_mats

#宏定义北斗GEO卫星PRN码, 用于D1/D2星历下的卫星位置计算(D1/D2不提供卫星类型, 无法判断GEO)
BDS_GEO_PRNS=['C01','C02','C03','C04','C05','C59','C60','C61','C62','C63','C64','C65']
#宏定义全局变量用于承接定位精度因子
Q_spp=np.zeros((6,6))

def reconstruct_obs_mat(obs_mat):
    #函数: 重整有效观测数据字典(根据Epoch_OK标识)
    #输入: 有效观测数据字典
    #输出: 重整后的观测数据字典
    r_obsmat=[]
    for i in range(len(obs_mat)):
        if(obs_mat[i][0]['Epoch_OK'])!=0:
            continue
        else:
            r_obsmat.append(obs_mat[i])
    #返回观测数据字典
    return r_obsmat

def check_obs_mats(obs_mats):
    #函数: 校验多系统观测值匹配性
    #输入: 多系统观测值列表
    #输出: 观测值数据
    lens=[]
    for i in range(len(obs_mats)):
        lens.append(len(obs_mats[i]))
    #列表仅单值
    if(lens==1):
        print("Only one system observed")
        return True
    #列表有多系统观测
    for i in range(len(lens)-1):
        if(lens[i]!=lens[i+1]):
            print("Observations among systems not equal")
            return False
    #长度校验通过, 开始时间校验
    for i in range(len(obs_mats[0])):
        GPS_week=obs_mats[0][i][0]["GPSweek"]
        GPS_sec=obs_mats[0][i][0]["GPSsec"]
        for j in range(len(obs_mats)):
            check_week=obs_mats[j][i][0]['GPSweek']
            check_sec=obs_mats[j][i][0]['GPSsec']
            if(check_week!=GPS_week or check_sec!=GPS_sec):
                print("Observations among systems not in the same time period")
                return False
    
    #全部校验通过, 返回True
    return True

#函数: GPS星历读取
def BRDC2GPSEPH(filename):
    # 函数: 读取星历的GPS部分并保存为字典
    # 输入: RINEX格式的GPS星历
    # 输出: GPSK模型参数数组; B2a星历数组
    ion_params_GPSK=[]
    GPS_eph=[]
    data_split={}
    nav_line=0
    with open(filename,"r",errors='ignore') as f:
        lines=f.readlines()
        header_in=1
        for line in lines:
            #读取电离层参数
            if("IONOSPHERIC CORR" in line):
                #GPSK_A 
                line=line.replace('D','E')  #兼容D
                line=line.replace('d','e')  #兼容d
                ls=line.split()
                if("GPSA" in line):
                    ion_params_GPSK.append(float(ls[1]))
                    ion_params_GPSK.append(float(ls[2]))
                    ion_params_GPSK.append(float(ls[3]))
                    ion_params_GPSK.append(float(ls[4]))
                #BDSK_B
                if("GPSB" in line):
                    ion_params_GPSK.append(float(ls[1]))
                    ion_params_GPSK.append(float(ls[2]))
                    ion_params_GPSK.append(float(ls[3]))
                    ion_params_GPSK.append(float(ls[4]))
            if("END OF HEADER" in line):
                header_in=0                 #文件头读取终止
            if(not header_in):
                line=line.replace('D','E')  #兼容D
                line=line.replace('d','e')  #兼容d
                #读到GPSS星历, 处理第一行
                if(line[0]=='G'):
                    #开始读取一个新的GPS星历数据块
                    nav_line=1
                    prn=line[0:3]#卫星PRN
                    toc_y=int(line[3:3+5])#Toc年
                    toc_m=int(line[3+5:3+5+3])#Toc月
                    toc_d=int(line[3+5+3:3+5+3+3])#Toc日
                    toc_h=int(line[3+5+3+3:3+5+3+3+3])#Toc时
                    toc_min=int(line[3+5+3+3+3:3+5+3+3+3+3])#Toc分
                    toc_sec=float(line[3+5+3+3+3+3:3+5+3+3+3+3+3])#Toc秒
                    toc=epoch2time(COMMTIME(toc_y,toc_m,toc_d,toc_h,toc_min,toc_sec))#卫星钟时间(BDT)
                    a0=float(line[3+5+3+3+3+3+3:3+5+3+3+3+3+3+19])
                    a1=float(line[3+5+3+3+3+3+3+19:3+5+3+3+3+3+3+19+19])
                    a2=float(line[3+5+3+3+3+3+3+19+19:3+5+3+3+3+3+3+19+19+19])
                    data_split['prn']=prn
                    data_split['toc']=toc
                    data_split['a0']=a0
                    data_split['a1']=a1
                    data_split['a2']=a2
                    nav_line+=1
                    continue
                if(nav_line==2):
                    #处理星历第二行
                    IODE=line[4:4+19]
                    crs=line[4+19:4+19+19]
                    delta_n0=line[4+19+19:4+19+19+19]
                    M0=line[4+19+19+19:]
                    data_split['IODE']=round(float(IODE))
                    data_split['crs']=float(crs)
                    data_split['delta_n0']=float(delta_n0)#临时
                    data_split['M0']=float(M0)#临时
                    nav_line+=1
                    continue
                if(nav_line==3):
                    cuc=line[4:4+19]
                    e=line[4+19:4+19+19]
                    cus=line[4+19+19:4+19+19+19]
                    sqrtA=line[4+19+19+19:]
                    data_split['cuc']=float(cuc)
                    data_split['e']=float(e)
                    data_split['cus']=float(cus)
                    data_split['sqrtA']=float(sqrtA)
                    nav_line+=1
                    continue
                if(nav_line==4):
                    toe=line[4:4+19]            #星历参考时刻(BDT)
                    cic=line[4+19:4+19+19]
                    OMEGA0=line[4+19+19:4+19+19+19]
                    cis=line[4+19+19+19:]
                    data_split['toe']=float(toe)
                    data_split['cic']=float(cic)
                    data_split['OMEGA0']=float(OMEGA0)#临时
                    data_split['cis']=float(cis)
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==5):
                    i0=line[4:4+19]
                    crc=line[4+19:4+19+19]
                    omega=line[4+19+19:4+19+19+19]
                    OMEGA_DOT=line[4+19+19+19:]
                    data_split['i0']=float(i0)#临时
                    data_split['crc']=float(crc)
                    data_split['omega']=float(omega)#临时
                    data_split['OMEGA_DOT']=float(OMEGA_DOT)#临时
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==6):
                    IDOT=line[4:4+19]
                    Codes_on_L2=line[4+19:4+19+19]#电文来源
                    GPSweek=line[4+19+19:4+19+19+19]#GPS周
                    A_DOT=line[4+19+19+19:]#长半轴变化率
                    data_split['IDOT']=float(IDOT)#临时
                    data_split['Codes_on_L2']=round(float(Codes_on_L2))#电文来源
                    data_split['week']=round(float(GPSweek))    
                    data_split['A_DOT']=float(A_DOT)
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==7):
                    SV_accuracy=line[4:4+19]
                    Health=line[4+19:4+19+19]
                    TGD=line[4+19+19:4+19+19+19]
                    IODC=line[4+19+19+19:]
                    data_split['SV_accuracy']=round(float(SV_accuracy))
                    data_split['Health']=round(float(Health))
                    data_split['TGD']=float(TGD)
                    data_split['IODC']=round(float(IODC))
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==8):
                    HOW=line[4:4+19]
                    BNK=line[4+19:4+19+19]
                    spare=line[4+19+19:4+19+19+19]
                    spare=line[4+19+19+19:]
                    data_split['HOW']=round(float(HOW))
                    data_split['BNK']=float(BNK)
                    #星历数据读取完毕, 开始计算卫星类型相关变量
                    data_split['A']=data_split['sqrtA']**2
                    #保存星历
                    GPS_eph.append(data_split.copy())
                    #重置星历块
                    data_split={}
                    nav_line=0
                    continue
    return ion_params_GPSK,GPS_eph

#函数: GPS星历读取
def BRDC2GALEPH(filename,target_Data=517):
    # 函数: 读取星历的GAL部分并保存为字典
    # 输入: RINEX格式的GAL星历, GAL星历类型(默认INAV/FNAV混合, E1/E5b组合钟差)
    # 输出: GAL电离层模型Nequick-G参数数组; GAL星历数组
    ion_params_GAL=[]
    GAL_eph=[]
    data_split={}
    nav_line=0
    with open(filename,"r",errors='ignore') as f:
        lines=f.readlines()
        header_in=1
        for line in lines:
            #读取电离层参数
            if("IONOSPHERIC CORR" in line):
                #GAL广播电离层模型发播系数
                line=line.replace('D','E')  #兼容D
                line=line.replace('d','e')  #兼容d
                ls=line.split()
                if("GAL" in line):
                    ion_params_GAL.append(float(ls[1]))
                    ion_params_GAL.append(float(ls[2]))
                    ion_params_GAL.append(float(ls[3]))
            if("END OF HEADER" in line):
                header_in=0                 #文件头读取终止
            if(not header_in):
                line=line.replace('D','E')  #兼容D
                line=line.replace('d','e')  #兼容d
                #读到GPSS星历, 处理第一行
                if(line[0]=='E'):
                    #开始读取一个新的GPS星历数据块
                    nav_line=1
                    prn=line[0:3]#卫星PRN
                    toc_y=int(line[3:3+5])#Toc年
                    toc_m=int(line[3+5:3+5+3])#Toc月
                    toc_d=int(line[3+5+3:3+5+3+3])#Toc日
                    toc_h=int(line[3+5+3+3:3+5+3+3+3])#Toc时
                    toc_min=int(line[3+5+3+3+3:3+5+3+3+3+3])#Toc分
                    toc_sec=float(line[3+5+3+3+3+3:3+5+3+3+3+3+3])#Toc秒
                    toc=epoch2time(COMMTIME(toc_y,toc_m,toc_d,toc_h,toc_min,toc_sec))#卫星钟时间(BDT)
                    a0=float(line[3+5+3+3+3+3+3:3+5+3+3+3+3+3+19])
                    a1=float(line[3+5+3+3+3+3+3+19:3+5+3+3+3+3+3+19+19])
                    a2=float(line[3+5+3+3+3+3+3+19+19:3+5+3+3+3+3+3+19+19+19])
                    data_split['prn']=prn
                    data_split['toc']=toc
                    data_split['a0']=a0
                    data_split['a1']=a1
                    data_split['a2']=a2
                    nav_line+=1
                    continue
                if(nav_line==2):
                    #处理星历第二行
                    IODnav=line[4:4+19]
                    crs=line[4+19:4+19+19]
                    delta_n0=line[4+19+19:4+19+19+19]
                    M0=line[4+19+19+19:]
                    data_split['IODnav']=round(float(IODnav))
                    data_split['crs']=float(crs)
                    data_split['delta_n0']=float(delta_n0)#临时
                    data_split['M0']=float(M0)#临时
                    nav_line+=1
                    continue
                if(nav_line==3):
                    cuc=line[4:4+19]
                    e=line[4+19:4+19+19]
                    cus=line[4+19+19:4+19+19+19]
                    sqrtA=line[4+19+19+19:]
                    data_split['cuc']=float(cuc)
                    data_split['e']=float(e)
                    data_split['cus']=float(cus)
                    data_split['sqrtA']=float(sqrtA)
                    nav_line+=1
                    continue
                if(nav_line==4):
                    toe=line[4:4+19]            #星历参考时刻(BDT)
                    cic=line[4+19:4+19+19]
                    OMEGA0=line[4+19+19:4+19+19+19]
                    cis=line[4+19+19+19:]
                    data_split['toe']=float(toe)
                    data_split['cic']=float(cic)
                    data_split['OMEGA0']=float(OMEGA0)#临时
                    data_split['cis']=float(cis)
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==5):
                    i0=line[4:4+19]
                    crc=line[4+19:4+19+19]
                    omega=line[4+19+19:4+19+19+19]
                    OMEGA_DOT=line[4+19+19+19:]
                    data_split['i0']=float(i0)#临时
                    data_split['crc']=float(crc)
                    data_split['omega']=float(omega)#临时
                    data_split['OMEGA_DOT']=float(OMEGA_DOT)#临时
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==6):
                    IDOT=line[4:4+19]
                    Data=line[4+19:4+19+19]#电文来源
                    GPSweek=line[4+19+19:4+19+19+19]#GPS周
                    spare=line[4+19+19+19:]#长半轴变化率
                    data_split['IDOT']=float(IDOT)#临时
                    data_split['Data']=round(float(Data))#电文来源
                    data_split['week']=round(float(GPSweek))    
                    data_split['A_DOT']=0.0
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==7):
                    SISA=line[4:4+19]
                    Health=line[4+19:4+19+19]
                    BGD1=line[4+19+19:4+19+19+19]
                    BGD2=line[4+19+19+19:]
                    data_split['SISA']=round(float(SISA))
                    data_split['Health']=round(float(Health))
                    data_split['BGD1']=float(BGD1)
                    data_split['BGD2']=float(BGD2)
                    #print(ls)
                    nav_line+=1
                    continue
                if(nav_line==8):
                    TOW=line[4:4+19]
                    spare=line[4+19:4+19+19]
                    spare=line[4+19+19:4+19+19+19]
                    spare=line[4+19+19+19:]
                    data_split['TOW']=float(TOW)
                    #星历数据读取完毕, 开始计算卫星类型相关变量
                    data_split['A']=data_split['sqrtA']**2
                    if(data_split['Data']!=target_Data):
                        #重置星历块, 星历类型不匹配
                        data_split={}
                        nav_line=0
                        continue
                    #保存星历
                    GAL_eph.append(data_split.copy())
                    #重置星历块
                    data_split={}
                    nav_line=0
                    continue
    return ion_params_GAL,GAL_eph

#函数: BDS星历(D1/D2)读取
def BRDC2BDSEPH(filename):
    # 函数: 读取D1/D2电文并保存为字典
    # 输入: RINEX格式的BDS D1/D2电文
    # 输出: BDSK模型参数数组; D1/D2星历数组
    ion_params_BDSK={}
    BDS_EPH=[]
    data_split={}
    cnav_line=0
    with open(filename,"r",errors='ignore') as f:
        lines=f.readlines()
        header_in=1
        for line in lines:
            
            #读取电离层参数
            if("IONOSPHERIC CORR" in line):
                #BDSK_A
                if("BDSA" in line):
                    line=line.replace('D','E')  #兼容D
                    line=line.replace('d','e')  #兼容d
                    ls=line.split()
                    ion_params_BDSK[ls[-4]]=[]
                    ion_params_BDSK[ls[-4]].append(float(ls[1]))
                    ion_params_BDSK[ls[-4]].append(float(ls[2]))
                    ion_params_BDSK[ls[-4]].append(float(ls[3]))
                    ion_params_BDSK[ls[-4]].append(float(ls[4]))
                #BDSK_B
                if("BDSB" in line):
                    line=line.replace('D','E')  #兼容D
                    line=line.replace('d','e')  #兼容d
                    ls=line.split()
                    ion_params_BDSK[ls[-4]].append(float(ls[1]))
                    ion_params_BDSK[ls[-4]].append(float(ls[2]))
                    ion_params_BDSK[ls[-4]].append(float(ls[3]))
                    ion_params_BDSK[ls[-4]].append(float(ls[4]))
            if("END OF HEADER" in line):
                header_in=0                 #文件头读取终止
            if(not header_in):
                #读到BDS星历, 处理第一行
                line=line.replace('D','E')  #兼容D
                line=line.replace('d','e')  #兼容d
                if(line[0]=='C'):
                    #开始读取一个新的北斗星历数据块
                    cnav_line=1
                    prn=line[0:3]#卫星PRN
                    toc_y=int(line[3:3+5])#Toc年
                    toc_m=int(line[3+5:3+5+3])#Toc月
                    toc_d=int(line[3+5+3:3+5+3+3])#Toc日
                    toc_h=int(line[3+5+3+3:3+5+3+3+3])#Toc时
                    toc_min=int(line[3+5+3+3+3:3+5+3+3+3+3])#Toc分
                    toc_sec=float(line[3+5+3+3+3+3:3+5+3+3+3+3+3])#Toc秒
                    toc=epoch2time(COMMTIME(toc_y,toc_m,toc_d,toc_h,toc_min,toc_sec))#卫星钟时间(BDT)
                    a0=float(line[3+5+3+3+3+3+3:3+5+3+3+3+3+3+19])
                    a1=float(line[3+5+3+3+3+3+3+19:3+5+3+3+3+3+3+19+19])
                    a2=float(line[3+5+3+3+3+3+3+19+19:3+5+3+3+3+3+3+19+19+19])
                    data_split['prn']=prn
                    data_split['toc']=toc
                    data_split['a0']=a0
                    data_split['a1']=a1
                    data_split['a2']=a2
                    cnav_line+=1
                    continue
                #处理星历第二行
                if(cnav_line==2):
                    AODE=line[4:4+19]
                    crs=line[4+19:4+19+19]
                    delta_n0=line[4+19+19:4+19+19+19]
                    M0=line[4+19+19+19:]
                    data_split['AODE']=float(AODE)
                    data_split['crs']=float(crs)
                    data_split['delta_n0']=float(delta_n0)#临时
                    data_split['M0']=float(M0)#临时
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第三行
                if(cnav_line==3):
                    cuc=line[4:4+19]
                    e=line[4+19:4+19+19]
                    cus=line[4+19+19:4+19+19+19]
                    sqrtA=line[4+19+19+19:]
                    data_split['cuc']=float(cuc)
                    data_split['e']=float(e)
                    data_split['cus']=float(cus)
                    data_split['sqrtA']=float(sqrtA)
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第四行
                if(cnav_line==4):
                    toe=line[4:4+19]            #星历参考时刻(BDT)
                    cic=line[4+19:4+19+19]
                    OMEGA0=line[4+19+19:4+19+19+19]
                    cis=line[4+19+19+19:]
                    data_split['toe']=float(toe)
                    data_split['cic']=float(cic)
                    data_split['OMEGA0']=float(OMEGA0)#临时
                    data_split['cis']=float(cis)
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第五行
                if(cnav_line==5):
                    i0=line[4:4+19]
                    crc=line[4+19:4+19+19]
                    omega=line[4+19+19:4+19+19+19]
                    OMEGA_DOT=line[4+19+19+19:]
                    data_split['i0']=float(i0)#临时
                    data_split['crc']=float(crc)
                    data_split['omega']=float(omega)#临时
                    data_split['OMEGA_DOT']=float(OMEGA_DOT)#临时
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第六行
                if(cnav_line==6):
                    IDOT=line[4:4+19]
                    Data=line[4+19:4+19+19]#电文来源
                    BDSweek=line[4+19+19:4+19+19+19]#北斗周
                    data_split['IDOT']=float(IDOT)#临时
                    data_split['Data']=round(float(Data))#电文来源1=B1C CNAV1; 2=B2b CNAV1; 4=B2a CNAV2
                    data_split['week']=round(float(BDSweek))+1356    #BDS周转GPS周, 适配旧有函数time2gpst
                    data_split['A_DOT']=float(0.0)#D1/D2电文无长半轴变化率
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第七行
                if(cnav_line==7):
                    SV_accuracy=line[4:4+19]
                    Health=line[4+19:4+19+19]
                    TGD1=line[4+19+19:4+19+19+19]
                    TGD2=line[4+19+19+19:]
                    data_split['SV_accuracy']=round(float(SV_accuracy))
                    data_split['Health']=round(float(Health))
                    data_split['TGD1']=float(TGD1)
                    data_split['TGD2']=float(TGD2)
                    #print(ls)
                    cnav_line+=1
                    continue
                #处理星历第八行
                if(cnav_line==8):
                    SOW=line[4:4+19]
                    AODC=line[4+19:4+19+19]
                    delta_n0_DOT=line[4+19+19:4+19+19+19]
                    Sat_Type=line[4+19+19+19:]
                    data_split['SOW']=float(SOW)
                    data_split['AODC']=round(float(AODC))
                    data_split['delta_n0_DOT']=float(0.0)
                    #星历数据读取完毕, 开始计算卫星类型相关变量
                    #GEO
                    data_split['A']=data_split['sqrtA']**2
                    #保存星历
                    BDS_EPH.append(data_split.copy())
                    #重置星历块
                    data_split={}
                    #print(ls)
                    cnav_line=0
                    continue    
    return ion_params_BDSK,BDS_EPH

#GPS卫星位置计算
def GPSEPH2Satpos(GPS_eph,rt,prn,iodc=None,rho=0, eph_pre_h=0):
    #函数: 根据GPS-LNAV星历计算卫星位置
    #输入: 导航电文, 观测时间(GPST), 卫星PRN, 可选参数: 星历编号
    #输出: 卫星位置, 速度, 钟差, 钟速
    f1,f2=1575.42e6,1227.60e6#GPS L1、L2频点中心频率
    GM_GPS=3.986004418e14   #WGS-84地球引力常数
    OMGE_GPS= 7.2921151467e-5  #WGS-84地球自转速度
    #首先确定目标星历
    lnav_prns=[t['prn'] for t in GPS_eph]
    id=np.where(np.array(lnav_prns)==prn)
    satellite={}
    #GPST
    rt_gps=rt
    min_dt=9999999
    for i in id[0]:
        if(GPS_eph[i]['toc']<=rt_gps+eph_pre_h and abs(rt_gps-GPS_eph[i]['toc'])<min_dt):
            satellite=GPS_eph[i]
            min_dt=abs(rt_gps-GPS_eph[i]['toc'])
    #星历存在性与有效性校验
    if(satellite=={} or satellite['Health']!=0):
        #print(rt,"No",prn)
        return False
    ##开始计算卫星位置
    toe=gpst2time(satellite['week'],satellite['toe'])
    toc=satellite['toc']
    #计算卫星轨道长半轴
    A=satellite['A']
    #计算参考时刻的平均运动角速度
    n0 =sqrt(GM_GPS/A/A/A)
    #计算相对星历参考历元的时间
    tk=rt-rho/clight-toe
    if(tk>302400):
        tk=tk-604800
    elif(tk<-302400):
        tk=tk+604800
    else:
        tk=tk
    #计算长半轴
    Ak=A+satellite['A_DOT']*tk
    #计算卫星平均运动角速度的偏差(GPS星历不提供二阶改正)
    dnA=satellite['delta_n0']#+0.5*satellite['delta_n0_DOT']*tk
    #计算改正后的卫星平均角速度
    n=n0+dnA
    #计算平近点角
    Mk=satellite['M0']+n*tk
    #迭代计算偏近点角
    Ek=Mk
    Ek1=0.0
    ek_count=0
    while(1):
        Ek1=Mk+satellite['e']*sin(Ek)
        ek_count+=1
        if(abs(Ek1-Ek)<1e-13 or ek_count>200):
            break  
        Ek=Ek1
    Ek=Ek1
    #计算真近点角
    vk = atan2(sqrt(1-satellite['e']*satellite['e'])*sin(Ek)/(1-satellite['e']*cos(Ek)), (cos(Ek)-satellite['e'])/(1-satellite['e']*cos(Ek)))
    #计算升交点角距Phik
    Phik=vk+satellite['omega']
    #计算二阶调和改正数
    deltauk = satellite['cus'] * sin(2 * Phik) + satellite['cuc'] * cos(2 * Phik)
    deltark = satellite['crs'] * sin(2 * Phik) + satellite['crc'] * cos(2 * Phik)
    deltaik = satellite['cis'] * sin(2 * Phik) + satellite['cic'] * cos(2 * Phik)
    #计算改进的升交角距
    uk=Phik+deltauk
    #计算改进的向径
    rk= Ak * (1- satellite['e'] * cos(Ek)) + deltark;   
    #计算改正的轨道倾角
    ik = satellite['i0'] + deltaik + satellite['IDOT'] * tk
    #计算卫星在轨道平面内的坐标
    xk1 = rk * cos(uk)
    yk1 = rk * sin(uk)
    #计算经过改正的升交点经度(MEO/IGSO)
    Omegak = satellite['OMEGA0'] + (satellite['OMEGA_DOT']- OMGE_GPS) * tk- OMGE_GPS *satellite['toe']
    #计算GPS-MEO卫星在WGS84-ECEF下的位置
    Xk=xk1*cos(Omegak)-yk1*cos(ik)*sin(Omegak)
    Yk=xk1*sin(Omegak)+yk1*cos(ik)*cos(Omegak)
    Zk=yk1*sin(ik)
    ##开始计算卫星速度
    #计算偏近点角的时间导数
    dotEk=n/(1-satellite['e']*cos(Ek))
    #计算真近点角的时间导数
    dotvk=sqrt((1+satellite['e'])/(1-satellite['e']))*cos(vk/2)*cos(vk/2)/cos(Ek/2)/cos(Ek/2)*dotEk
    #计算升交角距的时间导数
    dotuk=dotvk*(1+2*satellite['cus']*cos(2*Phik)-2*satellite['cuc']*sin(2*Phik))
    #计算向径的时间导数
    dotrk=dotEk*Ak*satellite['e']*sin(Ek)+2*dotvk*(satellite['crs']*cos(2*Phik)-satellite['crc']*sin(2*Phik))
    #计算轨道倾角的时间导数
    dotik=satellite['IDOT']+2*dotvk*(satellite['cis']*cos(2*Phik)-satellite['cic']*sin(2*Phik))
    #计算卫星轨道平面位置的时间导数
    dotxk = cos(uk) * dotrk- rk * sin(uk) * dotuk
    dotyk = sin(uk) * dotrk + rk * cos(uk) * dotuk
    #计算IGSO/MEO卫星在BDCS-ECEF下的速度
    dotR12=[[cos(Omegak),-sin(Omegak)*cos(ik),-xk1*sin(Omegak)-yk1*cos(Omegak)*cos(ik), yk1*sin(Omegak)*sin(ik)],
            [sin(Omegak), cos(Omegak)*cos(ik), xk1*cos(Omegak)-yk1*sin(Omegak)*cos(ik),-yk1*cos(Omegak)*sin(ik)],
            [0,sin(ik),0,Yk*cos(ik)]]
    tx=[[dotxk],[dotyk],[satellite['OMEGA_DOT']-OMGE_GPS],[dotik]]
    dxyz=np.array(dotR12).dot(np.array(tx))
    Fclk=-2*sqrt(GM_GPS)/clight/clight*satellite['e']*sqrt(Ak)*sin(Ek)#相对论效应
    tdts=satellite['a0']+satellite['a1']*(rt-toc)+satellite['a2']*((rt-toc)**2)+Fclk
    tdtss=satellite['a1']+2*satellite['a2']*(rt-toc)
    return [Xk,Yk,Zk,tdts,float(dxyz[0][0]),float(dxyz[1][0]),float(dxyz[2][0]),tdtss,satellite['TGD'],f1*f1/f2/f2*satellite['TGD'],satellite['IODC']]

#BDS卫星位置计算
def BDSEPH2SatPos(BDS_eph,rt,prn,iodc=None,rho=0):
    #函数: 根据BDS-D1/D2星历计算卫星位置
    #输入: 导航电文, 观测时间(GPST), 卫星PRN, 可选参数: 星历编号
    #输出: 卫星位置, 速度, 钟差, 钟速
    GM_BDS=3.986004418e14   #BDCS地球引力常数
    OMGE_BDS= 7.2921150e-5  #BDCS地球自转速度

    #首先确定目标星历
    cnav_prns=[t['prn'] for t in BDS_eph]
    id=np.where(np.array(cnav_prns)==prn)
    satellite={}
    #GPST转BDST
    rt_gps=rt
    rt=rt-14
    min_dt=9999999
    for i in id[0]:
        if(BDS_eph[i]['toc']<=rt_gps and rt_gps-BDS_eph[i]['toc']<min_dt):
            satellite=BDS_eph[i]
            min_dt=rt_gps-BDS_eph[i]['toc']
    #星历存在性与有效性校验
    if(satellite=={} or satellite['Health']!=0):
        #print(rt,"No",prn)
        return False
    ##开始计算卫星位置
    toe=gpst2time(satellite['week'],satellite['toe'])
    toc=satellite['toc']
    #计算卫星轨道长半轴
    A=satellite['A']
    #计算参考时刻的平均运动角速度
    n0 =sqrt(GM_BDS/A/A/A)
    #计算相对星历参考历元的时间
    tk=rt-rho/clight-toe#(satellite['a0']+satellite['a1']*(rt-toe)+satellite['a2']*((rt-toe)**2))+satellite['TGD']-toe
    if(tk>302400):
        tk=tk-604800
    elif(tk<-302400):
        tk=tk+604800
    else:
        tk=tk
    #计算长半轴
    Ak=A+satellite['A_DOT']*tk
    #计算卫星平均运动角速度的偏差
    dnA=satellite['delta_n0']+0.5*satellite['delta_n0_DOT']*tk
    #计算改正后的卫星平均角速度
    n=n0+dnA
    #计算平近点角
    Mk=satellite['M0']+n*tk
    #迭代计算偏近点角
    Ek=Mk
    Ek1=0.0
    ek_count=0
    while(1):
        Ek1=Mk+satellite['e']*sin(Ek)
        ek_count+=1
        if(abs(Ek1-Ek)<1e-13 or ek_count>200):
            break  
        Ek=Ek1
    Ek=Ek1
    #计算真近点角
    vk = atan2(sqrt(1-satellite['e']*satellite['e'])*sin(Ek)/(1-satellite['e']*cos(Ek)), (cos(Ek)-satellite['e'])/(1-satellite['e']*cos(Ek)))
    #计算升交点角距Phik
    Phik=vk+satellite['omega']
    #计算二阶调和改正数
    deltauk = satellite['cus'] * sin(2 * Phik) + satellite['cuc'] * cos(2 * Phik)
    deltark = satellite['crs'] * sin(2 * Phik) + satellite['crc'] * cos(2 * Phik)
    deltaik = satellite['cis'] * sin(2 * Phik) + satellite['cic'] * cos(2 * Phik)
    #计算改进的升交角距
    uk=Phik+deltauk
    #计算改进的向径
    rk= Ak * (1- satellite['e'] * cos(Ek)) + deltark;   
    #计算改正的轨道倾角
    ik = satellite['i0'] + deltaik + satellite['IDOT'] * tk
    #计算卫星在轨道平面内的坐标
    xk1 = rk * cos(uk)
    yk1 = rk * sin(uk)
    #计算经过改正的升交点经度(MEO/IGSO)
    Omegak = satellite['OMEGA0'] + (satellite['OMEGA_DOT']- OMGE_BDS) * tk- OMGE_BDS *satellite['toe']
    #计算经过改正的升交点经度(GEO)
    Omegak1= satellite['OMEGA0'] + (satellite['OMEGA_DOT']) * tk - OMGE_BDS * satellite['toe']
    #计算IGSO/MEO卫星在BDCS-ECEF下的位置
    Xk=xk1*cos(Omegak)-yk1*cos(ik)*sin(Omegak)
    Yk=xk1*sin(Omegak)+yk1*cos(ik)*cos(Omegak)
    Zk=yk1*sin(ik)
    #计算GEO卫星在BDCS-ECEF下的位置
    if(prn in BDS_GEO_PRNS):
        Ometk = OMGE_BDS * tk
        f = -5.0 / 180 * pi
        Rx = [[1,0,0],  
              [0,cos(f),sin(f)],  
              [0,-sin(f),cos(f)]]
        Rz = [[cos(Ometk),sin(Ometk),0],  
              [-sin(Ometk),cos(Ometk),0],
              [0,0,1]]
        xyzgk=[0.0,0.0,0.0]
        xyzgk[0] = xk1 * cos(Omegak1) - yk1 * cos(ik) * sin(Omegak1)
        xyzgk[1] = xk1 * sin(Omegak1) + yk1 * cos(ik) * cos(Omegak1)
        xyzgk[2] = yk1 * sin(ik)
        xyzr=(np.array(Rz).dot(np.array(Rx))).dot(np.array(xyzgk).reshape(3,1))
        Xk,Yk,Zk=xyzr[0][0],xyzr[1][0],xyzr[2][0]
    
    ##开始计算卫星速度
    #计算偏近点角的时间导数
    dotEk=n/(1-satellite['e']*cos(Ek))
    #计算真近点角的时间导数
    dotvk=sqrt((1+satellite['e'])/(1-satellite['e']))*cos(vk/2)*cos(vk/2)/cos(Ek/2)/cos(Ek/2)*dotEk
    #计算升交角距的时间导数
    dotuk=dotvk*(1+2*satellite['cus']*cos(2*Phik)-2*satellite['cuc']*sin(2*Phik))
    #计算向径的时间导数
    dotrk=dotEk*Ak*satellite['e']*sin(Ek)+2*dotvk*(satellite['crs']*cos(2*Phik)-satellite['crc']*sin(2*Phik))
    #计算轨道倾角的时间导数
    dotik=satellite['IDOT']+2*dotvk*(satellite['cis']*cos(2*Phik)-satellite['cic']*sin(2*Phik))
    #计算卫星轨道平面位置的时间导数
    dotxk = cos(uk) * dotrk- rk * sin(uk) * dotuk
    dotyk = sin(uk) * dotrk + rk * cos(uk) * dotuk
    #计算IGSO/MEO卫星在BDCS-ECEF下的速度
    dotR12=[[cos(Omegak),-sin(Omegak)*cos(ik),-xk1*sin(Omegak)-yk1*cos(Omegak)*cos(ik), yk1*sin(Omegak)*sin(ik)],
            [sin(Omegak), cos(Omegak)*cos(ik), xk1*cos(Omegak)-yk1*sin(Omegak)*cos(ik),-yk1*cos(Omegak)*sin(ik)],
            [0,sin(ik),0,Yk*cos(ik)]]
    tx=[[dotxk],[dotyk],[satellite['OMEGA_DOT']-OMGE_BDS],[dotik]]
    dxyz=np.array(dotR12).dot(np.array(tx))
    #计算GEO卫星在BDCS-ECEF下的速度
    if(prn in BDS_GEO_PRNS):
        Ometk = OMGE_BDS * tk
        f = -5.0 / 180 * pi
        dotR1=np.zeros((3,5),dtype=np.float64)
        tx=np.array([dotxk,dotyk,satellite['OMEGA_DOT'],dotik,OMGE_BDS]).reshape(5,1)
        dotR1[0][0] = cos(Ometk) * cos(Omegak1) + sin(Ometk) * cos(f) * sin(Omegak1)
        dotR1[0][1] = -cos(Ometk) * cos(ik) * sin(Omegak1) + sin(Ometk) * cos(f) * cos(ik) * cos(Omegak1) + sin(Ometk) * sin(f) * sin(ik)
        dotR1[0][2] = -xk1 * cos(Ometk) * sin(Omegak1) - yk1 * cos(ik) * cos(Ometk) * cos(Omegak1) + sin(Ometk) * cos(f) * xk1 * cos(Omegak1) - yk1 * cos(ik) * sin(Omegak1) * sin(Ometk) * cos(f)
        dotR1[0][3] = yk1 * sin(ik) * sin(Omegak1) * cos(Ometk) - yk1 * sin(Ometk) * cos(f) * sin(ik) * cos(Omegak1) + sin(Ometk) * sin(f) * yk1 * cos(ik)
        dotR1[0][4] = -sin(Ometk) * (xk1 * cos(Omegak1) - yk1 * cos(ik) * sin(Omegak1)) + cos(Ometk) * cos(f) * (xk1 * sin(Omegak1) + yk1 * cos(ik) * cos(Omegak1)) + cos(Ometk) * sin(f) * yk1 * sin(ik)
        
        dotR1[1][0] = -sin(Ometk) * cos(Omegak1) + cos(Ometk) * cos(f) * sin(Omegak1)
        dotR1[1][1] = sin(Ometk) * cos(ik) * sin(Omegak1) + cos(Ometk) * cos(f) * cos(ik) * cos(Omegak1) + cos(Ometk) * sin(f) * sin(ik)
        dotR1[1][2] = sin(Ometk) * xk1 * sin(Omegak1) + yk1 * sin(Ometk) * cos(ik) * cos(Omegak1) + cos(Ometk) * cos(f) * xk1 * cos(Omegak1) - yk1 * cos(Ometk) * cos(f) * cos(ik) * sin(Omegak1)
        dotR1[1][3] = -yk1 * sin(ik) * sin(Omegak1) * sin(Ometk) - cos(Ometk) * cos(f) * yk1 * sin(ik) * cos(Omegak1) + cos(Ometk) * sin(f) * yk1 * cos(ik)
        dotR1[1][4] = -cos(Ometk) * (xk1 * cos(Omegak1) - yk1 * cos(ik) * sin(Omegak1)) - sin(Ometk) * cos(f) * (xk1 * sin(Omegak1) + yk1 * cos(ik) * cos(Omegak1)) - sin(Ometk) * sin(f) * yk1 * sin(ik)
        
        dotR1[2][0] = -sin(f) * sin(Omegak1)
        dotR1[2][1] = -sin(f) * cos(ik) * cos(Omegak1) + cos(f) * sin(ik)
        dotR1[2][2] = -xk1 * sin(f) * cos(Omegak1) + yk1 * cos(ik) * sin(Omegak1) * sin(f)
        dotR1[2][3] = sin(f) * yk1 * sin(ik) * cos(Omegak1) + cos(f) * yk1 * cos(ik)
        dotR1[2][4] = 0

        dotxyzrr=dotR1.dot(tx)
        dxyz=dotxyzrr
    #
    #计算卫星钟差&钟速
    #tdts=satellite['a0']+satellite['a1']*(rt-toe)+satellite['a2']*((rt-toe)**2)
    Fclk=-2*sqrt(GM_BDS)/clight/clight*satellite['e']*sqrt(Ak)*sin(Ek)#相对论效应
    tdts=satellite['a0']+satellite['a1']*(rt-toc)+satellite['a2']*((rt-toc)**2)+Fclk
    tdtss=satellite['a1']+2*satellite['a2']*(rt-toc)
    return [Xk,Yk,Zk,tdts,float(dxyz[0][0]),float(dxyz[1][0]),float(dxyz[2][0]),tdtss,satellite['TGD1'],satellite['TGD2']]

#GAL卫星位置计算
def GALEPH2Satpos(GPS_eph,rt,prn,iodc=None,rho=0):
    #函数: 根据GPS-LNAV星历计算卫星位置
    #输入: 导航电文, 观测时间(GPST), 卫星PRN, 可选参数: 星历编号
    #输出: 卫星位置, 速度, 钟差, 钟速
    f1,f2=1575.42e6,1176.45e6#GAL E1、L5a频点中心频率
    GM_GAL=3.986004418e14   #WGS-84地球引力常数
    OMGE_GAL= 7.2921151467e-5  #WGS-84地球自转速度
    #首先确定目标星历
    lnav_prns=[t['prn'] for t in GPS_eph]
    id=np.where(np.array(lnav_prns)==prn)
    satellite={}
    #GPST
    rt_gps=rt
    min_dt=9999999
    for i in id[0]:
        if(GPS_eph[i]['toc']<=rt_gps and rt_gps-GPS_eph[i]['toc']<min_dt):
            satellite=GPS_eph[i]
            min_dt=rt_gps-GPS_eph[i]['toc']
    #星历存在性与有效性校验
    if(satellite=={} or satellite['Health']!=0):
        #print(rt,"No",prn)
        return False
    ##开始计算卫星位置
    toe=gpst2time(satellite['week'],satellite['toe'])
    toc=satellite['toc']
    #计算卫星轨道长半轴
    A=satellite['A']
    #计算参考时刻的平均运动角速度
    n0 =sqrt(GM_GAL/A/A/A)
    #计算相对星历参考历元的时间
    tk=rt-rho/clight-toe
    if(tk>302400):
        tk=tk-604800
    elif(tk<-302400):
        tk=tk+604800
    else:
        tk=tk
    #计算长半轴
    Ak=A+satellite['A_DOT']*tk
    #计算卫星平均运动角速度的偏差(GPS星历不提供二阶改正)
    dnA=satellite['delta_n0']#+0.5*satellite['delta_n0_DOT']*tk
    #计算改正后的卫星平均角速度
    n=n0+dnA
    #计算平近点角
    Mk=satellite['M0']+n*tk
    #迭代计算偏近点角
    Ek=Mk
    Ek1=0.0
    ek_count=0
    while(1):
        Ek1=Mk+satellite['e']*sin(Ek)
        ek_count+=1
        if(abs(Ek1-Ek)<1e-13 or ek_count>200):
            break  
        Ek=Ek1
    Ek=Ek1
    #计算真近点角
    vk = atan2(sqrt(1-satellite['e']*satellite['e'])*sin(Ek)/(1-satellite['e']*cos(Ek)), (cos(Ek)-satellite['e'])/(1-satellite['e']*cos(Ek)))
    #计算升交点角距Phik
    Phik=vk+satellite['omega']
    #计算二阶调和改正数
    deltauk = satellite['cus'] * sin(2 * Phik) + satellite['cuc'] * cos(2 * Phik)
    deltark = satellite['crs'] * sin(2 * Phik) + satellite['crc'] * cos(2 * Phik)
    deltaik = satellite['cis'] * sin(2 * Phik) + satellite['cic'] * cos(2 * Phik)
    #计算改进的升交角距
    uk=Phik+deltauk
    #计算改进的向径
    rk= Ak * (1- satellite['e'] * cos(Ek)) + deltark;   
    #计算改正的轨道倾角
    ik = satellite['i0'] + deltaik + satellite['IDOT'] * tk
    #计算卫星在轨道平面内的坐标
    xk1 = rk * cos(uk)
    yk1 = rk * sin(uk)
    #计算经过改正的升交点经度(MEO/IGSO)
    Omegak = satellite['OMEGA0'] + (satellite['OMEGA_DOT']- OMGE_GAL) * tk- OMGE_GAL *satellite['toe']
    #计算GPS-MEO卫星在WGS84-ECEF下的位置
    Xk=xk1*cos(Omegak)-yk1*cos(ik)*sin(Omegak)
    Yk=xk1*sin(Omegak)+yk1*cos(ik)*cos(Omegak)
    Zk=yk1*sin(ik)
    ##开始计算卫星速度
    #计算偏近点角的时间导数
    dotEk=n/(1-satellite['e']*cos(Ek))
    #计算真近点角的时间导数
    dotvk=sqrt((1+satellite['e'])/(1-satellite['e']))*cos(vk/2)*cos(vk/2)/cos(Ek/2)/cos(Ek/2)*dotEk
    #计算升交角距的时间导数
    dotuk=dotvk*(1+2*satellite['cus']*cos(2*Phik)-2*satellite['cuc']*sin(2*Phik))
    #计算向径的时间导数
    dotrk=dotEk*Ak*satellite['e']*sin(Ek)+2*dotvk*(satellite['crs']*cos(2*Phik)-satellite['crc']*sin(2*Phik))
    #计算轨道倾角的时间导数
    dotik=satellite['IDOT']+2*dotvk*(satellite['cis']*cos(2*Phik)-satellite['cic']*sin(2*Phik))
    #计算卫星轨道平面位置的时间导数
    dotxk = cos(uk) * dotrk- rk * sin(uk) * dotuk
    dotyk = sin(uk) * dotrk + rk * cos(uk) * dotuk
    #计算IGSO/MEO卫星在BDCS-ECEF下的速度
    dotR12=[[cos(Omegak),-sin(Omegak)*cos(ik),-xk1*sin(Omegak)-yk1*cos(Omegak)*cos(ik), yk1*sin(Omegak)*sin(ik)],
            [sin(Omegak), cos(Omegak)*cos(ik), xk1*cos(Omegak)-yk1*sin(Omegak)*cos(ik),-yk1*cos(Omegak)*sin(ik)],
            [0,sin(ik),0,Yk*cos(ik)]]
    tx=[[dotxk],[dotyk],[satellite['OMEGA_DOT']-OMGE_GAL],[dotik]]
    dxyz=np.array(dotR12).dot(np.array(tx))
    Fclk=-2*sqrt(GM_GAL)/clight/clight*satellite['e']*sqrt(Ak)*sin(Ek)#相对论效应
    tdts=satellite['a0']+satellite['a1']*(rt-toc)+satellite['a2']*((rt-toc)**2)+Fclk
    tdtss=satellite['a1']+2*satellite['a2']*(rt-toc)
    return [Xk,Yk,Zk,tdts,float(dxyz[0][0]),float(dxyz[1][0]),float(dxyz[2][0]),tdtss,satellite['BGD2'],f1*f1/f2/f2*satellite['BGD2']]


#函数, 计算卫星位置/广播星历SPP
def SPPM_form_BRDC(obs_mats,obs_index,GPS_eph=[],BDS_eph=[],GAL_eph=[],sat_out=[],ion_param=[],sol_mode='IF',freqs=[[1575.42e6,1227.60e6]],el_threthod=7.0,obslist=[],pre_rr=[],out_mode='Pos'):
    rr=[100,100,100]
    #观测值列表构建(异常值剔除选星)
    obs_sys_list=[]
    if(not len(obslist)):
        #双频观测值有效性校验
        obslist=[]
        for sys_obs in obs_mats:
            obs_mat=sys_obs#逐系统读取观测值
            for i in range(len(obs_mat[obs_index][1])):
                obsdata=obs_mat[obs_index][1][i]['OBS']
                obshealth=1
                if(sol_mode=='IF' and (obsdata[0]==0.0 or obsdata[1]==0.0 or obsdata[5]==0.0 or obsdata[6]==0.0)):
                    obshealth=0
                elif(sol_mode!='IF' and obsdata[0]==0.0):
                    obshealth=0
                #双频观测值有效性校验
                if(obshealth):
                    if obs_mat[obs_index][1][i]['PRN'] not in sat_out:
                        if(obs_mat[obs_index][1][i]['PRN'][0] not in obs_sys_list):
                            obs_sys_list.append(obs_mat[obs_index][1][i]['PRN'][0])#记录PRN顺序
                        obslist.append(obs_mat[obs_index][1][i])  
    obslist_new=obslist.copy()#高度角截至列表
    sat_num=len(obslist)
    ex_index=np.zeros(sat_num,dtype=int)#排除卫星标识
    #广播星历计算卫星位置
    peph_sat_pos={}
    for i in range(sat_num):
        #光速
        clight=2.99792458e8
        #观测时间&观测值
        rt_week=obs_mat[obs_index][0]['GPSweek']
        rt_sec=obs_mat[obs_index][0]['GPSsec']
        rt_unix=gpst2time(rt_week,rt_sec)
        
        #原始伪距
        p1=obslist[i]['OBS'][0]
        s1=obslist[i]['OBS'][4]
        p2=obslist[i]['OBS'][5]
        s2=obslist[i]['OBS'][9]
        #卫星位置内插
        si_PRN=obslist[i]['PRN']
        #逐系统计算卫星位置
        if(si_PRN[0]=='C'):
            rs_info=BDSEPH2SatPos(BDS_eph,rt_unix,si_PRN,None,rho=p1)
            if(not rs_info):
                ex_index[i]=1
                continue
        elif(si_PRN[0]=='G'):
            rs_info=GPSEPH2Satpos(GPS_eph,rt_unix,si_PRN,None,rho=p1)
            if(not rs_info):
                ex_index[i]=1
                continue
        elif(si_PRN[0]=='E'):
            rs_info=GALEPH2Satpos(GAL_eph,rt_unix,si_PRN,None,rho=p1)
            if(not rs_info):
                ex_index[i]=1
                continue
        #暂不支持的卫星系统
        else:
            ex_index[i]=1
            continue
        rs=[rs_info[0],rs_info[1],rs_info[2]]    #广播卫星轨道
        drs=[rs_info[4],rs_info[5],rs_info[6]]   #广播卫星速度
        dts=rs_info[3]                           #广播卫星钟差(含相对论效应)
        #码间偏差信息校正(对于BDS而言, 钟差基准建立在B3I单频上, 因此需要对无电离层组合进行码间偏差校正)
        dcb1,dcb2=0.0,0.0
        if(si_PRN[0]=='C'):
            if(obs_mats[obs_sys_list.index('C')][obs_index][0]['obstype'][0] in ['C2I','C1I']):
                dcb1=rs_info[8]
            if(obs_mats[obs_sys_list.index('C')][obs_index][0]['obstype'][4] in ['C7I']):
                dcb2=rs_info[9]
        if(si_PRN[0]=='G'):
        #对于GPS而言, 钟差基准建立在L1/L2组合观测值上, IF组合无需进行码偏差校正
            if(obs_mats[obs_sys_list.index('G')][obs_index][0]['obstype'][0] in ['C1C']):
                dcb1=rs_info[8]
            if(obs_mats[obs_sys_list.index('G')][obs_index][0]['obstype'][4] in ['C2W','C2X','C2P']):
                dcb2=rs_info[9]
        #对于GAL而言, 钟差基准建立在E1/E5b组合观测值上, IF组合无需进行码偏差校正
        if(si_PRN[0]=='E'):
            if(obs_mats[obs_sys_list.index('E')][obs_index][0]['obstype'][0] in ['C1C','C1X','C1P']):
                dcb1=rs_info[8]
            if(obs_mats[obs_sys_list.index('E')][obs_index][0]['obstype'][4] in ['C5Q','C5X']):
                dcb2=rs_info[9]
        peph_sat_pos[si_PRN]=[rs[0],rs[1],rs[2],dts,drs[0],drs[1],drs[2],rs_info[7],dcb1,dcb2]
    
    if(out_mode=="Sat only"):
        return obslist, peph_sat_pos
    #排除星历不支持的卫星
    for i in range(sat_num):
        if(ex_index[i]):
            obslist_new.remove(obslist[i])
    obslist=obslist_new.copy()#更新观测列表
    obslist_new=obslist.copy()#高度角截至列表
    sat_num=len(obslist)#重新计算卫星数量
    ex_index=np.zeros(sat_num,dtype=int)#排除卫星标识
    #记录各系统卫星数量
    sat_num_C,sat_num_G,sat_num_E=0,0,0
    for obs in obslist:
        if(obs['PRN'][0]=='C'):
            sat_num_C+=1
        if(obs['PRN'][0]=='G'):
            sat_num_G+=1
        if(obs['PRN'][0]=='E'):
            sat_num_E+=1
    #确定卫星系统数量
    sys_num=0
    if(sat_num_G!=0):
        sys_num+=1
    if(sat_num_C!=0):
        sys_num+=1
    if(sat_num_E!=0):
        sys_num+=1
    #初始化高度角排除列表
    ex_index=np.zeros(sat_num,dtype=int)
    
    if(sat_num<3+sys_num):
        print("The number of Satellites is not enough, pass epoch.")
        return [0,0,0,0,0,0],[],[]
    
    #卫星位置计算
    #伪距单点定位SPP
    if(len(pre_rr)):
        #有先验位置
        rr[0]=pre_rr[0]
        rr[1]=pre_rr[1]
        rr[2]=pre_rr[2]
    result=np.zeros((6),dtype=np.float64)
    result[0:3]=rr
    result[3]=1.0
    if(len(pre_rr)):
        result[3]=pre_rr[3]
    #最小二乘求解滤波初值
    ls_count=0
    while(1):
        #光速, GPS系统维持的地球自转角速度(弧度制)
        clight=2.99792458e8
        OMGE=7.2921151467E-5

        #观测值矩阵初始化
        Z=np.zeros(sat_num,dtype=np.float64)
        #设计矩阵初始化
        H=np.zeros((sat_num,4+2),dtype=np.float64)
        #单位权中误差矩阵初始化(虚拟观测存在, 随机模型向量长度不定)
        var=[]
        #观测值、设计矩阵构建
        for i in range(0,sat_num):        
            #观测时间&观测值
            rt_week=obs_mat[obs_index][0]['GPSweek']
            rt_sec=obs_mat[obs_index][0]['GPSsec']
            rt_unix=gpst2time(rt_week,rt_sec)
            si_PRN=obslist[i]['PRN']
            #print(rt_week,rt_sec,rt_unix)
        
            #伪距
            p1=obslist[i]['OBS'][0]
            s1=obslist[i]['OBS'][4]
            p2=obslist[i]['OBS'][5]
            s2=obslist[i]['OBS'][9]
            #print(p1,p2,phi1,phi2)
            #卫星PRN
            si_PRN=obslist[i]['PRN']
            #频率分发
            f1=freqs[obs_sys_list.index(si_PRN[0])][0]
            f2=freqs[obs_sys_list.index(si_PRN[0])][1]
            #卫星位置
            rs=[peph_sat_pos[si_PRN][0],peph_sat_pos[si_PRN][1],peph_sat_pos[si_PRN][2]]
            dts=peph_sat_pos[si_PRN][3]
            
            r0=sqrt( (rs[0]-rr[0])*(rs[0]-rr[0])+(rs[1]-rr[1])*(rs[1]-rr[1])+(rs[2]-rr[2])*(rs[2]-rr[2]) )
            #线性化的站星单位向量
            urs_x=(rr[0]-rs[0])/r0
            urs_y=(rr[1]-rs[1])/r0
            urs_z=(rr[2]-rs[2])/r0
            
            #单卫星设计矩阵赋值
            if(si_PRN[0]=='G'):#对每颗卫星所属系统钟差对应项系数进行赋值
                H[i]=[urs_x,urs_y,urs_z,1,0,0]
                rr_clk=result[3]
            if(si_PRN[0]=='C'):
                H[i]=[urs_x,urs_y,urs_z,0,1,0]
                rr_clk=result[4]
            if(si_PRN[0]=='E'):
                H[i]=[urs_x,urs_y,urs_z,0,0,1]
                rr_clk=result[5]

            #地球自转改正到卫地距上
            r0=r0+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/clight

            #观测矩阵
            if(sol_mode=='SF'):
                Z[i]=p1-r0-rr_clk-get_Tropdelay(rr,rs)-get_ion_GPS(rt_unix,rr,rs,ion_param)+clight*dts
            
            #双频无电离层延迟组合
            elif(sol_mode=='IF'):
                f12=f1*f1
                f22=f2*f2
                p_IF=f12/(f12-f22)*(p1-clight*peph_sat_pos[si_PRN][8])-f22/(f12-f22)*(p2-clight*peph_sat_pos[si_PRN][9])
                Z[i]=p_IF-r0-rr_clk-get_Tropdelay(rr,rs)+clight*dts
            elif(sol_mode=='BDGIM'):
                Z[i]=p1-r0-rr_clk-get_Tropdelay(rr,rs)-40.28*get_BDSGIM(rt_unix,ion_param,rr,rs,MF_mode=2)/(f1/1e8)/(f1/1e8)+clight*dts

            #随机模型
            #var[i][i]= 0.00224*10**(-s1 / 10) 
            _,el=getazel(rs,rr)
            var.append(0.3*0.3+0.3*0.3/sin(el)/sin(el))
            if(el*180.0/pi<el_threthod):
                var[i]=var[i]*100#低高度角拒止
                ex_index[i]=1
            if(ex_index[i]==1 and el*180.0/pi>=el_threthod):
                ex_index[i]=0
            if(sol_mode=='IF'):
                var[i]=var[i]*9
        #各系统钟差失配,添加虚拟零观测确保H满秩
        if(sat_num_G==0):
            H_sub_G=np.array([0,0,0,1,0,0]).reshape(1,6)
            var.append(0.01)#虚拟方差
            H=np.concatenate((H,H_sub_G)) #设计矩阵处理
            Z=np.append(Z,0.0)
        if(sat_num_C==0):
            H_sub_C=np.array([0,0,0,0,1,0]).reshape(1,6)
            var.append(0.01)#虚拟方差
            H=np.concatenate((H,H_sub_C)) #设计矩阵处理
            Z=np.append(Z,0.0)
        if(sat_num_E==0):
            H_sub_E=np.array([0,0,0,0,0,1]).reshape(1,6)
            var.append(0.01)#虚拟方差
            H=np.concatenate((H,H_sub_E)) #设计矩阵处理
            Z=np.append(Z,0.0)

        #权重矩阵
        W=np.zeros((len(var),len(var)),dtype=np.float64)
        for i in range(len(var)):
            W[i][i]=1.0/var[i]
        
        #最小二乘求解:
        dresult=getLSQ_solution(H,Z,W=W,weighting_mode='S')
        
        #迭代值更新
        result[0]+=dresult[0]
        result[1]+=dresult[1]
        result[2]+=dresult[2]
        result[3]+=dresult[3]
        result[4]+=dresult[4]
        result[5]+=dresult[5]

        #更新测站位置
        rr[0]=result[0]
        rr[1]=result[1]
        rr[2]=result[2]
        #print(dresult)
        ls_count+=1
        if(abs(dresult[0])<1e-4 and abs(dresult[1])<1e-4 and abs(dresult[2])<1e-4):
            #估计先验精度因子
            break
        if(ls_count>200):
            break
    #排除低高度角卫星
    for i in range(sat_num):
        if(ex_index[i]):
            obslist_new.remove(obslist[i])
    try:
        global Q_spp
        Q_spp=inv(H.T.dot(W).dot(H))#估计先验精度因子
    except:
        Q_spp=np.zeros((len(result),len(result)))
    return result,obslist_new,peph_sat_pos


if __name__ == "__main__":
    
    #广播星历读取
    eph_path= 'data/BRDC/brdc1320.24p' 
    ion_parmas_GPSK,GPS_eph=BRDC2GPSEPH(eph_path)
    ion_parmas_BDSK,BDS_eph=BRDC2BDSEPH(eph_path)
    ion_parmas_GAL,GAL_eph=BRDC2GALEPH(eph_path)

    #观测文件
    obs_path='data/OBS/WUH2/wuh21320.24o'
    #观测系统与频点(需对应)
    sys_indexes=['G','C','E']
    obs_types=[
               ['C1C','L1C','D1C','S1C','C2W','L2W','D2W','S2W'],
               ['C2I','L2I','D2I','S2I','C6I','L6I','D6I','S6I'],
               ['C1X','L1X','D1X','S1X','C5X','L5X','D5X','L5X']
               ]
    freqs=[
           [1575.420e6,1176.45e6],
           [1561.098e6,1268.52e6],
           [1575.420e6,1227.60e6]
           ]
    #排除卫星列表
    sat_out=[]

    #单频/双频(SF/IF)
    sol_mode='IF'
    
    #结果文件保存路径
    out_path='nav_result'
        
    
    #观测值校验
    obs_mats=[]
    for i in range(len(sys_indexes)):
        sys=sys_indexes[i]
        t_obs_mat=reconstruct_obs_mat(RINEX3_to_obsmat(obs_path,obs_types[i],sys,0,f1=freqs[i][0],f2=freqs[i][1]))
        obs_mats.append(t_obs_mat.copy())
    multi_GNSS_YES=check_obs_mats(obs_mats)
    if(not multi_GNSS_YES):
        ValueError("Muti-GNSS Observations Not Valid")
    
    
    print("Easy4SPP Initial Success, Configurations: ")
    print("GNSS Systems:", sys_indexes)
    print("GNSS Code Types:",obs_types)
    print("GNSS Frequency:",freqs)
    print("SPP Observation Combination Mode: ",sol_mode)
    print("SPP Results Output Path: ",out_path)

    #逐历元SPP
    sol_logs=[]
    
    for i in range(0,len(obs_mats[0])):
        
        #观测时间&观测值
        rt_week=obs_mats[0][i][0]['GPSweek']
        rt_sec=obs_mats[0][i][0]['GPSsec']
        rt_unix=gpst2time(rt_week,rt_sec)
        
        #标准单点定位SPP解算
        result,obslist_new,peph_sat_pos=SPPM_form_BRDC(obs_mats,i,sat_out=sat_out,BDS_eph=BDS_eph,GPS_eph=GPS_eph,GAL_eph=GAL_eph,freqs=freqs,ion_param=ion_parmas_GPSK,sol_mode=sol_mode)
        #后验精度校验
        if(result[0]!=0.0):
            q_neu=xyz2neu([result[0]+sqrt(Q_spp[0][0]),result[1]+sqrt(Q_spp[1][1]),result[2]+sqrt(Q_spp[2][2])],result)
        ct=time2COMMONTIME(rt_unix)
        print("[{}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}]".format(ct['year'],ct['month'],ct['day'],ct['hour'],ct['minute'],int(ct['second'])),len(obslist_new),abs(q_neu[0]),abs(q_neu[1]),abs(q_neu[2]),end='\r')
        sol_logs.append([rt_week,rt_sec]+list(result)+list(np.abs(q_neu))+[len(obslist_new)])
    try:
        np.save(out_path+'/{}.out'.format(os.path.basename(obs_path)),sol_logs,allow_pickle=True)
        print("\nNavigation results saved at ",out_path+'/{}.out'.format(os.path.basename(obs_path)))
    except:
        np.save('nav_result'+'/{}.out'.format(os.path.basename(obs_path)),sol_logs,allow_pickle=True)
        print("\nNavigation results saved at ",'nav_result'+'/{}.out'.format(os.path.basename(obs_path)))
    print("Easy4PPP Solving Success")
    #print(sqrt(sum(np.array(neus)[:,0]**2)/len(neus)),sqrt(sum(np.array(neus)[:,1]**2)/len(neus)),sqrt(sum(np.array(neus)[:,2]**2)/len(neus)))
    
