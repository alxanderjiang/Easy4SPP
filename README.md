# Easy4SPP

[[中]](./README-zh.md) &ensp; [[EN]](./README.md)

An easily ported multi-GNSS Standard Point Positioning (SPP) and common PNT model generation library toolbox coded in Python.

This is a part of open-source toolbox Easy4PNT. Other toolboxs of Easy4PNT is listed here (clicked to jump to the target): [[Easy4B2b]](https://github.com/alxanderjiang/Easy4B2b), [[Easy4RTK]](https://github.com/alxanderjiang/Easy4RTK), [[Easy4PPP]](https://github.com/alxanderjiang/Easy4PPP), [[Easy4PTK]](https://github.com/alxanderjiang/Easy4PTK).

## Quick Start
1. Running the main function of src/brdc_proc.py to get the example IF-SPP results of WUH2 on DOY 132, 2024. The results stored in the "nav_result" folder.
2. The result can be visualized by running the Jupyter Notebook file "nav_result.ipynb". The positioning errors, clock bias of different systems is shown after running all the blocks in "nav_result.ipynb"
3. Users can transfer the results log in '.npy' format into ASCII log in '.txt' format by running the last block of "nav_result.ipynb"

## Downloading and preperations
1. Download the **zip pakeage directly** or using git clone by running the following commend:
```bash
git clone https://github.com/alxanderjiang/Easy4SPP.git
```
2. Unzip the sample data folders: data.zip and nav_result.zip to the same path of Easy4SPP. If linux but no GUI, please run the following commends:

```bash
cd Easy4SPP
unzip data.zip
```
3. Ensure that the numpy, tqdm, ipykernel, numba are available in your Python environment. If not, please run the following commends to install:

```bash
pip install numpy
pip install tqdm
pip install ipykernel
pip install numba
```

  numpy and tqdm is used in the core codes while ipykernel is necessary to run Jupyter Notebook tutorials. numba is used to accelerate the computation (this can be ignored by change all the "numba_inv" function to simple "inv()" function). Easy4SPP only support running from a __main__ function with variables definition for post solving.
Some problems may happen when install or use numba because of laking the library scipy, please install it by running the following commends:

```bash
pip install scipy
```

## Configurations
Easy4SPP gives many configuration choices for different dataset collected by different GNSS devices. All configurations shold be set in the __main__ function of src/brdc_proc.py, from about line 1141 to 1167. The path of broadcast ephemeris in version of RINEX 3.x shold be set as:
```Python
    #广播星历读取
    eph_path= 'data/BRDC/brdc1320.24p' 
```
The path of obervation file in version of RINEX 3.x or RINEX 4.x shold be set as:
```Python
    #观测文件
    obs_path='data/OBS/WUH2/wuh21320.24o'
```
The systems, selected signals and frequencies shold be set as:
```Python
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
```
Make sure that all the system&signal&frequency selections are corresponding in above lists. The outlier satellites shold be set as:
```Python
    #排除卫星列表
    sat_out=[]
```
If sat_out is empty, that means no satellites will be outlier in Easy4SPP. Choose single-freqency (SF) or ion-free (IF) observation model as:
```Python
    #单频/双频(SF/IF)
    sol_mode='IF'
```
The SF mode means using the GPS Klobchar model parameters readed in broadcast ephemeris, the IF means that using IF obsevation model to eliminate ionospheric delay. Finally, the results saving path shold be set as:
```Python
    #结果文件保存路径
    out_path='nav_result'
```
