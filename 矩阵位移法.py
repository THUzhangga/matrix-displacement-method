#!/usr/bin/python3.5

# -*- coding: utf-8 -*-

#author:ZhangGa(Peter Zhang)

#date:2016.10.9

#copyright:Tsinghua University

import numpy as np  # 导入numpy模块并命名为np
import numpy.linalg as nplg  # 导入numpy模块的linalg子模块并命名为nplg
import matplotlib.pyplot as plt  # 导入matplotlib模块的pyplot子模块并命名为plt
import os  # 导入操作系统函数模块
import win32ui  # 导入win32api的win32ui模块
import re  #导入正则表达式模块
class Unit:
    num = 0  # 单元编号
    beginnode = []  # 始端
    endnode = []  # 末端
    vector = []  # 定位向量
    EA = 0
    EI = 0
    l = 0  # 长度
    alpha = 0  # 角度α（整体坐标系顺时针旋转到局部坐标系为正）
    q = 0  # 均布荷载
    _k = np.zeros((6, 6))  # 刚度矩阵（局部坐标系下）
    k = np.zeros((6, 6))  # 刚度矩阵（整体坐标系下）
    _Fp = np.zeros(6)  # 单元固端约束力向量（局部坐标系下）
    Fp = np.zeros(6)  # 单元固端约束力向量（整体坐标系下）
    _F = np.zeros(6)  # 杆端内力（局部坐标系下）
    F = np.zeros(6)  # 杆端内力（整体坐标系下）
    P = np.zeros(6)  # 单元等效结点荷载（整体坐标系下）
    T = np.zeros((6, 6))  # 坐标转换矩阵
    delta = np.zeros(6)  # 结点位移分量（局部坐标系下）
    _delta = np.zeros(6)  # 结点位移分量（整体坐标系下）
node = {}
unit = {}
shift = []#结点位移数
ground = []#基础
P_node = []#结点外荷载
filename = []#文件名
def SaveData(unit):
    print(">>>保存计算结果")
    data = open(filename[0].split('.')[0]+'计算结果.txt','w')
    for u in unit.values():
        data.write('单元编号：'+str(u.num)+'\n')
        data.write('局部坐标系下单元'+str(u.num)+'刚度矩阵：\n' + str(u._k) + '\n')
        data.write('整体坐标系下单元'+str(u.num)+'刚度矩阵：\n' + str(u.k) + '\n')
        data.write('局部坐标系下单元'+str(u.num)+'杆端内力：\n' + str(u._F)+'\n')
        data.write('整体坐标系下单元'+str(u.num)+'杆端内力：\n' + str(u.F) + '\n')
        data.write('局部坐标系下单元'+str(u.num)+'结点位移分量：\n' + str(u._delta) + '\n')
        data.write('整体坐标系下单元'+str(u.num)+'结点位移分量：\n' + str(u.delta) + '\n\n')

def ReadData(shfit,node,unit,P_node):
    dlg = win32ui.CreateFileDialog(1)  # 1表示打开文件对话框
    dlg.SetOFNInitialDir(os.getcwd())  # 设置打开文件对话框中的初始显示目录
    dlg.DoModal()#打开文件对话框
    filename.append(dlg.GetPathName()) # 获取选择的文件名称
    data = open(filename[0],'r')
    for line in data:
        text = line.split(',')
        if text[0] == "结点位移数":
            shfit.append(int(text[1]))
        if text[0] == "结点":
            num = int(text[1])
            node[num] = list(map(lambda x: float(x), text[2:4]))
        if text[0] == "基础":
            ground.append(int(text[1]))
            ground.append(int(text[2]))
        if text[0] == "单元":
            num = int(text[1])
            unit[num] = Unit()
            unit[num].num = num
            unit[num].beginnode = node[int(text[2])]
            unit[num].endnode = node[int(text[3])]
            unit[num].vector = list(map(lambda x: int(x), text[4:10]))
            unit[num].EA = float(text[10])
            unit[num].EI = float(text[11])
            unit[num].l = float(text[12])
            unit[num].alpha = (float(text[13]) / 180) * np.pi
            unit[num].q = (float(text[14]))
        if text[0] == "局部荷载":
            num = int(text[1])
            unit[num]._Fp = list(map(lambda x: float(x), text[2:8]))
        if text[0] == "结点荷载":
            for x in text[1:]:
                P_node.append(float(x))
    data.close()

def K0(EA,EI,L):#本函数用于生成局部坐标系下的单元刚度矩阵
    a = EA / L
    b = EI / (L ** 1)
    c = EI / (L ** 2)
    d = EI / (L ** 3)
    return np.array([
        [a,0,0,-a,0,0],
        [0,12*d,6*c,0,-12*d,6*c],
        [0,6*c,4*b,0,-6*b,2*b],
        [-a,0,0,a,0,0],
        [0,-12*d,-6*c,0,12*d,-6*c],
        [0,6*c,2*b,0,-6*c,4*b]
    ])

def T0(alpha):#本函数用于通过旋转角α(alpha)生成单元坐标转换矩阵
    return np.array([
        [np.cos(alpha),np.sin(alpha),0,0,0,0],
        [-np.sin(alpha), np.cos(alpha), 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0,0,0,np.cos(alpha), np.sin(alpha), 0],
        [0,0,0,-np.sin(alpha), np.cos(alpha), 0],
        [0,0,0,0, 0,  1],
    ])

def P2W(K,k,lamda):#Part to Whole，局部矩阵定位到整体矩阵中
                    #K为整体阶段矩阵，k为当前局部矩阵(整体坐标系下)，lamda为当前定位向量
    order = k.shape[0] #order为局部矩阵的阶数
    nzlamda = []  # nonzero lamda，定位向量中的非零项
    dict = {}#建立一个字典，字典的键值为定位向量中的非零项（即nzlamda中的元素），值为局部矩阵中对应该键值的数
    #提取定位向量中的非零项并与原矩阵对应编码建立字典
    for i in range(order):
        if lamda[i] != 0:
            dict[lamda[i]] = i + 1
            nzlamda.append(lamda[i])
    for i in nzlamda:
        for j in nzlamda:
            i0 = dict[i]
            j0 = dict[j]
            K[i-1][j-1] += k[i0-1][j0-1]#这里各个项都要减1是因为计算机是从0开始存储的
    return K

def P2W_Vector(P,p,lamda):#局部等效结点荷载定位到整体荷载
    length = len(p)#获得p的长度
    nzlamda = []  # nonzero lamda，定位向量中的非零项
    dict = {}  # 建立一个字典，字典的键值为定位向量中的非零项（即nzlamda中的元素），值为局部矩阵中对应该键值的数
    # 提取定位向量中的非零项并与原矩阵对应编码建立字典
    for i in range(length):
        if lamda[i] != 0:
            dict[lamda[i]] = i + 1
            nzlamda.append(lamda[i])
    for i in nzlamda:
        i0 = dict[i]
        P[i-1] += p[i0 - 1]  # 这里各个项都要减1是因为计算机是从0开始存储的
    return P

def W2P(D,d,lamda):#整体坐标系下，整体结点位移D反向定位到局部结点位移d
    length = len(d)  # 获得d的长度
    d = np.zeros(length)
    nzle = []  # nonzero local encoing，对应的总码不为零的局部码
    dict = {}  # 建立一个字典，字典的键值为对应的总码不为零的局部码（即nz中的元素），值为对应的总码
    # 提取定位向量中的非零项并与原矩阵对应编码建立字典
    for i in range(length):
        if lamda[i] != 0:
            dict[i] = lamda[i]
            nzle.append(i+1)
    for i in nzle:
        i0 = dict[i-1]
        d[i-1] += D[i0 - 1]  # 这里各个项都要减1是因为计算机是从0开始存储的
    return d

#以下为画图函数
def line(Point1,Point2,width=1):#端点为Point1与Point2的线段
    if Point2[0] != Point1[0]:
        k = (Point2[1] - Point1[1]) / ((Point2[0] - Point1[0]) * 1.0)
        for x in np.linspace(Point1[0],Point2[0],200):
            plt.scatter(x,k * (x - Point1[0]) + Point1[1], color='k',s=width,marker='o',label=str)
    elif Point2[1] != Point1[1]:
        m = (Point2[0] - Point1[0]) / ((Point2[1] - Point1[1]) * 1.0)
        for y in np.linspace(Point1[1], Point2[1], 200):
            plt.scatter( m * (y-Point1[1]) + Point1[0],y, color='k',s=width, marker='o', label=str)

def base(Point):#画出基础
    for x in np.linspace(Point[0] - 0.3, Point[0] + 0.3, 100):
        plt.scatter(x, Point[1], color='k', s=0.1, marker='o', label=str)
    for t in np.linspace(Point[0] - 0.2, Point[0] + 0.2, 3):
        P1 = [t, Point[1]]
        P2 = [t - 0.15, Point[1] - 0.25]
        line(P1, P2, width=0.1)

def BMD(u,q,scanfactor=1):#Bending Moment Diagram弯矩图
    #先在局部坐标系下求出M的方程
    Fp = u._F#注意应为局部坐标系下的局部荷载！
    alpha = -u.alpha
    L = u.l
    begin = [0,-Fp[2]]
    end = [L,Fp[5]]
    k = (end[1] - begin[1]) /((end[0] - begin[0]) * 1.0)
    t = np.array([
        [np.cos(alpha), -np.sin(alpha)],
        [np.sin(alpha), np.cos(alpha)]
    ])
    for x in np.linspace(0,L,200):
        y = (k * (x - begin[0]) + begin[1] + q / 2.0 * x * (x- L)) / scanfactor
        draw = np.dot(t,[x,y]) + u.beginnode
        plt.scatter(draw[0], draw[1], color='k', s=0.5, marker='o', label=str)
    for x in np.linspace(0, L, 10):
        y = (k * (x - begin[0]) + begin[1] +  q / 2.0 * x * (x- L)) / scanfactor
        p1 = np.dot(t,[x, 0]) + u.beginnode
        p2 = np.dot(t,[x, y]) + u.beginnode
        line(p1,p2)

#读入文件
print('>>>读入文件');
ReadData(shift, node, unit, P_node)
examplename = re.split(r'\\|\.',filename[0])[-2]
num_shift = shift[0]
print('>>>计算'+examplename)

K = np.zeros((num_shift,num_shift))#整体刚度矩阵初始化
print('>>>计算整体刚度矩阵')
for k in unit.keys():
    u = unit[k]
    u._k = K0(u.EA,u.EI,u.l)
    u.T = T0(u.alpha)
    u.k = np.dot(np.dot(u.T.T,u._k),u.T)
    K = P2W(K,u.k,u.vector)

#计算等效结点荷载
print('>>>计算等效结点荷载')
P = np.zeros(num_shift)
for u in unit.values():
    u.Fp = np.dot(u.T.T,u._Fp)
    u.P = -u.Fp
    P = P2W_Vector(P,u.P,u.vector)
if len(P_node)==0:
    P_node = np.zeros(6)
P = P + P_node

#计算结点位移
print('>>>计算结点位移')
delta = np.dot(nplg.inv(K),P)#解基本方程得结点位移向量△(delta)

#计算单元杆端内力
print('>>>计算单元杆端内力')
for u in unit.values():
    u.delta = W2P(delta,u.delta,u.vector)
    u._delta = np.dot(u.T,u.delta)
    u._F = np.dot(u._k,u._delta) + u._Fp
for k in node.keys():
    n = node[k]
    plt.text(n[0], n[1], str(k), color='r')

#画基础
for g in ground:
    base(node[g])

#下面研究决定画图的图幅以及弯矩图的放缩（太大了会超出图...）
x_max = 0
y_max = 0
F_max = 0
for n in node.values():
    x_max = max(x_max, n[0])
    y_max = max(y_max, n[1])
plt.xlim(-x_max, x_max * 2)
plt.ylim(-y_max, y_max * 2)
for u in unit.values():
    F_max = max(F_max, abs(u._F[2] / u.l),abs(u._F[5] / u.l))
scanfactor = F_max * 2

#绘图
print('>>>绘制弯矩图')
for u in unit.values():
    line(u.beginnode,u.endnode)
    BMD(u,u.q,scanfactor)

#设置图片标题并保存图片
'''
显示中文时需要用u'中文',
但是由于matplotlib.pyplot在显
示时无法找到合适的字体，会显示
乱码（我的显示为方框），因此需
要设定字体
'''
plt.title(examplename+'弯矩图',fontproperties='SimHei')
plt.savefig(examplename+"弯矩图.png",dpi=200 )
SaveData(unit)  #保存数据
plt.show()
