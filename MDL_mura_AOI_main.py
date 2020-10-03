from math import e
#import cv2
#import pylab as plt
from PIL import Image
import numpy as np
import sys
import os
import datetime
import tkinter
import time
import csv

window = tkinter.Tk()
window.title('AOI光學自動檢測系統')
window.geometry('400x100')

def action():
    def LOG(x, y, c):
        return (1 / (c ** 4)) * ((x ** 2 + y ** 2) / (c ** 2) - 2) * e ** (-(x ** 2 + y ** 2) / (2 * c ** 2))

    def FilterSize(c):
        return (6 * c + 1)

#[Wi權重分析]

#-----------------------------------------------Point 2D Filter

    #print('輸入2D濾波器西格瑪Size:')
    c = num1.get()
    
    print(int(FilterSize(c)), 'X', int(FilterSize(c)))

    n = (int(FilterSize(c)) - 1) / 2
    first = n
    end = - n
    LOGF = [] #連續曲線LOG濾波器
    for i in range(int(first), int(end) - 1, -1):
        y = i
        a = []
        for j in range(int(end), int(first) + 1):
            x = j
            a.append(- LOG(x, y, c))
        LOGF.append(a)

    T = 0
    F = 0
    for i in range(len(LOGF)):
        for j in range(len(LOGF)):
            if LOGF[i][j] < 0:
                T += LOGF[i][j]
            else:
                F += LOGF[i][j]

    N = 1 #倍數

    LOGFilter = [] #離散LOG濾波器
    for i in range(len(LOGF)):
        p = []
        for j in range(len(LOGF)):
            if LOGF[i][j] < 0:
                p.append(N * - LOGF[i][j] / T)
            else:
                p.append(N * LOGF[i][j] / F)
        LOGFilter.append(p)

    LOGFilter_array = np.array(LOGFilter) 

#-----------------------------------------------H Band 1D Filter

    #print('輸入1D濾波器西格瑪Size:')
    cb = num2.get()
    
    print(int(FilterSize(cb)), 'X', 1)

    n = (int(FilterSize(cb)) - 1) / 2
    first = n
    end = - n

    LOGBF = [] #連續曲線LOG濾波器
    y = 0
    for i in range(int(end), int(first) + 1):
        x = i
        LOGBF.append(- LOG(x, y, cb))
    T = 0
    F = 0
    for i in range(len(LOGBF)):
        if LOGBF[i] < 0:
            T += LOGBF[i]
        else:
            F += LOGBF[i]

    N = 1#倍數
    LOGBFilter1 = [] #離散LOG濾波器
    for i in range(len(LOGBF)):
        if LOGBF[i] < 0:
            LOGBFilter1.append(N * - LOGBF[i] / T)
        else:
            LOGBFilter1.append(N * LOGBF[i] / F)

    LOGBFilter = []
    LOGBFilter.append(LOGBFilter1)

    LOGBFilter_array = np.array(LOGBFilter)

#-----------------------------------------------V Band 1D Filter    

    #print('輸入1Dv濾波器西格瑪Size:')
    cv = num5.get()

    print(1, 'X', int(FilterSize(cv)))

    n = (int(FilterSize(cv)) - 1) / 2
    first = n
    end = - n
    LOGVF = [] #連續曲線LOG濾波器
    for i in range(int(first), int(end) - 1, -1):
        y = i
        x = 0
        v = []
        v.append(-LOG(x, y, cv))
        LOGVF.append(v)
    #print(LOGVF)


    T = 0
    F = 0
    for i in range(len(LOGVF)):
        if LOGVF[i][0] < 0:
            T += LOGVF[i][0]
        else:
            F += LOGVF[i][0]

    N = 1 #倍數

    LOGVFilter = [] #離散LOG濾波器
    for i in range(len(LOGVF)):
        pv = []
        if LOGVF[i][0] < 0:
            pv.append(N * - LOGVF[i][0] / T)
        else:
            pv.append(N * LOGVF[i][0] / F)
        LOGVFilter.append(pv)

    LOGVFilter_array = np.array(LOGVFilter)

#===================================================================================================

#[Xi輸入值 x Wi權重分析]
#-----------------------------------------------Point 2D Filter 影像卷積convolution處理

    with open('ID座標資料.csv', 'a', encoding='cp950', newline='') as ff:
        idx = 0
        while True:
            t1 = datetime.datetime.now()
            idx += 1
            filename = "picture {0}.jpg".format(idx)
            if os.path.exists(filename):
                img_pil = Image.open(filename)#開啟CCD擷取的不良照片
                X, Y = img_pil.size
                img_gray = img_pil.convert('L')#轉換成灰階圖片
                img_p = np.array(img_pil)#原圖片轉成陣列
                img = np.array(img_gray)#灰階圖片轉成陣列

                AfterImage = [[0] * X for _ in range(Y)]
                m = 0
                n = 0
                W = 0
                for j in range(Y - int(FilterSize(c)) // 2 * 2): #卷積陣列計算
                    for i in range(X - int(FilterSize(c)) // 2 * 2):
                        W = img[j:int(FilterSize(c)) + m, i:int(FilterSize(c)) + n] * LOGFilter_array
                        #print(sum(sum(W)),j, int(FilterSize(c)) + m, i, int(FilterSize(c)) + n)
                        AfterImage[j + int(FilterSize(c)) // 2][i + int(FilterSize(c)) // 2] = sum(sum(W))
                        n += 1
                    n = 0
                    m += 1
#-----------------------------------------------H Band 1D Filter 影像卷積convolution處理
                AfterImageB = [[0] * X for _ in range(Y)]
                n = 0
                W = 0
                for j in range(Y): #卷積陣列計算
                    for i in range(X - int(FilterSize(cb)) // 2 * 2):
                        W = img[j, i:int(FilterSize(cb)) + n] * LOGBFilter_array
                        AfterImageB[j][i + int(FilterSize(cb)) // 2] = sum(sum(W))
                        n += 1
                    n = 0        
                
#----------------------------------------------V Band 1D Filter 影像卷積convolution處理
                AfterImageV = [[0] * X for _ in range(Y)]
                W = 0
                m = 0
                for j in range(Y - int(FilterSize(cv)) // 2 * 2):
                    for i in range(X):
                        W = img[j:int(FilterSize(cv)) + m, i] * LOGVFilter_array[0:int(FilterSize(cv)), 0]
                        AfterImageV[j + int(FilterSize(cv)) // 2][i] = sum(W)
                    m += 1
  
#=============================================================================================
        
#[閥值]
                LimitValue = num3.get()#LOG門檻設定值
                LimitValueB = num4.get()#LOG門檻設定值
                LimitValueV = num6.get()#LOG門檻設定值
                
                
                for y in range(Y):
                    for x in range(X):
                        if abs(AfterImage[y][x]) > LimitValue:
                            img_p[y, x, 0] = 255
                            img_p[y, x, 1] = 0
                            img_p[y, x, 2] = 0
                        

                        if abs(AfterImageB[y][x]) > LimitValueB:
                            img_p[y, x, 0] = 0
                            img_p[y, x, 1] = 255
                            img_p[y, x, 2] = 0
                        
                       
                        if abs(AfterImageV[y][x]) > LimitValueV:
                            img_p[y, x, 0] = 0
                            img_p[y, x, 1] = 0
                            img_p[y, x, 2] = 255
                       
                        
                        
                
        
                im = Image.fromarray(img_p)
                im.show()
                T = "Tes2tpicture {0}.jpg".format(idx)
                im.save(T)
                
            
             
        
                k = [[0 for i in range(X)] for j in range(Y)]
                kb = [[0 for i in range(X)] for j in range(Y)]
                kv = [[0 for i in range(X)] for j in range(Y)]
                
                u = [[-1, 0], [-1, 1], [0, 1], [1, 1], [1, 0], [1, -1], [0, -1], [-1, -1]]
        
                def limit():#不良紅框圈製作
                    tempxmax = 0
                    tempymax = 0
                    tempxmin = temp[0][0]
                    tempymin = temp[0][1]
                    for y in range(len(temp)):
                        if temp[y][1] > tempymax:
                            tempymax = temp[y][1]
                        if temp[y][0] > tempxmax:
                            tempxmax = temp[y][0]
                        if temp[y][1] < tempymin:
                            tempymin = temp[y][1]
                        if temp[y][0] <tempxmin:
                            tempxmin = temp[y][0]
                        jmax, imax = tempxmax, tempymax
                        jmin, imin = tempxmin, tempymin
                        jmiddle = (jmax - jmin) // 2 + jmin
                        imiddle = (imax - imin) // 2 + imin
        
                    for i in range(imax - imin + 1):
                        img_p[jmin, imin + i, 0] = 255
                        img_p[jmin, imin + i, 1] = 0
                        img_p[jmin, imin + i, 2] = 0
                        img_p[jmax, imin + i, 0] = 255
                        img_p[jmax, imin + i, 1] = 0
                        img_p[jmax, imin + i, 2] = 0
                    for j in range(jmax - jmin + 1):
                        img_p[jmin + j, imax, 0] = 255
                        img_p[jmin + j, imax, 1] = 0
                        img_p[jmin + j, imax, 2] = 0
                        img_p[jmin + j, imin, 0] = 255
                        img_p[jmin + j, imin, 1] = 0
                        img_p[jmin + j, imin, 2] = 0
        
                    print()
                    print("不良位置座標:", imiddle, jmiddle)
                    return imiddle, jmiddle
                
                    
                def check(A, j, i, t, L, k):
                    k[j][i] = 1
                    stack = [[j, i]]
                    while len(stack) > 0:
                        j, i = stack.pop()
                        t.append([j, i])
                        for e in u:
                            jj, ii = j + e[0], i + e[1]
                            if 0 <= jj < Y and 0 <= ii < X:
                                x = A[jj][ii]
                                if abs(x) > L and k[jj][ii] != 1:
                                    t.append([jj, ii])
                                    stack.append([jj, ii])
                                    k[jj][ii] = 1
        
        
                result = []
                Defect = []
                for j in range(Y):
                    for i in range(X):
                        temp = []
                        if abs(AfterImage[j][i]) > LimitValue and k[j][i] != 1:
                            check(AfterImage, j, i, temp, LimitValue, k)
                            
                        if abs(AfterImageB[j][i]) > LimitValueB and kb[j][i] != 1:
                            check(AfterImageB, j, i, temp, LimitValueB, kb)
                            
                        if abs(AfterImageV[j][i]) > LimitValueV and kv[j][i] != 1:
                            check(AfterImageV, j, i, temp, LimitValueV, kv)
                            
                        if len(temp) > 0:
                            #limit()
                            imi, jmi = limit()
                            Defect.append([imi, jmi])
                            result.append(temp)
#============================================================================        
                im = Image.fromarray(img_p)
                im.show()
                T = "Testpicture {0}.jpg".format(idx)
                im.save(T)
                print()
                print("不良項目數量:", len(result))
                t2 = datetime.datetime.now()
                print()
                print("檢查時間:", t2 - t1)
                t3 = time.ctime()

                data = [[T, Defect, t3, len(result)]]
                csv_writer = csv.writer(ff)
                csv_writer.writerows(data)
            else:
                break
 
        finish.set('E~N~D')   
    
   
#操作介面======================================================================    

num1 = tkinter.IntVar()
num2 = tkinter.IntVar()
num3 = tkinter.DoubleVar()
num4 = tkinter.DoubleVar()
num5 = tkinter.IntVar()
num6 = tkinter.DoubleVar()
finish = tkinter.StringVar()

#--------------------------------------------------------------------
TwoD_f_item = tkinter.Entry(window, width=10, textvariable=num1)
OneDB_f_item = tkinter.Entry(window, width=10, textvariable=num2)
OneDV_f_item = tkinter.Entry(window, width=10, textvariable=num5)
TwoD_f_item.grid(row=1, column=1) 
OneDB_f_item.grid(row=2, column=1)
OneDV_f_item.grid(row=3, column=1)

TwoD_limitvalue = tkinter.Entry(window, width=10, textvariable=num3)
OneDB_limitvalue = tkinter.Entry(window, width=10, textvariable=num4)
OneDV_limitvalue = tkinter.Entry(window, width=10, textvariable=num6)
TwoD_limitvalue.grid(row=1, column=2)
OneDB_limitvalue.grid(row=2, column=2)
OneDV_limitvalue.grid(row=3, column=2)

#-------------------------------------------------------------------

btn1 = tkinter.Button(window, width=15, text='Action', command=action)
btn1.grid(row=0, column=3)

label = tkinter.Label(window, width=15, textvariable=finish)
label.grid(row=1, column=3)

#--------------------------------------------------------------------

TwoD_f_label = tkinter.Label(window, width=10, text='2D濾波器')
TwoD_f_label.grid(row=1, column=0)

OneDB_f_label = tkinter.Label(window, width=10, text='1DH濾波器')
OneDB_f_label.grid(row=2, column=0)

OneDV_f_label = tkinter.Label(window, width=10, text='1DV濾波器')
OneDV_f_label.grid(row=3, column=0)

filtervalue_label = tkinter.Label(window, width=10, text='西格瑪參數')
filtervalue_label.grid(row=0, column=1)

limitvalue_label = tkinter.Label(window, width=10, text='閥值')
limitvalue_label.grid(row=0, column=2)



window.mainloop()

