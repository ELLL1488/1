import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import integrate
import pylab
def faza(z):
    a = z.real
    b = z.imag
    pi1 = math.acos(-1)
    if a == 0:
        if b > 0:
            faz = pi1 / 2
        else:
            faz = (-1) * pi1 / 2
    else:
        faz = math.atan(b / a)
        if a < 0:
            faz += pi1
    return faz
def trapezoid(ff, a, b, n):
    global f
    if b < a:
        a, b = b, a
    h = (b - a) / n
    g = [(0.5 * h * (ff(a + (i * h)) + ff(a + ((i + 1) * h)))) for i in range(0, n)]
    return math.sqrt(sum(g) * f)
iz_Q = 0
n = 40
sch = 0
QQQ1 = []
QQQ2 = []
QQQ3 = []
QQQ4 = []
QQQ5 = []
while sch <= n + 1:
    R1=50
    R2=4
    R4=50
    R5=2000
    C1=2000*10**(-6)
    C2=700*10**(-6)
    C4=300*10**(-6)
    L1=2*10**(-3)
    L2=25*10**(-3)
    f=30
    if iz_Q == 0:
        U=380
    pi=3.14
    Q1=380
    Q2=38
    hQ=(Q2-Q1)/n
    w=2*pi*f
    XL1=w*L1
    XL2=w*L2
    XC1=1/w/C1
    XC2=1/w/C2
    XC4=1/w/C4
    UAD = complex(U, 0)
    Z1=complex(0,XL1)
    Z2=R4*complex(R2,-XC1)
    Z3=complex(0,XL2)
    Z4=Z1*Z2*Z3/(Z1*Z2+Z1*Z3+Z2*Z3)
    ZAD=Z4+complex(R1+R4+R5,-XC4-XC2)
    IAD=UAD/ZAD
    I1=IAD
    I5=IAD
    UAB=I1*R1
    UCD=I1*complex(R5,-XC2)
    UEC=I1*complex(R4,-XC4)
    UBE=I1*Z4
    I2=UBE/Z1
    I3=UBE/Z2
    I4=UEC/Z3
    I5=UEC/R4
    UL2=I4*complex(0,XL2)
    UC2=I1*complex(0,-XC2)
    I1MAX=abs(I1)
    I2MAX=abs(I2)
    I3MAX=abs(I3)
    I4MAX=abs(I4)
    I5MAX=abs(I5)
    UADMAX=abs(UAD)
    UABMAX=abs(UAB)
    UC2MAX=abs(UC2)
    UBEMAX=abs(UBE)
    UECMAX=abs(UEC)
    UCDMAX=abs(UCD)
    UL2MAX=abs(UL2)
    FI1=faza(I1)
    FI2=faza(I2)
    FI3=faza(I3)
    FI4=faza(I4)
    FI5=faza(I5)
    FUAD=faza(UAD)
    FUAB=faza(UAB)
    FUC2=faza(UC2)
    FUCD=faza(UCD)
    FUBE=faza(UBE)
    FUEC=faza(UEC)
    FUL2=faza(UL2)
    T=1/f
    t=0
    h=1/(f*n)
    Img1,Img2,Img3,Img4,Img5=[],[],[],[],[]
    UmgAD,UmgAB,UmgC2,UmgBE,UmgEC,UmgCD,UmgL2=[],[],[],[],[],[],[]
    Idei = []
    for i in range(1, n + 2):
        Img1.append(I1MAX * math.sin(w * t + FI1))
        Img2.append(I2MAX * math.sin(w * t + FI2))
        Img3.append(I3MAX * math.sin(w * t + FI3))
        Img4.append(I4MAX * math.sin(w * t + FI4))
        Img5.append(I5MAX * math.sin(w * t + FI5))
        UmgAD.append(UADMAX * math.sin(w * t + FUAD))
        UmgAB.append(UABMAX * math.sin(w * t + FUAB))
        UmgC2.append(UC2MAX * math.sin(w * t + FUC2))
        UmgBE.append(UBEMAX * math.sin(w * t + FUBE))
        UmgEC.append(UECMAX * math.sin(w * t + FUEC))
        UmgCD.append(UCDMAX * math.sin(w * t + FUCD))
        UmgL2.append(UL2MAX * math.sin(w * t + FUL2))
        t+=h
    F1 = lambda t: (I1MAX * math.sin(w * t + FI1)) ** 2
    F2 = lambda t: (I2MAX * math.sin(w * t + FI2)) ** 2
    F3 = lambda t: (I3MAX * math.sin(w * t + FI3)) ** 2
    F4 = lambda t: (I4MAX * math.sin(w * t + FI4)) ** 2
    F5 = lambda t: (I5MAX * math.sin(w * t + FI5)) ** 2
    F8 = lambda t: (UADMAX * math.sin(w * t + FUAD)) ** 2
    F9 = lambda t: (UABMAX * math.sin(w * t + FUAB)) ** 2
    F10 = lambda t: (UC2MAX * math.sin(w * t + FUC2)) ** 2
    F11 = lambda t: (UBEMAX * math.sin(w * t + FUBE)) ** 2
    F12 = lambda t: (UECMAX * math.sin(w * t + FUEC)) ** 2
    F13 = lambda t: (UCDMAX * math.sin(w * t + FUCD)) ** 2
    F14 = lambda t: (UL2MAX * math.sin(w * t + FUL2)) ** 2
    Id1 = trapezoid(F1, 0, T, n)
    Id2 = trapezoid(F2, 0, T, n)
    Id3 = trapezoid(F3, 0, T, n)
    Id4 = trapezoid(F4, 0, T, n)
    Id5 = trapezoid(F5, 0, T, n)
    UdAD = trapezoid(F8, 0, T, n)
    UdAB = trapezoid(F9, 0, T, n)
    UdC2=trapezoid(F10, 0, T, n)
    UdBE = trapezoid(F11, 0, T, n)
    UdEC = trapezoid(F12, 0, T, n)
    UdCD = trapezoid(F13, 0, T, n)
    UdL2 = trapezoid(F14, 0, T, n)
    Idei.append(( Id1, Id2, Id3, Id4, Id5, 
                  UdAD, UdAB,UdBE, UdEC, UdCD,
                  UdAD, UdL2,UdC2))
    from prettytable import PrettyTable
    mytable = PrettyTable()
    mytable.field_names = ["параметр", "МАХ", "Действующие", "Угол"]
    mytable.add_row(["I1", I1MAX, Id1, FI1])
    mytable.add_row(["I2", I2MAX, Id2, FI2])
    mytable.add_row(["I3", I3MAX, Id3, FI3])
    mytable.add_row(["I4", I4MAX, Id4, FI4])
    mytable.add_row(["I5", I5MAX, Id5, FI5])
    mytable.add_row(["UAD", UADMAX, UdAD, FUAD])
    mytable.add_row(["UAB", UABMAX, UdAB, FUAB])
    mytable.add_row(["UBE", UBEMAX, UdBE, FUBE])
    mytable.add_row(["UEC", UECMAX, UdEC, FUEC])
    mytable.add_row(["UCD", UCDMAX, UdCD, FUCD])
    mytable.add_row(["UL2", UL2MAX, UdL2, FUL2])
    mytable.add_row(["UC2", UC2MAX, UdC2, FUC2])
    print(mytable)
    if iz_Q != 0:
        QQQ1.append(Id1)
        QQQ2.append(Id4)
        QQQ3.append(UdC2)
        QQQ4.append(UdL2)
        QQQ5.append(UdAD)
            
        t = np.linspace(0, 40, 41)
        Toki = {'Img1': Img1, 'Img2': Img2, 'Img3': Img3, 'Img4': Img4,'Img5': Img5}
        Napr1 = {'UmgAD': UmgAD, 'UmgAB': UmgAB, 'UmgBE': UmgBE, 'UmgEC': UmgEC, 'UmgCD': UmgCD}
        Napr2 = {'UmgAD': UmgAD, 'UmgL2': UmgL2,'UmgC2': UmgC2, }
        data1 = pd.DataFrame(data=Toki)
        data2 = pd.DataFrame(data=Napr1)
        data3 = pd.DataFrame(data=Napr2)
        print(data1)
        print(data2)
        print(data3)
        plt.figure(figsize=(12, 12))
        plt.subplot(2, 2, 1)
        plt.plot(t, Img1, t, Img2, t, Img3, t, Img4, t, Img5)
        plt.title("Осцилограмма токов")
        plt.xlabel("t")
        plt.ylabel("I")
        plt.grid(True)
        plt.subplot(2, 2, 2)
        plt.plot(t, UmgAD, t, UmgAB,  t, UmgBE, t, UmgEC, t, UmgCD)
        plt.title("Осцилограмма напряжений")
        plt.xlabel("t")
        plt.ylabel("U")
        plt.grid(True)
        plt.subplot(2, 2, 3)
        plt.plot(t, UmgAD, t, UmgL2, t, UmgC2)
        plt.title("Осцилограмма напряжений на элементах")
        plt.xlabel("t")
        plt.ylabel("U")
        plt.grid(True)
        plt.show()
        
    if iz_Q == 0:
           print(' ')
    else:
        for i in Idei:
            print('  |  '.join(['{:.5f}'.format(float(x)) for x in i]))
    iz_Q = 1
    U = (Q1 + hQ * sch) 
    sch += 1
t = np.linspace(0, 40, 41)
plt.figure(2)
plt.plot(t, QQQ1, t, QQQ2)
plt.figure(3)
plt.plot(t, QQQ4, t, QQQ3, t, QQQ5)
plt.show()
