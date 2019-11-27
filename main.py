import math
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

f1 = open("dane/data_spf1.txt", "r")
f2 = open("dane/data_spf2.txt", "r")
gx1 = open("dane/sp2_F_GX3_1.jpg", 'r')
gx2 = open("dane/sp2_F_GX3_2.jpg", 'r')

f1_data = [[],[],[],[],[],[],[]]
f2_data = [[],[],[],[],[],[],[]]
gx1_data = [[],[],[],[],[],[],[]]
gx2_data = [[],[],[],[],[],[],[]]


def quat_init(phi, theta, psi):
    q0 = math.cos(psi/2)*math.cos(theta/2)*math.cos(phi/2)+math.sin(psi/2)*math.sin(theta/2)*math.sin(phi/2)
    q1 = math.cos(psi/2)*math.cos(theta/2)*math.sin(phi/2)-math.sin(psi/2)*math.sin(theta/2)*math.cos(phi/2)
    q2 = math.cos(psi/2)*math.sin(theta/2)*math.cos(phi/2)+math.sin(psi/2)*math.cos(theta/2)*math.sin(phi/2)
    q3 = math.sin(psi/2)*math.cos(theta/2)*math.cos(phi/2)-math.cos(psi/2)*math.sin(theta/2)*math.sin(phi/2)
    return [q0, q1, q2, q3]
    
def Vlength(tab):
    sum = 0
    for i in range(len(tab)):
        sum += tab[i]*tab[i]
    return math.sqrt(sum)

def Qb(w, dt):
    w_ = Vlength(w)
    phi = w_ * dt
    q0 = math.cos(phi/2)
    q1 = w[0]/w_ * math.sin(phi/2)
    q2 = w[1]/w_ * math.sin(phi/2)
    q3 = w[2]/w_ * math.sin(phi/2)
    return [q0, q1, q2, q3]
    
def mulQ(Q1, Q2):
    q0 = Q1[0]*Q2[0]-Q1[1]*Q2[1]-Q1[2]*Q2[2]-Q1[3]*Q2[3]
    q1 = Q1[0]*Q2[1]+Q1[1]*Q2[0]+Q1[2]*Q2[3]-Q1[3]*Q2[2]
    q2 = Q1[0]*Q2[2]-Q1[1]*Q2[3]+Q1[2]*Q2[0]+Q1[3]*Q2[1]
    q3 = Q1[0]*Q2[3]+Q1[1]*Q2[2]-Q1[2]*Q2[1]+Q1[3]*Q2[0]
    return [q0, q1, q2, q3]

def Cmatrix(Q):
    return [[2*(Q[0]*Q[0]+Q[1]*Q[1])-1, 2*(Q[1]*Q[2]-Q[0]*Q[3]), 2*(Q[1]*Q[3]+Q[0]*Q[2])],[2*(Q[1]*Q[2]+Q[0]*Q[3]), 2*(Q[0]*Q[0]+Q[2]*Q[2])-1, 2*(Q[2]*Q[3]-Q[0]*Q[1])],[2*(Q[1]*Q[3]-Q[0]*Q[2]), 2*(Q[2]*Q[3]+Q[0]*Q[1]), 2*(Q[0]*Q[0]+Q[3]*Q[3])-1]]

def Qangles(Q):
    Mat = Cmatrix(Q)
    theta = math.asin(-Mat[0][2])
    phi = -2*math.atan(Mat[1][2]/(Mat[2][2]+math.cos(theta)))
    psi = -2*math.atan(Mat[0][1]/(Mat[0][0]+math.cos(theta)))
    return [phi, theta, psi]

def mulKin(angles, pqr):
    dphi = pqr[0]+math.sin(angles[0])*math.tan(angles[1])*pqr[1]+math.cos(angles[0])*math.tan(angles[1])*pqr[2]
    dtheta = math.cos(angles[0])*pqr[1]-math.sin(angles[0])*pqr[2]
    dpsi = math.sin(angles[0])/math.cos(angles[1])*pqr[1]+math.cos(angles[0])/math.cos(angles[1])*pqr[2]
    return (dphi, dtheta, dpsi)

def GenAB(tab):
    #przyjmuje U jako tab[0] i a jako tab[1]
    A = tab[1][0]/(tab[0][0]-tab[0][1])+tab[1][1]/(tab[0][1]-tab[0][0])
    B = tab[1][0]*tab[0][1]/(tab[0][1]-tab[0][0])+tab[1][1]*tab[0][0]/(tab[0][0]-tab[0][1])
    return [A, B]
    
def FindIndex(tab, row, val):
    indx = 0
    for i in tab[row]:
        if(i<val):
            indx+=1
        else:
            break
    return int(indx)

def GenAverage(tab, start, end, row):
    st = FindIndex(tab, row[0], start)
    ed = FindIndex(tab, row[0], end)
    return np.average(tab[row[1]][st:ed])
    
def Integral(tab, start, end, row):
    st = FindIndex(tab, row[0], start)
    ed = FindIndex(tab, row[0], end)
    sum = 0
    for i in range(st, ed):
        sum+=0.5*(tab[row[1]][i]+tab[row[1]][i+1])*(tab[row[0]][i+1]-tab[row[0]][i])
    return sum    

for line in f1:
	ll = line.split("\t")
	for i in range(len(ll)):
		f1_data[i].append(float(ll[i]))

for line in f2:
	ll = line.split("\t")
	for i in range(len(ll)):
		f2_data[i].append(float(ll[i]))

for line in gx1:
	ll = line.split("\t")
	for i in range(len(ll)):
		gx1_data[i].append(float(ll[i]))
		
for line in gx2:
	ll = line.split("\t")
	for i in range(len(ll)):
		gx2_data[i].append(float(ll[i]))
		
f1.close()
f2.close()
gx1.close()
gx2.close()

t_scale_wait = 2.37#czas oczekiwania na pomiar skali

#skalowanie przyspieszenia
f1_file_scale = open("output/f1_skale.txt", 'w')
f1_file_scale.write("oś\tA\tB\tU_avg_min\tU_avg_max\n")
f1_a_x_scale = [[GenAverage(f1_data, 0, 2.37, [0,1]), GenAverage(f1_data, 4.5, 7.21, [0,1])], [-9.81, 9.81]]
ab = GenAB(f1_a_x_scale)
f1_file_scale.write("x\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 0, 2.37, [0,1]), GenAverage(f1_data, 4.5, 7.21, [0,1])))
f1_a_y_scale = [[GenAverage(f1_data, 8.67, 10.40, [0,2]), GenAverage(f1_data, 12.84, 16.43, [0,2])], [-9.81, 9.81]]
ab = GenAB(f1_a_y_scale)
f1_file_scale.write("y\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 8.67, 10.40, [0,2]), GenAverage(f1_data, 12.84, 16.43, [0,2])))
f1_a_z_scale = [[GenAverage(f1_data, 18.08, 19.81, [0,3]), GenAverage(f1_data, 23.24, 27.0, [0,3])], [-9.81, 9.81]]
ab = GenAB(f1_a_z_scale)
f1_file_scale.write("z\t{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\n".format(ab[0], ab[1], GenAverage(f1_data, 18.08, 19.81, [0,3]), GenAverage(f1_data, 23.24, 27.0, [0,3])))
f1_file_scale.close()

#skalowanie giroskopów
f2_file_scale = open("output/f2_skala.txt", 'w')
f2_times = [[[0.0, 11.35], [0.0, 11.35], [14.78, 31.46]], [[39.17, 45.91], [19.0, 28.27], [1.88, 9.42]]]
f2_Uoi_w = []#[GenAverage(f2_data, f2_times[0][0], 11.35, [0,4]), GenAverage(f2_data, 0.0, 11.35, [0,5]), GenAverage(f2_data, 14.78, 31.46, [0,6])]
for i in range(len(f2_times[0])):
    f2_Uoi_w.append(GenAverage(f2_data, f2_times[0][i][0], f2_times[0][i][1], [0, i+4]))
f2_file_scale.write("U_ox: {0:.4f}\nU_oy: {1:.4f}\nU_oz: {2:.4f}\n".format(f2_Uoi_w[0], f2_Uoi_w[1], f2_Uoi_w[2]))
f2_rescaled_vals = []
for i in range(4):
    f2_rescaled_vals.append([])
    for o in range(len(f2_data[4])):
        if(i>0):
            f2_rescaled_vals[i].append(f2_data[i+3][o]-f2_Uoi_w[i-1])
        else:
            f2_rescaled_vals[i].append(f2_data[i][o])
        #print(f2_rescaled_vals[i][o])
f2_ki = []
for i in range(1, 4):
    print("{0}, {1}\n".format(f2_times[1][i-1][0], f2_times[1][i-1][1]))
    intg = Integral(f2_rescaled_vals, f2_times[1][i-1][0], f2_times[1][i-1][1], [0, i])
    print(intg)
    f2_ki.append(90.0/intg)
    f2_file_scale.write("k{0}: {1:.4f}\n".format(i, f2_ki[i-1]))
f2_file_scale.close()


#obliczenia 3DM-GX2
gx1_file = open("output/gx2_wyniki.txt", 'w')
gx1_0_time = 3.22
gx1_ang_theta0 = math.asin(GenAverage(gx1_data, 0, 3.22, [0, 1]))*180.0/math.pi
gx1_ang_phi0 = math.asin(-GenAverage(gx1_data, 0, 3.22, [0, 2])/math.cos(gx1_ang_theta0))*180.0/math.pi
gx1_file.write("theta_0 = {0:.4f}\tphi_0 = {0:.4f}".format(gx1_ang_theta0, gx1_ang_phi0))
gx1_file.close()

#IMU-ZAIOL-01
#przyspieszenie
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(f1_data[0], f1_data[1], label="x")
plt.plot(f1_data[0], f1_data[2], label="y")
plt.plot(f1_data[0], f1_data[3], label="z")
#plt.plot([f1_data[0][indx], f1_data[0][indx]], [-300, 300], label = "granica pomiaru do skalowania")
plt.grid(True)
plt.title("Wykres odczytanego napięcia na IMU-ZAIOL-01 w funkcji czasu (składowa \"$a$\").")
plt.ylabel("U [V]")
plt.xlabel("t [s]")
plt.ylim(1.85, 3)
plt.legend(loc = "lower right")
plt.savefig("output/IMU-ZAIOL-01_1_a_U_t.png")

#obrót
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(f2_data[0], f2_data[4], label="oś x")
plt.plot(f2_data[0], f2_data[5], label="oś y")
plt.plot(f2_data[0], f2_data[6], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanego napięcia na IMU-ZAIOL-01 w funkcji czasu (składowa \"$\omega$\").")
plt.ylabel(r"U $\left [ V \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/IMU-ZAIOL-01_2_U_w_t.png")

#3DM-GX2
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx1_data[0], gx1_data[1], label="oś x")
plt.plot(gx1_data[0], gx1_data[2], label="oś y")
plt.plot(gx1_data[0], gx1_data[3], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanych wartości na 3DM-GX2 w funkcji czasu (składowa \"$a$\").")
plt.ylabel(r"$\frac{a}{g} \left [ - \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_1_a.png")

plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx1_data[0], gx1_data[4], label="oś x")
plt.plot(gx1_data[0], gx1_data[5], label="oś y")
plt.plot(gx1_data[0], gx1_data[6], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanych wartości na 3DM-GX2 w funkcji czasu (składowa \"$\omega$\").")
plt.ylabel(r"$\omega \left [ \frac{m}{s^2} \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_1_w.png")
#plt.show()

plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx2_data[0], gx2_data[1], label="oś x")
plt.plot(gx2_data[0], gx2_data[2], label="oś y")
plt.plot(gx2_data[0], gx2_data[3], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanych wartości na 3DM-GX2 w funkcji czasu (składowa \"$a$\").")
plt.ylabel(r"$\frac{a}{g} \left [ - \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_2_a.png")

plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx2_data[0], gx2_data[4], label="oś x")
plt.plot(gx2_data[0], gx2_data[5], label="oś y")
plt.plot(gx2_data[0], gx2_data[6], label="oś z")
plt.grid(True)
plt.title("Wykres odczytanych wartości na 3DM-GX2 w funkcji czasu (składowa \"$\omega$\").")
plt.ylabel(r"$\omega \left [ \frac{m}{s^2} \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_2_w.png")
#plt.show()

#orientacja przestrzenna
gx1_file_quat = open("output/gx1_quat.txt", 'w')
gx1_ang_theta0 = math.asin(GenAverage(gx1_data, 0, 3.22, [0, 1]))
gx1_ang_phi0 = math.asin(-GenAverage(gx1_data, 0, 3.22, [0, 2])/math.cos(gx1_ang_theta0))
qa = quat_init(gx1_ang_phi0, gx1_ang_theta0, 0)
gx1_angles = [[],[],[]]
gx1_time_passed = []

for i in range(len(gx1_data[0])-1):
    w = [gx1_data[4][i], gx1_data[5][i], gx1_data[6][i]]
    dt = gx1_data[0][i+1]-gx1_data[0][i]
    qb = Qb(w, dt)
    qq = mulQ(qa, qb)
    gx1_angles_a = Qangles(qq)
    gx1_angles[0].append(gx1_angles_a[0]*180.0/math.pi)
    gx1_angles[1].append(gx1_angles_a[1]*180.0/math.pi)
    gx1_angles[2].append(gx1_angles_a[2]*180.0/math.pi)
    qa=qq
    if(i>1):
        gx1_time_passed.append(gx1_time_passed[i-1]+dt)
    else:
        gx1_time_passed.append(dt)
    gx1_file_quat.write("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(gx1_time_passed[i], gx1_angles[0][i], gx1_angles[1][i], gx1_angles[2][i]))

gx1_file_quat.close()

#kinematyka
gx1_file_kin = open("output/gx1_kin.txt", 'w')

gx1_time_passed_kin = [0]
gx1_angles_kin = [[gx1_ang_phi0],[gx1_ang_theta0],[0]]

for i in range(1, len(gx1_data[0])):
    gx1_angles_a = [gx1_angles_kin[0][i-1], gx1_angles_kin[1][i-1], gx1_angles_kin[2][i-1]]
    gx1_pqr_a = [gx1_data[4][i], gx1_data[5][i], gx1_data[6][i]]
    gx1_matrix_kin = mulKin(gx1_angles_a, gx1_pqr_a)
    dt = gx1_data[0][i]-gx1_data[0][i-1]
    for o in range(3):
        gx1_angles_kin[o].append(gx1_angles_kin[o][i-1]+gx1_matrix_kin[o]*dt)
    gx1_time_passed_kin.append(gx1_time_passed[i-1]+dt)
    gx1_file_kin.write("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(gx1_time_passed_kin[i], gx1_angles_kin[0][i], gx1_angles_kin[1][i], gx1_angles_kin[2][i]))

for i in range(len(gx1_angles_kin[0])):
    for o in range(3):
        gx1_angles_kin[o][i] *= 180.0/math.pi

gx1_file_kin.close()

#rysowanie wykresu 3d kwaternionów
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx1_time_passed, gx1_angles[0], label="$\phi$")
plt.plot(gx1_time_passed, gx1_angles[1], label=r"$\theta$")
plt.plot(gx1_time_passed, gx1_angles[2], label="$\psi$")
plt.grid(True)
plt.title("Wykres kątów orientacji obliczonych kwaternionami w funkcji czasu.")
plt.ylabel(r"kąt $\left [ ^\circ \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_2_ang.png")

plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx1_time_passed_kin, gx1_angles_kin[0], label="$\phi$")
plt.plot(gx1_time_passed_kin, gx1_angles_kin[1], label=r"$\theta$")
plt.plot(gx1_time_passed_kin, gx1_angles_kin[2], label="$\psi$")
plt.grid(True)
plt.title("Wykres kątów orientacji obliczonych równaniami kinematycznymi w funkcji czasu.")
plt.ylabel(r"kąt $\left [ ^\circ \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_2_ang_kin.png")

#ruch swobodny
gx2_file_quat = open("output/gx2_quat.txt", 'w')
gx2_ang_theta0 = math.asin(GenAverage(gx2_data, 0, 3.22, [0, 1]))
gx2_ang_phi0 = math.asin(-GenAverage(gx2_data, 0, 3.22, [0, 2])/math.cos(gx2_ang_theta0))
qa = quat_init(gx2_ang_phi0, gx2_ang_theta0, 0)
gx2_angles = [[],[],[]]
gx2_time_passed = []

for i in range(len(gx2_data[0])-1):
    w = [gx2_data[4][i], gx2_data[5][i], gx2_data[6][i]]
    dt = gx2_data[0][i+1]-gx2_data[0][i]
    qb = Qb(w, dt)
    qq = mulQ(qa, qb)
    gx2_angles_a = Qangles(qq)
    gx2_angles[0].append(gx2_angles_a[0]*180.0/math.pi)
    gx2_angles[1].append(gx2_angles_a[1]*180.0/math.pi)
    gx2_angles[2].append(gx2_angles_a[2]*180.0/math.pi)
    qa=qq
    if(i>1):
        gx2_time_passed.append(gx2_time_passed[i-1]+dt)
    else:
        gx2_time_passed.append(dt)
    gx2_file_quat.write("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(gx2_time_passed[i], gx2_angles[0][i], gx2_angles[1][i], gx2_angles[2][i]))

gx2_file_quat.close()

#kinematyka
gx2_file_kin = open("output/gx2_kin.txt", 'w')

gx2_time_passed_kin = [0]
gx2_angles_kin = [[gx2_ang_phi0],[gx2_ang_theta0],[0]]

for i in range(1, len(gx2_data[0])):
    gx2_angles_a = [gx2_angles_kin[0][i-1], gx2_angles_kin[1][i-1], gx2_angles_kin[2][i-1]]
    gx2_pqr_a = [gx2_data[4][i], gx2_data[5][i], gx2_data[6][i]]
    gx2_matrix_kin = mulKin(gx2_angles_a, gx2_pqr_a)
    dt = gx2_data[0][i]-gx2_data[0][i-1]
    for o in range(3):
        gx2_angles_kin[o].append(gx2_angles_kin[o][i-1]+gx2_matrix_kin[o]*dt)
    gx2_time_passed_kin.append(gx2_time_passed[i-1]+dt)
    gx2_file_kin.write("{0:.3f}\t{1:.3f}\t{2:.3f}\t{3:.3f}\n".format(gx2_time_passed_kin[i], gx2_angles_kin[0][i], gx2_angles_kin[1][i], gx2_angles_kin[2][i]))

for i in range(len(gx2_angles_kin[0])):
    for o in range(3):
        gx2_angles_kin[o][i] *= 180.0/math.pi

gx2_file_kin.close()

#przesunięcie

#rysowanie wykresu 3d kwaternionów
plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx2_time_passed, gx2_angles[0], label="$\phi$")
plt.plot(gx2_time_passed, gx2_angles[1], label=r"$\theta$")
plt.plot(gx2_time_passed, gx2_angles[2], label="$\psi$")
plt.grid(True)
plt.title("Wykres kątów orientacji obliczonych kwaternionami w funkcji czasu.")
plt.ylabel(r"kąt $\left [ ^\circ \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_3_ang.png")

plt.figure(figsize = (25/2.54, 20/2.54))
plt.subplot(111)
plt.plot(gx2_time_passed_kin, gx2_angles_kin[0], label="$\phi$")
plt.plot(gx2_time_passed_kin, gx2_angles_kin[1], label=r"$\theta$")
plt.plot(gx2_time_passed_kin, gx2_angles_kin[2], label="$\psi$")
plt.grid(True)
plt.title("Wykres kątów orientacji obliczonych równaniami kinematycznymi w funkcji czasu.")
plt.ylabel(r"kąt $\left [ ^\circ \right ]$")
plt.xlabel("t [s]")
plt.legend(loc = "lower right")
plt.savefig("output/3DM-GX2_3_ang_kin.png")
plt.show()