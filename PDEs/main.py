import numpy as np
import sys
import os

def f(u,v,mu,sigma,lamba,kk):
    return lamba*np.multiply(u,v);

def g(u,v,mu,sigma,lamba,kk):
    return -(sigma/kk)*np.multiply(v,v) - lamba*np.multiply(u,v);

#This function will be used to compute the values of the coeffcients for the ETD2 method
def compute_values(c, dt, x_points):
    exp_c_dt = np.exp(c * dt)
    ones_array = np.ones((x_points, int(x_points / 2) + 1))

    ff2 = np.divide(exp_c_dt - 1, c, out=dt * ones_array, where=c != 0)
    ff3 = np.divide((1 + c * dt) * exp_c_dt - 1 - 2 * c * dt, dt * c**2, out=1.5 * dt * ones_array, where=c != 0)
    ff4 = np.divide(-exp_c_dt + 1 + c * dt, dt * c**2, out=-0.5 * dt * ones_array, where=c != 0)
    return exp_c_dt, ff2, ff3, ff4

if len(sys.argv) != 13:
    print("Usage: python3 main.py t_start t_end dt x_start x_end x_points mu sigma lambda kk dx dy")
    sys.exit()

np.seterr(over='raise')

t_start = float(sys.argv[1])
t_end = float(sys.argv[2])
dt = float(sys.argv[3])
x_start = float(sys.argv[4])
x_end = float(sys.argv[5])
x_points = int(sys.argv[6])
mu = float(sys.argv[7])
sigma = float(sys.argv[8])
lamba = float(sys.argv[9])
kk = float(sys.argv[10])
ddx = float(sys.argv[11])
ddy = float(sys.argv[12])

os.system("rm output/*")

# Initialiazing variables
t_points = int(np.ceil((t_end-t_start)/dt))
dx = (x_end-x_start)/(x_points)
shape_x=(x_points,x_points)
shape_k=(x_points,int(x_points/2+1))
u_x , v_x, f_x, g_x = [np.zeros(shape_x) for _ in range(4)]
u_k , v_k, f_k, g_k, f_k_prev, g_k_prev = [np.zeros(shape_k,dtype=complex) for _ in range(6)]
c_u, c_v = [np.zeros(shape_k) for _ in range(2)]
x, y = np.meshgrid(np.linspace(x_start, x_end, x_points, endpoint=False), np.linspace(x_start, x_end, x_points, endpoint=False))

#Initial conditions
# u_x[:x_points//2, :x_points//2] = 1
# u_x[x_points//2:, x_points//2:] = 1
# v_x[x_points//2:, :x_points//2] = 1
# v_x[:x_points//2, x_points//2:] = 1
u_x=np.square(np.sin(2*np.pi*x/x_end))+0.1
v_x=np.square(np.cos(2*np.pi*x/x_end)+np.sin(2*np.pi*y/x_end))+0.1
f_x=f(u_x,v_x,mu,sigma,lamba,kk)
g_x=g(u_x,v_x,mu,sigma,lamba,kk)

#Initial Fourier transform
u_k = np.fft.rfft2(u_x,norm="ortho")
f_k=np.fft.rfft2(f_x,norm="ortho")
v_k = np.fft.rfft2(v_x,norm="ortho")
g_k=np.fft.rfft2(g_x,norm="ortho")

#Precomputing wave vectors and constants
kx=np.fft.fftfreq(x_points,d=dx) * 2 * np.pi
ky=np.fft.rfftfreq(x_points,d=dx) * 2 * np.pi
kx, ky = np.meshgrid(kx, ky, indexing='ij')
c_u=-ddx*np.square(kx)-ddy*np.square(ky)-mu
c_v=-ddx*np.square(kx)-ddy*np.square(ky)+sigma
f1, f2, f3, f4 = compute_values(c_u, dt, x_points)
g1, g2, g3, g4 = compute_values(c_v, dt, x_points)

#Outputting data
data = np.column_stack([x.ravel(), y.ravel(), u_x.ravel(), v_x.ravel()])
np.savetxt("output/0_x.dat", data, fmt='%f')
data = np.column_stack([kx.ravel(), ky.ravel(), np.abs(u_k).ravel(), np.abs(v_k).ravel()])
np.savetxt("output/0_k.dat", data, fmt='%f')

print(f"t = 0, RTEy={np.abs(u_k[0,int(x_points/2)])/np.abs(u_k[0,0]):.1e}, RTEx={np.abs(u_k[int(x_points/2),0])/np.abs(u_k[0,0]):.1e},RTExy={np.abs(u_k[x_points-1,int(x_points/2)])/np.abs(u_k[0,0]):.1e}")
counter=1
for i in range(1,t_points+1):
    t = t_start + i*dt
    print(f"t = {t:.2f}, RTEy={np.abs(u_k[0,int(x_points/2)])/np.abs(u_k[0,0]):.1e}, RTEx={np.abs(u_k[int(x_points/2),0])/np.abs(u_k[0,0]):.1e},RTExy={np.abs(u_k[x_points-1,int(x_points/2)])/np.abs(u_k[0,0]):.1e}",np.max(u_x),np.max(v_x),np.min(u_x),np.min(v_x))
    if i==1:
       u_k=u_k*f1 + f2*f_k
       v_k=v_k*g1 + g2*g_k
    else:
        u_k=u_k*f1 +f3*f_k + f4*f_k_prev
        v_k=v_k*g1 +g3*g_k + g4*g_k_prev
    f_k_prev=f_k.copy()
    u_x=np.fft.irfft2(u_k,(x_points,x_points),norm="ortho")
    v_x=np.fft.irfft2(v_k,(x_points,x_points),norm="ortho")
    f_x=f(u_x,v_x,mu,sigma,lamba,kk)
    f_k=np.fft.rfft2(f_x,norm="ortho")
    g_k_prev=g_k.copy()
    g_x=g(u_x,v_x,mu,sigma,lamba,kk)
    g_k=np.fft.rfft2(g_x,norm="ortho")
    if t >= counter:
        data = np.column_stack([x.ravel(), y.ravel(), u_x.ravel(), v_x.ravel()])
        np.savetxt("output/"+str(counter)+"_x.dat", data, fmt='%f')
        data = np.column_stack([kx.ravel(), ky.ravel(), np.abs(u_k).ravel(), np.abs(v_k).ravel()])
        np.savetxt("output/"+str(counter)+"_k.dat", data, fmt='%f')
        counter+=1
