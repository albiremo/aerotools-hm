
#import
#-----------------------------------------------------
import numpy             as np
from   scipy.linalg      import solve,solve_banded
import matplotlib        as mp
mp.use("Qt4Agg")
import scipy             as sp
#import matplotlib.pyplot as plt
from   numpy   import pi
#-----------------------------------------------------
#-----------------------------------------------------
#dictionary definition for NACA 5 digit

def zero():
    q=0.0580
    k=361.4
    return q,k
def one():
    q=0.1260
    k=361.4
    return q,k
def two():
    q=0.2025
    k=51.64
    return q,k
def three():
    q=0.2900
    k=15.957
    return(q,k)
def four():
    q=0.3910
    k=3.230
    return q,k
#----------------------------------------------------
#-------------------------------------------------
def NACA4camberline(xc,mc,pos_mc):
    m=pos_mc/100
    p=pos_mc/10
    yc = ((m/(p**2))*(2*p-xc)*xc)*(xc<p)+((m/(1-p)**2)*(1- 2*p + (2*p-xc)*xc))*(xc>=p)
    Dyc = ((m/p**2)*2*(p-xc))*(xc<p)+((m/(1-p)**2) * 2*(p-xc))*(xc>=p)
    return yc,Dyc

def NACA5camberline(xc,tipologia):
    options={210: zero,
    220: one,
    230: two,
    240: three,
    250: four,}
    try:
        q,k=options[tipologia]()
    except KeyError:
        print('code not supported')
    yc = np.zeros([len(xc)])
    Dyc = np.zeros([len(xc)])
    yc = ((k/6)*(xc**3-3*q*xc**2+q**2*(3-q)*xc))*(xc<q) +((k/6)*q**3*(1-xc))*(xc>=q)
    Dyc = ((k/6)*(3*xc**2 -6*q*xc+ q**2*(3-q)*np.ones(len(xc))))*(xc<q)-((k/6)*np.ones(len(xc))*q**3)*(xc>=q)
    return yc,Dyc

def NACAthick(xc,SS,xtec):
    s=0.01*SS
    #classical NACA thickness
    tk = 5*s*(0.29690*np.sqrt(xc) -0.12600*xc -0.35160*(xc**2) + 0.28430*(xc**3) -0.10150*xc**4)
    # is possible to evaluate correction of this coefficients, due to the fact
    # that we need 0 thickness exit
    if xtec<1:
        tkte=tk[-1]
        A1=np.zeros([4,4],np.float64)
        A1 = 5*s*np.matrix([[1.0, 1.0, 1.0, 1.0],[xtec, xtec**2, xtec**3, xtec**4],[1, 2*xtec, 3*xtec**2, 4*xtec**3],[0, 2, 6*xtec, 12*xtec**2]])
        # remember the tonda
        rhs=np.zeros([4],np.float64)
        rhs=np.array([-tkte,0,0,0])
        b=solve(A1,rhs)
        tk = 5*s*(0.29690*np.sqrt(xc) -0.12600*xc -0.35160*xc**2 +\
         0.28430*xc**3 -0.10150*xc**4)*(xc<xtec) +\
        5*s*(0.29690*np.sqrt(xc)+(-0.12600+b[0])*xc +(-0.35160+b[1])*xc**2 \
        + (0.28430+b[2])*xc**3 +(-0.10150+b[3])*xc**4)*(xc>=xtec)
    tk[-1]=0
    return tk



def NACA(code,nc):
    nc=nc+1
    #I've to make 1 more every time, because the last one of arange
    #is not complained
#    xc = 0.5*(1-np.cos(np.arange(0,pi,pi/(nc-1))))
    xc = np.linspace(0,1,num = nc-1,endpoint = True)
    nv = 2*nc-1 #number of panel vertices, must be double of number of nodes but
    #minus 1, because we have 1 more node on chord than number of vertices
    xv = np.zeros([nv],dtype=np.float64)
    yv = np.zeros([nv],dtype=np.float64)
    nn = len(code)
    if nn>5 or nn<4:
        print('error enter a NACA 4 or 5 digit code')
        return
    else:
        if nn==4:
            #4 digit case
            A=list(code)
            mc=np.int(A[0])        #maximum camber
            pos_mc=np.int(A[1])     #position of maximum camber
            SS=np.int(A[2])*10+np.int(A[3])
            print('max_camber:',mc,'pos_max_camber:',pos_mc,'max_thick:',SS)
            xtec=np.float64(1)
             #maximum thickness
            # camber line construction
            yc,Dyc=NACA4camberline(xc,mc,pos_mc)

            #thickness construction
            tk=NACAthick(xc,SS,xtec)
            #print(xv[0:nc-1].shape,xv[nc:nv-1].shape)
            theta=np.arctan(Dyc)
            #print(tk.shape, theta.shape,xc.shape,nv)
            #xv=np.zeros([nv],np.float64)
            #yv=np.zeros([nv],np.float64)
            xv[0:nc-1]=xc[0:nc-1]-tk*np.sin(theta[0:nc-1])
            yv[0:nc-1]=yc[0:nc-1]+tk*np.cos(theta[0:nc-1])
            xv[0:nc-1]=xv[nc-1:0:-1]
            yv[0:nc-1]=yv[nc-1:0:-1]
            xv[nc:nv-1]=xc[1:nc-1]+tk[1:nc-1]*np.sin(theta[1:nc-1])
            yv[nc:nv-1]=yc[1:nc-1]-tk[1:nc-1]*np.cos(theta[1:nc-1])
            xvnew=np.zeros([nv-3],np.float64)
            yvnew=np.zeros([nv-3],np.float64)
            xvnew=xv[1:nv-1]
            yvnew=yv[1:nv-1]
            with open("naca_airfoil.txt","w") as air:
                for i in range(0,nv-2):
                    print(xvnew[i],yvnew[i], file=air, sep=' ')
            return(xvnew,yvnew)
        else:
            #5 digit case
            A=list(code)
            tipologia=np.int(A[0])*100+np.int(A[1])*10+np.int(A[2])
            SS=np.int(A[3])*10+np.int(A[4])
            #camberline construction
            yc,Dyc=NACA5camberline(xc,tipologia)
            #thickness construction
            xtec=1.0
            tk=NACAthick(xc,SS,xtec)
            #print(xv[0:nc-1].shape,xv[nc:nv-1].shape)
            theta=np.arctan(Dyc)
            #print(tk.shape, theta.shape,xc.shape,nv)
            #xv=np.zeros([nv],np.float64)
            #yv=np.zeros([nv],np.float64)
            xv[0:nc-1]=xc[0:nc-1]-tk*np.sin(theta[0:nc-1])
            yv[0:nc-1]=yc[0:nc-1]+tk*np.cos(theta[0:nc-1])
            xv[0:nc-1]=xv[nc-1:0:-1]
            yv[0:nc-1]=yv[nc-1:0:-1]
            xv[nc:nv-1]=xc[1:nc-1]+tk[1:nc-1]*np.sin(theta[1:nc-1])
            yv[nc:nv-1]=yc[1:nc-1]-tk[1:nc-1]*np.cos(theta[1:nc-1])
            xvnew=np.zeros([nv-3],np.float64)
            yvnew=np.zeros([nv-3],np.float64)
            xvnew=xv[1:nv-1]
            yvnew=yv[1:nv-1]
            with open("naca_airfoil.txt","w") as air:
                for i in range(0,nv-2):
                    print(xvnew[i],yvnew[i], file=air, sep=' ')
            return(xvnew,yvnew)

            return(yv)
