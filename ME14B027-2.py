
# coding: utf-8

# In[1]:

X=[[]]
g=[]


# In[2]:

#input via text file. It is tedius to enter values when prompted.
with open('data2.txt') as f:
    Data = []
    for line in f:
        line = line.split()
        if line:
            line = [float(x) for x in line]
            Data.append(line)


# In[3]:

n=Data[0][0]
nActive=Data[0][1]
g=Data[1]
i=2
while i<n+2:
    X[0].append(Data[i])
    i+=1


# In[4]:

#List of velocities of all point vortices
def U(g,X,n):
    pi=3.14159
    i=0
    f=[]
    while i<n:
        j=0
        ux=0
        uy=0
        while j<n:
            if j==i:
                j+=1
                continue
            else:
                r2=(X[i][1]-X[j][1])*(X[i][1]-X[j][1])+(X[i][0]-X[j][0])*(X[i][0]-X[j][0])
                ux+=-g[j]*(X[i][1]-X[j][1])/(2*pi*r2)
                uy+=g[j]*(X[i][0]-X[j][0])/(2*pi*r2)
            j+=1
        f.append([ux,uy])
        i+=1
        
    return f


# In[5]:

#Elementwise add two lists
def AddLists(a,b,n):
    c=[]
    i=0
    while i<n:
        c.append([0.0,0.0])
        i+=1

    i=0
    while i<n:
        j=0
        while j<2:
            c[i][j]=a[i][j]+b[i][j]
            j+=1
        i+=1
        
    return c


# In[6]:

#RK4 method
def Rk4Scheme(g,n,X0,U,dt):
    Xf=[]
    i=0
    while i<n:
        Xf.append([0.0,0.0])
        i+=1

    k1=U(g,X0,n)
    X1=AddLists(X0,k1,n)
    k2=U(g,X1,n)
    X2=AddLists(X1,k2,n)
    k3=U(g,X2,n)
    X3=AddLists(X2,k3,n)
    k4=U(g,X3,n)
    i=0
    while i<n:
        j=0
        while j<2:
            Xf[i][j]=X0[i][j]+(dt/6.0)*(k1[i][j]+2*k2[i][j]+2*k3[i][j]+k4[i][j])
            j+=1
        i+=1
        
    return Xf


# In[7]:

#Euler method
def EulerScheme(g,n,X0,U,dt):
    Xf=[]
    i=0
    while i<n:
        Xf.append([0.0,0.0])
        i+=1

    k1=U(g,X0,n)
    i=0
    while i<n:
        j=0
        while j<2:
            Xf[i][j]=X0[i][j]+dt*k1[i][j]
            j+=1
        i+=1
        
    return Xf


# In[8]:

#Finding trajectory of the point vortices
i=1
dt=0.001
N=1000000
while i<N:
    #X.append(Rk4Scheme(g,n,X[-1],U,dt))
    X.append(EulerScheme(g,n,X[-1],U,dt))
    i+=1


# In[9]:

xt=[]
yt=[]
i=0
while i<n:
    j=0
    xt.append([])
    yt.append([])
    while j<N:
        xt[i].append(X[j][i][0])
        yt[i].append(X[j][i][1])
        j+=1
    i+=1


# In[10]:

import matplotlib.pyplot as plt


# In[11]:

#plotting trajectories and path lines
i=0
while i<n:
    plt.plot(xt[i],yt[i])
    plt.show()
    i+=1


# In[12]:

#velocity at a point at an instant
def Velocity(x,y,X,g,n):
    i=0
    ux=0
    uy=0
    pi=3.14159
    while i<n:
        r2=(x-X[i][0])*(x-X[i][0])+(y-X[i][1])*(y-X[i][1])
        if r2>0.01:
            ux+=-g[i]*(y-X[i][1])/(2*pi*r2)
            uy+=g[i]*(x-X[i][0])/(2*pi*r2)
        i+=1
    return [ux,uy]


# In[13]:

#linear combination of positions of point vortices.
def linComb(X,n,ns,j):
    i=0
    a=[]
    k=0
    while k<n:
        a.append(1)
        k+=1
    a[int((j-1)%n)]+=1+(j-1)//n
    xTempor=[0,0]
    norm=0
    while i<n:
        xTempor[0]+=a[i]*X[i][0]
        xTempor[1]+=a[i]*X[i][1]
        norm+=a[i]
        i+=1
    xTempor[0]/=norm
    xTempor[1]/=norm
    return xTempor


# In[14]:

#Draw stream lines originating from a point a specific time
def streamLines(g,X,n,t,xi,ns,Ns,dt):
    x=[]
    y=[]
    i=0
    while i<ns:
        x.append([xi[i][0]])
        y.append([xi[i][1]])
        j=0
        while j<Ns:
            u=Velocity(x[i][-1],y[i][-1],X[t],g,n)
            x[i].append(x[i][-1]+u[0]*dt)
            y[i].append(y[i][-1]+u[1]*dt)
            j+=1
        plt.plot(x[i],y[i])
        i+=1
    plt.show()


# In[15]:

xi=[]
t=750000
j=0
ns=7
while j<ns:
    temp=linComb(X[t],nActive,ns,j)
    xi.append(temp) 
    j+=1
Ns=1000000
dt=0.001
streamLines(g,X,nActive,t,xi,ns,Ns,dt)


# In[ ]:



