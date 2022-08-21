from feon.sa import *
import matplotlib.pyplot as plt

def verify_value(value):
    def is_int(n):
        try:
            int(n)
            return True
        except ValueError:
            return False
        else:
            return int(n).is_integer()
    
    while is_int(value)==False:
        value = input("Enter a correct value: ")
    return int(value)

def verify_value_0(value):
    def is_float(n):
        try:
            float(n)
            return True
        except ValueError:
            return False
        else:
            return float(n).is_float()
    
    while is_float(value)==False or float(value)<=0:
        value = input("Enter a correct value: ")
    return float(value)

def verify_value_00(value):
    
    def is_float(n):
        try:
            float(n)
            return True
        except ValueError:
            return False
        else:
            return float(n).is_integer()
    
    while is_float(value)==False:
        value = input("Enter a correct value: ")
    return float(value)

def verify_value_01(value):
    def is_float(n):
        try:
            float(n)
            return True
        except ValueError:
            return False
        else:
            return float(n).is_float()
    
    while is_float(value)==False:
        value = input("Enter a correct value: ")
    return float(value)

def verify_value_1(value):
    def is_int(n):
        try:
            int(n)
            return True
        except ValueError:
            return False
        else:
            return int(n).is_integer()
    
    while is_int(value)==False or int(value)<=1:
        value = input("Enter a correct value: ")
    return int(value)

def verify_value_2(value):
    def is_int(n):
        try:
            int(n)
            return True
        except ValueError:
            return False
        else:
            return int(n).is_integer()
    
    while is_int(value)==False or int(value)<0:
        value = input("Enter a correct value: ")
    return int(value)

def verify_value_3(value):
    def is_int(n):
        try:
            int(n)
            return True
        except ValueError:
            return False
        else:
            return int(n).is_integer()
    
    while is_int(value)==False:
        value = input("Enter a correct value: ")
    return int(value)

def nnod(fill, col):
    from feon.tools import pair_wise
    import feon as  fn
    lista = []
    for i in range(fill):
        for j in range(col):
            lista.append(Node(j, i))
    lista = lista[col:]
    return len(lista)

def k_extract(M, coord):
    L = []
    for i in coord:
        x, y = i
        L.append([int(x-1), int(y)])
    valores = []
    for x, y in L:
        valor = M[x, y]
        valores.append(valor)
    return valores

def lista_nodos(fill, col):
    y = fill-1
    ty = (y*col)
    lista = []
    for i in range(fill):
        for j in range(col):
            lista.append(Node(j, i))
    return lista

def get_values():
    print("Enter the parameters:\n")
    v_I = input("I: ")
    v_I = verify_value_00(v_I)
    v_R = input("R: ")
    v_R = verify_value(v_R)

    v_Z = input("Z: ")
    v_Z = verify_value(v_Z)
    v_n = input("n: ")
    v_n = verify_value(v_n)
    v_Ct = input("Ct: ")
    v_Ct = verify_value(v_Ct)
    v_hn = input("hn: ")
    v_hn = verify_value(v_hn)
    v_a = input("a: ")
    v_a = verify_value(v_a)
    
    v_Fa = input("Fa: ")
    v_Fa = verify_value(v_Fa)
    v_Fd = input("Fd: ")
    v_Fd = verify_value(v_Fd)
    v_Fs = input("Fs: ")
    v_Fs = verify_value(v_Fs)
    
    v_r = input("r: ")
    v_r = verify_value(v_r)
    v_Øp = input("Øp: ")
    v_Øp = verify_value(v_Øp)
    v_Øe = input("Øe: ")
    v_Øe = verify_value(v_Øe)
    
    v_g = input("g: ")
    v_g = verify_value(v_g)
    
    T1 = v_Ct*(v_hn**v_a)
    T2 = 1.3*T1
    t0 = 0.1*v_Fs*(v_Fd//v_Fa)
    tc = 0.55*v_Fs*(v_Fd//v_Fa)
    tl = 2.4*(v_Fd)
    Sa = v_Z*v_Fa
    Sa0 = v_n*v_Z*v_Fa
    f = 1/(v_R*v_Øp*v_Øe)
    
    m_T = [0]
    t_0 = t0
    while t_0<=3.1:
        m_T.append(t_0)
        t_0+=0.01
    
    m_Sa = []
    for i in m_T:
        if i<tc:
            m_Sa.append(Sa0)
        else:
            val = Sa0*(tc/i)**(v_r)
            m_Sa.append(val)
    plt.title("ESPECTRO DEL DISEÑO")
    plt.xlabel('T')
    plt.ylabel('Sa')
    plt.plot(m_T, m_Sa)
    
def mmr(AA, nniv, ncol):

    ID = np.arange(ncol*3,ncol*nniv*3)
    print(AA.shape)
    print()
    npiso = nniv-1
    x, y = AA.shape
    
    first = [ID[0]*i for i in range(1, nniv)]
    
    ID = first + [i for i in ID if i not in first]
    print(ID)
    print()

    for i in range(len(ID)):
        for j in range(len(ID)):
            valor1 = ID[i]
            valor2 = ID[j]
            AA[i,j] = AA[valor1,valor2]
            
    kxx, kxy = AA[0:ncol,0:ncol], AA[ncol:ncol*nniv*3,0:ncol]
    kyx, kyy = kxy.T, AA[ncol:ncol*nniv*3,ncol:ncol*nniv*3]
    k1 = np.dot(np.linalg.inv(kyy), kxy)
    k2 = np.dot(kxy.T, k1)
    k3 = kxx - k2
    return k3

def drawing_estructure(filas, columnas):
    print("Estructura ({}x{}):".format(filas, columnas))
    piso = "╢═══"*(columnas-1)+ "╢"
    edificio = [piso for i in range(filas-1)]
    floor = [j[1:-1].replace("╢", "╦") for j in edificio]
    print("{}{}{}".format("╔",floor[0],"╗"))
    for e in range(len(edificio)-1):
        print(edificio[e])
    floor = [j[1:-1].replace("╢", "╩") for j in edificio]
    print("{}{}{}".format("╚",floor[0],"╝"))

def total_elementos(fill, col):
    
    h = col-1
    th = h*(fill-1)
    
    y = fill-1
    ty = (y*col)
    
    total = th+ty
    
    print("Elementos en x: {} en {} niveles total: {}".format(h, y, th))
    print("Elementos en y: {} en {} niveles total: {}".format(col, y, ty))
    print()
    return total

def k_(M, filas, columnas):
    rango = np.arange(columnas*3, columnas*filas*3)
    ID1 = [rango[i] for i in range(0,len(rango),3)]
    ID2 = [i for i in rango if i not in ID1]
    ID3 = []
            
    for i in range(1,len(rango),3):
        ID3.append(rango[i])
        ID3.append(rango[i+1])

    kxx_coord = []
    for i in ID1:
        for j in ID1:
            kxx_coord.append([i, j])
    
    kxy_coord = []
    for i in ID1:
        for k in ID2:
            kxy_coord.append([i, k])
            
    kyy_coord = []
    for i in ID3:
        for j in ID3:
            kyy_coord.append([i, j])
            
    kxx = np.array(k_extract(M, kxx_coord))
    kxy = np.array(k_extract(M, kxy_coord))
    kyy = np.array(k_extract(M, kyy_coord))
    
    nnodos = nnod(filas, columnas)
    kxx = kxx.reshape((nnodos, nnodos))
    kxy = kxy.reshape((nnodos, nnodos*2))
    kyx = kxy.T
    kyy = kyy.reshape((nnodos*2, nnodos*2))
    
    k1 = np.dot(np.linalg.inv(kyy), kxy.T)
    k2 = np.dot(kxy, k1)
    k3 = kxx - k2
    return k3

def m_r(E, A, I, Ncol, Nniv, SC, dic):

    Nele = (Nniv-1)*(Ncol)+(Nniv-1)*(Ncol-1)  
    H = [i for i in range(0,Nniv)]     
    L = [j for j in range(0,Ncol)]

    Nluz = (Ncol -1) 
    Npiso = (Nniv -1)
    Nnod = (Ncol*Nniv)
    Nele_viga = Nluz*Npiso
    Nele_col = Ncol*Npiso
    
    ID=[]                      
    contador = 0
    for i in range(1,Npiso+1):
        for j in range(Nluz):
            a = (i)*(Ncol*3)+(3*j)
            ID.append([a+1,a+2,a+3,a+4,a+5,a+6])
            contador = contador + 1

    for j in range(1,Npiso+1):
        for i in range (1,Ncol+1):
            a = (3*i-2)+(j-1)*(Ncol*3)
            b = (3*i-2)+(Ncol*3*j)
            ID.append([a+1,a,a+2,b+1,b,b+2])
            
    n = []
    for i in range(0,Nniv):
        for j in range(0,Ncol):
            nnod= j+Ncol*i
            x =L[j]
            y =H[i]
            n.append(Node(x,y))

    v = []
    for j in range(1,Nniv):
        for i in range(Ncol-1):
            v.append(Beam2D11((n[j*Ncol+i],n[j*Ncol+i+1]),E,A,I))    
    c = []
    for i in range(1,Nniv):
        for j in range(Ncol):
            c.append(Beam2D11((n[(i-1)*Ncol+j],n[Ncol*(i-1)+j+Ncol]),E,A,I))

    s = System()
    s.add_nodes(n)
    s.add_elements(v,c)
    
    for Nc in dic:
        fx, fy = dic[Nc]
        s.add_node_force(Nc,Fx = fx,Fy = fy)

    s.add_fixed_sup(SC)
    matriz_de_rigidez = s.KG
    s.solve()
    
    return matriz_de_rigidez, n

def ujns_drifts(FX, KOND, MM, lista_h):
    Lambda= np.linalg.solve(KOND,MM)
    Lambda = np.absolute(Lambda) #Periodos 
    T1 = np.pi*2/Lambda[0] 
    T2 = np.pi*2/Lambda[1] 
    s = 0.02 
    a0 = 2*s*(T1**0.5*T2**0.5)/(T1**0.5+T2**0.5) 
    a1 = 2*s*1/(T1**0.5+T2**0.5) 
    C = a0*MM + al*kconde 
    Ca = a0/(Lambda*2) 
    m1 = Lambda.T*m*np.linalg.inv(Lambda) 
    Lnl = Lambda.T*MM
    Mn = Lambda.T*m1 
    Tn = np.linalg.solve(Lnl,Mn) 
    Dn = fx/Lambda**2 
    Ujn = Tn*Lambda*Dn
    
    drifts = []
    for i in range(len(h)):
        diff = Ujn[i+1, i] - Ujn[i, i]
        drift = diff/alturas[i]
        drifts.append(drift)
    
    return drifts