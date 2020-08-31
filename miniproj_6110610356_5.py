#CN101 mini project
#นาย อภิวิชญ์ ภาคีสุข 6110610356
#==========================LU-Factorization==================================
def LU():#ฟังก์ชันที่ใช้แก้ระบบสมการเชิงเส้นโดยวิธีการแยกLU
    print('Fine X with LU-Factorization\n          AX = B')   
    n = int(input('Enter scale(rows and colums) of Matrix A : '))
    A = []
    for i in range(0,n):
        A.append([])
    for i in range(0,n):
        for j in range(0,n):
            A[i].append(float(input('Enter value of rows '+str(i+1)+' colums '+str(j+1)+' of Matrix A : ')))
    B = []
    for i in range(0,n):
        B.append(float(input('Enter value of rows '+str(i+1)+' of Matrix B : ')))
    L = []
    for i in range(0,n):
        L.append([])
    for i in range(0,n):
        for j in range(0,n):
            if i != j:
                L[i].append(0)
            else:
                L[i].append(1)
    U = []
    for i in range(0,n):
        U.append([])
    for i in range(0,n):
        for j in range(0,n):
            U[i].append(0)
    Y = []
    for i in range(0,n):
        Y.append(0)
    X = []
    for i in range(0,n):
        X.append(0)
    for j in range(0,n):
        U[0][j] = A[0][j]
    for i in range(0,n):
        L[i][0] = A[i][0]/U[0][0]
    for i in range(1,n):
        for j in range(1,n):
            Sum = 0
            if j<i:
                for k in range(0,j):
                    Sum += L[i][k]*U[k][j]
                if U[j][j] == 0:
                    print('\nNo solution of system of linear equations')
                    print('------------------------------------------------------------\n')
                    return main()
                L[i][j] = (A[i][j]-Sum)/U[j][j]
            else:
                for k in range(0,i):
                    Sum += L[i][k]*U[k][j]
                U[i][j] = A[i][j]-Sum
    Y[0] = B[0]
    for i in range(1,n):
        Sum = 0
        for j in range(0,i):
            Sum += Y[j]*L[i][j]
        Y[i] = B[i]-Sum

    if U[n-1][n-1] == 0:
        print('\nNo solution of system of linear equations')
        print('------------------------------------------------------------\n')
        return main()
    X[n-1] = Y[n-1]/U[n-1][n-1]
    for i in range(n-2,-1,-1):
        Sum = 0
        j = n-1
        while j > i:
            Sum += X[j]*U[i][j]
            j=j-1
        X[i] = (Y[i]-Sum)/U[i][j]
    print()
    
    expression2('L',L,n)
    expression2('U',U,n)
    expression1('Y',Y,n)
    expression1('X',X,n)

def expression1(s,x,n):#ฟังก์ชันที่ใช้แสดงค่าในเมทริกซ์ที่มี 1 แถวหรือ 1 หลัก
    print(s+' = ',end = '')
    for i in range(0,n):
        print(str(format(x[i],'20.15f'))+'\n    ',end = '')
    print('')

def expression2(s,x,n):#ฟังก์ชันที่ใช้แสดงค่าในเมทริกซืจตุรัส
    print(s+' = ',end = '')
    for i in range(0,n):
        if i > 0:
            print('    ',end='')
        for j in range(0,n):
            print(format(x[i][j],'20.15f'),'  ',end = '')
        print()
    print()
#==========================Eigenvalues Problems==============================            
def Trn(Z):#ฟังก์ชันที่ใช้ทรานสโพสเมทริกซ์
    R = len(Z)
    C = len(Z[0])
    ZT = []
    for i in range(0,C):
        ZT.append([])
    for i in range(0,C):
        for j in range(0,R):
            ZT[i].append(0)
    
    for j in range(0,R):
        for i in range(0,C):
            ZT[i][j] = Z[j][i]
    return ZT

def multiply_Matrix(x,y):#ฟังก์ชันที่ใช้คูณเมทริกซ์
    z = []
    R1 = len(x)
    R2 = len(y)
    C2 = len(y[0])
    for i in range(0,R1):
        z.append([])
    for i in range(0,R1):
        for j in range(0,C2):
            z[i].append(0)
    
    for i in range(0,R1):
        for j in range(0,C2):
            for k in range(0,R2):
                z[i][j] += x[i][k]*y[k][j]
    return z

def define_A():#ฟังก์ชันที่ใช้กำหนดค่าในเมทริกซ์ A
    n = int(input('Enter scale(rows and colums) of Matrix A : '))
    A = []
    for i in range(0,n):
        A.append([])
    for i in range(0,n):
        for j in range(0,n):
            A[i].append(float(input('Enter value of rows '+str(i+1)+' colums '+str(j+1)+' of Matrix A : ')))
    return A

def Schur(A):#ฟังก์ชันที่ใช้ประมาณค่า Eigenvalues โดย Schur's Theorem
    n = len(A)
    L = []
    for j in range(0,n):
        Sum = 0
        for k in range(0,n):
            Sum += A[j][k]**2
        L.append(Sum)
    Sum = 0
    print('Eigenvalues = ',end='')
    for i in range(0,n):
        Sum += L[i]
        L[i] = L[i]**0.5
        if i == n-1:
            print(L[i],end=' ')
        else:
            print(L[i],end=' , ')
    Lm = Sum**0.5
    print()

def Collatz_Inclusion(A):#ฟังก์ที่ใช้ประมาณค่า Eigenvalues โดย Collatz_Inclusion
    n = len(A)
    x = []
    for i in range(0,n):
        x.append([])
        x[i].append(float(input('Enter the components of initial Eigenvector(Can\'t be zero vector) : ')))
    count = 0
    for k in range(0,n):
        if x[k][0] == 0:
            count += 1
        if count == n:
            print('\nInitial Eigenvector must not be zero vector')
            print('------------------------------------------------------------\n')
            return Collatz_Inclusion(A)
    
    iteration  = int(input('Enter iteration : '))
    y = x
    for i in range(0,iteration):
        x = y
        y = multiply_Matrix(A,x)
        q = []
        q.append(y[0][0]/x[0][0])
        Max = Min = q[0]
        for j in range(1,n):
            q.append(y[j][0]/x[j][0])
            if Max <= q[j]:
                Max = q[j]
            elif Min >= q[j]:
                Min = q[j]
        print('\niteration =',i+1,'  ',Min,'<= Eigenvalue <=',Max)
        print('                 Lrngth =',Max-Min,end='\n')

def Power_Method(A):#ฟังก์ที่ใช้ประมาณค่า Eigenvalues โดย Power Method
    n = len(A)
    x = []
    for i in range(0,n):
        x.append([])
        x[i].append(float(input('Enter the components of initial Eigenvector(Can\'t be zero vector) : ')))
    count = 0
    for k in range(0,n):
        if x[k][0] == 0:
            count += 1
        if count == n:
            print('\nInitial Eigenvector must not be zero vector')
            print('------------------------------------------------------------\n')
            return Power_Method(A)
    
    iteration  = int(input('Enter iteration : '))
    q = []
    delta = []
    y = x
    print('iteration              Quotient            bound for error            Error')
    for i in range(0,iteration):
        x = y
        y = multiply_Matrix(A,x)
        m0 = multiply_Matrix(Trn(x),x)[0][0]
        m1 = multiply_Matrix(Trn(x),y)[0][0]
        m2 = multiply_Matrix(Trn(y),y)[0][0]
        q.append(m1/m0)
        delta.append((m2/m0-q[i]**2)**0.5)

    for j in range(0,iteration):
        print('   j =',j+1,end='      ')
        print(str(format(q[j],'20.15f')),end='     ')
        print(str(format(delta[j],'20.15f')),end='     ')
        print(str(format(q[i]-q[j],'20.15f')))
    print('Eigenvalue is approximately',q[i],end='\n')

def Find_E(A):#ฟังก์ชันที่เรียกใช้ฟังก์ชันหา Eigenvalues แบบต่างๆ
        print('\nFind eigenvalue with')
        print('Press 1 : Schur\'s theorem')
        print('Press 2 : Collatz_Inclusion(CAUTION : elements of matrix A must not all negaitive values)')
        print('Press 3 : Power Method')
        select = int(input('You choose a topic : '))
        print()
        if select == 1:
            Schur(A)
        elif select == 2:
            Collatz_Inclusion(A)
        elif select == 3:
            Power_Method(A)
            
def main():#ฟังก์ชันหลักที่ใช้ควบคุมการทำงานของโปรแกรม
    print('Press 1 : Find the answer to a system of linear equations')
    print('Press 2 : Find Eigenvalues')
    print('Press 0 : Finish')
    select = int(input('You choose a topic : '))
    if select == 1:
        i = 1
        print()
        LU()
        while 1:
            if i>1:
                print('Press 1 : Try again')
                print('Press 9 : Back')
                print('Press 0 : Finish')
                select = int(input('You choose a topic : '))
                if select == 1:
                    print()
                    LU()
                elif select == 9:
                    print('\n============================================================\n')
                    return main()
                elif select == 0:
                    print('Done!')
                    print('============================================================\n')
                    return
            print('------------------------------------------------------------\n')
            i+=1
    elif select == 2:
        i = 1
        print('\nDefine Matrix A')
        A = define_A()
        Find_E(A)
        while 1:
            if i>1:
                print('Press 1 : Try again')
                print('Press 2 : Define new Matrix A and try again')
                print('Press 9 : Back')
                print('Press 0 : Finish')
                select = int(input('You choose a topic : '))
                if select == 1:
                    Find_E(A)
                elif select == 2:
                    print('\n------------------------------------------------------------\n')
                    print('Define Matrix A')
                    A = define_A()
                    Find_E(A)
                elif select == 9:
                    print('\n============================================================\n')
                    return main()
                elif select == 0:
                    print('Done!')
                    print('============================================================\n')
                    return
            print('\n------------------------------------------------------------\n')
            i+=1
    elif select == 0:
        print('Done!')
        print('============================================================\n')
        return
    return main()
main()
