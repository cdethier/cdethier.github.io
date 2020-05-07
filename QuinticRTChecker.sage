import os
import ast

os.chdir("/home/sage/Dropbox/ResearchProjects/ReductionTheory/ComputationalExamples")
f = open("QuinticForms.tex","w")

f.write("\\documentclass[11pt]{article}\r\n")
f.write("\\usepackage{longtable}\r\n")
f.write("\\usepackage{amsmath}\r\n")
f.write("\\textwidth = 7 in\r\n")
f.write("\\oddsidemargin = -0.5 in \r\n")
f.write("\\evensidemargin = 0.0 in \r\n")
f.write("\\headsep = 0.3 in\r\n")
f.write("\\begin{document}\r\n")
f.write("\\title{Quintic Diagonalizable Forms with $m < 255137$}\r\n")
f.write("\\date{\\today}\r\n")
f.write("\\maketitle\r\n")
f.write("\\begin{longtable}{p{6.5cm}|p{1cm}|p{3.5cm}|p{3.5cm}}\r\n")
f.write("Form & $m$ & $F(x,y) = 1$ & $F(x,y) = -1$ \\\\ \\hline \\endhead \r\n")



def mSolver(n,disc):

    if is_even(n):
        if disc > 0:
            m = (disc/n)^(1/(n-1))
        else:
            m = -(-disc/n)^(1/(n-1))
    else:
        m = (disc^2/n^2)^(1/(n-1))

    return m



def jnPossibilities(n, m):

    if is_even(n):
        jnList = [m/((-1)^((n+2)/2)*n)]

    else:
        jnList = [m/(n^2),-m/(n^2)]

    return jnList



def DChinPossibilities(n, jn):

    DChinList = []

    nmod4 = mod(n,4)

    if nmod4 == 0:
        DBound = Integer(n^2*(n-1)^2*jn)
        DInitPoss = divisors(DBound)
        DInitPoss = [D for D in DInitPoss if D^(2/(n-2)) in DInitPoss]
        DInitPoss = DInitPoss + [-D for D in DInitPoss]
        Dmod4 = mod(D,4)
        DPoss = [D for D in DInitPoss if (mod(D,4) == 0 or mod(D,4) == 1)]

    if nmod4 == 1:
        DBound = Integer(n^4*(n-1)^4*jn)
        DInitPoss = divisors(DBound)
        DInitPoss = [D for D in DInitPoss if D^(n-2) in DInitPoss]
        if sign(jn) == -1:
            DInitPoss = [-D for D in DInitPoss]
        DPoss = [D for D in DInitPoss if (mod(D,4) == 0 or mod(D,4) == 1)]

    if nmod4 == 2:
        DBound = Integer(n^2*(n-1)^2*jn)
        DInitPoss = divisors(DBound)
        DInitPoss = [D for D in DInitPoss if D^(2/(n-2)) in DInitPoss]
        DInitPoss = DInitPoss + [-D for D in DInitPoss]
        Dmod4 = mod(D,4)
        DPoss = [D for D in DInitPoss if (mod(D,4) == 0 or mod(D,4) == 1)]

    if nmod4 == 3:
        DBound = Integer(n^4*(n-1)^4*jn)
        DInitPoss = divisors(DBound)
        DInitPoss = [D for D in DInitPoss if D^(n-2) in DInitPoss]
        if sign(jn) == 1:
            DInitPoss = [-D for D in DInitPoss]
        Dmod4 = mod(D,4)
        DPoss = [D for D in DInitPoss if (mod(D,4) == 0 or mod(D,4) == 1)]

    if 1 in DPoss:
        DPoss.remove(1)
    if 4 in DPoss:
        DPoss.remove(4)
    if -4 in DPoss:
        DPoss.remove(-4)
    if 8 in DPoss:
        if mod(m,12800) != 0:
            DPoss.remove(8)
    if -8 in DPoss:
        if mod(m,12800) != 0:
            DPoss.remove(-8)


    for D in DPoss:
        if is_even(n):
            Chin = jn/(D^(n/2))
            DChinList = DChinList + [[D,Chin,m]]

        else:
            if is_square(jn/D^n):
                Chin = sqrt(jn/(D^n))
                DChinList = DChinList + [[D,Chin]]
                DChinList = DChinList + [[D,-Chin]]

    DChinList = list(Set(DChinList))

    return DChinList



def QuadraticPossibilities(D):

    if D.is_square():
        QuadList = []
        b = sqrt(D)
        for a in range(1,b):
            if gcd(a,b) == 1:
                QuadList = QuadList + [BinaryQF([a,b,0])]

    else:
        QuadList = list(set(BinaryQF_reduced_representatives(D)))
        QuadList = [quad for quad in QuadList if gcd(gcd(quad[0],quad[1]),quad[2]) == 1]

    return QuadList



def nEvenDiagonalReduction(n,m):
    #print 'check1.5'
    K.<x,y> = PolynomialRing(QQ)
    Forms = set([])
    aPossibilities = divisors(m/n)
    aPossibilities = aPossibilities + [-a for a in aPossibilities]
    for a in aPossibilities:
        b = (m*(-1)^(n/2))/(n*a)
        diagForm = a*x^n - b*y^n
        Forms = Forms.union(Set([diagForm]))

    return list(Forms)
    


def nOddDiagonalReduction(n,m):
    #print 'check1'
    K.<x,y> = PolynomialRing(QQ)
    Forms = set([])
    aSquarePossibilities = divisors(m/n^2)
    aPossibilities = [a for a in aSquarePossibilities if a^2 in aSquarePossibilities]
    for a in aPossibilities:
        bSq = m/(n^2*a^2)
        if bSq > 0:
            if is_square(bSq):
                b = sqrt(bSq)
                diagForm = a*x^n - b*y^n
                Forms = Forms.union(set([diagForm]))
                b = -sqrt(bSq)
                diagForm = a*x^n - b*y^n
                Forms = Forms.union(set([diagForm]))

    return list(Forms)



def D4Reduction(n,m):
    #print "check6"
    K.<x,y> = PolynomialRing(QQ)
    Forms = set([])
    aSquarePossibilities = divisors(m/25600)
    aPossibilities = [a for a in aSquarePossibilities if a^2 in aSquarePossibilities]
    aPossibilities = aPossibilities + [-a for a in aPossibilities]
    for a in aPossibilities:
        bSq = m/(25600*a^2)
        if bSq > 0:
            if is_square(bSq):
                b = sqrt(bSq)
                if mod(a,2) == mod(b,2):
                    diagForm = expand(a*x^5 + b*(x+2)^5)
                    Forms = Forms.union(set([diagForm]))
                    b = -b
                    diagForm = expand(a*x^5 + b*(x+2)^5)
                    Forms = Forms.union(set([diagForm]))

    return list(Forms)



def Dneg4Reduction(n,m):
    #print "check5"
    K.<x,y> = PolynomialRing(QQ[I])
    Forms = set([])
    zSquarePossibilities = divisors(m/1600)
    zPossibilities = [z for z in zSquarePossibilities if z^2 in zSquarePossibilities]
    for z in zPossibilities:
        for aSq in range(0,z+1):
            if is_square(aSq):
                bSq = z - aSq
                if is_square(bSq):
                    for a in [-sqrt(aSq),sqrt(aSq)]:
                        for b in [-sqrt(bSq),sqrt(bSq)]:
                            diagForm = expand((1/2)*(a+b*I)*(x+I)^5 + (1/2)*(a-b*I)*(x-I)^5)
                            Forms = Forms.union(set([diagForm]))

    return list(Forms)



def DSquareReduction(n,disc,Chin,A,B):
    K.<x,y> = PolynomialRing(QQ)
    D = B^2
    #print 'check2'

    Forms = set([])

    beta1 = 0
    delta1 = -B/A

    nalpha1bound = n^2*(n-1)^2*B^2*Chin*A^n

    if nalpha1bound.is_integer():
        nalpha1List = divisors(nalpha1bound)
        if is_even(n):
            nalpha1List = nalpha1List + [-nalpha1 for nalpha1 in nalpha1List]
        for nalpha1 in nalpha1List:
            ngamma1 = nalpha1bound/nalpha1
            alpha1 = nalpha1/(sqrt(D)*n*(n-1))
            gamma1 = ngamma1/(sqrt(D)*n*(n-1))
            if (alpha1 - gamma1).is_integer():
                diagForm = expand(alpha1*(x - beta1*y)^n - gamma1*(x - delta1*y)^n)
                coefficients = diagForm.coefficients()
                coefficientsBool = [coefficient.is_integer() for coefficient in coefficients]
                if all(coefficientsBool):
                    disctest = Integer((diagForm.discriminant(x))(y=1))
                    if disctest == disc:
                        Forms = Forms.union(set([diagForm]))
                        #print diagForm

    return list(Forms)



def DNegativeReduction(n,disc,Chin,D,A,B,C):
    #print 'check3'
    
    Forms = set([])

    d = squarefree_part(D)
    L.<a> = QuadraticField(d)
    K.<x,y> = PolynomialRing(L)
    f = a.coordinates_in_terms_of_powers()
    sqrtD = a*sqrt(D/d)
    beta1 = (-B + sqrtD)/(2*A)
    delta1 = (-B - sqrtD)/(2*A)
    alpha1Norm = n^2*(n-1)^2*D*(Chin)*A^n
    alpha1List = list(gp("bnfisintnorm(bnfinit(a^2 + " + str(-d) + ")," + str(alpha1Norm) + ")"))
    alpha1List = alpha1List + [-alpha1 for alpha1 in alpha1List]
    #print alpha1List
    for alpha1 in alpha1List:
         alpha1 = eval(str(alpha1))
         alpha1 = alpha1/(n*(n-1)*sqrtD)
         gamma1 = -f(alpha1)[0] + f(alpha1)[1]*a
         diagForm = expand(alpha1*(x-beta1*y)^n - gamma1*(x-delta1*y)^n)
         coefficients = diagForm.coefficients()
         coefficientsBool = [coefficient.is_integer() for coefficient in coefficients]
         if all(coefficientsBool):
             disctest = Integer((diagForm.discriminant(x))(y=1))
             if disctest == disc:
                 Forms = Forms.union(Set([diagForm]))
                 #print diagForm
    #print 'check3.5'
    return list(Forms)



def DPositiveReduction(n,disc,Chin,D,A,B,C):
    #print 'check4'

    Forms = set([])

    d = squarefree_part(D)
    L.<a> = QuadraticField(d)
    K.<x,y> = PolynomialRing(L)
    sqrtD = a*sqrt(D/d)
    beta1 = (-B + sqrtD)/(2*A)
    delta1 = (-B - sqrtD)/(2*A)
    UK = UnitGroup(L)
    u = (UK.gens_values())[1]
    f = a.coordinates_in_terms_of_powers()
    check = False
    s = 1
    mu = f(beta1)[0]
    nu = f(beta1)[1]
    while check == False:
        c = f(u^s)[0]
        b = f(u^s)[1]
        coordinates = [c - (b*mu)/nu, (b*mu^2)/nu - b*nu*d, -b/nu, (b*mu)/nu + c]
        coordinatesBool = [coordinate.is_integer() for coordinate in coordinates]
        if all(coordinatesBool):
            check = True
        else:
            s = s + 1
    s = s*n
    N = n^2*(n-1)^2*D*(Chin)*A^n
    FundamentalSolutions = list(gp("bnfisintnorm(bnfinit(a^2-" + str(d) + ")," + str(N) + ")"))
    FS = [eval(str(fu))/(n*(n-1)*sqrtD) for fu in FundamentalSolutions]
    for fu in FS:
        for stemp in range(0,s):
            fuus = fu*u^(stemp)
            fuusbarneg = -f(fuus)[0]+f(fuus)[1]*a
            diagForm = expand(fuus*(x-beta1*y)^n - fuusbarneg*(x-delta1*y)^n)
            coefficients = diagForm.coefficients()
            coefficientsBool = [coefficient.is_integer() for coefficient in coefficients]
            if all(coefficientsBool):
                disctest = Integer((diagForm.discriminant(x))(y=1))
                if disctest == disc:
                    Forms = Forms.union(Set([diagForm]))
                    #print diagForm

    return list(Forms)



def ReductionAlgorithm(n,disc): 

    Forms = []
    m = mSolver(n,disc)

    if is_even(n):
        if mod(m,n) == 0:
            Forms = Forms + nEvenDiagonalReduction(n,m)
    else:
        if mod(m,n^2) == 0:
            Forms = Forms + nOddDiagonalReduction(n,m)

    if mod(m,25600) == 0:
        Forms = Forms + D4Reduction(n,m)

    if mod(m,1600) == 0:
        Forms = Forms + Dneg4Reduction(n,m)

    DChinList = []
    jnList = jnPossibilities(n,m)
    for jn in jnList:
        DChinList = DChinList + DChinPossibilities(n,jn)

    for DChin in DChinList:

        D = DChin[0]
        Chin = DChin[1]
        #print 'D = ' + str(D)
        #print 'Chin = ' + str(Chin)
        #print ''

        for quad in QuadraticPossibilities(D):
            if D > 0:
                if is_square(D):
                    Forms = Forms + DSquareReduction(n,disc,Chin,quad[0],quad[1])
                else:
                    Forms = Forms + DPositiveReduction(n,disc,Chin,D,quad[0],quad[1],quad[2])
            else:
                Forms = Forms + DNegativeReduction(n,disc,Chin,D,quad[0],quad[1],quad[2])

    return Forms

for m in range(1,255137):
    print "m = " + str(m)
    for form in ReductionAlgorithm(5,5*m^2):
        poly = form(y = 1)
        
        gp("tnf = thueinit(" + str(poly) + ")")
        solposString = str(gp("thue(tnf,1)"))
        solnegString = str(gp("thue(tnf,-1)"))

        
        solpos = ast.literal_eval(solposString)
        solneg = ast.literal_eval(solnegString)


        j = 1
        while j < len(solpos):
            k = 0
            while k < j:
                if solpos[j][0] == -solpos[k][0] and solpos[j][1] == -solpos[k][1]:
                    del solpos[j]
                    k = k - 1
                    j = j - 1
                k = k + 1
            j = j + 1
        j = 1
        while j < len(solneg):
            k = 0
            while k < j:
                if solneg[j][0] == -solneg[k][0] and solneg[j][1] == -solneg[k][1]:
                    del solneg[j]
                    k = k - 1
                    j = j - 1
                k = k + 1
            j = j + 1


        f.write("$" + latex(form) + "$ &" + str(m) + "&")
        if len(solpos) != 0:
            f.write("$")
            for j in range(0,len(solpos)):
                f.write("(" + str(-solpos[j][0]) + "," + str(-solpos[j][1]) + ") \\text{ }")
            f.write("$")
        f.write(" & ")
        if len(solneg) != 0:
            f.write("$")
            for j in range(0,len(solneg)):
                f.write("(" + str(-solneg[j][0]) + "," + str(-solneg[j][1]) + ") \\text{ }")
            f.write("$")
        f.write(" \\\\ \r\n")

f.write("\\end{longtable}\r\n")
f.write('\\end{document}')
f.close()


os.system("pdflatex QuinticForms.tex")