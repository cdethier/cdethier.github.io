import math
import ast
import os
R.<x,y> = PolynomialRing(QQ)



#####

#To run this code, you should insert a directory of your choosing here:

os.chdir("")

#####



f = open("FormsAA.tex","w")

#f is the file Forms.tex where we will record our forms and their solutions.
#The following sets up the tex formatting for Forms.tex:

f.write("\\documentclass[11pt]{article}\r\n")
f.write("\\usepackage{longtable}\r\n")
f.write("\\usepackage{amsmath}\r\n")
f.write("\\textwidth = 7 in\r\n")
f.write("\\oddsidemargin = -0.0 in \r\n")
f.write("\\evensidemargin = 0.0 in \r\n")
f.write("\\headsep = 0.3 in\r\n")
f.write("\\begin{document}\r\n")
f.write("\\title{Forms with $J_F = 0$ and $0 > I_F \\geq -3000$}\r\n")
f.write("\\date{\\today}\r\n")
f.write("\\maketitle\r\n")
f.write("\\begin{longtable}{p{6.5cm}|p{1cm}|p{3.5cm}|p{3.5cm}}\r\n")
f.write("Form & $I_F$ & $F(x,y) = 1$ & $F(x,y) = -1$ \\\\ \\hline \\endhead \r\n")

#Statistics is a matrix whose first coordinate represents the number of solutions to F(x,y) = 1
#and whose second coordinate represents the number of solutions to F(x,y) = -1.
#Each entry is the number of forms with that many solutions to F(x,y) = 1 and F(x,y) = -1.

Statistics = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

#formCount - the total number of forms generated.

formCount = 0

#The following loops generate our forms. For a given I, we go through all possibilities of a,b,c.
#The bounds on a,b,c are those given in Cremona.
#From a,b,c we can solve for H,R,d,e checking for integrality at each step.

for Itemp in range(1,1001):
    I = -3*Itemp

    #The bounds on a are from Cremona, equation (59) 
    
    aBound = (2*math.sqrt(-I))/(3*math.sqrt(3))
    L = range(int(math.ceil(-aBound)),int(math.floor(aBound)+1))
    L.remove(0)

    for a in L:

        Ba = (4/3)*math.sqrt(-I)*math.sqrt(-4*I - 27*a^2)

        Hmin = max(4*I/3,-Ba)
        Hmax = 0
            
        #These H bounds translate to c bounds under the identity H = 8ac - 3b^2 given in (19) of Cremona.
        #As Cremona notes, finding c first removes the integrality check on c, as H is always integral
        #provided that a,b,c are.

        #The bounds on b are -2|a| < b leq 2|a|, they can be found in section 4.6 of Cremona

        for b in range(-2*abs(a) + 1, 2*abs(a)+1):

            #The bounds on H are from (60) in Cremona
            
            cminInitial = (Hmin+3*b**2)/(8*a)
            cmaxInitial = (Hmax+3*b**2)/(8*a)

            cmin = min(cminInitial,cmaxInitial)
            cmax = max(cminInitial,cmaxInitial)

            for c in range(int(math.ceil(cmin)),int(math.floor(cmax))+1):

                H = 8*a*c - 3*b**2
                
                #Compute R using the seminvariant syzygy (22) from Cremona:

                if (H**3-48*I*a**2*H)/(-27) >= 0:

                    Rpossibilities = [-math.sqrt((H**3-48*I*a**2*H)/(-27)),math.sqrt((H**3-48*I*a**2*H)/(-27))]
 
                    if (H**3-48*I*a**2*H) == 0:
                        Rpossibilities = [0]

                    for R in Rpossibilities:

                        if R.is_integer() == True:
                  
                            #Then d using (20), the definition of R:
                
                            d = (R - b**3 + 4*a*b*c)/(8*a**2)

                            if d.is_integer() == True:
                    
                                #and finally e using (15), the definition of I:
                    
                                e = (I+3*b*d-c**2)/(12*a)

                                if e.is_integer() == True:
                                    d = int(d)
                                    e = int(e)

                                    #Checking whether a,b,c,d,e have a common factor reduces
                                    #repeated forms up to SL2(Z) action. (although repeats may still occur)

                                    if gcd(a,gcd(b,gcd(c,gcd(d,e)))) == 1:

                                        formCount = formCount + 1
                                        form = a*x^4 + b*x^3*y + c*x^2*y^2 + d*x*y^3 +e*y^4
                                
                                        #poly is the unhomogonized form. It is needed because PARI prefers 
                                        #solving Thue equations with polynomials isntead of forms.
                                 
                                        poly = form(y = 1)
                                  
                                        #The following calls PARI (specifically GP) and finds the solutions to F(x,y) = 1 and -1
                                
                                        gp("tnf = thueinit(" + str(poly) + ")")
                                        solposString = str(gp("thue(tnf,1)"))
                                        solnegString = str(gp("thue(tnf,-1)"))
                                
                                        #Calling GP returns a string formatted like a list, this reads the string as a list:
                                  
                                        solpos = ast.literal_eval(solposString)
                                        solneg = ast.literal_eval(solnegString)
                                
                                        #Next we get rid of repeated solutions, as we consider (x,y) and (-x,-y) the same solution.
                                        #The j and kth solutions are compared, the jth is removed if they are identical.
                                
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
                                   
                                        #We increase the entry of statistics corresponding to the number of solutions found:
                                
                                        Statistics[len(solpos)][len(solneg)] = Statistics[len(solpos)][len(solneg)] + 1
                                
                                        #To check for errors, we compute I and J again (J should be 0) then compare Itest with expected I                              
                                        Itest = 12*a*e - 3*b*d + c**2
                                        Jtest = 72*a*c*e + 9*b*c*d - 27*a*d**2 - 27*e*b**2 - 2*c**3
                                        if I != Itest or Jtest != 0:
                                            print "Uh oh"
                                    
                                        #Recording the form in Forms.tex. Also serves as a check if F(x,y) = 1 or F(x,y) = -1
                                        #has more than three solutions, as the indexes of Statistics only go up to 3.
                                        #This only activates when F(x,y) has more than 
                                        #zero solutions to avoid extra $'s in tex formatting.   
                                    
                                        f.write("$" + latex(form) + "$ &" + str(Itest) + "&")
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
                                  
                                        #Print I to record progress.
                                
                                        print I

f.write("\\end{longtable}\r\n")
f.write(str(formCount) + ' forms generated\r\n')
f.write('\\end{document}')
f.close()

#g is the file FormStatistics.tex where we record statistics about the forms generated.
#This was not included in Forms.tex because Forms.tex can get unwieldy to compile and browse through,
#particularly if the range of I is increased substantially.

g = open("FormStatisticsAA.tex","w")

g.write("\\documentclass[11pt]{article}\r\n")
g.write("\\usepackage{longtable}\r\n")
g.write("\\usepackage{amsmath}\r\n")
g.write("\\textwidth = 7 in\r\n")
g.write("\\oddsidemargin = -0.0 in \r\n")
g.write("\\evensidemargin = 0.0 in \r\n")
g.write("\\headsep = 0.3 in\r\n")
g.write("\\begin{document}\r\n")
g.write("\\title{Forms with $J_F = 0$ and $0 > I_F \\geq -3000$ Statistics}\r\n")
g.write("\\date{\\today}\r\n")
g.write("\\maketitle\r\n")

g.write("\\begin{tabular}{l|l|l}  \r\n")
g.write("$F(x,y) = 1$ & $F(x,y) = -1$ & Number of Forms \\\\ \hline \hline \r\n")
for pos in range(0,len(Statistics)):
    for neg in range(0,len(Statistics[pos])):
        g.write(str(pos) + " & " + str(neg) + " & " + str(Statistics[pos][neg]) + "\\\\ \r\n")
g.write("\\end{tabular}\r\n")
g.write("\\end{document}")
g.close()

os.system("pdflatex FormsAA.tex")
os.system("pdflatex FormStatisticsAA.tex")

