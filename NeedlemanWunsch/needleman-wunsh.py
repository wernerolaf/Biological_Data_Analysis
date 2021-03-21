import argparse

def Needleman_Wunsh(CONFIG=None,SEQ1="",SEQ2="",GP=-2,DIFF=-5,SAME=5,MAX_SEQ_LENGTH=-1,MAX_NUMBER_OF_PATH=1,OUTPUT="solutions.txt"):

    def read_config(CONFIG=None):
        variables={"SEQ1","SEQ2","GP","DIFF","SAME","MAX_SEQ_LENGTH","MAX_NUMBER_OF_PATH","OUTPUT"}
        toUpdate={}
        if CONFIG==None:
            return
        f=open(CONFIG, "r")
        for variable in f:
            temp=variable.split("=",1)
            if(len(temp)==2):
                name=temp[0]
                value=temp[1].strip("\n")
                if name in variables and len(value)>0:
                    if name[0:3]=="SEQ":
                        toUpdate[name]=value
                    elif name=="OUTPUT":
                        toUpdate[name]=value
                    else:
                        toUpdate[name]=int(value)
        f.close()
        return(toUpdate)
           
    
    
    def create_paths(MatrixOfDirections,MAX_NUMBER_OF_PATH,SEQ1,SEQ2):
        rows=len(MatrixOfDirections)
        columns=len(MatrixOfDirections[0])
        solutions=[]
        foundSolutions=0
        def _traverse(x,y,text1,text2):
            nonlocal foundSolutions
            if foundSolutions>=MAX_NUMBER_OF_PATH:
                return
            elif y==-1:
                solutions.append([(x+1)*"_"+text1,SEQ2[:x+1]+text2])
                foundSolutions+=1
            elif x==-1:
                solutions.append([SEQ1[:y+1]+text1,(y+1)*"_"+text2])
                foundSolutions+=1
            else:
                if(MatrixOfDirections[x][y] % 2 == 0):
                    _traverse(x,y-1,SEQ1[y]+text1,"_"+text2)

                if(MatrixOfDirections[x][y] % 5 == 0):
                    _traverse(x-1,y,"_"+text1,SEQ2[x]+text2)

                if(MatrixOfDirections[x][y] % 3 == 0):
                    _traverse(x-1,y-1,SEQ1[y]+text1,SEQ2[x]+text2)

        _traverse(rows-1,columns-1,"","")
        return(solutions)
    
    def write_solutions(solutions,OUTPUT):
        f=open(OUTPUT, "w")
        f.write("Score:{}\n".format(MatrixOfValues[len(SEQ2)][len(SEQ1)]))
        for solution in solutions:
            f.write(solution[0])
            f.write("\n")
            f.write(solution[1])
            f.write("\n\n")

        f.close()
    
    def calculate_values_directions(MatrixOfValues,MatrixOfDirections,SEQ1,SEQ2,GP,DIFF,SAME):
        for x in range(len(SEQ2)):
            for y in range(len(SEQ1)):
                left=MatrixOfValues[x+1][y]+GP
                up=MatrixOfValues[x][y+1]+GP
                if(SEQ2[x]==SEQ1[y]):
                    slant=MatrixOfValues[x][y]+SAME
                else:
                    slant=MatrixOfValues[x][y]+DIFF

                newValue=max(left,up,slant)
                MatrixOfValues[x+1][y+1]=newValue
                if(newValue==left):
                    MatrixOfDirections[x][y]*=2
                if(newValue==slant):
                    MatrixOfDirections[x][y]*=3
                if(newValue==up):
                    MatrixOfDirections[x][y]*=5
    
    #checking arguments
    conf=None
    try:
        conf=read_config(CONFIG)
    except:
        print("Reading config failed")
    if conf!=None:
        for variable in conf.keys():
            if variable=="SEQ1":
                SEQ1=conf[variable]
            elif variable=="SEQ2":
                SEQ2=conf[variable]
            elif variable=="GP":
                GP=conf[variable]
            elif variable=="DIFF":
                DIFF=conf[variable]
            elif variable=="SAME":
                SAME=conf[variable]
            elif variable=="MAX_SEQ_LENGTH":
                MAX_SEQ_LENGTH=conf[variable]
            elif variable=="MAX_NUMBER_OF_PATH":
                MAX_NUMBER_OF_PATH=conf[variable]
            elif variable=="OUTPUT":
                OUTPUT=conf[variable]
    
    if(MAX_SEQ_LENGTH<1):
        print(MAX_SEQ_LENGTH)
        print("Set MAX_SEQ_LENGTH to positive integer")
        return
    elif(MAX_SEQ_LENGTH<max(len(SEQ1),len(SEQ2))):
        print("MAX_SEQ_LENGTH reached")
        return
    
    #setting up needed matrices
    MatrixOfValues = [[0 for x in range(len(SEQ1)+1)] for y in range(len(SEQ2)+1)] 
    for i in range(len(SEQ1)+1):
        MatrixOfValues[0][i]=i*GP
    
    for i in range(len(SEQ2)+1):
        MatrixOfValues[i][0]=i*GP
    
    MatrixOfDirections=[[1 for x in range(len(SEQ1))] for y in range(len(SEQ2))]
    
    #calculating values and updating directions
    calculate_values_directions(MatrixOfValues,MatrixOfDirections,SEQ1,SEQ2,GP,DIFF,SAME)
    #creating solution based on MatrixOfDirections
    solutions=create_paths(MatrixOfDirections,MAX_NUMBER_OF_PATH,SEQ1,SEQ2)
    #writing solutions
    write_solutions(solutions,OUTPUT)
    print("Done")
    
    
parser = argparse.ArgumentParser()
parser.add_argument("-a","--SEQ1", default="",type=str, nargs='?')
parser.add_argument("-b","--SEQ2", default="",type=str, nargs='?')
parser.add_argument("-c","--CONFIG", default=None,type=str, nargs='?')
parser.add_argument("-o","--OUTPUT", default="solutions.txt",type=str, nargs='?')
parser.add_argument("-g","--GP", default=-2,type=int, nargs='?')
parser.add_argument("-d","--DIFF",default=-5,type=int, nargs='?')
parser.add_argument("-s","--SAME", default=5,type=int, nargs='?')
parser.add_argument("-m","--MAX_SEQ_LENGTH", default=-1,type=int, nargs='?')
parser.add_argument("-p","--MAX_NUMBER_OF_PATH",default=1,type=int, nargs='?')
args=parser.parse_args()

Needleman_Wunsh(CONFIG=args.CONFIG,SEQ1=args.SEQ1,SEQ2=args.SEQ2,GP=args.GP,DIFF=args.DIFF,SAME=args.SAME,MAX_SEQ_LENGTH=args.MAX_SEQ_LENGTH,MAX_NUMBER_OF_PATH=args.MAX_NUMBER_OF_PATH,OUTPUT=args.OUTPUT)

